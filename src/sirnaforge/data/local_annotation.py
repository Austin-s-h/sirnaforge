"""Offline transcript structure from pinned GENCODE GTF or NCBI RefSeq GTF/GFF3 files."""

from __future__ import annotations

import asyncio
import gzip
import hashlib
import os
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO
from urllib.parse import unquote

from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.data.base import AbstractTranscriptAnnotationClient
from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcript_annotation import EnsemblTranscriptModelClient
from sirnaforge.models.transcript_annotation import Interval, TranscriptAnnotation, TranscriptAnnotationBundle

_GTF_ATTRIBUTE = re.compile(r'([^\s;]+)\s+"([^"]*)"\s*;?')
_TRANSCRIPT_FEATURES = {"transcript", "mRNA", "ncRNA", "lnc_RNA", "rRNA", "tRNA", "pseudogenic_transcript"}


def stable_id(identifier: str) -> str:
    """Remove a numeric accession version while retaining PAR and remapping suffixes."""
    return re.sub(r"\.\d+(?=_|$)", "", identifier)


def identifier_matches(query: str, identifier: str) -> bool:
    """Match explicit versions exactly, allowing a versionless accession query."""
    return query == identifier or (stable_id(query) == query and query == stable_id(identifier))


def open_annotation(path: Path) -> TextIO:
    """Open plain or gzip-compressed annotation and FASTA files as text."""
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


@dataclass(frozen=True)
class AnnotationFeature:
    """One feature with 1-based inclusive coordinates and parsed attributes."""

    chrom: str
    kind: str
    start: int
    end: int
    strand: int
    attributes: dict[str, list[str]]

    def value(self, *names: str) -> str:
        """Return the first nonempty value for the requested attribute names."""
        return next((self.attributes[name][0] for name in names if self.attributes.get(name)), "")

    def interval(self) -> Interval:
        """Retain the source coordinate convention without conversion."""
        return Interval(seq_region_name=self.chrom, start=self.start, end=self.end, strand=self.strand)


def iter_features(path: Path) -> Iterator[AnnotationFeature]:
    """Read GTF or GFF3, retaining repeated GTF tags and decoded GFF attributes."""
    with open_annotation(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"{path}:{number}: expected nine GTF/GFF3 columns")
            chrom, _, kind, start, end, _, strand, _, raw = columns
            if kind not in _TRANSCRIPT_FEATURES | {"gene", "exon", "CDS"}:
                continue
            if strand not in {"+", "-"}:
                raise ValueError(f"{path}:{number}: transcript features require a strand")
            attributes: dict[str, list[str]] = {}
            if re.match(r"[^\s=;]+=", raw):
                for item in raw.split(";"):
                    key, separator, value = item.partition("=")
                    if separator:
                        attributes[key] = [unquote(v) for v in value.split(",")]
            else:
                for key, value in _GTF_ATTRIBUTE.findall(raw):
                    attributes.setdefault(key, []).append(value)
            feature = AnnotationFeature(chrom, kind, int(start), int(end), 1 if strand == "+" else -1, attributes)
            if not 1 <= feature.start <= feature.end:
                raise ValueError(f"{path}:{number}: invalid 1-based interval")
            yield feature


def _gene_identifier(feature: AnnotationFeature) -> str:
    for cross_reference in feature.attributes.get("Dbxref", []):
        if cross_reference.startswith("GeneID:"):
            return cross_reference.partition(":")[2]
    return feature.value("gene_id", "ID")


def _overlaps(feature: AnnotationFeature, region: tuple[str, int, int]) -> bool:
    chrom, start, end = region
    # RefSeq contig accession names remain distinct; only a literal chr prefix is normalized.
    return (
        feature.chrom.removeprefix("chr") == chrom.removeprefix("chr") and feature.start <= end and feature.end >= start
    )


class LocalTranscriptModelClient(AbstractTranscriptAnnotationClient):
    """Query local GENCODE or RefSeq annotation without contacting any REST endpoint.

    Queries stream the file in up to three passes, storing genes and selected transcript models only.
    Files must be from the declared species and assembly; explicit requests must match exactly.
    """

    def __init__(self, annotation_path: str | Path, *, provider: str, species: str, assembly: str):
        """Bind a provider to a user-selected annotation file and reference identity."""
        super().__init__()
        if provider not in {"gencode", "refseq"}:
            raise ValueError("Local annotation provider must be gencode or refseq")
        if not assembly.strip():
            raise ValueError("An explicit assembly is required for local annotations")
        self.path = Path(annotation_path).expanduser().resolve()
        if not self.path.is_file():
            raise ValueError(f"Annotation file does not exist: {self.path}")
        self.provider = f"{provider}_file"
        self.species = normalize_species_name(species)
        self.assembly = assembly
        self.base_url = self.path.as_uri()
        self.annotation_version: str | None = None
        self._validate_source_header()

    def _validate_source_header(self) -> None:
        declared_builds = set()
        declared_versions = set()
        with open_annotation(self.path) as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                declared = re.match(r"#!genome-build(?:-accession)?\s+(.+)", line.strip())
                if declared:
                    declared_builds.add(declared[1].removeprefix("NCBI_Assembly:"))
                version = re.match(r"#!annotation-version\s+(.+)", line.strip())
                if version:
                    declared_versions.add(version[1].strip())
        if declared_builds and self.assembly not in declared_builds:
            raise ValueError(f"File declares assembly {sorted(declared_builds)}, configured {self.assembly}")
        if len(declared_versions) > 1:
            raise ValueError(f"File declares conflicting annotation versions {sorted(declared_versions)}")
        if declared_versions:
            self.annotation_version = declared_versions.pop()

    def source_metadata(self) -> dict[str, str]:
        """Return the exact local source identity for workflow provenance."""
        digest = hashlib.sha256()
        with self.path.open("rb") as handle:
            while block := handle.read(1024 * 1024):
                digest.update(block)
        metadata = {
            "provider": self.provider,
            "endpoint": self.base_url,
            "assembly": self.assembly,
            "species": self.species,
            "sha256": digest.hexdigest(),
        }
        if self.annotation_version:
            metadata["annotation_version"] = self.annotation_version
        return metadata

    def _validate_reference(self, species: str, reference: ReferenceChoice) -> None:
        if normalize_species_name(species) != self.species:
            raise ValueError(f"Annotation species is {self.species}, requested {species}")
        if not reference.enabled or reference.value != self.assembly:
            raise ValueError(f"Annotation assembly is {self.assembly}, requested {reference.value}")

    def _select(self, ids: list[str], regions: list[tuple[str, int, int]]) -> dict[str, TranscriptAnnotation]:
        if not ids and not regions:
            return {}
        genes: dict[str, AnnotationFeature] = {}
        selected: dict[str, AnnotationFeature] = {}
        for feature in iter_features(self.path):
            if feature.kind == "gene":
                for key in {feature.value("ID"), feature.value("gene_id"), _gene_identifier(feature)} - {""}:
                    genes[key] = feature
            elif feature.kind in _TRANSCRIPT_FEATURES:
                aliases = [
                    feature.value("transcript_id", "Name", "ID"),
                    feature.value("gene_id"),
                    feature.value("gene_name", "gene"),
                    *feature.attributes.get("Parent", []),
                ]
                if any(identifier_matches(query, alias) for query in ids for alias in aliases if alias) or any(
                    _overlaps(feature, region) for region in regions
                ):
                    self._select_feature(selected, feature)
        # Gene IDs/symbols in GFF3 may be present only on the parent feature.
        wanted_genes = {
            key
            for key, gene in genes.items()
            if any(
                identifier_matches(query, alias)
                for query in ids
                for alias in [key, gene.value("gene_name", "gene", "Name"), _gene_identifier(gene)]
                if alias
            )
        }
        if wanted_genes:
            for feature in iter_features(self.path):
                if feature.kind in _TRANSCRIPT_FEATURES and (
                    feature.value("gene_id") in wanted_genes
                    or wanted_genes.intersection(feature.attributes.get("Parent", []))
                ):
                    self._select_feature(selected, feature)
        models = self._models(selected, genes)
        return self._attach_structure(models)

    @staticmethod
    def _select_feature(selected: dict[str, AnnotationFeature], feature: AnnotationFeature) -> None:
        identifier = feature.value("ID", "transcript_id")
        if not identifier:
            raise ValueError("Transcript feature lacks an identifier")
        if identifier in selected and selected[identifier] != feature:
            raise ValueError(f"Conflicting transcript features: {identifier}")
        selected[identifier] = feature

    def _attach_structure(self, models: dict[str, TranscriptAnnotation]) -> dict[str, TranscriptAnnotation]:
        if not models:
            return {}
        for feature in iter_features(self.path):
            if feature.kind not in {"exon", "CDS"}:
                continue
            for parent in feature.attributes.get("Parent", []) or [feature.value("transcript_id")]:
                if parent in models:
                    model = models[parent]
                    if feature.chrom != model.seq_region_name or feature.strand != model.strand:
                        raise ValueError(f"Feature locus disagrees with transcript {model.transcript_id}")
                    if feature.start < model.start or feature.end > model.end:
                        raise ValueError(f"Feature extends beyond transcript {model.transcript_id}")
                    (model.exons if feature.kind == "exon" else model.cds).append(feature.interval())
        output = {}
        for model in models.values():
            model.exons = sorted(set(model.exons), key=lambda item: item.start, reverse=model.strand == -1)
            model.cds = sorted(set(model.cds), key=lambda item: item.start, reverse=model.strand == -1)
            if model.transcript_id in output:
                raise ValueError(f"Ambiguous transcript accession in annotation: {model.transcript_id}")
            output[model.transcript_id] = model
        return output

    def _models(
        self, selected: dict[str, AnnotationFeature], genes: dict[str, AnnotationFeature]
    ) -> dict[str, TranscriptAnnotation]:
        models = {}
        for parent, feature in selected.items():
            gene_key = feature.value("gene_id", "Parent")
            gene = genes.get(gene_key)
            transcript_id = feature.value("transcript_id", "Name", "ID")
            if not transcript_id or not gene_key:
                raise ValueError(f"Transcript lacks a stable identifier or parent gene in {self.path}")
            models[parent] = TranscriptAnnotation(
                transcript_id=transcript_id,
                transcript_name=feature.value("transcript_name") or None,
                is_canonical=bool({"Ensembl_canonical", "MANE_Select"}.intersection(feature.attributes.get("tag", []))),
                gene_id=_gene_identifier(gene) if gene else gene_key,
                symbol=feature.value("gene_name", "gene")
                or (gene.value("gene_name", "gene", "Name") if gene else None),
                biotype=feature.value("transcript_type", "transcript_biotype", "gene_biotype")
                or (gene.value("gene_type", "gene_biotype") if gene else None)
                or None,
                seq_region_name=feature.chrom,
                start=feature.start,
                end=feature.end,
                strand=feature.strand,
                gene_interval=gene.interval() if gene else None,
                provider=self.provider,
                endpoint=self.base_url,
                reference_choice=self.assembly,
            )
        return models

    def _bundle(
        self, ids: list[str], regions: list[tuple[str, int, int]], reference: ReferenceChoice
    ) -> TranscriptAnnotationBundle:
        models = self._select(ids, regions)
        unresolved = [
            query
            for query in ids
            if not any(
                identifier_matches(query, model.transcript_id)
                or identifier_matches(query, model.gene_id)
                or query == model.symbol
                for model in models.values()
            )
        ]
        return TranscriptAnnotationBundle(transcripts=models, unresolved=unresolved, reference_choice=reference)

    async def fetch_by_ids(
        self, ids: list[str], *, species: str, reference: ReferenceChoice
    ) -> TranscriptAnnotationBundle:
        """Resolve transcript accessions, gene IDs, or symbols from a local source."""
        self._validate_reference(species, reference)
        return await asyncio.to_thread(self._bundle, list(dict.fromkeys(ids)), [], reference)

    async def fetch_by_regions(
        self, regions: list[str], *, species: str, reference: ReferenceChoice
    ) -> TranscriptAnnotationBundle:
        """Resolve transcript spans overlapping 1-based inclusive source-coordinate regions."""
        self._validate_reference(species, reference)
        parsed = []
        for region in regions:
            match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
            if not match or not 1 <= int(match[2]) <= int(match[3]):
                raise ValueError(f"Invalid 1-based region: {region}")
            parsed.append((match[1], int(match[2]), int(match[3])))
        bundle = await asyncio.to_thread(self._bundle, [], parsed, reference)
        bundle.unresolved = [
            region
            for region, coordinates in zip(regions, parsed, strict=True)
            if not any(
                model.seq_region_name.removeprefix("chr") == coordinates[0].removeprefix("chr")
                and model.start <= coordinates[2]
                and model.end >= coordinates[1]
                for model in bundle.transcripts.values()
            )
        ]
        return bundle


def annotation_client_from_environment(database: str, species: str) -> AbstractTranscriptAnnotationClient | None:
    """Select an explicit local provider, or retain the existing Ensembl default."""
    explicit = os.environ.get("SIRNAFORGE_ANNOTATION_PROVIDER")
    provider = explicit or ("gencode" if database == "gencode" else "ensembl")
    if provider == "none":
        return None
    if provider == "ensembl":
        return EnsemblTranscriptModelClient()
    if provider not in {"gencode", "refseq"}:
        raise ValueError("SIRNAFORGE_ANNOTATION_PROVIDER must be ensembl, gencode, refseq, or none")
    prefix = f"SIRNAFORGE_{provider.upper()}"
    path = os.environ.get(f"{prefix}_GTF") or os.environ.get(f"{prefix}_ANNOTATION")
    assembly = os.environ.get(f"{prefix}_ASSEMBLY")
    if not path or not assembly:
        raise ValueError(f"{provider} annotation requires {prefix}_ANNOTATION (or _GTF) and {prefix}_ASSEMBLY")
    declared_species = os.environ.get(f"{prefix}_SPECIES", "human")
    if normalize_species_name(species) != normalize_species_name(declared_species):
        raise ValueError(f"{provider} annotation is configured for {declared_species}, requested {species}")
    return LocalTranscriptModelClient(path, provider=provider, species=declared_species, assembly=assembly)
