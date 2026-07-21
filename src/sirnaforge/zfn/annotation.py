"""GTF/GFF-backed annotation provider for ZFN off-target sites."""

from __future__ import annotations

import gzip
import re
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, TypeVar

from sirnaforge.models.zfn import GenomicAnnotationConfig, ZFNOffTargetSite
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

_ATTR_PATTERN = re.compile(r"\s*([^\s=;]+)\s*[=\s]\s*\"?([^\";]+)\"?\s*$")


class _IntervalLike(Protocol):
    start_1based: int
    end_1based: int


_TInterval = TypeVar("_TInterval", bound=_IntervalLike)


@dataclass(slots=True)
class _IntervalRecord:
    start_1based: int
    end_1based: int
    gene: str | None


@dataclass(slots=True)
class _GeneRecord:
    start_1based: int
    end_1based: int
    gene: str | None
    strand: str

    @property
    def tss_1based(self) -> int:
        return self.start_1based if self.strand != "-" else self.end_1based


@dataclass(slots=True)
class _ChromFeatureIndex:
    exons: list[_IntervalRecord]
    exon_starts: list[int]
    genes: list[_GeneRecord]
    gene_starts: list[int]
    tss_positions: list[int]
    tss_records: list[_GeneRecord]


class GTFZFNAnnotationProvider:
    """Annotate predicted ZFN sites using local GTF/GFF intervals."""

    def __init__(self) -> None:
        """Initialize provider-level cache of parsed annotation indexes."""
        self._cache: dict[Path, dict[str, _ChromFeatureIndex]] = {}

    def annotate(self, sites: list[ZFNOffTargetSite], config: GenomicAnnotationConfig) -> list[ZFNOffTargetSite]:
        """Attach region and nearest gene labels to each site."""
        annotation_path = config.resolved_annotation_path()
        if annotation_path is None:
            logger.warning("ZFN annotation skipped: no resolved annotation path")
            return sites

        feature_index = self._cache.get(annotation_path)
        if feature_index is None:
            feature_index = self._parse_annotation_file(annotation_path)
            self._cache[annotation_path] = feature_index

        if not feature_index:
            logger.warning("ZFN annotation skipped: parsed no usable features from %s", annotation_path)
            return sites

        updated: list[ZFNOffTargetSite] = []
        for site in sites:
            chrom_index = self._find_chrom_index(feature_index, site.chrom)
            if chrom_index is None:
                updated.append(site.model_copy(update={"region": "unknown", "nearest_gene": None}))
                continue

            region, nearest_gene = self._annotate_one_site(site, chrom_index, config)
            updated.append(site.model_copy(update={"region": region, "nearest_gene": nearest_gene}))

        return updated

    @staticmethod
    def _normalize_chrom(chrom: str) -> str:
        token = chrom.strip().lower()
        if token.startswith("chr"):
            token = token[3:]
        return token

    def _find_chrom_index(self, feature_index: dict[str, _ChromFeatureIndex], chrom: str) -> _ChromFeatureIndex | None:
        if chrom in feature_index:
            return feature_index[chrom]

        normalized = self._normalize_chrom(chrom)
        for key, value in feature_index.items():
            if self._normalize_chrom(key) == normalized:
                return value
        return None

    def _annotate_one_site(
        self,
        site: ZFNOffTargetSite,
        chrom_index: _ChromFeatureIndex,
        config: GenomicAnnotationConfig,
    ) -> tuple[str, str | None]:
        start_1based = site.start_1based
        end_1based = site.end_1based

        overlapping_exons = self._find_overlaps(start_1based, end_1based, chrom_index.exons, chrom_index.exon_starts)
        exon_genes = [record.gene for record in overlapping_exons if record.gene]
        if overlapping_exons:
            return "exon", exon_genes[0] if exon_genes else self._nearest_gene_for_interval(
                start_1based, end_1based, chrom_index
            )

        promoter_gene = self._find_promoter_gene(start_1based, end_1based, chrom_index, config)
        if promoter_gene is not None:
            return "promoter", promoter_gene

        overlapping_genes = self._find_overlaps(start_1based, end_1based, chrom_index.genes, chrom_index.gene_starts)
        gene_names = [record.gene for record in overlapping_genes if record.gene]
        if overlapping_genes:
            return "intron", gene_names[0] if gene_names else self._nearest_gene_for_interval(
                start_1based, end_1based, chrom_index
            )

        return "intergenic", self._nearest_gene_for_interval(start_1based, end_1based, chrom_index)

    @staticmethod
    def _find_overlaps(
        query_start_1based: int,
        query_end_1based: int,
        records: list[_TInterval],
        starts: list[int],
    ) -> list[_TInterval]:
        if not records:
            return []

        right = bisect_right(starts, query_end_1based)
        overlaps: list[_TInterval] = []
        idx = right - 1
        while idx >= 0:
            record = records[idx]
            if record.end_1based < query_start_1based:
                break
            if record.start_1based <= query_end_1based:
                overlaps.append(record)
            idx -= 1
        return overlaps

    def _find_promoter_gene(
        self,
        query_start_1based: int,
        query_end_1based: int,
        chrom_index: _ChromFeatureIndex,
        config: GenomicAnnotationConfig,
    ) -> str | None:
        for gene in chrom_index.genes:
            tss = gene.tss_1based
            if gene.strand == "-":
                promoter_start = max(1, tss - config.promoter_downstream_bp)
                promoter_end = tss + config.promoter_upstream_bp
            else:
                promoter_start = max(1, tss - config.promoter_upstream_bp)
                promoter_end = tss + config.promoter_downstream_bp

            if promoter_start <= query_end_1based and promoter_end >= query_start_1based:
                return gene.gene
        return None

    @staticmethod
    def _nearest_gene_for_interval(
        query_start_1based: int,
        query_end_1based: int,
        chrom_index: _ChromFeatureIndex,
    ) -> str | None:
        if not chrom_index.tss_positions:
            return None

        center = (query_start_1based + query_end_1based) // 2
        insert_at = bisect_right(chrom_index.tss_positions, center)

        candidates: list[_GeneRecord] = []
        if insert_at < len(chrom_index.tss_records):
            candidates.append(chrom_index.tss_records[insert_at])
        if insert_at > 0:
            candidates.append(chrom_index.tss_records[insert_at - 1])

        if not candidates:
            return None

        nearest = min(candidates, key=lambda gene: abs(gene.tss_1based - center))
        return nearest.gene

    def _parse_annotation_file(self, annotation_path: Path) -> dict[str, _ChromFeatureIndex]:
        opener = gzip.open if annotation_path.suffix == ".gz" else open
        exons_by_chrom: dict[str, list[_IntervalRecord]] = {}
        genes_by_chrom: dict[str, list[_GeneRecord]] = {}

        with opener(annotation_path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parsed = self._parse_feature_line(line)
                if parsed is None:
                    continue

                chrom, feature_type, start_1based, end_1based, strand, gene_name = parsed
                if feature_type == "exon":
                    exons_by_chrom.setdefault(chrom, []).append(
                        _IntervalRecord(start_1based=start_1based, end_1based=end_1based, gene=gene_name)
                    )
                else:
                    genes_by_chrom.setdefault(chrom, []).append(
                        _GeneRecord(
                            start_1based=start_1based,
                            end_1based=end_1based,
                            gene=gene_name,
                            strand=strand,
                        )
                    )

        # If the GTF/GFF has no explicit gene features, infer coarse gene spans from exons.
        if not genes_by_chrom and exons_by_chrom:
            for chrom, exons in exons_by_chrom.items():
                by_gene: dict[str, tuple[int, int]] = {}
                for exon in exons:
                    gene = exon.gene or "unknown"
                    if gene not in by_gene:
                        by_gene[gene] = (exon.start_1based, exon.end_1based)
                    else:
                        prev_start, prev_end = by_gene[gene]
                        by_gene[gene] = (min(prev_start, exon.start_1based), max(prev_end, exon.end_1based))
                genes_by_chrom[chrom] = [
                    _GeneRecord(start_1based=start, end_1based=end, gene=gene, strand="+")
                    for gene, (start, end) in by_gene.items()
                ]

        output: dict[str, _ChromFeatureIndex] = {}
        for chrom in set(exons_by_chrom) | set(genes_by_chrom):
            exons = sorted(exons_by_chrom.get(chrom, []), key=lambda r: (r.start_1based, r.end_1based))
            genes = sorted(genes_by_chrom.get(chrom, []), key=lambda r: (r.start_1based, r.end_1based))
            tss_records = sorted(genes, key=lambda r: r.tss_1based)
            output[chrom] = _ChromFeatureIndex(
                exons=exons,
                exon_starts=[r.start_1based for r in exons],
                genes=genes,
                gene_starts=[r.start_1based for r in genes],
                tss_positions=[r.tss_1based for r in tss_records],
                tss_records=tss_records,
            )

        logger.info(
            "Parsed ZFN annotation file %s with %s contigs",
            annotation_path,
            len(output),
        )
        return output

    def _parse_feature_line(self, line: str) -> tuple[str, str, int, int, str, str | None] | None:
        """Parse one GTF/GFF feature line into normalized fields."""
        if not line or line.startswith("#"):
            return None

        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            return None

        chrom, _source, feature_type, start_raw, end_raw, _score, strand, _phase, attr_raw = cols
        if feature_type not in {"gene", "transcript", "exon"}:
            return None

        try:
            start_1based = int(start_raw)
            end_1based = int(end_raw)
        except ValueError:
            return None

        if start_1based <= 0 or end_1based < start_1based:
            return None

        attrs = self._parse_attributes(attr_raw)
        gene_name = self._resolve_gene_name(attrs)
        normalized_strand = strand if strand in {"+", "-"} else "+"
        return chrom, feature_type, start_1based, end_1based, normalized_strand, gene_name

    @staticmethod
    def _parse_attributes(raw_attributes: str) -> dict[str, str]:
        attrs: dict[str, str] = {}
        for chunk in raw_attributes.split(";"):
            token = chunk.strip()
            if not token:
                continue
            match = _ATTR_PATTERN.match(token)
            if not match:
                continue
            key = match.group(1).strip()
            value = match.group(2).strip()
            if key and value:
                attrs[key] = value
        return attrs

    @staticmethod
    def _resolve_gene_name(attrs: dict[str, str]) -> str | None:
        for key in ("gene_name", "Name", "gene_id", "ID", "Parent"):
            value = attrs.get(key)
            if value:
                return value
        return None
