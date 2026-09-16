#!/usr/bin/env python3
"""Shared Ensembl reference table for transcriptome (cDNA) and genome (DNA) sources.

Both :class:`~sirnaforge.data.transcriptome_manager.TranscriptomeManager` and
:class:`~sirnaforge.data.genome_manager.GenomeManager` build their ``SOURCES`` dicts
from the single :data:`ENSEMBL_ASSEMBLIES` table defined here, so that adding a new
species is a one-line change rather than duplicated, hand-written URLs.

Notes on Ensembl FASTA layout:
- cDNA files are always published as ``<Species>.<Assembly>.cdna.all.fa.gz``.
- Genomic DNA is published as ``dna.primary_assembly`` **only** for assemblies that
  separate the primary assembly from haplotype/patch sequences (human, mouse). Other
  assemblies (e.g. rat GRCr8, macaque Mmul_10) publish only ``dna.toplevel``. A wrong
  ``dna_kind`` yields a 404 at download time, so it is pinned per species and verified
  against ``ftp.ensembl.org`` when a species is added.
- Assembly names are embedded in the filename. ``current_fasta`` tracks the latest
  Ensembl release, so a new release can rename an assembly and break a URL; keeping the
  assembly in this one table makes that a single-line fix.
- The assembly name is also inside each cDNA header (``chromosome:<Assembly>:...``), which is
  what :func:`infer_species_from_cdna_headers` reads to label a custom reference with the
  species it actually contains rather than the one the parameter claimed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .transcriptome_manager import TranscriptomeSource

ENSEMBL_FTP_BASE = "https://ftp.ensembl.org/pub/current_fasta"


@dataclass(frozen=True, slots=True)
class EnsemblAssembly:
    """One Ensembl species/assembly and how to address its cDNA and DNA FASTA files."""

    species: str
    """Canonical species key used across siRNAforge (e.g. ``human``, ``mouse``)."""

    ensembl_slug: str
    """Ensembl species directory / filename stem (e.g. ``homo_sapiens``)."""

    assembly: str
    """Assembly name embedded in the FASTA filename (e.g. ``GRCh38``)."""

    dna_kind: str
    """Genomic DNA file kind: ``primary_assembly`` or ``toplevel`` (verified per species)."""

    common_name: str
    """Human-readable name used in source descriptions (e.g. ``rhesus macaque``)."""

    genome_source_key: str
    """Stable key for the genome ``SOURCES`` dict (human preserves its legacy key)."""

    def _species_stem(self) -> str:
        """Return the capitalized ``<Species>`` stem used in Ensembl filenames."""
        return self.ensembl_slug.capitalize()

    def cdna_url(self) -> str:
        """Full URL for the all-transcript cDNA FASTA."""
        filename = f"{self._species_stem()}.{self.assembly}.cdna.all.fa.gz"
        return f"{ENSEMBL_FTP_BASE}/{self.ensembl_slug}/cdna/{filename}"

    def dna_url(self) -> str:
        """Full URL for the genomic DNA FASTA (kind pinned per species)."""
        filename = f"{self._species_stem()}.{self.assembly}.dna.{self.dna_kind}.fa.gz"
        return f"{ENSEMBL_FTP_BASE}/{self.ensembl_slug}/dna/{filename}"

    @property
    def transcriptome_source_key(self) -> str:
        """Stable key for the transcriptome (cDNA) ``SOURCES`` dict."""
        return f"ensembl_{self.species}_cdna"

    def to_transcriptome_source(self) -> TranscriptomeSource:
        """Build the cDNA :class:`TranscriptomeSource` for this species."""
        # Lazy import avoids a circular import: transcriptome_manager builds its SOURCES
        # from this module at class-definition time.
        from .transcriptome_manager import TranscriptomeSource  # noqa: PLC0415

        return TranscriptomeSource(
            name="ensembl_cdna",
            url=self.cdna_url(),
            species=self.species,
            format="fasta",
            compressed=True,
            description=f"Ensembl {self.common_name} cDNA sequences ({self.assembly})",
        )

    def to_genome_source(self) -> TranscriptomeSource:
        """Build the genomic DNA :class:`TranscriptomeSource` for this species."""
        from .transcriptome_manager import TranscriptomeSource  # noqa: PLC0415

        assembly_label = "primary assembly" if self.dna_kind == "primary_assembly" else "toplevel assembly"
        return TranscriptomeSource(
            name=self.genome_source_key,
            url=self.dna_url(),
            species=self.species,
            format="fasta",
            compressed=True,
            description=f"Ensembl {self.common_name} {self.assembly} {assembly_label} genomic reference",
        )


# Single source of truth for Ensembl references. Ordering matches the historical
# transcriptome SOURCES ordering (human, mouse, rat, macaque) so default-source
# selection is unchanged.
#
# dna_kind verified against ftp.ensembl.org/pub/current_fasta/<slug>/dna/:
#   human, mouse  -> primary_assembly (200)
#   rat, macaque  -> toplevel         (primary_assembly 404s)
ENSEMBL_ASSEMBLIES: tuple[EnsemblAssembly, ...] = (
    EnsemblAssembly(
        species="human",
        ensembl_slug="homo_sapiens",
        assembly="GRCh38",
        dna_kind="primary_assembly",
        common_name="human",
        # Preserve the pre-existing key referenced by the CLI, ZFN model, and tests.
        genome_source_key="ensembl_human_hg38_primary",
    ),
    EnsemblAssembly(
        species="mouse",
        ensembl_slug="mus_musculus",
        assembly="GRCm39",
        dna_kind="primary_assembly",
        common_name="mouse",
        genome_source_key="ensembl_mouse_grcm39_primary",
    ),
    EnsemblAssembly(
        species="rat",
        ensembl_slug="rattus_norvegicus",
        assembly="GRCr8",
        dna_kind="toplevel",
        common_name="rat",
        genome_source_key="ensembl_rat_grcr8_toplevel",
    ),
    EnsemblAssembly(
        species="macaque",
        ensembl_slug="macaca_mulatta",
        assembly="Mmul_10",
        dna_kind="toplevel",
        common_name="rhesus macaque",
        genome_source_key="ensembl_macaque_mmul10_toplevel",
    ),
)


def infer_species_from_cdna_headers(fasta_path: Path, max_records: int = 200) -> str | None:
    """Infer the canonical species of a cDNA FASTA from its own Ensembl headers.

    A custom transcriptome path used to be labelled with the literal species ``"transcriptome"``,
    which is not a species: it never matches the query species, so every hit is cross-species with
    no resolvable orthologue, and a screen against mouse cDNA published 10,217 hits as unqualified
    ``off_target`` with zero orthologs. Issue #99 requires the species come from the *resolved
    reference* rather than from the parameter it arrived through, and Ensembl cDNA headers carry it::

        >ENSMUST00000108658.11 cdna chromosome:GRCm39:11:69471109:69482701:1 gene:ENSMUSG00000059552.15

    The assembly token is matched against :data:`ENSEMBL_ASSEMBLIES`, so this needs no new table and
    stays correct as species are added. Only the first ``max_records`` headers are read: a cDNA file
    is single-species, and a mixed file should not be silently assigned one label.

    Args:
        fasta_path: Path to a (plain-text) cDNA FASTA.
        max_records: How many headers to sample before giving up.

    Returns:
        Canonical species name (e.g. ``mouse``), or None when the headers do not identify one --
        an honest unknown, so the caller can say so rather than invent a species.
    """
    assembly_to_species = {assembly.assembly: assembly.species for assembly in ENSEMBL_ASSEMBLIES}
    seen: set[str] = set()
    try:
        with fasta_path.open() as handle:
            records = 0
            for line in handle:
                if not line.startswith(">"):
                    continue
                records += 1
                if records > max_records:
                    break
                for token in line.split():
                    if token.startswith("chromosome:") or token.startswith("scaffold:"):
                        parts = token.split(":")
                        if len(parts) > 1 and parts[1] in assembly_to_species:
                            seen.add(assembly_to_species[parts[1]])
    except (OSError, UnicodeDecodeError):
        # A compressed or binary file at this path infers nothing rather than aborting the run: the
        # caller may have named an index prefix whose neighbouring file is not plain-text cDNA.
        return None

    if len(seen) == 1:
        return seen.pop()
    return None


def build_transcriptome_sources() -> dict[str, TranscriptomeSource]:
    """Build the transcriptome (cDNA) ``SOURCES`` dict from the assembly table."""
    return {assembly.transcriptome_source_key: assembly.to_transcriptome_source() for assembly in ENSEMBL_ASSEMBLIES}


def build_genome_sources() -> dict[str, TranscriptomeSource]:
    """Build the genome (DNA) ``SOURCES`` dict from the assembly table."""
    return {assembly.genome_source_key: assembly.to_genome_source() for assembly in ENSEMBL_ASSEMBLIES}
