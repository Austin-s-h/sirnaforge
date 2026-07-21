"""Integration tests for ZFN annotation and reference-resolution behavior."""

from __future__ import annotations

from pathlib import Path

import pytest

from sirnaforge.data.genome_manager import GenomeManager
from sirnaforge.models.zfn import (
    GenomicAnnotationConfig,
    ZFNDesignParameters,
    ZFNHalfSiteConstraints,
    ZFNSpacerConstraints,
)
from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher

REAL_CHR3_SLICE = Path(__file__).resolve().parent / "data" / "zfn" / "ensembl_chr3_grch38_slice.fa"
REAL_CHR3_SEQUENCE = (
    "AGTGTAAATGAGCTCCAGCATCTTTATTTTTTTTTATTTTTTATTTGTTTTATTATACTTTAAGTTTTAGGGTACACGTG"
    "CACAAAGTGCAGGTTTGTTACATATGTATACATGTGCCATGTTGGTGTGCTGCACCCATTA"
)
LEFT = "TAAGTTTTAGGG"
RIGHT = "CACTTTGTGCAC"


def _fixture_fasta() -> Path:
    return REAL_CHR3_SLICE


def _write_gtf(tmp_path: Path, chrom: str, gene_name: str) -> Path:
    gtf = tmp_path / f"annotation_{chrom}.gtf"
    gtf.write_text(
        (
            f'{chrom}\ttest\tgene\t55\t100\t.\t+\t.\tgene_id "{gene_name}"; gene_name "{gene_name}";\n'
            f'{chrom}\ttest\texon\t55\t100\t.\t+\t.\tgene_id "{gene_name}"; transcript_id "TX1"; gene_name "{gene_name}";\n'
        ),
        encoding="utf-8",
    )
    return gtf


def _search_sites(
    *,
    fasta: Path | None = None,
    search_space_reference: str | None = None,
    annotation: GenomicAnnotationConfig,
) -> list:
    params = ZFNDesignParameters(
        search_space_fasta=str(fasta) if fasta is not None else None,
        search_space_reference=search_space_reference,
        left_half_site=LEFT,
        right_half_site=RIGHT,
        top_n_sites=10,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )
    return ExhaustiveZFNOffTargetSearcher().search(params, annotation=annotation)


@pytest.mark.integration
def test_zfn_real_ensembl_chr3_slice_annotation_alias_matches_between_fasta_and_gtf(tmp_path: Path) -> None:
    """A real Ensembl chr3 slice should annotate when FASTA and GTF differ only by chr prefix."""
    gtf = _write_gtf(tmp_path, chrom="chr3", gene_name="GENE_CHR3")

    sites = _search_sites(
        fasta=_fixture_fasta(),
        annotation=GenomicAnnotationConfig(annotation_path=str(gtf)),
    )

    assert sites
    assert any(site.chrom == "chr3" for site in sites)
    assert any(site.region == "exon" for site in sites)
    assert any(site.nearest_gene == "GENE_CHR3" for site in sites)


@pytest.mark.integration
def test_zfn_hg38_named_reference_annotation_alias_works_end_to_end(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Named hg38 reference searches should preserve chr alias annotation behavior through reference resolution."""
    gtf = _write_gtf(tmp_path, chrom="3", gene_name="HG38_GENE")

    def _fake_get_genome(
        self: GenomeManager,
        source_name: str,
        force_refresh: bool = False,
        build_index: bool = True,
    ) -> dict[str, Path] | None:
        del self, force_refresh, build_index
        assert source_name == "ensembl_human_hg38_primary"
        return {"fasta": _fixture_fasta()}

    monkeypatch.setattr(GenomeManager, "get_genome", _fake_get_genome)

    sites = _search_sites(
        search_space_reference="ensembl_human_hg38_primary",
        annotation=GenomicAnnotationConfig(annotation_path=str(gtf)),
    )

    assert sites
    assert any(site.chrom == "chr3" for site in sites)
    assert any(site.region == "exon" for site in sites)
    assert any(site.nearest_gene == "HG38_GENE" for site in sites)
