"""Unit tests for ZFN cache-resolved search space and design support."""

from pathlib import Path

import pytest
from Bio.Seq import Seq

from sirnaforge.data.genome_manager import GenomeManager
from sirnaforge.models.zfn import (
    DimerMode,
    IUPACMode,
    ZFNAlgorithm,
    ZFNDesignParameters,
    ZFNHalfSiteConstraints,
    ZFNSpacerConstraints,
)
from sirnaforge.zfn.design import ZFNDesigner
from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher

LEFT = "GCGTACGTA"
RIGHT = "TACGGCATA"


def _write_fasta(tmp_path: Path, sequence: str, name: str = "chr1") -> Path:
    fasta = tmp_path / "search_space.fa"
    fasta.write_text(f">{name}\n{sequence}\n", encoding="utf-8")
    return fasta


def _canonical_site(left: str, right: str, spacer: str = "AAAAA") -> str:
    return f"{left}{spacer}{str(Seq(right).reverse_complement())}"


def test_explicit_search_space_fasta_finds_heterodimer_sites(tmp_path: Path) -> None:
    """Explicit FASTA search spaces should produce heterodimer-compatible sites."""
    sequence = f"TTTT{_canonical_site(LEFT, RIGHT)}CCCC{_canonical_site(LEFT, RIGHT)}"
    fasta = _write_fasta(tmp_path, sequence)

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )

    sites = ExhaustiveZFNOffTargetSearcher().search(params)
    assert len(sites) >= 2
    assert all(site.spacer_len == 5 for site in sites)


def test_iupac_allow_mode_matches_and_none_mode_rejects(tmp_path: Path) -> None:
    """IUPAC-aware mode should match ambiguous queries while NONE mode does not."""
    left_iupac = "ACGTNRYWS"
    right_iupac = "TGCANRYWS"
    left_concrete = "ACGTAGTAC"
    right_concrete = "TGCAAGTAC"
    sequence = f"AAA{_canonical_site(left_concrete, right_concrete)}TTT"
    fasta = _write_fasta(tmp_path, sequence)

    allow_params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=left_iupac,
        right_half_site=right_iupac,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0, iupac_mode=IUPACMode.ALLOW_IUPAC),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )
    none_params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=left_iupac,
        right_half_site=right_iupac,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0, iupac_mode=IUPACMode.NONE),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )

    searcher = ExhaustiveZFNOffTargetSearcher()
    assert len(searcher.search(allow_params)) >= 1
    assert len(searcher.search(none_params)) == 0


def test_include_homodimers_enumerates_ll_or_rr_sites(tmp_path: Path) -> None:
    """Homodimer mode should enumerate LL/RR orientation sites when present."""
    ll_site = f"{LEFT}AAAAA{LEFT}"
    rr_site = f"{str(Seq(RIGHT).reverse_complement())}AAAAA{str(Seq(RIGHT).reverse_complement())}"
    sequence = f"GGG{ll_site}TTT{rr_site}CCC"
    fasta = _write_fasta(tmp_path, sequence)

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        dimer_mode=DimerMode.INCLUDE_HOMODIMERS,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5], require_opposite_strands=False),
    )

    sites = ExhaustiveZFNOffTargetSearcher().search(params)
    orientations = {site.orientation.value for site in sites}
    assert "L...L" in orientations or "R...R" in orientations


def test_search_space_reference_uses_genome_cache_manager(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Named search-space references should resolve through GenomeManager."""
    sequence = f"TTT{_canonical_site(LEFT, RIGHT)}AAA"
    fasta = _write_fasta(tmp_path, sequence)

    def _fake_get_genome(self: GenomeManager, source_name: str, force_refresh: bool = False, build_index: bool = True):
        assert source_name == "ensembl_human_hg38_primary"
        return {"fasta": fasta}

    monkeypatch.setattr(GenomeManager, "get_genome", _fake_get_genome)

    params = ZFNDesignParameters(
        search_space_reference="ensembl_human_hg38_primary",
        search_space_fasta=None,
        left_half_site=LEFT,
        right_half_site=RIGHT,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )

    sites = ExhaustiveZFNOffTargetSearcher().search(params)
    assert len(sites) >= 1


def test_sliding_window_stride_controls_scan_density(tmp_path: Path) -> None:
    """Sliding-window stride should tune genomic scan density deterministically."""
    sequence = f"A{_canonical_site(LEFT, RIGHT)}AA{_canonical_site(LEFT, RIGHT)}A"
    fasta = _write_fasta(tmp_path, sequence)

    dense_params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0, window_stride=1),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )
    sparse_params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0, window_stride=2),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )

    searcher = ExhaustiveZFNOffTargetSearcher()
    dense_sites = searcher.search(dense_params)
    sparse_sites = searcher.search(sparse_params)

    assert len(dense_sites) >= len(sparse_sites)


@pytest.mark.parametrize("algorithm", [ZFNAlgorithm.HOMOLOGY, ZFNAlgorithm.CONSERVED_G, ZFNAlgorithm.ZFN_V2])
def test_designer_supports_multiple_algorithms(tmp_path: Path, algorithm: ZFNAlgorithm) -> None:
    """ZFN designer should support and score all configured algorithms."""
    sequence = f"AAA{_canonical_site(LEFT, RIGHT)}CCC"
    fasta = _write_fasta(tmp_path, sequence)

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        algorithm=algorithm,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )

    result = ZFNDesigner().evaluate_pair(params)
    assert len(result.candidates) == 1
    candidate = result.candidates[0]
    assert 0.0 <= candidate.composite_score <= 100.0
    assert isinstance(candidate.passes_offtarget_filters, (bool, str))
