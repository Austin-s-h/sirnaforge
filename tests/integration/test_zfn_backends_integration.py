"""Integration tests validating parity across ZFN search backends."""

from pathlib import Path

import pytest
from Bio.Seq import Seq

from sirnaforge.models.zfn import (
    ZFNDesignParameters,
    ZFNHalfSiteConstraints,
    ZFNSearchBackend,
    ZFNSpacerConstraints,
)
from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher, build_zfn_search_index

pytest.importorskip("ahocorasick")
pytest.importorskip("fm_index")

LEFT = "GCGTACGTA"
RIGHT = "TACGGCATA"


def _canonical_site(left: str, right: str, spacer: str = "AAAAA") -> str:
    """Build one canonical ZFN site in genomic L...R orientation."""
    return f"{left}{spacer}{str(Seq(right).reverse_complement())}"


def _write_multicontig_fasta(tmp_path: Path) -> Path:
    """Create one small two-contig FASTA containing deterministic ZFN sites."""
    chr2 = f"TTTT{_canonical_site(LEFT, RIGHT)}CCCC{_canonical_site(LEFT, RIGHT, spacer='AAAAAA')}"
    chr3 = f"GGGG{_canonical_site(LEFT, RIGHT)}AAAA"

    fasta = tmp_path / "zfn_backends.fa"
    fasta.write_text(f">chr2\n{chr2}\n>chr3\n{chr3}\n", encoding="utf-8")
    return fasta


def _run_search(
    *,
    backend: ZFNSearchBackend,
    fasta: Path,
    search_space_index: str | None = None,
) -> set[tuple[str, int, int, str, int, int]]:
    """Run one backend and return stable off-target identity keys."""
    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        search_space_index=search_space_index,
        left_half_site=LEFT,
        right_half_site=RIGHT,
        search_backend=backend,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5, 6]),
    )
    sites = ExhaustiveZFNOffTargetSearcher().search(params)
    return {
        (
            s.chrom,
            s.start_1based,
            s.end_1based,
            s.orientation.value,
            s.left_mismatches,
            s.right_mismatches,
        )
        for s in sites
    }


@pytest.fixture
def zfn_multicontig_fasta(tmp_path: Path) -> Path:
    """Small deterministic FASTA shared by all backend case variants."""
    return _write_multicontig_fasta(tmp_path)


@pytest.fixture
def exhaustive_baseline(zfn_multicontig_fasta: Path) -> set[tuple[str, int, int, str, int, int]]:
    """Canonical result set used as the backend parity baseline."""
    return _run_search(backend=ZFNSearchBackend.EXHAUSTIVE_PYTHON, fasta=zfn_multicontig_fasta)


@pytest.fixture(
    params=[
        pytest.param((ZFNSearchBackend.PYAHOCORASICK, False), id="pyahocorasick"),
        pytest.param((ZFNSearchBackend.FM_INDEX, False), id="fm_index_live"),
        pytest.param((ZFNSearchBackend.FM_INDEX, True), id="fm_index_persisted"),
    ]
)
def backend_case(
    request: pytest.FixtureRequest,
    tmp_path: Path,
    zfn_multicontig_fasta: Path,
) -> tuple[ZFNSearchBackend, str | None]:
    """Provide backend + optional persisted index bundle path for parity checks."""
    backend, use_persisted_index = request.param
    search_space_index: str | None = None
    if backend == ZFNSearchBackend.FM_INDEX and use_persisted_index:
        bundle_dir = tmp_path / "fm_bundle"
        summary = build_zfn_search_index(
            backend=ZFNSearchBackend.FM_INDEX,
            genome_fasta=zfn_multicontig_fasta,
            output_dir=bundle_dir,
        )
        assert summary["backend"] == ZFNSearchBackend.FM_INDEX.value
        search_space_index = str(bundle_dir)
    return backend, search_space_index


@pytest.mark.integration
def test_zfn_backends_match_on_real_optional_dependencies(
    zfn_multicontig_fasta: Path,
    exhaustive_baseline: set[tuple[str, int, int, str, int, int]],
    backend_case: tuple[ZFNSearchBackend, str | None],
) -> None:
    """Each optional backend variant should match exhaustive_python on identical inputs."""
    backend, search_space_index = backend_case
    backend_result = _run_search(
        backend=backend,
        fasta=zfn_multicontig_fasta,
        search_space_index=search_space_index,
    )
    assert backend_result == exhaustive_baseline
