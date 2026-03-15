"""Unit tests for ZFN cache-resolved search space and design support."""

from pathlib import Path

import pytest
from Bio.Seq import Seq

from sirnaforge.data.annotation_manager import AnnotationManager
from sirnaforge.data.genome_manager import GenomeManager
from sirnaforge.models.zfn import (
    DimerMode,
    GenomicAnnotationConfig,
    IUPACMode,
    MatchOrientation,
    Strand,
    ZFNAlgorithm,
    ZFNDesignParameters,
    ZFNHalfSiteConstraints,
    ZFNOffTargetFilterCriteria,
    ZFNOffTargetSite,
    ZFNShardingConfig,
    ZFNSpacerConstraints,
)
from sirnaforge.zfn.design import ZFNDesigner
from sirnaforge.zfn.rank import rank_sites
from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher, _HalfSiteHit

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
        algorithm=ZFNAlgorithm.HOMOLOGY,
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


def test_exhaustive_stride_1_finds_all_sites_in_toy_fasta(tmp_path: Path) -> None:
    """Stride=1 scanning should recover all inserted canonical sites in deterministic FASTA."""
    right_rc = str(Seq(RIGHT).reverse_complement())
    site = f"{LEFT}AAAAA{right_rc}"
    sequence = f"GGG{site}TTTTTTTT{site}CCC"
    fasta = _write_fasta(tmp_path, sequence)

    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        search_space_fasta=str(fasta),
        half_site_constraints=ZFNHalfSiteConstraints(window_stride=1, max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )

    sites = ExhaustiveZFNOffTargetSearcher().search(params)
    starts = {site.start_1based for site in sites}
    assert len(sites) == 2
    assert starts == {4, 35}


def test_iupac_allow_mode_matches_ambiguity_codes(tmp_path: Path) -> None:
    """ALLOW_IUPAC should match ambiguous right half-sites (e.g. W=A/T) better than NONE."""
    left = "ACGTACGTACG"
    right_concrete = "TTGGAATTCCA"
    right_ambiguous = "TTGGAAWWCCA"
    sequence = f"CCC{_canonical_site(left, right_concrete, spacer='AAAAAA')}GGG"
    fasta = _write_fasta(tmp_path, sequence)

    allow_params = ZFNDesignParameters(
        left_half_site=left,
        right_half_site=right_ambiguous,
        search_space_fasta=str(fasta),
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0, iupac_mode=IUPACMode.ALLOW_IUPAC),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[6]),
    )
    none_params = ZFNDesignParameters(
        left_half_site=left,
        right_half_site=right_ambiguous,
        search_space_fasta=str(fasta),
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0, iupac_mode=IUPACMode.NONE),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[6]),
    )

    searcher = ExhaustiveZFNOffTargetSearcher()
    assert len(searcher.search(allow_params)) >= 1
    assert len(searcher.search(none_params)) == 0


def test_spacer_length_filters_sites(tmp_path: Path) -> None:
    """Allowed spacer lengths should strictly filter predicted site sets."""
    sequence = f"AAA{_canonical_site(LEFT, RIGHT, spacer='AAAAA')}TTT{_canonical_site(LEFT, RIGHT, spacer='AAAAAA')}CCC"
    fasta = _write_fasta(tmp_path, sequence)

    base = {
        "left_half_site": LEFT,
        "right_half_site": RIGHT,
        "search_space_fasta": str(fasta),
        "half_site_constraints": ZFNHalfSiteConstraints(max_mismatches=0),
    }

    sites_5 = ExhaustiveZFNOffTargetSearcher().search(
        ZFNDesignParameters(**base, spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]))
    )
    sites_6 = ExhaustiveZFNOffTargetSearcher().search(
        ZFNDesignParameters(**base, spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[6]))
    )

    assert sites_5
    assert sites_6
    assert {s.spacer_len for s in sites_5} == {5}
    assert {s.spacer_len for s in sites_6} == {6}
    assert {(s.chrom, s.start_1based, s.end_1based) for s in sites_5} != {
        (s.chrom, s.start_1based, s.end_1based) for s in sites_6
    }


def test_sharded_search_dedupes_overlap_boundary_hits(tmp_path: Path) -> None:
    """Chunk overlap should avoid boundary misses and collapse duplicate shard hits."""
    canonical = _canonical_site(LEFT, RIGHT, spacer="AAAAA")
    sequence = "A" * 40 + canonical + "C" * 60
    fasta = _write_fasta(tmp_path, sequence, name="chr3")

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
        sharding=ZFNShardingConfig(
            enabled=True,
            chunk_size_bp=50,
            overlap_bp=0,
            chromosomes=["chr3"],
            max_workers=2,
        ),
    )

    sites = ExhaustiveZFNOffTargetSearcher().search(params)
    assert len(sites) == 1
    site = sites[0]
    assert site.chrom == "chr3"
    assert site.start_1based == 41


def test_sharded_search_matches_unsharded_coordinates(tmp_path: Path) -> None:
    """Sharded search should match unsharded genomic site set for deterministic toy FASTA."""
    sequence = (
        "TT"
        + _canonical_site(LEFT, RIGHT, spacer="AAAAA")
        + "GGGGG"
        + _canonical_site(LEFT, RIGHT, spacer="AAAAAA")
        + "CC"
    )
    fasta = _write_fasta(tmp_path, sequence, name="chr3")

    base_kwargs = {
        "search_space_fasta": str(fasta),
        "left_half_site": LEFT,
        "right_half_site": RIGHT,
        "half_site_constraints": ZFNHalfSiteConstraints(max_mismatches=0),
        "spacer_constraints": ZFNSpacerConstraints(allowed_spacer_lengths=[5, 6]),
    }

    unsharded = ExhaustiveZFNOffTargetSearcher().search(ZFNDesignParameters(**base_kwargs))
    sharded = ExhaustiveZFNOffTargetSearcher().search(
        ZFNDesignParameters(
            **base_kwargs,
            sharding=ZFNShardingConfig(
                enabled=True,
                chunk_size_bp=40,
                overlap_bp=0,
                chromosomes=["3"],
                max_workers=2,
            ),
        )
    )

    unsharded_key = {(s.chrom, s.start_1based, s.end_1based, s.orientation.value) for s in unsharded}
    sharded_key = {(s.chrom, s.start_1based, s.end_1based, s.orientation.value) for s in sharded}
    assert sharded_key == unsharded_key


def test_single_contig_sharding_chunks_large_contig(tmp_path: Path) -> None:
    """Single-contig inputs should still chunk when sharding is enabled."""
    sequence = "A" * 200
    fasta = _write_fasta(tmp_path, sequence, name="chr3")

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        sharding=ZFNShardingConfig(enabled=True, chunk_size_bp=25, overlap_bp=0),
    )
    searcher = ExhaustiveZFNOffTargetSearcher()
    chrom_sequences = searcher._load_fasta(fasta)

    shards = searcher._build_shard_specs(chrom_sequences, params)
    assert len(shards) == 8
    assert shards[0].core_start0 == 0
    assert shards[0].core_end0 == 25
    assert shards[0].scan_start0 == 0
    assert shards[0].scan_end0 == 75
    assert shards[-1].core_end0 == len(sequence)


def test_memory_based_worker_cap_allows_two_workers_on_stable_mid_memory_host(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The live-memory cap should reflect the current default reserve on a host with about 5 GiB free."""
    searcher = ExhaustiveZFNOffTargetSearcher()

    monkeypatch.setattr(searcher, "_available_memory_gb", lambda: 4.9)

    expected = int((4.9 - searcher._DEFAULT_MEMORY_RESERVE_GB) / searcher._DEFAULT_PER_WORKER_GB)
    assert searcher._memory_based_worker_cap(8) == max(1, expected)


def test_memory_based_worker_cap_honors_internal_reserve_override(monkeypatch: pytest.MonkeyPatch) -> None:
    """Internal env override should allow restoring a more conservative reserve when profiling."""
    searcher = ExhaustiveZFNOffTargetSearcher()

    monkeypatch.setattr(searcher, "_available_memory_gb", lambda: 4.9)
    monkeypatch.setenv("SIRNAFORGE_ZFN_MEMORY_RESERVE_GB", "4.0")

    assert searcher._memory_based_worker_cap(8) == 1


def test_sharded_search_matches_unsharded_on_large_complex_synthetic_genome(tmp_path: Path) -> None:
    """Sharded and unsharded scans should agree on a larger, multi-contig synthetic genome."""
    site_5 = _canonical_site(LEFT, RIGHT, spacer="AAAAA")
    site_6 = _canonical_site(LEFT, RIGHT, spacer="AAAAAA")
    decoy = "G" * len(site_5)

    chr3 = ("A" * 30_000) + site_5 + ("C" * 31_000) + decoy + ("T" * 31_000) + site_6 + ("A" * 30_000)
    chr5 = ("C" * 25_000) + site_6 + ("G" * 40_000) + site_5 + ("T" * 25_000)

    fasta = tmp_path / "complex_synthetic.fa"
    fasta.write_text(f">chr3\n{chr3}\n>chr5\n{chr5}\n", encoding="utf-8")

    base_kwargs = {
        "search_space_fasta": str(fasta),
        "left_half_site": LEFT,
        "right_half_site": RIGHT,
        "half_site_constraints": ZFNHalfSiteConstraints(max_mismatches=0),
        "spacer_constraints": ZFNSpacerConstraints(allowed_spacer_lengths=[5, 6]),
    }

    searcher = ExhaustiveZFNOffTargetSearcher()
    unsharded = searcher.search(ZFNDesignParameters(**base_kwargs, sharding=ZFNShardingConfig(enabled=False)))
    sharded = searcher.search(
        ZFNDesignParameters(
            **base_kwargs,
            sharding=ZFNShardingConfig(
                enabled=True,
                chunk_size_bp=20_000,
                overlap_bp=0,
                chromosomes=["chr3", "chr5"],
                max_workers=4,
            ),
        )
    )

    unsharded_key = {(s.chrom, s.start_1based, s.end_1based, s.orientation.value) for s in unsharded}
    sharded_key = {(s.chrom, s.start_1based, s.end_1based, s.orientation.value) for s in sharded}
    assert sharded_key == unsharded_key
    assert len(sharded_key) >= 4


def test_chromosome_filter_tokens_support_groups_ranges_and_globs() -> None:
    """Chromosome filters should support aliases, groups, numeric ranges, and glob patterns."""
    searcher = ExhaustiveZFNOffTargetSearcher()
    chrom_sequences = {
        "chr1": "A" * 10,
        "chr2": "A" * 10,
        "chrX": "A" * 10,
        "chrM": "A" * 10,
        "chrUn_gl000220": "A" * 10,
    }

    assert searcher._resolve_target_contigs(chrom_sequences, ["autosomes"]) == ["chr1", "chr2"]
    assert searcher._resolve_target_contigs(chrom_sequences, ["sex"]) == ["chrX"]
    assert searcher._resolve_target_contigs(chrom_sequences, ["mito"]) == ["chrM"]
    assert searcher._resolve_target_contigs(chrom_sequences, ["1-2"]) == ["chr1", "chr2"]
    assert searcher._resolve_target_contigs(chrom_sequences, ["chrUn_*"]) == ["chrUn_gl000220"]
    assert searcher._resolve_target_contigs(chrom_sequences, ["3", "chr1"]) == ["chr1"]


def test_ranking_tiebreak_region_exon_promoter_intron_intergenic() -> None:
    """Ranking API should tie-break equal scores by region priority."""
    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        algorithm=ZFNAlgorithm.HOMOLOGY,
    )

    sites = [
        ZFNOffTargetSite(
            site_id="s1",
            chrom="chr1",
            start_1based=1,
            end_1based=23,
            strand=Strand.PLUS,
            orientation=MatchOrientation.LR,
            spacer_len=5,
            sequence="N" * 23,
            left_mismatches=1,
            right_mismatches=1,
            total_mismatches=2,
            score=80.0,
            region="intergenic",
            nearest_gene=None,
            left_aligned="x",
            right_aligned="y",
        ),
        ZFNOffTargetSite(
            site_id="s2",
            chrom="chr1",
            start_1based=2,
            end_1based=24,
            strand=Strand.PLUS,
            orientation=MatchOrientation.LR,
            spacer_len=5,
            sequence="N" * 23,
            left_mismatches=1,
            right_mismatches=1,
            total_mismatches=2,
            score=80.0,
            region="exon",
            nearest_gene="GENE",
            left_aligned="x",
            right_aligned="y",
        ),
        ZFNOffTargetSite(
            site_id="s3",
            chrom="chr1",
            start_1based=3,
            end_1based=25,
            strand=Strand.PLUS,
            orientation=MatchOrientation.LR,
            spacer_len=5,
            sequence="N" * 23,
            left_mismatches=1,
            right_mismatches=1,
            total_mismatches=2,
            score=80.0,
            region="promoter",
            nearest_gene="GENE",
            left_aligned="x",
            right_aligned="y",
        ),
    ]

    ranked = rank_sites(sites, params)
    assert [s.site_id for s in ranked[:3]] == ["s2", "s3", "s1"]


def test_ranking_tiebreak_prefers_lower_mismatch_burden() -> None:
    """For equal score, ranking should prefer lower mismatch burden deterministically."""
    params = ZFNDesignParameters(
        left_half_site=LEFT, right_half_site=RIGHT, spacer_constraints={"allowed_spacer_lengths": [5]}
    )
    base = {
        "chrom": "chr1",
        "start_1based": 10,
        "end_1based": 32,
        "strand": Strand.PLUS,
        "orientation": MatchOrientation.LR,
        "spacer_len": 5,
        "sequence": "N" * 23,
        "score": 90.0,
        "region": "unknown",
        "nearest_gene": None,
        "left_aligned": "x",
        "right_aligned": "y",
    }
    low_mm = ZFNOffTargetSite(
        site_id="low_mm",
        left_mismatches=0,
        right_mismatches=0,
        total_mismatches=0,
        **base,
    )
    high_mm = ZFNOffTargetSite(
        site_id="high_mm",
        left_mismatches=1,
        right_mismatches=1,
        total_mismatches=2,
        **base,
    )

    ranked = rank_sites([high_mm, low_mm], params)
    assert [site.site_id for site in ranked] == ["low_mm", "high_mm"]


def test_resolve_search_space_fasta_raises_when_custom_resolution_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    """Custom FASTA path should raise when cache manager cannot resolve a fasta path."""

    def _fake_get_custom_genome(
        self: GenomeManager,
        fasta_path: str,
        build_index: bool = True,
        cache_name: str | None = None,
    ) -> dict[str, Path] | None:
        return None

    monkeypatch.setattr(GenomeManager, "get_custom_genome", _fake_get_custom_genome)
    params = ZFNDesignParameters(left_half_site=LEFT, right_half_site=RIGHT, search_space_fasta="/tmp/missing.fa")

    with pytest.raises(ValueError, match="Unable to resolve search_space_fasta"):
        ExhaustiveZFNOffTargetSearcher()._resolve_search_space_fasta(params)


def test_resolve_search_space_fasta_raises_when_named_reference_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    """Named source should raise when cache manager cannot resolve a fasta path."""

    def _fake_get_genome(
        self: GenomeManager,
        source_name: str,
        force_refresh: bool = False,
        build_index: bool = True,
    ) -> dict[str, Path] | None:
        return None

    monkeypatch.setattr(GenomeManager, "get_genome", _fake_get_genome)
    params = ZFNDesignParameters(left_half_site=LEFT, right_half_site=RIGHT, search_space_fasta=None)

    with pytest.raises(ValueError, match="Unable to resolve search_space_reference"):
        ExhaustiveZFNOffTargetSearcher()._resolve_search_space_fasta(params)


def test_resolve_annotation_path_handles_missing_file(tmp_path: Path) -> None:
    """Annotation resolver should return None when configured path does not exist."""
    annotation = GenomicAnnotationConfig(annotation_path=str(tmp_path / "missing.gtf"))
    resolved = ExhaustiveZFNOffTargetSearcher()._resolve_annotation_path(annotation)
    assert resolved is None


def test_resolve_annotation_path_uses_existing_reference_path(tmp_path: Path) -> None:
    """Annotation resolver should accept annotation_reference when it points to a local file."""
    annotation_file = tmp_path / "annotations.gtf"
    annotation_file.write_text("# gtf\n", encoding="utf-8")

    annotation = GenomicAnnotationConfig(annotation_reference=str(annotation_file))
    resolved = ExhaustiveZFNOffTargetSearcher()._resolve_annotation_path(annotation)

    assert resolved is not None
    assert resolved.exists()
    assert resolved.name.endswith(".fa")
    assert annotation.annotation_path == str(resolved)


def test_resolve_annotation_path_downloads_reference_url(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Annotation resolver should download URL-like references into a local cache path."""
    downloaded = tmp_path / "downloaded.gtf.gz"
    downloaded.write_text("# gtf\n", encoding="utf-8")

    def _fake_get_custom_annotation(
        self: AnnotationManager,
        annotation_path_or_url: str | Path,
        cache_name: str | None = None,
    ) -> Path:
        assert input_location == "https://example.org/annotations.gtf.gz"
        assert cache_name is None
        return downloaded

    input_location = "https://example.org/annotations.gtf.gz"
    monkeypatch.setattr(AnnotationManager, "get_custom_annotation", _fake_get_custom_annotation)

    annotation = GenomicAnnotationConfig(
        annotation_reference=input_location,
        cache_dir=str(tmp_path),
    )
    resolved = ExhaustiveZFNOffTargetSearcher()._resolve_annotation_path(annotation)

    assert resolved == downloaded
    assert annotation.annotation_path == str(downloaded)


def test_load_fasta_raises_on_missing_and_empty(tmp_path: Path) -> None:
    """FASTA loader should fail on missing file and empty fasta content."""
    searcher = ExhaustiveZFNOffTargetSearcher()

    with pytest.raises(FileNotFoundError):
        searcher._load_fasta(tmp_path / "absent.fa")

    empty_fasta = tmp_path / "empty.fa"
    empty_fasta.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="No sequences found"):
        searcher._load_fasta(empty_fasta)


def test_seed_positions_and_evaluate_window_edge_branches() -> None:
    """Seed-position and mismatch-threshold edge branches should behave deterministically."""
    searcher = ExhaustiveZFNOffTargetSearcher()

    assert searcher._seed_positions("L", seq_len=10, seed_len=None) == set()
    assert searcher._seed_positions("R", seq_len=10, seed_len=0) == set()
    assert searcher._seed_positions("L", seq_len=6, seed_len=10) == {0, 1, 2, 3, 4, 5}
    assert searcher._seed_positions("R", seq_len=6, seed_len=2) == {0, 1}

    assert (
        searcher._evaluate_window(
            kind="L",
            query="AAAA",
            observed="TTTT",
            mode=IUPACMode.NONE,
            seed_len=2,
            max_mm=1,
            max_seed_mm=1,
        )
        is None
    )
    assert (
        searcher._evaluate_window(
            kind="L",
            query="AAAA",
            observed="AATT",
            mode=IUPACMode.NONE,
            seed_len=2,
            max_mm=4,
            max_seed_mm=0,
        )
        is None
    )


def test_scan_half_site_skips_short_contigs() -> None:
    """Scanning should skip contigs shorter than query length."""
    params = ZFNDesignParameters(left_half_site=LEFT, right_half_site=RIGHT)
    hits = ExhaustiveZFNOffTargetSearcher()._scan_half_site("L", LEFT, {"chrShort": "AAAA"}, params)
    assert hits == []


def test_build_site_rejects_invalid_pairing_modes_and_spacer() -> None:
    """Pair builder should reject same-strand, invalid spacer, and LL in heterodimer-only mode."""
    searcher = ExhaustiveZFNOffTargetSearcher()
    params = ZFNDesignParameters(left_half_site=LEFT, right_half_site=RIGHT)
    chrom_seq = "A" * 200

    first = _HalfSiteHit("L", "chr1", 20, 29, Strand.PLUS, LEFT, LEFT, 0, 0, [], LEFT)
    same_strand = _HalfSiteHit("R", "chr1", 34, 43, Strand.PLUS, RIGHT, RIGHT, 0, 0, [], RIGHT)
    assert searcher._build_site(first, same_strand, chrom_seq, {5}, params) is None

    opposite_strand = _HalfSiteHit("R", "chr1", 34, 43, Strand.MINUS, RIGHT, RIGHT, 0, 0, [], RIGHT)
    assert searcher._build_site(first, opposite_strand, chrom_seq, {6}, params) is None

    left_left = _HalfSiteHit("L", "chr1", 34, 43, Strand.MINUS, LEFT, LEFT, 0, 0, [], LEFT)
    assert searcher._build_site(first, left_left, chrom_seq, {5}, params) is None


def test_pair_hits_deduplicates_higher_score_for_same_site_key() -> None:
    """Pairing should retain the highest-score site among duplicates with same key."""
    searcher = ExhaustiveZFNOffTargetSearcher()
    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        dimer_mode=DimerMode.INCLUDE_HOMODIMERS,
        algorithm=ZFNAlgorithm.CONSERVED_G,
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5], require_opposite_strands=False),
    )
    chrom_seq = {"chr1": "A" * 400}

    left_hits = [
        _HalfSiteHit("L", "chr1", 10, 19, Strand.PLUS, "GGGGGGGGG", "GGGGGGGGG", 0, 0, [], "GGGGGGGGG"),
        _HalfSiteHit("L", "chr1", 10, 19, Strand.PLUS, "AAAAAAAAA", "AAAAAAAAA", 0, 0, [], "AAAAAAAAA"),
    ]
    right_hits = [_HalfSiteHit("R", "chr1", 24, 33, Strand.MINUS, "GGGGGGGGG", "GGGGGGGGG", 0, 0, [], "GGGGGGGGG")]

    deduped = searcher._pair_hits(left_hits, right_hits, chrom_seq, params)
    assert len(deduped) == 1
    assert deduped[0].score == pytest.approx(100.0)


def test_score_helpers_cover_zero_g_and_non_empty_penalties() -> None:
    """Scoring helper branches should handle zero-G and non-empty mismatch lists."""
    searcher = ExhaustiveZFNOffTargetSearcher()

    left = _HalfSiteHit("L", "chr1", 0, 9, Strand.PLUS, "AAAAAAAAA", "AAAAAAAAA", 0, 0, [], "AAAAAAAAA")
    right = _HalfSiteHit("R", "chr1", 14, 23, Strand.MINUS, "TTTTTTTTT", "TTTTTTTTT", 0, 0, [], "TTTTTTTTT")
    assert searcher._score_conserved_g(left, right) == searcher._score_homology(left, right)

    penalized = _HalfSiteHit(
        "L",
        "chr1",
        0,
        9,
        Strand.PLUS,
        "ACGTACGTA",
        "ACGTACGTA",
        2,
        1,
        [1, 7],
        "AcGTACGtA",
    )
    penalty = searcher._half_site_v2_penalty(penalized, is_left=True)
    assert penalty > 0.0


def test_paired_sites_include_score_component_breakdown(tmp_path: Path) -> None:
    """Predicted sites should expose explainable score components for deterministic ranking."""
    sequence = f"AAA{_canonical_site(LEFT, RIGHT)}CCC"
    fasta = _write_fasta(tmp_path, sequence)
    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        algorithm=ZFNAlgorithm.HOMOLOGY,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )

    sites = ExhaustiveZFNOffTargetSearcher().search(params)
    assert sites
    components = sites[0].score_components
    assert components
    assert "algorithm_weighted_penalty" in components
    assert "mismatch_penalty" in components
    assert "seed_penalty" in components
    assert components["algorithm_weighted_penalty"] == pytest.approx(0.0)


def test_designer_filter_and_empty_candidate_branches() -> None:
    """Designer should exercise filter-failure strings and empty-site candidate paths."""
    designer = ZFNDesigner()
    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        off_target_filters=ZFNOffTargetFilterCriteria(max_total_sites=0, min_site_score_to_count=70.0),
    )

    empty_candidate = designer._build_candidate(params, [])
    assert empty_candidate.predicted_sites_total == 0
    assert empty_candidate.worst_site_score is None
    assert empty_candidate.best_offtarget_score is None

    site_exon = ZFNOffTargetSite(
        site_id="site_exon",
        chrom="chr1",
        start_1based=1,
        end_1based=23,
        strand=Strand.PLUS,
        orientation=MatchOrientation.LR,
        spacer_len=5,
        sequence="N" * 23,
        left_mismatches=1,
        right_mismatches=1,
        total_mismatches=2,
        score=90.0,
        region="exon",
        nearest_gene="GENE",
        left_aligned="x",
        right_aligned="y",
    )

    assert "TOTAL_SITES_EXCEEDED" in designer._apply_filters(
        ZFNOffTargetFilterCriteria(max_total_sites=0), [site_exon], [site_exon], []
    )
    assert "EXONIC_SITES_EXCEEDED" in designer._apply_filters(
        ZFNOffTargetFilterCriteria(max_exonic_sites=0), [site_exon], [site_exon], []
    )
    assert "PROMOTER_SITES_EXCEEDED" in designer._apply_filters(
        ZFNOffTargetFilterCriteria(max_promoter_sites=0), [site_exon], [], [site_exon]
    )


def test_designer_tool_versions_contains_expected_keys() -> None:
    """Tool-version map should include sirnaforge and python version keys."""
    versions = ZFNDesigner()._tool_versions()
    assert "sirnaforge" in versions
    assert "python" in versions


def test_search_path_with_annotation_provider_and_top_n(tmp_path: Path) -> None:
    """Searcher should annotate results and enforce top_n truncation."""

    class _StubAnnotationProvider:
        def annotate(
            self,
            sites: list[ZFNOffTargetSite],
            _annotation: GenomicAnnotationConfig,
        ) -> list[ZFNOffTargetSite]:
            updated: list[ZFNOffTargetSite] = []
            for site in sites:
                updated.append(site.model_copy(update={"region": "promoter"}))
            return updated

    sequence = f"AAA{_canonical_site(LEFT, RIGHT, spacer='AAAAA')}TTT{_canonical_site(LEFT, RIGHT, spacer='AAAAA')}CCC"
    fasta = _write_fasta(tmp_path, sequence)
    annotation_path = tmp_path / "annotations.gtf"
    annotation_path.write_text("# gtf\n", encoding="utf-8")

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        top_n_sites=1,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )
    annotation = GenomicAnnotationConfig(annotation_path=str(annotation_path))

    sites = ExhaustiveZFNOffTargetSearcher(annotation_provider=_StubAnnotationProvider()).search(
        params, annotation=annotation
    )
    assert len(sites) == 1
    assert all(site.region == "promoter" for site in sites)


def test_designer_default_annotation_provider_applies_gtf_annotations(tmp_path: Path) -> None:
    """Default designer/searcher wiring should annotate sites from a valid GTF."""
    site_seq = _canonical_site(LEFT, RIGHT, spacer="AAAAA")
    fasta = _write_fasta(tmp_path, f"TTTT{site_seq}CCCC")

    gtf_path = tmp_path / "annotation.gtf"
    gtf_path.write_text(
        (
            'chr1\ttest\tgene\t1\t200\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";\n'
            'chr1\ttest\texon\t1\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";\n'
        ),
        encoding="utf-8",
    )

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        top_n_sites=10,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )
    annotation = GenomicAnnotationConfig(annotation_path=str(gtf_path))

    result = ZFNDesigner().evaluate_pair(params=params, annotation=annotation)

    assert result.off_target_sites
    assert any(site.region == "exon" for site in result.off_target_sites)
    assert any(site.nearest_gene == "GENE1" for site in result.off_target_sites)


def test_designer_annotation_matches_chr_aliases_between_fasta_and_gtf(tmp_path: Path) -> None:
    """Annotation should work when FASTA and GTF differ only by chr prefix."""
    site_seq = _canonical_site(LEFT, RIGHT, spacer="AAAAA")
    fasta = _write_fasta(tmp_path, f"TTTT{site_seq}CCCC", name="3")

    gtf_path = tmp_path / "annotation_chr3.gtf"
    gtf_path.write_text(
        (
            'chr3\ttest\tgene\t1\t200\t.\t+\t.\tgene_id "GENE3"; gene_name "GENE3";\n'
            'chr3\ttest\texon\t1\t200\t.\t+\t.\tgene_id "GENE3"; transcript_id "TX3"; gene_name "GENE3";\n'
        ),
        encoding="utf-8",
    )

    params = ZFNDesignParameters(
        search_space_fasta=str(fasta),
        left_half_site=LEFT,
        right_half_site=RIGHT,
        top_n_sites=10,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0),
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
    )
    annotation = GenomicAnnotationConfig(annotation_path=str(gtf_path))

    result = ZFNDesigner().evaluate_pair(params=params, annotation=annotation)

    assert result.off_target_sites
    assert any(site.region == "exon" for site in result.off_target_sites)
    assert any(site.nearest_gene == "GENE3" for site in result.off_target_sites)


def test_resolve_search_space_fasta_raises_when_source_and_fasta_missing() -> None:
    """Resolver should raise when both explicit fasta and named source are absent."""
    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        search_space_fasta=None,
        search_space_reference=None,
    )
    with pytest.raises(ValueError, match="Provide either search_space_fasta or search_space_reference"):
        ExhaustiveZFNOffTargetSearcher()._resolve_search_space_fasta(params)


def test_resolve_annotation_path_returns_none_when_unset() -> None:
    """Annotation resolver should return None when no annotation is configured."""
    assert ExhaustiveZFNOffTargetSearcher()._resolve_annotation_path(GenomicAnnotationConfig()) is None


def test_build_site_supports_rr_orientation_branch() -> None:
    """RR pairing in homodimer mode should build a site through non-L branch selection."""
    searcher = ExhaustiveZFNOffTargetSearcher()
    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        dimer_mode=DimerMode.INCLUDE_HOMODIMERS,
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5], require_opposite_strands=False),
    )
    chrom_seq = "A" * 200
    first = _HalfSiteHit("R", "chr1", 20, 29, Strand.PLUS, RIGHT, RIGHT, 0, 0, [], RIGHT)
    second = _HalfSiteHit("R", "chr1", 34, 43, Strand.PLUS, RIGHT, RIGHT, 0, 0, [], RIGHT)

    site = searcher._build_site(first, second, chrom_seq, {5}, params)
    assert site is not None
    assert site.orientation == MatchOrientation.RR


def test_pair_hits_homodimer_skips_cross_chromosome_pairs() -> None:
    """Homodimer loops should ignore cross-chromosome pairings."""
    searcher = ExhaustiveZFNOffTargetSearcher()
    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        dimer_mode=DimerMode.INCLUDE_HOMODIMERS,
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5], require_opposite_strands=False),
    )
    left_hits = [
        _HalfSiteHit("L", "chr1", 10, 19, Strand.PLUS, LEFT, LEFT, 0, 0, [], LEFT),
        _HalfSiteHit("L", "chr2", 20, 29, Strand.PLUS, LEFT, LEFT, 0, 0, [], LEFT),
    ]
    right_hits = [
        _HalfSiteHit("R", "chr1", 24, 33, Strand.MINUS, RIGHT, RIGHT, 0, 0, [], RIGHT),
        _HalfSiteHit("R", "chr2", 34, 43, Strand.MINUS, RIGHT, RIGHT, 0, 0, [], RIGHT),
    ]
    chrom_sequences = {"chr1": "A" * 200, "chr2": "C" * 200}

    sites = searcher._pair_hits(left_hits, right_hits, chrom_sequences, params)
    assert sites
    assert all(site.chrom in {"chr1", "chr2"} for site in sites)


def test_build_site_second_kind_l_branch_is_reached() -> None:
    """Build-site should map left/right hits from second arg when second.kind is L."""
    searcher = ExhaustiveZFNOffTargetSearcher()
    params = ZFNDesignParameters(
        left_half_site=LEFT,
        right_half_site=RIGHT,
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=[5], require_opposite_strands=True),
    )
    chrom_seq = "A" * 200

    first = _HalfSiteHit("R", "chr1", 34, 43, Strand.MINUS, RIGHT, RIGHT, 0, 0, [], RIGHT)
    second = _HalfSiteHit("L", "chr1", 20, 29, Strand.PLUS, LEFT, LEFT, 0, 0, [], LEFT)
    site = searcher._build_site(first, second, chrom_seq, {5}, params)

    assert site is not None
    assert site.left_mismatches == second.mismatches
    assert site.right_mismatches == first.mismatches


def test_score_manufacturability_applies_homopolymer_repeat_penalty() -> None:
    """Manufacturability should apply repeat penalty when AAAA/TTTT/CCCC/GGGG is present."""
    designer = ZFNDesigner()
    params = ZFNDesignParameters(
        left_half_site="AAAACCCCC",
        right_half_site="GGGGTTTTT",
    )
    score = designer._score_manufacturability(params)
    assert score == pytest.approx(85.0)
