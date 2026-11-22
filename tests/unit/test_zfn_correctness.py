"""Independent geometry, interval, export and backend regressions for issue #82."""

import csv
from pathlib import Path

import pandas as pd
import pytest

from sirnaforge.models.zfn import (
    ZFN_OFFTARGET_CSV_COLUMNS,
    GenomicAnnotationConfig,
    IUPACMode,
    MatchOrientation,
    Strand,
    ZFNDesignParameters,
    ZFNDesignResult,
    ZFNHalfSiteConstraints,
    ZFNOffTargetFilterCriteria,
    ZFNOffTargetSite,
    ZFNSearchBackend,
    ZFNShardingConfig,
    ZFNSpacerConstraints,
)
from sirnaforge.zfn.annotation import GTFZFNAnnotationProvider
from sirnaforge.zfn.design import ZFNDesigner
from sirnaforge.zfn.nextflow_bridge import _site_from_row, aggregate_zfn_shard_results, run_zfn_shard_search
from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher, _HalfSiteHit

# Literal genomic plus-strand sequence: no reverse-complement helper is used to
# manufacture a sequence satisfying the implementation under test.
GENOMIC = "GTCATCCTCATCTACACAAACTGCAAAAG"
LEFT = "GTCATCCTCATC"
RIGHT = "CTTTTGCAGTTT"


@pytest.fixture(autouse=True)
def _isolated_reference_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Keep regression inputs independent of the user's genome cache."""
    monkeypatch.setenv("SIRNAFORGE_CACHE_DIR", str(tmp_path / "cache"))


def _params(tmp_path: Path, sequence: str, **kwargs: object) -> ZFNDesignParameters:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr3\n{sequence}\n")
    defaults = {
        "left_half_site": LEFT,
        "right_half_site": RIGHT,
        "search_space_fasta": str(fasta),
        "search_backend": ZFNSearchBackend.EXHAUSTIVE_PYTHON,
        "half_site_constraints": ZFNHalfSiteConstraints(max_mismatches=1, seed_len_from_fokI=6, seed_max_mismatches=1),
        "spacer_constraints": ZFNSpacerConstraints(allowed_spacer_lengths=[5]),
        "sharding": ZFNShardingConfig(enabled=False),
    }
    defaults.update(kwargs)
    return ZFNDesignParameters(**defaults)


@pytest.mark.parametrize("backend", list(ZFNSearchBackend))
def test_right_foki_seed_follows_genomic_spacer_geometry(tmp_path: Path, backend: ZFNSearchBackend) -> None:
    """The spacer-adjacent right genomic base is LAST in the 5′→3′ right query."""
    # Genomic position 17 is first base after the 5bp spacer; position 28 is distal.
    proximal = GENOMIC[:17] + "C" + GENOMIC[18:]
    distal = GENOMIC[:-1] + "C"
    searcher = ExhaustiveZFNOffTargetSearcher()
    proximal_params = _params(tmp_path, proximal, search_backend=backend)
    prox_sites = searcher.search(proximal_params)
    assert len(prox_sites) == 1
    assert prox_sites[0].right_mismatch_positions == [11]
    assert prox_sites[0].right_seed_mismatches == 1
    assert prox_sites[0].orientation == MatchOrientation.LR
    distal_params = _params(tmp_path, distal, search_backend=backend)
    dist_sites = searcher.search(distal_params)
    assert len(dist_sites) == 1
    assert dist_sites[0].right_mismatch_positions == [0]
    assert dist_sites[0].right_seed_mismatches == 0
    assert prox_sites[0].score < dist_sites[0].score
    strict = proximal_params.model_copy(
        update={
            "half_site_constraints": proximal_params.half_site_constraints.model_copy(update={"seed_max_mismatches": 0})
        }
    )
    _params(tmp_path, proximal)
    assert searcher.search(strict) == []


def test_polarity_is_symmetric_in_query_coordinates() -> None:
    """Both arrays' 3′ bases sit nearest FokI despite opposite genomic strands."""
    searcher = ExhaustiveZFNOffTargetSearcher()
    penalties = {}
    for side, strand in (("L", Strand.PLUS), ("R", Strand.MINUS)):
        for position in (0, 11):
            hit = _HalfSiteHit(side, "chr3", 0, 12, strand, "A" * 12, "A" * 12, 1, 0, [position], "A" * 12)
            penalties[side, position] = searcher._half_site_v2_penalty(hit, is_left=side == "L")
    assert penalties["L", 11] == penalties["R", 11]
    assert penalties["L", 0] == penalties["R", 0]
    assert penalties["R", 11] > penalties["R", 0]


def test_mirrored_pair_is_inward_and_outward_pair_is_rejected(tmp_path: Path) -> None:
    """Swapping genomic orientation changes L...R to R...L while keeping 3′ ends inward."""
    # Literal reverse complement of GENOMIC, independently specified.
    mirrored = "CTTTTGCAGTTTGTGTAGATGAGGATGAC"
    searcher = ExhaustiveZFNOffTargetSearcher()
    params = _params(tmp_path, mirrored, half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0))
    sites = searcher.search(params)
    assert len(sites) == 1
    assert sites[0].orientation == MatchOrientation.RL
    left = _HalfSiteHit("L", "chr3", 0, 12, Strand.MINUS, LEFT, LEFT, 0, 0, [], LEFT)
    right = _HalfSiteHit("R", "chr3", 17, 29, Strand.PLUS, RIGHT, RIGHT, 0, 0, [], RIGHT)
    assert searcher._build_site(left, right, GENOMIC, {5}, params) is None


def test_live_fm_cache_uses_sequence_content(tmp_path: Path) -> None:
    """Equal snapshots reuse an index; a same-length different genome cannot reuse it."""
    params = _params(tmp_path, GENOMIC, search_backend=ZFNSearchBackend.FM_INDEX)
    engine = ExhaustiveZFNOffTargetSearcher()._scan_engine_for(params)
    first = engine._live_multi_index({"chr3": "ACGTACGTACGT"})
    identical = engine._live_multi_index({"chr3": "".join(["ACGT", "ACGT", "ACGT"])})
    changed = engine._live_multi_index({"chr3": "ACGTACGTACGA"})
    assert first is identical
    assert first is not changed


def _site(**updates: object) -> ZFNOffTargetSite:
    values = {
        "site_id": "chr3:5000-5028",
        "chrom": "chr3",
        "start_1based": 5000,
        "end_1based": 5028,
        "strand": Strand.PLUS,
        "orientation": MatchOrientation.LR,
        "spacer_len": 5,
        "sequence": GENOMIC,
        "left_mismatches": 1,
        "right_mismatches": 2,
        "total_mismatches": 3,
        "score": 80,
        "left_aligned": LEFT,
        "right_aligned": RIGHT,
    }
    values.update(updates)
    return ZFNOffTargetSite(**values)


@pytest.mark.parametrize("with_exon,expected", [(True, "exon"), (False, "intron")])
def test_nested_features_keep_containing_interval(tmp_path: Path, with_exon: bool, expected: str) -> None:
    """A short intervening feature cannot hide a large gene/exon that contains the cut."""
    path = tmp_path / "nested.gtf"
    rows = [
        'chr3\ttest\tgene\t1000\t9000\t.\t+\t.\tgene_id "BIG";',
        'chr3\ttest\tgene\t2000\t2100\t.\t+\t.\tgene_id "SMALL";',
        'chr3\ttest\texon\t2000\t2100\t.\t+\t.\tgene_id "SMALL";',
    ]
    if with_exon:
        rows.append('chr3\ttest\texon\t1000\t9000\t.\t+\t.\tgene_id "BIG";')
    path.write_text("\n".join(rows) + "\n")
    sites = GTFZFNAnnotationProvider().annotate([_site()], GenomicAnnotationConfig(annotation_path=str(path)))
    assert sites[0].region == expected
    assert sites[0].nearest_gene == "BIG"
    if with_exon:
        params = _params(tmp_path, GENOMIC, off_target_filters=ZFNOffTargetFilterCriteria(max_exonic_sites=0))
        assert "EXONIC_SITES_EXCEEDED" in str(ZFNDesigner().build_candidate(params, sites).passes_offtarget_filters)


def test_empty_csv_has_readable_canonical_header(tmp_path: Path) -> None:
    """No-hit output is a valid empty table with the same columns as a populated run."""
    result = ZFNDesignResult(parameters=_params(tmp_path, GENOMIC), processing_time_s=0)
    path = tmp_path / "empty.csv"
    result.save_offtargets_csv(str(path))
    frame = pd.read_csv(path)
    assert frame.empty
    assert tuple(frame.columns) == ZFN_OFFTARGET_CSV_COLUMNS


def test_candidate_score_names_follow_offtarget_risk(tmp_path: Path) -> None:
    """The worst concerning site is the strongest predicted off-target."""
    params = _params(tmp_path, GENOMIC)
    candidate = ZFNDesigner().build_candidate(params, [_site(score=80), _site(score=95)])
    assert candidate.worst_site_score == 95
    assert candidate.best_offtarget_score == 80


def test_bridge_roundtrip_preserves_ranking_evidence(tmp_path: Path) -> None:
    """Merging CSV shards must retain fields used by rank_sites and match direct ordering."""
    params = _params(tmp_path, GENOMIC)
    sites = [
        _site(
            site_id="low-priority",
            start_1based=2,
            end_1based=30,
            dimer_compatible=False,
            score_components={"algorithm_weighted_penalty": 2.5},
        ),
        _site(
            site_id="high-priority",
            left_seed_mismatches=1,
            right_seed_mismatches=2,
            left_mismatch_positions=[11],
            right_mismatch_positions=[0, 11],
            score_components={"algorithm_weighted_penalty": 7.5},
        ),
    ]
    result = ZFNDesignResult(parameters=params, processing_time_s=0, off_target_sites=sites)
    shard = tmp_path / "shard.csv"
    result.save_offtargets_csv(str(shard))
    with shard.open() as handle:
        restored = [_site_from_row(row) for row in csv.DictReader(handle)]
    assert [site.model_dump() for site in restored] == [site.model_dump() for site in sites]
    aggregate_zfn_shard_results(
        shard_csv_glob=str(shard),
        output_sites_csv=tmp_path / "aggregate.csv",
        output_summary_json=tmp_path / "aggregate.json",
        params=params,
    )
    with (tmp_path / "aggregate.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert tuple(rows[0]) == ZFN_OFFTARGET_CSV_COLUMNS
    assert rows[0]["site_id"] == "high-priority"
    assert _site_from_row(rows[0]).model_dump() == sites[1].model_dump()


def test_boundary_crossing_site_has_exactly_one_owner(tmp_path: Path) -> None:
    """The site starts in core 1 and ends in core 2, and must survive once."""
    params = _params(
        tmp_path, "A" * 80 + GENOMIC + "A" * 80, half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=0)
    )
    searcher = ExhaustiveZFNOffTargetSearcher()
    full = searcher.search(params)
    assert len(full) == 1
    first = searcher.search_region(params, "chr3", 0, 140, 0, 90, top_n_sites=None)
    second = searcher.search_region(params, "chr3", 40, 189, 90, 189, top_n_sites=None)
    assert [site.model_dump() for site in first + second] == [site.model_dump() for site in full]


def test_actual_bridge_shards_match_direct_result_with_original_parameters(tmp_path: Path) -> None:
    """Typed bridge calls retain custom seeds, full evidence and boundary ownership."""
    params = _params(
        tmp_path,
        "A" * 80 + GENOMIC + "A" * 80,
        half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=1, seed_len_from_fokI=3, seed_max_mismatches=0),
    )
    direct = ExhaustiveZFNOffTargetSearcher().search(params)
    for number, scan_start, scan_end, core_start, core_end in ((0, 1, 140, 1, 90), (1, 41, 189, 91, 189)):
        run_zfn_shard_search(
            shard_id=str(number),
            shard_chrom="chr3",
            scan_start_1=scan_start,
            scan_end_1=scan_end,
            core_start_1=core_start,
            core_end_1=core_end,
            shard_max_mismatches=1,
            left_half_site=LEFT,
            right_half_site=RIGHT,
            genome_fasta=Path(params.search_space_fasta),
            search_backend=params.search_backend,
            search_space_index=None,
            algorithm=params.algorithm,
            dimer_mode=params.dimer_mode,
            spacer_lengths="5",
            annotation_file=None,
            output_sites_csv=tmp_path / f"shard{number}.csv",
            output_summary_json=tmp_path / f"shard{number}.json",
            params=params,
        )
    aggregate_zfn_shard_results(
        shard_csv_glob=str(tmp_path / "shard*.csv"),
        output_sites_csv=tmp_path / "merged.csv",
        output_summary_json=tmp_path / "merged.json",
        params=params,
    )
    with (tmp_path / "merged.csv").open() as handle:
        actual = [_site_from_row(row).model_dump() for row in csv.DictReader(handle)]
    assert actual == [site.model_dump() for site in direct]


@pytest.mark.parametrize("backend", [ZFNSearchBackend.PYAHOCORASICK, ZFNSearchBackend.FM_INDEX])
def test_four_mismatch_budget_and_ambiguous_reference_parity(tmp_path: Path, backend: ZFNSearchBackend) -> None:
    """Canonical pattern expansion remains bounded without losing N-overlapping matches."""
    params = _params(
        tmp_path,
        GENOMIC[:3] + "N" + GENOMIC[4:],
        half_site_constraints=ZFNHalfSiteConstraints(
            max_mismatches=4, seed_len_from_fokI=None, seed_max_mismatches=None
        ),
    )
    searcher = ExhaustiveZFNOffTargetSearcher()
    engine = searcher._scan_engine_for(params.model_copy(update={"search_backend": backend}))
    assert engine._candidate_pattern_count("L", LEFT, params) == 46666
    baseline = searcher.search(params)
    indexed = searcher.search(params.model_copy(update={"search_backend": backend}))
    assert baseline
    assert [site.model_dump() for site in indexed] == [site.model_dump() for site in baseline]


@pytest.mark.parametrize("mode", list(IUPACMode))
def test_ambiguous_runs_keep_scan_backend_parity(tmp_path: Path, mode: IUPACMode) -> None:
    """Long masked runs, IUPAC observed bases and partial regions match the baseline."""
    reference = "NNNNNNNNNNACGTNACGTACGTYACGTNNNNNNNNNNNACGTACGTACGTNNN"
    query = "ACGTNACGT"
    params = _params(
        tmp_path,
        reference,
        half_site_constraints=ZFNHalfSiteConstraints(
            max_mismatches=2, iupac_mode=mode, seed_len_from_fokI=3, seed_max_mismatches=1
        ),
    )
    searcher = ExhaustiveZFNOffTargetSearcher()
    results = []
    for backend in ZFNSearchBackend:
        engine = searcher._scan_engine_for(params.model_copy(update={"search_backend": backend}))
        hits = engine.scan_half_site(
            "R", query, {"chr3": reference}, params, target_chrom="chr3", region_start0=7, region_end0=50
        )
        results.append(
            sorted(
                (hit.start0, hit.strand.value, hit.mismatches, hit.seed_mismatches, hit.mismatch_positions)
                for hit in hits
            )
        )
    assert results[0]
    assert results[0] == results[1] == results[2]
