"""Pin the invariants of the tracked `baseline_0_7_1` fixture.

The fixture is a 292-row / 2,173-hit slice of the frozen 0.7.1 canonical public run — **0.8% of
its candidate rows**. It is therefore a *correctness* fixture, not a *scale* fixture: it pins
relations and column contracts that downstream 0.7.1 work packages consume, and it cannot be used
to exercise a payload budget or any absolute magnitude. Its provenance and verification limits
live in `tests/unit/data/baseline_0_7_1/README.md`; until this module existed, those claims were
prose only.

Thresholds are asserted as the **literals the run used**, never as live model defaults, because
0.7.1 deliberately moves several of those defaults (D3/D9). A frozen artifact must not start
disagreeing with itself when a default changes.
"""

from pathlib import Path

import pandas as pd
import pytest

DATA_DIR = Path(__file__).parent / "data" / "baseline_0_7_1"

# Row counts as cut. A regeneration against a different run changes these, and the README requires
# that to be a deliberate edit rather than a silent refresh.
EXPECTED_CANDIDATE_ROWS = 292
EXPECTED_OFFTARGET_ROWS = 2173
EXPECTED_MIRNA_ROWS = 268
EXPECTED_GUIDES = 28
EXPECTED_TRANSCRIPTS = 34
EXPECTED_MAX_ROWS_PER_GUIDE = 33

# Shape of the run this was cut from, for the scale disclaimer only.
FULL_RUN_CANDIDATE_ROWS = 34_863
FULL_RUN_OFFTARGET_ROWS = 80_520

# Thresholds as configured in the frozen run (MANIFEST.md §5), not as currently defaulted.
AS_RUN_MAX_OFF_TARGET_COUNT = 15
AS_RUN_MIN_ASYMMETRY_SCORE = 0.65

CANDIDATE_COLUMNS = (
    "id",
    "transcript_id",
    "position",
    "guide_sequence",
    "passenger_sequence",
    "gc_content",
    "asymmetry_score",
    "structure",
    "mfe",
    "paired_fraction",
    "target_accessibility_p",
    "target_accessibility_p_17mer",
    "target_accessibility_p_site",
    "duplex_stability_dg",
    "duplex_stability_score",
    "dg_5p",
    "dg_3p",
    "delta_dg_end",
    "melting_temp_c",
    "off_target_screened",
    "off_target_count",
    "off_target_penalty",
    "on_target_hits",
    "ortholog_hits",
    "repeat_hits",
    "ortholog_species",
    "repeat_flagged",
    "repeat_transcript_fraction",
    "transcriptome_hits_total",
    "transcriptome_hits_0mm",
    "transcriptome_hits_1mm",
    "transcriptome_hits_2mm",
    "transcriptome_hits_seed_0mm",
    "on_target_confirmed",
    "mirna_hits_total",
    "mirna_hits_0mm_seed",
    "mirna_hits_1mm_seed",
    "mirna_hits_high_risk",
    "transcript_hit_count",
    "transcript_hit_fraction",
    "isoform_coverage",
    "conservation_score",
    "design_score",
    "composite_score",
    "score_off_target",
    "score_target_accessibility",
    "score_asymmetry",
    "score_gc_content",
    "score_ago_start",
    "score_pos1_mismatch",
    "score_supp_13_16",
    "empirical_score",
    "scored_after_screening",
    "weight_set_version",
    "weight_vector",
    "passes_filters",
    "guide_overhang",
    "guide_modifications",
    "passenger_overhang",
    "passenger_modifications",
    "variant_mode",
    "allele_specific",
    "targeted_alleles",
    "overlapped_variants",
)

HIT_COLUMNS = (
    "qname",
    "qseq",
    "species",
    "rname",
    "coord",
    "strand",
    "cigar",
    "mapq",
    "as_score",
    "nm",
    "seed_mismatches",
    "offtarget_score",
)

MIRNA_COLUMNS = (
    "qname",
    "qseq",
    "species",
    "database",
    "mirna_id",
    "coord",
    "strand",
    "cigar",
    "mapq",
    "as_score",
    "nm",
    "seed_mismatches",
    "offtarget_score",
)

# The four scored terms of `postscreen_sirna_v4`, and the three miRNA-mode terms that must stay
# empty on a siRNA-mode run (see test_score_contributions_decompose_composite_exactly).
SIRNA_SCORE_TERMS = (
    "score_off_target",
    "score_target_accessibility",
    "score_asymmetry",
    "score_gc_content",
)
MIRNA_ONLY_SCORE_TERMS = ("score_ago_start", "score_pos1_mismatch", "score_supp_13_16")

REPEAT_COLUMNS = {"repeat_hits", "repeat_flagged", "repeat_transcript_fraction"}

pytestmark = [pytest.mark.unit]


@pytest.fixture(scope="module")
def candidates() -> pd.DataFrame:
    """The tracked candidate slice, read exactly as a consumer would."""
    return pd.read_csv(DATA_DIR / "candidates_sample.csv")


@pytest.fixture(scope="module")
def offtargets() -> pd.DataFrame:
    """The tracked per-alignment slice."""
    return pd.read_csv(DATA_DIR / "offtarget_hits_sample.tsv", sep="\t")


@pytest.fixture(scope="module")
def mirna_hits() -> pd.DataFrame:
    """The tracked per-miRNA-hit slice."""
    return pd.read_csv(DATA_DIR / "mirna_hits_sample.tsv", sep="\t")


def test_row_counts_match_the_readme(candidates, offtargets, mirna_hits):
    """The counts the README quotes are the counts on disk."""
    assert len(candidates) == EXPECTED_CANDIDATE_ROWS
    assert len(offtargets) == EXPECTED_OFFTARGET_ROWS
    assert len(mirna_hits) == EXPECTED_MIRNA_ROWS


def test_fixture_is_a_correctness_slice_not_a_scale_slice(candidates, offtargets):
    """0.8% of the run: nothing here can stand in for a full-scale payload or magnitude."""
    assert len(candidates) / FULL_RUN_CANDIDATE_ROWS < 0.01
    assert len(offtargets) / FULL_RUN_OFFTARGET_ROWS < 0.03


def test_candidate_column_layout_and_dtypes(candidates):
    """Column set, order and broad dtype class are what `candidates_all.csv` consumers expect."""
    assert tuple(candidates.columns) == CANDIDATE_COLUMNS

    for column in ("id", "transcript_id", "guide_sequence", "passenger_sequence", "passes_filters"):
        assert pd.api.types.is_string_dtype(candidates[column]), column
        assert candidates[column].notna().all(), column

    for column in ("off_target_count", "transcriptome_hits_total", "mirna_hits_0mm_seed", "mirna_hits_high_risk"):
        assert pd.api.types.is_integer_dtype(candidates[column]), column

    for column in ("off_target_screened", "repeat_flagged", "on_target_confirmed", "scored_after_screening"):
        assert pd.api.types.is_bool_dtype(candidates[column]), column

    for column in ("composite_score", "design_score", *SIRNA_SCORE_TERMS):
        assert pd.api.types.is_float_dtype(candidates[column]), column
        assert candidates[column].notna().all(), column


def test_hit_table_column_layout_and_dtypes(offtargets, mirna_hits):
    """Per-hit layouts match `combined_offtargets.tsv` / `combined_mirna_hits.tsv`."""
    assert tuple(offtargets.columns) == HIT_COLUMNS
    assert tuple(mirna_hits.columns) == MIRNA_COLUMNS

    for frame in (offtargets, mirna_hits):
        for column in ("coord", "mapq", "as_score", "nm", "seed_mismatches"):
            assert pd.api.types.is_integer_dtype(frame[column]), column
        for column in ("qname", "qseq", "species", "strand", "cigar"):
            assert pd.api.types.is_string_dtype(frame[column]), column
        assert pd.api.types.is_float_dtype(frame["offtarget_score"])

    # WP0c's columns are deliberately absent at 0ca8b88 (README limit 1); recomputing per-hit
    # classes is the consumer's job until it lands.
    assert not {"hit_class", "matched_symbol", "symbol_lookup_missing"} & set(offtargets.columns)


def test_transcriptome_hits_total_equals_off_target_count(candidates):
    """#101: one of the two columns carries no information, on every row."""
    assert (candidates["transcriptome_hits_total"] == candidates["off_target_count"]).all()


def test_mirna_double_gate_stays_demonstrable(candidates):
    """#101's redundant double gate: the two columns agree everywhere, on a non-zero population.

    `max_mirna_perfect_seed` reads `mirna_hits_0mm_seed` and `fail_on_high_risk_mirna` reads
    `mirna_hits_high_risk`, so relaxing either alone changes nothing. If this fixture ever loses
    its non-zero rows the trap stops being demonstrable and #101 loses its evidence.
    """
    perfect_seed = candidates["mirna_hits_0mm_seed"]
    high_risk = candidates["mirna_hits_high_risk"]

    assert (perfect_seed == high_risk).all()
    non_zero = int((perfect_seed > 0).sum())
    assert non_zero == 73
    # The gate populations, not merely the columns, are identical.
    assert set(candidates.index[perfect_seed > 0]) == set(candidates.index[high_risk > 0])


def test_guide_to_row_multiplicity_keeps_the_join_hazard_testable(candidates, offtargets):
    """#103: hits are keyed on one representative `id` per guide, so an `id` join drops the rest."""
    assert candidates["guide_sequence"].nunique() == EXPECTED_GUIDES
    assert candidates["transcript_id"].nunique() == EXPECTED_TRANSCRIPTS

    rows_per_guide = candidates["guide_sequence"].value_counts()
    assert rows_per_guide.max() == EXPECTED_MAX_ROWS_PER_GUIDE
    assert rows_per_guide.min() == 1

    hit_qnames = set(offtargets["qname"].astype(str))
    candidate_ids = set(candidates["id"].astype(str))
    assert hit_qnames <= candidate_ids
    # Exactly one id per guide was submitted to the aligner...
    assert len(hit_qnames) == EXPECTED_GUIDES
    carries_hits = candidates["id"].astype(str).isin(hit_qnames)
    # ...so joining evidence on `id` reads 264 of 292 rows (90.4%) as zero-off-target.
    assert int((~carries_hits).sum()) == 264

    # The safe join key is the guide sequence: every hit row's qseq is byte-identical to the
    # representative candidate's guide_sequence (#103/P8).
    guide_by_id = dict(zip(candidates["id"].astype(str), candidates["guide_sequence"], strict=True))
    hit_pairs = zip(offtargets["qname"].astype(str), offtargets["qseq"], strict=True)
    assert all(guide_by_id[qname] == qseq for qname, qseq in hit_pairs)


def test_zero_repeat_columns_are_not_readable_as_no_repeats(candidates):
    """README limit 3: the run recorded `repeat_summary: skipped / reference_unavailable`.

    These zeros are absent evidence, not a negative result — an offline `RepeatDetector` scan of
    the full run flags 424 of its 2,400 guides. Nothing in a row distinguishes "scanned, found
    none" from "never scanned", which is the defect #100 exists to close. If a repeat-evidence
    status column is ever added here, this assertion fails and the README claim must be revisited
    rather than the test relaxed.
    """
    assert (candidates["repeat_hits"] == 0).all()
    assert not candidates["repeat_flagged"].any()
    assert (candidates["repeat_transcript_fraction"] == 0.0).all()
    assert {column for column in candidates.columns if "repeat" in column} == REPEAT_COLUMNS


def test_ortholog_and_conservation_are_structurally_empty(candidates, offtargets, mirna_hits):
    """README limit 2: human-only run, so these files cannot test ortholog recognition."""
    assert (candidates["ortholog_hits"] == 0).all()
    assert candidates["conservation_score"].isna().all()
    assert candidates["ortholog_species"].isna().all()
    assert set(offtargets["species"].unique()) == {"human"}
    assert set(mirna_hits["species"].unique()) == {"hsa"}
    assert set(mirna_hits["database"].unique()) == {"mirgenedb"}


def test_score_contributions_decompose_composite_exactly(candidates):
    """The four `postscreen_sirna_v4` terms sum to `composite_score`; the miRNA terms stay null.

    The layout carries seven `score_*` columns because the writer is shared with miRNA mode. The
    exact decomposition — which section A's variance shares depend on — only holds while the three
    miRNA-mode terms are unpopulated.
    """
    assert set(candidates["weight_set_version"].unique()) == {"4.0.0"}
    assert set(candidates["weight_vector"].unique()) == {"postscreen_sirna_v4"}

    for column in MIRNA_ONLY_SCORE_TERMS:
        assert candidates[column].isna().all(), column

    residual = (candidates[list(SIRNA_SCORE_TERMS)].sum(axis=1) - candidates["composite_score"]).abs()
    assert residual.max() < 1e-9


def test_passes_filters_masks_gates_a_candidate_also_fails(candidates):
    """#103: the single status column reports only the first gate to fire.

    33 rows labelled `TRANSCRIPTOME_PERFECT_MATCH` independently fail the asymmetry floor as well,
    and none of them says so.
    """
    labels = candidates["passes_filters"].value_counts().to_dict()
    assert labels == {"TRANSCRIPTOME_PERFECT_MATCH": 254, "PASS": 28, "MIRNA_PERFECT_SEED": 10}

    low_asymmetry = candidates["asymmetry_score"] < AS_RUN_MIN_ASYMMETRY_SCORE
    assert int(low_asymmetry.sum()) == 33
    masked = candidates.loc[low_asymmetry, "passes_filters"]
    assert set(masked.unique()) == {"TRANSCRIPTOME_PERFECT_MATCH"}
    # Nothing labelled PASS silently fails the floor — the masking is one-directional here.
    assert not ((candidates["passes_filters"] == "PASS") & low_asymmetry).any()


def test_excess_offtarget_population_and_nm_tail(candidates, offtargets):
    """#101/D3: the cap population and the clipped-partial `nm` tail are both represented."""
    over_cap = candidates["off_target_count"] > AS_RUN_MAX_OFF_TARGET_COUNT
    assert int(over_cap.sum()) == 42
    assert candidates.loc[over_cap, "guide_sequence"].nunique() == 9

    nm_counts = offtargets["nm"].value_counts().reindex(range(7), fill_value=0).to_dict()
    assert nm_counts == {0: 1700, 1: 425, 2: 0, 3: 1, 4: 1, 5: 2, 6: 44}
    # nm >= 3 exists only as clipped partials: no full-length 21M record carries nm >= 3.
    assert not ((offtargets["cigar"] == "21M") & (offtargets["nm"] >= 3)).any()


def test_isoform_coverage_is_populated_but_never_gated(candidates):
    """README limit 5: the column is real, the gate was off — not exercised, merely inert."""
    assert candidates["isoform_coverage"].notna().all()
    assert candidates["isoform_coverage"].between(0.0, 1.0).all()
