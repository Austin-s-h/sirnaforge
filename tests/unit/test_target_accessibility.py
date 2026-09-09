"""Target-site accessibility: geometry, sign, missing-value policy and measured-efficacy guard.

Issue #95. The `accessibility` term used to fold the **guide against itself** and report the result
as target-site accessibility. It is now a real RNAplfold local-opening probability on the
transcript, and the thing most easily got wrong is which end of the target site is scored.

Geometry, once: target site T[1..L] 5'->3' and guide G[1..L] 5'->3' are antiparallel, so G[i] pairs
T[L+1-i]. The guide seed G[2..8] therefore pairs the **3' end** of the target site. Taking the
8-mer at the other end is a measured near-null (Spearman rho +0.05 vs +0.27 against knockdown), so
a geometry inversion leaves the pipeline working and only makes the score meaningless. The control
in `test_non_seed_end_control_stays_near_null` is what catches it.

Fixture data and its limitations (full provenance in `tests/unit/data/README.md`):

- Measured efficacy: Huesken D, Lange J, Mickanin C, Weiler J, Asselbergs F, Warner J, Meloon B,
  Engel S, Rosenberg A, Cohen D, Labow M, Reinhardt M, Natt F, Hall J. "Design of a genome-wide
  siRNA library using an artificial neural network." Nat Biotechnol. 2005 Aug;23(8):995-1001.
  doi:10.1038/nbt1118. PMID: 16025102. Obtained via the third-party redistribution
  https://github.com/apkrfi/unMod-siRNA-Pred (`raw_data/All_Dataset.csv`). ⚠️ The primary paper is
  not open access, so the values were **not** verified against it: these tests reproduce a
  redistribution, which bounds what they prove.
- Transcripts: NCBI Nucleotide records (`efetch`, `db=nuccore`), not part of the Huesken paper.
- TP53-201 (`ENST00000269305.9`): Ensembl.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.core.scoring import (
    SCORING_WEIGHT_SET_VERSION,
    ScoringError,
    compute_composite,
    target_accessibility_sub_score,
)
from sirnaforge.core.thermodynamics import (
    SEED_ANCHORED_WINDOW_NT,
    SEED_END_WINDOW_NT,
    TargetAccessibilityProfile,
    ThermodynamicCalculator,
)
from sirnaforge.models.sirna import (
    COMPOSITE_TERM_NAMES,
    DEFAULT_ACCESSIBILITY_LOG_FLOOR,
    DEFAULT_PLFOLD_MAX_BP_SPAN,
    DEFAULT_PLFOLD_WINDOW,
    DesignParameters,
    ScoringWeights,
    SiRNACandidate,
    TargetAccessibilityConfig,
)

DATA_DIR = Path(__file__).parent / "data"

# Measured on the vendored 180-row Huesken subset at W=150/L=100 (see tests/unit/data/README.md):
# scored seed-end rho +0.2437, non-seed control +0.0639. The full 2,779-row benchmark gives +0.267
# and +0.067. The bounds below leave headroom for the subset's noise; they are not the measurements,
# and the subset was not tuned to hit them.
SUBSET_MIN_SCORED_RHO = 0.18
SUBSET_MAX_CONTROL_RHO = 0.12
SUBSET_MIN_RHO_SEPARATION = 0.10

# pSUPER/pSuper-Retro anti-p53 target site: 78 Europe PMC full-text hits, the field's de facto
# positive control for p53 knockdown. Given as the target site (sense), unique in TP53-201.
TP53_PSUPER_SITE = "GACTCCAGTGGTAATCTAC"
TP53_PSUPER_OFFSET = 916


def _read_fasta(path) -> dict[str, str]:
    """Minimal FASTA reader keyed on the unversioned accession."""
    sequences: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks)
            name = line[1:].split()[0].split(".")[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def _spearman(a: pd.Series, b: pd.Series) -> float:
    """Spearman rho as Pearson on ranks; keeps scipy out of the test dependencies."""
    return float(np.corrcoef(a.rank().to_numpy(), b.rank().to_numpy())[0, 1])


@pytest.fixture(scope="module")
def benchmark_openings() -> pd.DataFrame:
    """Score the vendored measured-efficacy subset once per module.

    One RNAplfold pass per transcript (8 transcripts, 9,051 nt) covers all 180 rows.
    """
    table = pd.read_csv(DATA_DIR / "sirna_efficacy_subset.csv")
    transcripts = _read_fasta(DATA_DIR / "sirna_efficacy_subset_transcripts.fa")

    rows: list[dict[str, float | str]] = []
    for accession, group in table.groupby("accession"):
        sequence = transcripts[str(accession)].upper()
        profile = TargetAccessibilityProfile.fold(
            sequence,
            window_size=DEFAULT_PLFOLD_WINDOW,
            max_bp_span=DEFAULT_PLFOLD_MAX_BP_SPAN,
            u_max=int(group.guide_sequence.str.len().max()),
        )
        for row in group.itertuples():
            site = str(Seq(row.guide_sequence).reverse_complement())
            start = int(row.site_start)
            # The vendored offsets must still describe the vendored transcripts.
            assert sequence[start : start + len(site)] == site, f"{accession} row moved"
            rows.append(
                {
                    "accession": str(accession),
                    "efficacy": float(row.efficacy),
                    # Scored: the 8 nt at the site's 3' end, which pair guide positions 1-8.
                    "p_seed_end": profile.site_accessibility(start, len(site)).seed_end_8mer,
                    # Control: the 8 nt at the site's 5' end, same length, same fold, wrong end.
                    "p_non_seed_end": profile.site_accessibility(start, SEED_END_WINDOW_NT).seed_end_8mer,
                }
            )

    frame = pd.DataFrame(rows)
    assert len(frame) == len(table), "every vendored row must be scorable"
    return frame


@pytest.mark.unit
def test_scored_term_tracks_measured_knockdown(benchmark_openings):
    """The scored statistic must correlate with measured inhibition, or it is decoration.

    Modest by construction: rho ~+0.27 on the full benchmark is ~7% of rank variance, which is
    consistent with the term's 0.13 weight and is not an argument for weighting it higher.
    """
    rho = _spearman(benchmark_openings.efficacy, benchmark_openings.p_seed_end)
    assert rho >= SUBSET_MIN_SCORED_RHO, f"scored seed-end rho fell to {rho:+.3f}"


@pytest.mark.unit
def test_non_seed_end_control_stays_near_null(benchmark_openings):
    """Geometry guard: the wrong end of the target site must remain a near-null.

    Both windows are 8 nt on the same transcript from the same fold, differing only in which end of
    the site they cover. If this control rises to meet the scored metric, the site is being indexed
    backwards -- the bug this test exists for. The separation assertion is the noise-robust half:
    absolute rho on 180 rows moves around, the ordering does not.
    """
    scored = _spearman(benchmark_openings.efficacy, benchmark_openings.p_seed_end)
    control = _spearman(benchmark_openings.efficacy, benchmark_openings.p_non_seed_end)

    assert control <= SUBSET_MAX_CONTROL_RHO, f"non-seed control rose to {control:+.3f}"
    assert scored - control >= SUBSET_MIN_RHO_SEPARATION, (
        f"scored {scored:+.3f} is no longer clearly ahead of the non-seed control {control:+.3f}; "
        "the target site is probably indexed from the wrong end"
    )


@pytest.mark.unit
def test_published_tp53_construct_is_not_in_the_bottom_decile():
    """A construct this well established must not be scored as one of the least accessible sites.

    This is the check that failed for the whole-site statistic (1.8th percentile) and is what
    forced the seed-end anchoring; it sits at the 93rd percentile now.
    """
    sequence = next(iter(_read_fasta(DATA_DIR / "tp53_201.fa").values())).upper()
    assert sequence.count(TP53_PSUPER_SITE) == 1
    assert sequence.find(TP53_PSUPER_SITE) == TP53_PSUPER_OFFSET

    length = len(TP53_PSUPER_SITE)
    profile = TargetAccessibilityProfile.fold(
        sequence,
        window_size=DEFAULT_PLFOLD_WINDOW,
        max_bp_span=DEFAULT_PLFOLD_MAX_BP_SPAN,
        u_max=length,
    )
    windows = [
        opening
        for start in range(len(sequence) - length + 1)
        if (opening := profile.site_accessibility(start, length).seed_end_8mer) is not None
    ]
    site_value = profile.site_accessibility(TP53_PSUPER_OFFSET, length).seed_end_8mer
    percentile = 100.0 * sum(1 for value in windows if value <= site_value) / len(windows)

    assert percentile > 10.0, f"pSUPER site fell to the {percentile:.1f}th percentile"


@pytest.mark.unit
def test_scored_window_is_the_eight_mer_ending_at_the_site_three_prime_end():
    """Pin the indexing algebraically, independent of any folding outcome.

    The 21-mer site's scored window must be *identical* to the 8-mer site starting 13 nt further
    along -- i.e. both end at the same base. Off-by-one or a 5'-anchored window breaks this.
    """
    sequence = next(iter(_read_fasta(DATA_DIR / "tp53_201.fa").values())).upper()
    profile = TargetAccessibilityProfile.fold(sequence, u_max=21)

    start, length = 500, 21
    site = profile.site_accessibility(start, length)
    shifted = profile.site_accessibility(start + length - SEED_END_WINDOW_NT, SEED_END_WINDOW_NT)
    assert site.seed_end_8mer == shifted.seed_end_8mer

    # Same for the reported 17-mer: it is the 17 nt ending at the site's 3' end.
    seventeen = profile.site_accessibility(start + length - SEED_ANCHORED_WINDOW_NT, SEED_ANCHORED_WINDOW_NT)
    assert site.seed_anchored_17mer == seventeen.whole_site

    # Longer windows are subsets of shorter ones, so opening probability must be monotone decreasing.
    assert site.whole_site <= site.seed_anchored_17mer <= site.seed_end_8mer


def _hairpin_and_loop_transcript() -> tuple[str, int, int]:
    """A transcript with one target site buried in a GC stem and one in an unpairable A run.

    Returns (sequence, buried_site_start, open_site_start). The open region is pure A: A cannot
    pair with A, so there is no partner anywhere for it, which makes the expectation exact rather
    than merely likely.
    """
    stem_arm = "GCGGCCGCGGCCGCGGCCGCG"  # 21 nt, GC-rich
    partner = str(Seq(stem_arm).reverse_complement())

    open_start = 30
    prefix = "A" * open_start + "A" * 21 + "A" * 30 + "A" * 10
    buried_start = len(prefix)
    sequence = prefix + stem_arm + "AAAA" + partner + "A" * 30
    return sequence, buried_start, open_start


@pytest.mark.unit
def test_buried_site_scores_below_an_open_site():
    """Sign test: a site inside a hairpin stem must be less accessible than one in a loop.

    Pins the direction of the term. An inverted sign would make the pipeline prefer exactly the
    sites RISC cannot reach.
    """
    sequence, buried_start, open_start = _hairpin_and_loop_transcript()
    profile = TargetAccessibilityProfile.fold(sequence, u_max=21)

    buried = profile.site_accessibility(buried_start, 21).seed_end_8mer
    opened = profile.site_accessibility(open_start, 21).seed_end_8mer

    assert opened > 0.9, f"an unpairable poly-A site should be wide open, got {opened}"
    assert buried < 0.1, f"a site inside a GC stem should be closed, got {buried}"

    buried_score = target_accessibility_sub_score(buried)
    open_score = target_accessibility_sub_score(opened)
    assert open_score > buried_score
    assert open_score > 0.9
    assert 0.0 <= buried_score < 0.6


@pytest.mark.unit
def test_seed_end_and_non_seed_end_measure_different_ends():
    """A site open at its 3' end and closed at its 5' end must score high, not low.

    Built so the two ends disagree by construction: the site's first 8 nt are one arm of a stem,
    the remaining 13 are unpairable A.
    """
    arm = "GGCCGCGG"  # 8 nt
    partner = str(Seq(arm).reverse_complement())
    start = 30
    sequence = "A" * start + arm + "A" * 33 + partner + "A" * 30
    profile = TargetAccessibilityProfile.fold(sequence, u_max=21)

    seed_end = profile.site_accessibility(start, 21).seed_end_8mer
    non_seed_end = profile.site_accessibility(start, SEED_END_WINDOW_NT).seed_end_8mer

    assert seed_end > non_seed_end, "the open 3' end must score above the paired 5' end"
    assert seed_end > 0.9
    assert non_seed_end < seed_end / 2


@pytest.mark.unit
def test_missing_accessibility_is_none_not_a_default():
    """A value that could not be computed must be None all the way through."""
    assert target_accessibility_sub_score(None) is None

    profile = TargetAccessibilityProfile.fold("A" * 60, u_max=21)
    # A site running past the transcript end, and one too close to the 5' end for the window.
    assert profile.site_accessibility(50, 21).seed_end_8mer is None
    assert profile.site_accessibility(0, 5).seed_end_8mer is None


@pytest.mark.unit
def test_missing_accessibility_never_yields_a_score_at_all():
    """The None policy, as issue #96 leaves it: a missing input yields no score, not a rescaled one.

    Under #95 the term was dropped and the remaining weights renormalised. Issue #96 deleted every
    runtime weight operation, so there is nothing left to renormalise onto: a candidate without
    accessibility evidence has no design_score. That is stricter than the old behaviour and keeps
    the guarantee this test was written for -- a missing input must never score like a good one.
    """
    vector = ScoringWeights().design
    present = {"asymmetry": 0.8, "gc_content": 0.6}

    with pytest.raises(ScoringError, match="target_accessibility"):
        compute_composite(dict(present), vector)

    with_best = compute_composite({**present, "target_accessibility": 1.0}, vector)
    with_worst = compute_composite({**present, "target_accessibility": 0.0}, vector)
    assert with_worst.score < with_best.score

    # The vector sums to 1.0, so a perfect candidate scores 100 by construction, not by rescaling.
    perfect = compute_composite(dict.fromkeys(vector.terms, 1.0), vector)
    assert perfect.score == pytest.approx(100.0)


@pytest.mark.unit
def test_scoring_without_transcript_context_leaves_the_term_inactive():
    """Scoring a candidate with no transcript must omit the term, not invent a value."""
    candidate = SiRNACandidate(
        id="SIRNAF_TP53_1_21",
        transcript_id="ENST00000269305",
        position=1,
        guide_sequence="GTAGATTACCACTGGAGTCA",
        passenger_sequence="TGACTCCAGTGGTAATCTAC",
        gc_content=45.0,
        length=20,
        asymmetry_score=0.0,
    )
    SiRNADesigner(DesignParameters())._score_candidates([candidate])

    assert candidate.target_accessibility_p is None
    assert candidate.target_accessibility_p_17mer is None
    assert candidate.target_accessibility_p_site is None
    assert candidate.score_target_accessibility is None
    assert "target_accessibility" not in candidate.component_scores
    # design_v4 declares the term, so with no evidence for it there is no design_score either.
    assert candidate.design_score is None
    assert candidate.composite_score is None
    # Guide self-structure is still recorded: it is the EXCESS_PAIRING gate input.
    assert candidate.structure is not None
    assert candidate.mfe is not None


@pytest.mark.unit
def test_design_populates_accessibility_from_the_transcript():
    """The design path must fold the real transcript and record all three windows."""
    sequence = next(iter(_read_fasta(DATA_DIR / "tp53_201.fa").values())).upper()
    result = SiRNADesigner(DesignParameters(top_n=25)).design_from_sequence(sequence, "ENST00000269305.9")

    scored = [c for c in result.candidates if c.target_accessibility_p is not None]
    assert scored, "no candidate received a target accessibility value"

    for candidate in scored[:25]:
        assert 0.0 <= candidate.target_accessibility_p <= 1.0
        # Longer windows are subsets, so the probabilities must be ordered.
        assert candidate.target_accessibility_p_site <= candidate.target_accessibility_p_17mer
        assert candidate.target_accessibility_p_17mer <= candidate.target_accessibility_p
        assert candidate.score_target_accessibility is not None
        assert candidate.component_scores["target_accessibility"] == pytest.approx(
            target_accessibility_sub_score(candidate.target_accessibility_p)
        )


@pytest.mark.unit
def test_open_chain_guide_mfe_is_not_treated_as_a_failure():
    """An mfe of 0.0 with an all-dots structure is the open chain, the physical floor.

    It is a correct answer for an unstructured guide, not a sentinel and not a failed fold. A
    candidate reporting it must still be scored normally.
    """
    structure, mfe, paired_fraction = ThermodynamicCalculator().calculate_secondary_structure("A" * 21)

    assert structure == "." * 21
    assert mfe == pytest.approx(0.0)
    assert paired_fraction == 0.0


@pytest.mark.unit
def test_sub_score_normalisation_is_bounded_and_monotone():
    """The feature must land in [0, 1] -- the scorer rejects anything else -- and rise with P."""
    floor = DEFAULT_ACCESSIBILITY_LOG_FLOOR

    assert target_accessibility_sub_score(1.0) == pytest.approx(1.0)
    assert target_accessibility_sub_score(10.0**floor) == pytest.approx(0.0)
    assert target_accessibility_sub_score(10.0 ** (floor / 2)) == pytest.approx(0.5)
    # Below the floor clamps rather than going negative; exact zero (underflow) is the same case.
    assert target_accessibility_sub_score(10.0 ** (floor - 2)) == pytest.approx(0.0)
    assert target_accessibility_sub_score(0.0) == 0.0

    values = [target_accessibility_sub_score(10.0**-e) for e in range(0, 6)]
    assert values == sorted(values, reverse=True)
    assert all(0.0 <= v <= 1.0 for v in values)

    # A configurable floor changes the scale, which is why it belongs in the run manifest.
    assert target_accessibility_sub_score(1e-3, log_floor=-3.0) == pytest.approx(0.0)
    assert target_accessibility_sub_score(1e-3, log_floor=-6.0) == pytest.approx(0.5)


@pytest.mark.unit
def test_sub_score_rejects_impossible_inputs():
    """Out-of-range inputs are upstream bugs and must raise, not be silently clamped."""
    with pytest.raises(ValueError, match="probability"):
        target_accessibility_sub_score(1.5)
    with pytest.raises(ValueError, match="log_floor"):
        target_accessibility_sub_score(0.5, log_floor=0.0)


@pytest.mark.unit
def test_composite_term_set_names_the_quantity_it_computes():
    """Issue #95 renamed the term; the weight-set version must record the break.

    The weight and the version moved again in issue #96 (0.13 renormalised to 0.26 in practice ->
    0.35 declared at the design stage, 0.30 post-screen; 3.0.0 -> 4.0.0). What #95 pinned and this
    still pins is that the term is named for the quantity it computes and that `accessibility`,
    which named the wrong molecule, is gone from every vector.
    """
    assert "target_accessibility" in COMPOSITE_TERM_NAMES
    assert "accessibility" not in COMPOSITE_TERM_NAMES
    for vector in ScoringWeights().all_vectors():
        assert "target_accessibility" in vector.terms
        assert not hasattr(vector, "accessibility")
    assert ScoringWeights().design.target_accessibility == pytest.approx(0.35)
    assert ScoringWeights().postscreen_sirna.target_accessibility == pytest.approx(0.30)
    assert SCORING_WEIGHT_SET_VERSION == "4.0.0"

    candidate_fields = set(SiRNACandidate.model_fields)
    assert "score_target_accessibility" in candidate_fields
    assert "score_accessibility" not in candidate_fields
    assert {"target_accessibility_p", "target_accessibility_p_17mer", "target_accessibility_p_site"} <= candidate_fields


@pytest.mark.unit
def test_window_parameters_are_configurable_with_documented_defaults():
    """W and L are config, not constants: they move a site's accessibility percentile."""
    config = TargetAccessibilityConfig()
    assert (config.window_size, config.max_bp_span) == (DEFAULT_PLFOLD_WINDOW, DEFAULT_PLFOLD_MAX_BP_SPAN)
    assert config.log_floor == DEFAULT_ACCESSIBILITY_LOG_FLOOR
    assert DesignParameters().target_accessibility == config

    custom = DesignParameters(target_accessibility=TargetAccessibilityConfig(window_size=80, max_bp_span=40))
    assert custom.target_accessibility.window_size == 80

    # RNAplfold silently clamps L to W, so the combination is rejected instead of misreported.
    with pytest.raises(ValueError, match="max_bp_span"):
        TargetAccessibilityConfig(window_size=40, max_bp_span=100)


@pytest.mark.unit
def test_window_settings_change_the_stored_probability():
    """A configured window must actually reach the fold, or the knob is decoration."""
    sequence = next(iter(_read_fasta(DATA_DIR / "tp53_201.fa").values())).upper()
    narrow = TargetAccessibilityProfile.fold(sequence, window_size=40, max_bp_span=30, u_max=21)
    wide = TargetAccessibilityProfile.fold(sequence, window_size=240, max_bp_span=160, u_max=21)

    assert narrow.site_accessibility(916, 21).seed_end_8mer != wide.site_accessibility(916, 21).seed_end_8mer

    params = DesignParameters(target_accessibility=TargetAccessibilityConfig(window_size=40, max_bp_span=30))
    default_run = SiRNADesigner(DesignParameters()).design_from_sequence(sequence[:600], "t")
    narrow_run = SiRNADesigner(params).design_from_sequence(sequence[:600], "t")
    by_id = {c.id: c.target_accessibility_p for c in default_run.candidates}
    assert any(by_id[c.id] != c.target_accessibility_p for c in narrow_run.candidates)
