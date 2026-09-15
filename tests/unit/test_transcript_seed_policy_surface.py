"""The policy surface the transcript-seed channel and the isoform intent gates are declared through (#101).

Three separate contracts are pinned here, and they are separate on purpose:

* **What the five new gates ARE** -- id, setting, column, comparator, stage, scope -- so a client
  re-applying a descriptor reads the same number the run read. The scope matters more than usual: an
  intent gate that failed to declare the query species would fall back to every screened species and
  report ``unknown`` on any run where a secondary species failed to align.
* **What they DO by default**, which for the three transcript-seed ceilings is nothing at all. #102's
  rule is that an uncalibrated gate reports rather than rejects, and there is no calibration relating
  a 6/7/8mer site count to knockdown. ``max_excluded_isoform_hits`` is the deliberate exception: it
  restates a prohibition the caller stated, and a required limit uses ``fail``.
* **That the blast radius stayed clean** -- the seven design-stage gates are still the only
  unswitchable ones, every new gate's column is on the candidate row, and the profile content hash
  moved without anything pinning a literal.

The site table's own shape is pinned here too, including the one thing that would make the channel's
null fake if it were wrong: ``site_end == anchor_position``, i.e. the antiparallel map applied
forwards.
"""

import re
import subprocess
from typing import Any

import pandas as pd
import pandera.errors
import pytest

from sirnaforge.config.run_policy import (
    FILTER_SPEC_BY_ID,
    SETTING_BY_KEY,
    EntryPoint,
    ResolvedRunPolicy,
    RunPolicyError,
    declared_filter_ids,
    resolve_run_policy,
    switchable_filter_ids,
)
from sirnaforge.core.filtering import (
    Comparator,
    GateSpec,
    derive_passes_filters,
    evaluate_gates,
)
from sirnaforge.core.target_intent import CoverageStatus
from sirnaforge.core.transcript_seed import (
    COORDINATE_SYSTEM,
    GuideStrand,
    SeedClass,
    SiteRegion,
    SiteStrand,
    TranscriptSeedSite,
)
from sirnaforge.models.policy import (
    DECLARED_FILTER_IDS,
    FilterAction,
    FilterComparator,
    FilterEvaluation,
    FilterStage,
)
from sirnaforge.models.schemas import (
    COVERAGE_STATUS_VALUES,
    GUIDE_STRAND_VALUES,
    SEED_CLASS_VALUES,
    SITE_REGION_VALUES,
    SITE_STRAND_VALUES,
    TRANSCRIPT_SEED_COORDINATE_SYSTEM,
    GenomeAlignmentSchema,
    MiRNAAlignmentSchema,
    SiRNACandidateSchema,
    TranscriptSeedSiteSchema,
)
from sirnaforge.models.sirna import (
    OffTargetFilterCriteria,
    SiRNACandidate,
    build_candidate_row,
)

#: The three uncalibrated transcript-seed ceilings, which ship off, and the two intent gates, which
#: do not share their argument. Kept apart because the two groups have opposite defaults for opposite
#: reasons, and a single list of "the new gates" would hide that.
SEED_CEILING_IDS = (
    "max_transcript_seed_sites",
    "max_transcript_seed_transcripts",
    "max_transcript_seed_genes",
)
INTENT_GATE_IDS = ("max_excluded_isoform_hits", "max_unintended_isoform_hits")
NEW_FILTER_IDS = (*SEED_CEILING_IDS, *INTENT_GATE_IDS)

#: The candidate columns this slice adds. Five transcript-seed counters, two intent counters and the
#: reason an absent coverage fraction is absent.
NEW_CANDIDATE_COLUMNS = (
    "transcript_seed_sites_query",
    "transcript_seed_transcripts_query",
    "transcript_seed_genes_query",
    "transcript_seed_unresolved_gene_sites",
    "transcript_seed_sites_total",
    "excluded_isoform_hits",
    "unintended_isoform_hits",
    "intent_coverage_status",
)

#: The gates that could not be switched off before this slice and still cannot be: seven design-stage
#: thresholds that are plain floats with no absent state.
DESIGN_STAGE_UNSWITCHABLE = {
    "gc_content_min",
    "gc_content_max",
    "max_poly_runs",
    "max_repeat_transcript_fraction",
    "max_paired_fraction",
    "min_asymmetry_score",
    "min_empirical_score",
}


def _candidate(**overrides: Any) -> SiRNACandidate:
    """A minimal, otherwise-default candidate. Defaults are the point of most tests here."""
    fields: dict[str, Any] = {
        "id": "cand-1",
        "transcript_id": "ENST00000000001",
        "position": 1,
        "guide_sequence": "TGAGGTAGTAGGTTGTATAGT",
        "passenger_sequence": "ACTATACAACCTACTACCTCA",
        "gc_content": 42.9,
        "length": 21,
        "asymmetry_score": 0.7,
    }
    fields.update(overrides)
    return SiRNACandidate(**fields)


def _validated(candidate: SiRNACandidate) -> pd.DataFrame:
    """One candidate through ``build_candidate_row`` and the schema.

    The ``astype`` is the same pre-existing workaround ``test_selection_state_vocabulary`` and
    ``test_schemas`` apply, and it is unrelated to this slice: a single-row frame built from the row
    dict infers an all-``None`` column as ``object``, which pandera does not resolve to the declared
    nullable dtype. The workflow repairs exactly this before validating a real table, which is why
    this slice's new counters are declared as nullable floats rather than ``Int64`` -- that repair
    covers float columns generically and needed two hard-coded names for the ``Int64`` pair.
    """
    frame = pd.DataFrame([build_candidate_row(candidate)])
    all_null_objects = [name for name in frame.columns if frame[name].dtype == object and frame[name].isna().all()]
    declared = SiRNACandidateSchema.to_schema().columns
    repairs = {
        name: str(declared[name].dtype)
        for name in all_null_objects
        if name in declared and str(declared[name].dtype).startswith(("float", "Int64"))
    }
    return SiRNACandidateSchema.validate(frame.astype(repairs))


def _gate_specs(policy: ResolvedRunPolicy, filter_ids: tuple[str, ...]) -> tuple[GateSpec, ...]:
    """Build the pure evaluator's gates straight from the resolved descriptors.

    Deliberately not a hand-written threshold anywhere: this is what a client re-applying the
    published descriptors does, so a test that passes here says the *run's* gates behave this way and
    not merely that some gates would.
    """
    comparators = {FilterComparator.LE: Comparator.AT_MOST, FilterComparator.GE: Comparator.AT_LEAST}
    specs = []
    for filter_id in filter_ids:
        descriptor = policy.descriptor(filter_id)
        specs.append(
            GateSpec(
                filter_id=filter_id,
                threshold=descriptor.threshold,
                action=descriptor.action,
                comparator=comparators[descriptor.comparator],
            )
        )
    return tuple(specs)


def _observed_from_row(policy: ResolvedRunPolicy, candidate: SiRNACandidate) -> dict[str, float | None]:
    """Each gate's observation, read from the candidate ROW by the column the descriptor declares."""
    row = build_candidate_row(candidate)
    return {
        resolved.filter_id: row[resolved.descriptor.column]
        for resolved in policy.filters
        if resolved.descriptor.column in row
    }


@pytest.mark.unit
def test_every_new_gate_is_declared_in_all_three_places_it_has_to_be():
    """A gate is a registry entry, a setting name and a model field; two out of three is a defect.

    ``DECLARED_FILTER_IDS`` is the one that bites silently: it lives in ``models.policy`` because the
    candidate row needs it and ``models`` cannot import the registry, so a gate added only to
    ``FILTER_SPECS`` gets no verdict and no observed column on any row.
    """
    for filter_id in NEW_FILTER_IDS:
        spec = FILTER_SPEC_BY_ID[filter_id]
        assert spec.setting_key == filter_id, "the public setting name and the filter id are the same word here"
        assert spec.stage is FilterStage.POST_SCREEN, "every new gate reads screening evidence"
        assert spec.comparator is FilterComparator.LE, "all five are ceilings"
        assert filter_id in DECLARED_FILTER_IDS
        setting = SETTING_BY_KEY[filter_id]
        assert setting.model == "offtarget_filters"
        assert setting.field in OffTargetFilterCriteria.model_fields


@pytest.mark.unit
def test_an_uncalibrated_seed_ceiling_reports_rather_than_rejects():
    """#102's rule, applied to the channel with the least calibration in the build.

    No threshold and an ``off`` action, so the sites are counted, exported and acted on by nothing --
    the same shape ``max_transcriptome_seed_perfect`` ships, and for the same reason: nothing here
    relates a seed-site count to knockdown, and the count scales with the size of the reference, so a
    fixed integer would encode the reference rather than the biology. The definition has to carry that
    argument, because a reader who cannot see why a gate is off will "fix" it.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)

    for filter_id in SEED_CEILING_IDS:
        resolved = policy.filter(filter_id)
        assert resolved.descriptor.action is FilterAction.OFF, filter_id
        assert resolved.descriptor.threshold is None, filter_id
        assert not resolved.is_evaluated, filter_id
        assert "Reported, not enforced" in resolved.definition, filter_id
        assert "16 kb" in resolved.definition, f"{filter_id} must say why a fixed integer encodes the reference"


@pytest.mark.unit
def test_a_default_run_rejects_nothing_new():
    """The whole point of shipping three ceilings off: default behaviour does not change.

    A candidate carrying an enormous number of transcript-seed sites -- which a real full-cDNA scan
    produces, order 10^3-10^4 per guide -- survives a default run untouched, and each seed gate
    records ``not_evaluated`` rather than a pass it did not earn. The two intent counters are 0
    because nothing was declared excluded, so their gates pass trivially.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    candidate = _candidate(
        transcript_seed_sites_query=8412,
        transcript_seed_transcripts_query=3907,
        transcript_seed_genes_query=2611,
        transcript_seed_unresolved_gene_sites=44,
        transcript_seed_sites_total=15003,
    )

    outcomes = evaluate_gates(_gate_specs(policy, NEW_FILTER_IDS), _observed_from_row(policy, candidate))
    by_id = {outcome.filter_id: outcome for outcome in outcomes}

    for filter_id in SEED_CEILING_IDS:
        assert by_id[filter_id].evaluation is FilterEvaluation.NOT_EVALUATED, filter_id
        assert not by_id[filter_id].rejects, filter_id
    for filter_id in INTENT_GATE_IDS:
        assert not by_id[filter_id].rejects, filter_id

    label = derive_passes_filters(outcomes, {filter_id: filter_id for filter_id in NEW_FILTER_IDS})
    assert label is True, "no new gate may demote a default-configured candidate"


@pytest.mark.unit
def test_a_declared_exclusion_is_a_required_limit_and_ships_fail():
    """The one new gate that rejects, and why it is allowed to: the caller stated it.

    ``max_excluded_isoform_hits`` is not an uncalibrated prior. The transcripts came from the caller,
    so a hit on one restates a prohibition rather than guessing a threshold, and #101's own rule is
    that a required limit uses ``fail``. Being on-target by gene taxonomy does not make such a hit
    acceptable, which is why the counter is fed by the intent verdict and not by ``hit_class``.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    resolved = policy.filter("max_excluded_isoform_hits")

    assert resolved.descriptor.action is FilterAction.FAIL
    assert resolved.descriptor.threshold == 0
    assert "a required limit uses fail" in resolved.definition
    assert "on-target by gene taxonomy does not imply acceptable" in resolved.definition
    # The gate is about on-target rows only: no other class can be an isoform of the target gene.
    assert resolved.descriptor.scope.hit_classes == frozenset({"on_target", "undetermined"})

    offender = _candidate(excluded_isoform_hits=1)
    outcomes = evaluate_gates(_gate_specs(policy, INTENT_GATE_IDS), _observed_from_row(policy, offender))
    label = derive_passes_filters(outcomes, {"max_excluded_isoform_hits": SiRNACandidate.FilterStatus.EXCLUDED_ISOFORM})
    assert label is SiRNACandidate.FilterStatus.EXCLUDED_ISOFORM

    # ...and it is inert until an exclusion exists to be about, which is what keeps the default safe.
    clean = evaluate_gates(_gate_specs(policy, INTENT_GATE_IDS), _observed_from_row(policy, _candidate()))
    assert all(outcome.evaluation is not FilterEvaluation.FAIL for outcome in clean)


@pytest.mark.unit
def test_an_unrequired_isoform_ceiling_stays_off_until_a_caller_names_one():
    """'Same gene, not required' is a consequence of asking for selectivity, not a prohibition.

    Off by default, and promoted to a real limit the moment a caller states one -- the same promotion
    ``max_transcriptome_seed_perfect`` gets, so a user who names a ceiling gets a gate that rejects.
    """
    default = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    assert default.descriptor("max_unintended_isoform_hits").action is FilterAction.OFF
    assert default.descriptor("max_unintended_isoform_hits").threshold is None

    opted_in = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"max_unintended_isoform_hits": 2})
    descriptor = opted_in.descriptor("max_unintended_isoform_hits")
    assert descriptor.threshold == 2
    assert descriptor.action is FilterAction.FAIL, "a stated ceiling is a limit the caller chose"

    outcomes = evaluate_gates(
        _gate_specs(opted_in, ("max_unintended_isoform_hits",)),
        _observed_from_row(opted_in, _candidate(unintended_isoform_hits=3)),
    )
    assert outcomes[0].rejects


@pytest.mark.unit
@pytest.mark.parametrize("query_species", ["human", "mouse"])
def test_an_intent_gate_requires_only_the_query_species_pair(query_species: str):
    """The wiring trap this scope exists to defuse, pinned on a non-human query too.

    ``_gate_evidence_pairs`` falls back to every screened species when a gate declares no species, so
    an unscoped intent gate would demand transcriptome evidence for all of them and report ``unknown``
    on any run where a secondary species failed to align -- turning one failed alignment into a
    withheld shortlist. The scope is resolved per run rather than fixed, because the stratification is
    the run's query species and not literally human.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, query_species=query_species)

    for filter_id in NEW_FILTER_IDS:
        scope = policy.descriptor(filter_id).scope
        assert scope.species == frozenset({query_species}), filter_id


@pytest.mark.unit
def test_every_new_gate_can_be_switched_off():
    """Switchable because each threshold is ``int | None``, which is a state the model already has.

    This is not a formality: ``_switch_off_value`` refuses any field with no absent state rather than
    faking an inert extreme, so a ``bool`` or plain-``float`` threshold here would make the gate
    unswitchable and the resolver would say so by raising.
    """
    switchable = set(switchable_filter_ids())
    for filter_id in NEW_FILTER_IDS:
        assert filter_id in switchable, filter_id
        policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={filter_id: "off"})
        assert policy.descriptor(filter_id).action is FilterAction.OFF
        assert policy.value_of(filter_id) is None, "switching a gate off clears the number it compared"


@pytest.mark.unit
def test_every_new_gate_column_is_exported():
    """A gate whose column is not on the row declares evidence nobody can read.

    Both halves are checked: each descriptor's declared column is a key of ``build_candidate_row``,
    and every column this slice adds is there whatever a gate reads -- the two all-species/diagnostic
    counters are not gate inputs and would otherwise be droppable without a failure.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    row = build_candidate_row(_candidate())

    for filter_id in NEW_FILTER_IDS:
        column = policy.descriptor(filter_id).column
        assert column in row, f"{filter_id} declares column {column!r}, which no candidate row carries"
        assert policy.filter(filter_id).evidence_exported is True, filter_id

    for column in NEW_CANDIDATE_COLUMNS:
        assert column in row, column

    schema_columns = SiRNACandidateSchema.to_schema().columns
    for column in NEW_CANDIDATE_COLUMNS:
        assert column in schema_columns, f"{column} is emitted but never validated"
    for filter_id in NEW_FILTER_IDS:
        assert f"{filter_id}_observed" in schema_columns, filter_id


@pytest.mark.unit
def test_a_channel_that_did_not_run_is_unknown_and_not_a_clean_zero():
    """The reason the five seed counters are nullable while every counter beside them is not.

    The channel is opt-in. A 0 on a run that never scanned would read as "screened, no sites", which
    is the fabricated zero #101 is about; ``None`` makes the gate ``unknown``, which withholds from a
    qualified shortlist without rejecting. So a stated ceiling on an unrun channel decides nothing.
    """
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, stated={"max_transcript_seed_sites": 50})
    never_scanned = _candidate()
    assert never_scanned.transcript_seed_sites_query is None
    assert build_candidate_row(never_scanned)["transcript_seed_sites_query"] is None

    unrun = evaluate_gates(
        _gate_specs(policy, ("max_transcript_seed_sites",)), _observed_from_row(policy, never_scanned)
    )
    assert unrun[0].evaluation is FilterEvaluation.UNKNOWN
    assert not unrun[0].rejects

    scanned_clean = evaluate_gates(
        _gate_specs(policy, ("max_transcript_seed_sites",)),
        _observed_from_row(policy, _candidate(transcript_seed_sites_query=0)),
    )
    assert scanned_clean[0].evaluation is FilterEvaluation.PASS, "0 is a screen that found nothing, which is a claim"


@pytest.mark.unit
def test_a_gene_ceiling_says_that_it_is_a_lower_bound():
    """An unresolved gene is reported beside the gene count, not dropped from it.

    A gene-level cap that silently discarded the sites whose transcript resolved to no gene would
    read lower the worse the annotation got. The remainder therefore has its own column, and the
    gate's definition says the count is a declared lower bound while that column is above 0.
    """
    definition = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW).filter("max_transcript_seed_genes")
    assert "RESOLVED" in definition.definition
    assert "transcript_seed_unresolved_gene_sites" in definition.definition
    assert "lower bound" in definition.definition

    row = build_candidate_row(_candidate(transcript_seed_genes_query=12, transcript_seed_unresolved_gene_sites=5))
    assert row["transcript_seed_genes_query"] == 12
    assert row["transcript_seed_unresolved_gene_sites"] == 5, "the remainder survives to the row"


@pytest.mark.unit
def test_no_transcript_seed_count_is_folded_into_the_total_offtarget_sum():
    """Three channels, and the combined ceiling stays the two it documents.

    ``max_total_offtarget_hits`` sums transcriptome and miRNA hits. Adding seed sites to it would
    compare sites with alignments in one number and move a gate that ships ``fail``, so the third
    channel gets its own columns and its own ceilings instead.
    """
    candidate = _candidate(
        transcriptome_hits_0mm_query=1,
        mirna_hits_0mm_seed_query=2,
        total_offtarget_hits_query=3,
        transcript_seed_sites_query=9000,
        transcript_seed_sites_total=12000,
    )
    row = build_candidate_row(candidate)
    assert row["total_offtarget_hits_query"] == 3

    definition = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW).filter("max_transcript_seed_sites")
    assert "never folded into max_total_offtarget_hits" in definition.definition


@pytest.mark.unit
def test_the_coverage_status_column_only_speaks_the_intent_vocabulary():
    """An absent coverage fraction has to say WHY, in one of three declared words.

    ``core.target_intent.CoverageStatus`` owns the vocabulary; the schema restates it because
    ``models`` cannot import ``core``, so the two are pinned against each other here rather than
    trusted to stay aligned. Free text would drift into three spellings of "unknown", and a reader
    who cannot tell "no denominator" from "sequence missing" reads a null coverage as a low one.
    """
    assert set(COVERAGE_STATUS_VALUES) == {status.value for status in CoverageStatus}

    for status in CoverageStatus:
        candidate = _candidate(intent_coverage_status=status.value)
        assert build_candidate_row(candidate)["intent_coverage_status"] == status.value
        assert _validated(candidate).loc[0, "intent_coverage_status"] == status.value

    # "unknown" is the word a reader would reach for, and it is exactly the one that loses the
    # distinction between an absent denominator and an absent sequence, so it must be rejected.
    with pytest.raises(pandera.errors.SchemaError):
        _validated(_candidate(intent_coverage_status="unknown"))


@pytest.mark.unit
def test_a_default_candidate_row_still_validates():
    """Every column here is additive: the default row, with all five counters null, must validate."""
    validated = _validated(_candidate())
    for column in NEW_CANDIDATE_COLUMNS:
        assert pd.isna(validated.loc[0, column]) or validated.loc[0, column] == 0


def _site(**overrides: Any) -> TranscriptSeedSite:
    """The let-7a-like 7mer-m8 site from the geometry module's worked example, anchored at 107."""
    fields: dict[str, Any] = {
        "guide_id": "guide-1",
        "guide_sequence": "TGAGGTAGTAGGTTGTATAGT",
        "queried_strand": GuideStrand.GUIDE,
        "species": "human",
        "transcript_id": "ENST00000000001",
        "transcript_version": "3",
        "gene_id": "ENSG00000000001",
        "site_start": 101,
        "site_end": 107,
        "anchor_position": 107,
        "site_class": SeedClass.SEVENMER_M8,
        "site_strand": SiteStrand.TRANSCRIPT_SENSE,
        "region": SiteRegion.FULL_CDNA,
        "coordinate_system": COORDINATE_SYSTEM,
        "annotation_provenance": "ensembl_cdna:test-reference",
    }
    fields.update(overrides)
    return TranscriptSeedSite(**fields)


def _site_frame(*sites: TranscriptSeedSite) -> pd.DataFrame:
    """Sites as the published table writes them: enum members flattened to their values."""
    return pd.DataFrame(
        [
            {
                field: value.value if hasattr(value, "value") else value
                for field, value in ((name, getattr(site, name)) for name in TranscriptSeedSite.__slots__)
            }
            for site in sites
        ]
    )


@pytest.mark.unit
def test_the_site_schema_covers_the_scanner_s_own_field_set_and_nothing_else():
    """One declared shape for ``<species>_transcript_seed_sites.tsv``, derived from the dataclass.

    ``strict=True``, so an undeclared column means a producer published a site property this contract
    never described. Pinned against ``TranscriptSeedSite.__slots__`` rather than a hand-written list,
    because the failure mode is a field added to the scanner and forgotten here.
    """
    assert set(TranscriptSeedSiteSchema.to_schema().columns) == set(TranscriptSeedSite.__slots__)
    assert TranscriptSeedSiteSchema.to_schema().strict is True

    validated = TranscriptSeedSiteSchema.validate(_site_frame(_site()))
    assert validated.loc[0, "anchor_position"] == 107
    assert validated.loc[0, "coordinate_system"] == TRANSCRIPT_SEED_COORDINATE_SYSTEM


@pytest.mark.unit
def test_the_site_schema_is_not_an_alignment_schema():
    """Deliberately not derived from either aligner table, and it has to stay that way.

    A seed site is an exact string match found by scanning a cDNA reference: there is no CIGAR, no
    mapping quality and no edit distance to report, and its coordinates are 1-based transcript cDNA
    rather than the aligner's 0-based ``coord``. Inheriting would have told a reader of both tables
    that one coordinate convention covered both.
    """
    assert not issubclass(TranscriptSeedSiteSchema, MiRNAAlignmentSchema)
    assert not issubclass(TranscriptSeedSiteSchema, GenomeAlignmentSchema)

    columns = set(TranscriptSeedSiteSchema.to_schema().columns)
    assert columns.isdisjoint({"cigar", "mapq", "nm", "seed_mismatches", "coord", "offtarget_score"})
    assert "coordinate_system" in columns, "the convention is stated on every row, not implied"


@pytest.mark.unit
def test_a_site_whose_anchor_is_not_its_end_is_rejected():
    """The one geometry error that would make the channel's null fake, caught at the artifact boundary.

    Pairing is antiparallel -- ``p(i) = anchor + 2 - i`` -- so the anchor is the site's 3'-most
    complementary base and therefore its last transcript position. A row where the two disagree has
    had the map applied backwards, which finds guide-IDENTICAL (passenger-orientation) windows at the
    same expected rate as real sites: the counts look right and only the identity is wrong.
    """
    with pytest.raises(pandera.errors.SchemaError):
        TranscriptSeedSiteSchema.validate(_site_frame(_site(site_start=107, site_end=113, anchor_position=101)))

    with pytest.raises(pandera.errors.SchemaError):
        TranscriptSeedSiteSchema.validate(_site_frame(_site(site_start=108, site_end=107, anchor_position=107)))


@pytest.mark.unit
def test_the_site_table_vocabularies_match_the_scanner_s_enums():
    """The schema restates four enums, so the restatements are pinned to their originals.

    Same reason as the coverage-status vocabulary: ``models`` cannot import ``core``, so drift is
    prevented by a test rather than by an import. ``site_strand`` has exactly one member on purpose --
    a minus-strand spelling would only mean something for a genome-coordinate scan, and declaring it
    would let a reader believe this scan considered both.
    """
    assert set(GUIDE_STRAND_VALUES) == {member.value for member in GuideStrand}
    assert set(SITE_STRAND_VALUES) == {member.value for member in SiteStrand}
    assert set(SEED_CLASS_VALUES) == {member.value for member in SeedClass}
    assert set(SITE_REGION_VALUES) == {member.value for member in SiteRegion}
    assert TRANSCRIPT_SEED_COORDINATE_SYSTEM == COORDINATE_SYSTEM
    assert SITE_STRAND_VALUES == ("transcript_sense",)


@pytest.mark.unit
def test_only_the_seven_design_stage_gates_are_still_unswitchable():
    """The blast radius: five nullable post-screen thresholds join no unswitchable set.

    ``test_a_design_stage_gate_cannot_be_switched_off_in_0_7_1_and_says_so`` enumerates exactly this
    set, so it is restated here to say that the enumeration was not supposed to move.
    """
    assert set(declared_filter_ids()) - set(switchable_filter_ids()) == DESIGN_STAGE_UNSWITCHABLE


@pytest.mark.unit
def test_the_profile_hash_moves_and_nothing_pins_a_literal_one():
    """Adding gates changes ``content_hash``, which is correct and must stay unpinned.

    ``identity()`` hashes the per-filter default ACTIONS as well as the baseline numbers, exactly so
    that two runs reporting one hash mean one policy. Five new ids with declared actions therefore
    move it -- and a literal hash written into a test, a doc or a fixture would have to be edited by
    hand every time a gate is declared, which is how a hash stops meaning anything.
    """
    identity = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW).profile
    assert identity.content_hash.startswith("sha256:")

    tracked = subprocess.run(  # noqa: S603 - fixed argv, no shell
        ["git", "ls-files", "-z"],  # noqa: S607 - git resolved from PATH, as everywhere else in the suite
        capture_output=True,
        check=True,
        text=True,
    ).stdout.split("\0")
    literal = re.compile(r"sha256:[0-9a-f]{16,}")
    offenders = []
    for path in tracked:
        if not path or not path.endswith((".py", ".md", ".rst", ".json", ".toml", ".yml", ".yaml", ".nf")):
            continue
        try:
            text = open(path, encoding="utf-8").read()  # noqa: SIM115, PTH123
        except (OSError, UnicodeDecodeError):
            continue
        if literal.search(text):
            offenders.append(path)
    assert offenders == [], f"a literal content hash is pinned in {offenders}"


@pytest.mark.unit
def test_a_switched_off_gate_cannot_be_reported_as_enforced():
    """An off-by-default gate may be turned off, and may not be turned on by an action alone.

    The three seed ceilings resolve to off because they have no threshold, and the resolver refuses a
    ``fail`` override rather than publishing a manifest that claims a limit with no number in it.
    """
    for filter_id in SEED_CEILING_IDS:
        with pytest.raises(RunPolicyError, match="declares no threshold"):
            resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, filter_actions={filter_id: "fail"})
