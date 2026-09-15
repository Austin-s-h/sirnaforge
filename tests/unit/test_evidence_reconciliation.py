"""``core/screening_evidence.py``: envelope IO, plan building, and plan/observed reconciliation.

Reconciliation is what turns "no envelope was written" into a FAILED record instead of a silent
absence, and what stops a censored (truncated) search from being read as complete evidence. These
tests pin the join key, the four statuses' completion semantics, and the schema-version guard that
makes an unrecognised payload unobserved rather than a crash.
"""

import hashlib
import json
from pathlib import Path

import pytest

from sirnaforge.core.screening_evidence import (
    MISSING_EVIDENCE_DETAIL,
    EvidenceEnvelope,
    EvidenceProducer,
    EvidenceSource,
    build_plan,
    censored_counts,
    collect_evidence,
    completed_pairs,
    completed_pairs_without_plan,
    evidence_filename,
    failed_entry,
    guide_set_digest,
    join_key,
    not_requested_entry,
    parse_reconciliation_payload,
    read_evidence,
    reconcile,
    reconciliation_payload,
    write_evidence,
    write_reconciliation,
)
from sirnaforge.models.evidence import (
    EVIDENCE_SCHEMA_VERSION,
    EvidenceStatus,
    ObservedCount,
    ObservedCounts,
    ScreeningEvidence,
    ScreeningEvidenceEntry,
    ScreeningPlan,
    ScreeningPlanEntry,
)
from sirnaforge.models.policy import ScreeningChannel

DIGEST = "3f9c1a2b4d5e6f70"
OTHER_DIGEST = "0011223344556677"


def _complete_envelope(
    *,
    species: str = "human",
    channel: ScreeningChannel = ScreeningChannel.TRANSCRIPTOME,
    sites: int = 3,
    digest: str = DIGEST,
    source: EvidenceSource = EvidenceSource.ENVELOPE,
) -> EvidenceEnvelope:
    return EvidenceEnvelope(
        producer=EvidenceProducer.OFFTARGET_ANALYSIS,
        source=source,
        entry=ScreeningEvidenceEntry(
            channel=channel,
            species=species,
            guide_set_digest=digest,
            status=EvidenceStatus.COMPLETE,
            counts=ObservedCounts(sites=ObservedCount(value=sites)),
        ),
    )


@pytest.mark.unit
def test_a_species_with_no_envelope_reconciles_to_failed_with_a_reason():
    """The rat entry was planned but never published; it must fail, not vanish."""
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None), ("rat", None)], mirna_species=[])

    result = reconcile(plan, observed=[_complete_envelope(species="human")])

    rat_entries = [e for e in result.evidence.entries if e.species == "rat"]
    assert len(rat_entries) == 1
    assert rat_entries[0].status is EvidenceStatus.FAILED
    assert rat_entries[0].detail
    assert completed_pairs(result) == frozenset({("transcriptome", "human")})


@pytest.mark.unit
def test_missing_entries_carry_the_default_detail_and_it_is_exported():
    """The reconciler's synthesized detail is the documented constant, not an ad-hoc string."""
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("rat", None)], mirna_species=[])

    result = reconcile(plan, observed=[])

    assert result.evidence.entries[0].detail == MISSING_EVIDENCE_DETAIL
    assert result.sources[join_key(plan.entries[0])] is EvidenceSource.SYNTHESIZED


@pytest.mark.unit
def test_a_zero_hit_complete_envelope_round_trips_and_counts_as_completed(tmp_path: Path):
    """A zero is a measurement, not an absence: it must survive write/read and count as completion."""
    entry = ScreeningEvidenceEntry(
        channel=ScreeningChannel.TRANSCRIPTOME,
        species="human",
        guide_set_digest=DIGEST,
        status=EvidenceStatus.COMPLETE,
        counts=ObservedCounts(sites=ObservedCount(value=0)),
    )
    path = write_evidence(tmp_path, producer=EvidenceProducer.OFFTARGET_ANALYSIS, entry=entry)

    envelope = read_evidence(path)

    assert envelope is not None
    assert envelope.entry.counts.sites.value == 0
    assert envelope.entry.counts.sites.is_lower_bound is False
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None)], mirna_species=[])
    assert completed_pairs(reconcile(plan, observed=[envelope])) == frozenset({("transcriptome", "human")})


@pytest.mark.unit
def test_a_censored_envelope_is_excluded_from_completed_pairs():
    """A lower bound cannot show a ceiling was respected, so censored evidence never completes."""
    entry = ScreeningEvidenceEntry(
        channel=ScreeningChannel.TRANSCRIPTOME,
        species="human",
        guide_set_digest=DIGEST,
        status=EvidenceStatus.CENSORED,
        counts=censored_counts(retained=500, cap=500, pre_cap=None),
        detail="retained 500 of an unknown total: max_hits=500 truncated the list",
    )
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None)], mirna_species=[])
    envelope = EvidenceEnvelope(producer=EvidenceProducer.OFFTARGET_ANALYSIS, entry=entry)

    assert completed_pairs(reconcile(plan, observed=[envelope])) == frozenset()
    assert completed_pairs_without_plan(ScreeningEvidence(entries=(entry,))) == frozenset()
    assert entry.counts.sites.is_lower_bound is True
    assert entry.counts.sites.truncated is True


@pytest.mark.unit
def test_read_evidence_returns_none_for_an_unrecognised_schema_version(tmp_path: Path):
    """A version this reader does not know must never be read as complete."""
    path = tmp_path / "transcriptome_human_evidence.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": "1",
                "producer": "offtarget_analysis",
                "source": "envelope",
                "entry": {
                    "channel": "transcriptome",
                    "species": "human",
                    "guide_set_digest": DIGEST,
                    "status": "complete",
                },
            }
        )
    )

    assert read_evidence(path) is None


@pytest.mark.unit
def test_read_evidence_returns_none_for_missing_and_malformed_files(tmp_path: Path):
    """Missing, unreadable and malformed are all "unobserved", never a raised exception."""
    assert read_evidence(tmp_path / "does_not_exist_evidence.json") is None

    malformed = tmp_path / "transcriptome_mouse_evidence.json"
    malformed.write_text("{not json")
    assert read_evidence(malformed) is None

    wrong_shape = tmp_path / "transcriptome_dog_evidence.json"
    wrong_shape.write_text(json.dumps({"schema_version": EVIDENCE_SCHEMA_VERSION, "producer": "offtarget_analysis"}))
    assert read_evidence(wrong_shape) is None


@pytest.mark.unit
def test_an_unplanned_envelope_is_kept_and_listed_not_dropped():
    """The miRNA module screens more species than the plan requires; nothing it produced disappears."""
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None)], mirna_species=[])
    unplanned_envelope = _complete_envelope(species="chicken", channel=ScreeningChannel.MIRNA_SEED)

    result = reconcile(plan, observed=[_complete_envelope(species="human"), unplanned_envelope])

    assert join_key(unplanned_envelope.entry) in result.unplanned
    assert any(e.species == "chicken" for e in result.evidence.entries)


@pytest.mark.unit
def test_an_envelope_for_another_guide_set_cannot_complete_the_pair_it_reconciled_failed_against():
    """#100: the digest is in the join identity so one screen's counts cannot be read as another's.

    ``reconcile`` already honours that -- the plan entry fails and the envelope lands in
    ``unplanned`` -- but projecting to (channel, species) used to discard the digest, so the same
    envelope satisfied the very pair whose plan entry had just failed over the mismatch.
    """
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None)], mirna_species=[])
    foreign = _complete_envelope(species="human", digest=OTHER_DIGEST)

    result = reconcile(plan, observed=[foreign])

    assert result.keys_with_status(EvidenceStatus.FAILED) == (("transcriptome", "human", DIGEST),)
    assert join_key(foreign.entry) in result.unplanned
    assert completed_pairs(result) == frozenset()


@pytest.mark.unit
def test_an_envelope_for_a_pair_the_plan_never_names_still_completes_when_the_guide_set_matches():
    """Plan silence is not denial: the restated fallback plan carries transcriptome entries only.

    Gating on planned join keys alone would report a legitimately screened miRNA channel as never
    having run, so an envelope for this run's own guide set counts for a pair the plan is silent
    about -- while a foreign guide set never does, whatever pair it claims (#100).
    """
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None)], mirna_species=[])
    same_guides = _complete_envelope(species="human", channel=ScreeningChannel.MIRNA_SEED)
    other_guides = _complete_envelope(species="mouse", channel=ScreeningChannel.MIRNA_SEED, digest=OTHER_DIGEST)

    result = reconcile(plan, observed=[_complete_envelope(species="human"), same_guides, other_guides])

    assert completed_pairs(result) == frozenset({("transcriptome", "human"), ("mirna_seed", "human")})


@pytest.mark.unit
def test_completed_pairs_without_plan_stays_lenient_so_a_planless_caller_reports_nothing_unscreened():
    """A caller with no plan has no planned digest, and must not read its whole screen as absent."""
    entry = _complete_envelope(species="human", digest=OTHER_DIGEST).entry

    assert completed_pairs_without_plan(ScreeningEvidence(entries=(entry,))) == frozenset({("transcriptome", "human")})


@pytest.mark.unit
def test_strict_drops_a_legacy_summary_entry_through_the_reconciliation_it_is_recorded_in():
    """``strict`` reads the reconciliation's own sources now that it is no longer passed separately."""
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None), ("mouse", None)], mirna_species=[])
    inferred = _complete_envelope(species="mouse", source=EvidenceSource.LEGACY_SUMMARY)

    result = reconcile(plan, observed=[_complete_envelope(species="human"), inferred])

    assert completed_pairs(result) == frozenset({("transcriptome", "human"), ("transcriptome", "mouse")})
    assert completed_pairs(result, strict=True) == frozenset({("transcriptome", "human")})


@pytest.mark.unit
def test_reconciliation_adopts_the_plans_reference_id():
    """An in-container emitter knows only species and an index path; the plan's identity wins."""
    plan_entry = ScreeningPlanEntry(
        channel=ScreeningChannel.TRANSCRIPTOME,
        species="human",
        reference_id="ensembl_114_cdna",
        guide_set_digest=DIGEST,
    )
    envelope = EvidenceEnvelope(
        producer=EvidenceProducer.OFFTARGET_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.TRANSCRIPTOME,
            species="human",
            reference_id=None,
            guide_set_digest=DIGEST,
            status=EvidenceStatus.COMPLETE,
        ),
    )

    result = reconcile(plan=_plan_of(plan_entry), observed=[envelope])

    assert result.evidence.entries[0].reference_id == "ensembl_114_cdna"


def _plan_of(entry: ScreeningPlanEntry) -> "ScreeningPlan":
    return ScreeningPlan(entries=(entry,))


@pytest.mark.unit
def test_join_key_and_filename_are_the_documented_shape():
    """The 3-tuple and filename spelling other modules depend on."""
    entry = ScreeningPlanEntry(channel=ScreeningChannel.TRANSCRIPTOME, species="human", guide_set_digest=DIGEST)

    assert join_key(entry) == ("transcriptome", "human", DIGEST)
    assert evidence_filename(channel=ScreeningChannel.TRANSCRIPTOME, species="human") == (
        "transcriptome_human_evidence.json"
    )


@pytest.mark.unit
def test_guide_set_digest_matches_workflow_pys_own_hash(tmp_path: Path):
    """Two computations of one file's digest that disagreed would silently split one guide set in two."""
    fasta = tmp_path / "candidates.fasta"
    fasta.write_text(">g1\nACGUACGUACGU\n")

    expected = hashlib.sha256(fasta.read_bytes()).hexdigest()[:16]
    assert guide_set_digest(fasta) == expected


@pytest.mark.unit
def test_build_plan_emits_one_entry_per_requested_channel_and_species():
    """An unrequested channel gets no plan entry at all, so NOT_REQUESTED can be told from FAILED."""
    plan = build_plan(
        guide_set_digest=DIGEST,
        transcriptome=[("human", "ensembl_114_cdna"), ("rat", None)],
        mirna_species=["human"],
        search_settings={"max_hits": 500},
    )

    keys = {e.key for e in plan.entries}
    assert (ScreeningChannel.TRANSCRIPTOME.value, "human", "ensembl_114_cdna", DIGEST) in keys
    assert (ScreeningChannel.MIRNA_SEED.value, "human", None, DIGEST) in keys
    assert not any(e.channel is ScreeningChannel.MIRNA_SEED and e.species == "rat" for e in plan.entries)
    assert all(e.search_settings == {"max_hits": 500} for e in plan.entries)


@pytest.mark.unit
def test_censored_counts_rejects_a_pre_cap_below_the_retained_count():
    """A cap cannot have retained more than the search actually produced."""
    with pytest.raises(ValueError, match="pre-cap"):
        censored_counts(retained=10, cap=5, pre_cap=3)


@pytest.mark.unit
def test_failed_and_not_requested_entry_helpers_produce_valid_entries():
    """Convenience constructors used by callers that never see a plan/envelope pair directly."""
    plan_entry = ScreeningPlanEntry(
        channel=ScreeningChannel.TRANSCRIPTOME, species="human", reference_id="idx", guide_set_digest=DIGEST
    )

    failed = failed_entry(plan_entry, "index build ran out of memory")
    assert failed.status is EvidenceStatus.FAILED
    assert failed.reference_id == "idx"

    not_requested = not_requested_entry(ScreeningChannel.MIRNA_SEED, "mouse", DIGEST)
    assert not_requested.status is EvidenceStatus.NOT_REQUESTED
    assert not_requested.counts.observed_units == ()


@pytest.mark.unit
def test_reconciliation_round_trips_through_the_written_file(tmp_path: Path):
    """``write_reconciliation``/``parse_reconciliation_payload`` must agree with each other."""
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None), ("rat", None)], mirna_species=[])
    result = reconcile(plan, observed=[_complete_envelope(species="human")])

    path = write_reconciliation(tmp_path, result)
    assert path.name == "evidence.json"

    payload = json.loads(path.read_text())
    parsed = parse_reconciliation_payload(payload)

    assert parsed is not None
    assert parsed.plan == result.plan
    assert parsed.evidence == result.evidence
    assert parsed.sources == dict(result.sources)
    assert parsed.unplanned == result.unplanned


@pytest.mark.unit
def test_parse_reconciliation_payload_rejects_none_and_wrong_schema_version():
    """A missing or mismatched schema version is unusable, not a partial parse."""
    assert parse_reconciliation_payload(None) is None
    assert parse_reconciliation_payload({"schema_version": "1"}) is None


@pytest.mark.unit
def test_reconciliation_payload_serialises_sources_and_unplanned_in_the_documented_shape():
    """Sources as ``"channel|species|digest"`` strings, unplanned as 3-element arrays."""
    plan = build_plan(guide_set_digest=DIGEST, transcriptome=[("human", None)], mirna_species=[])
    unplanned_envelope = _complete_envelope(species="chicken", channel=ScreeningChannel.MIRNA_SEED)
    result = reconcile(plan, observed=[_complete_envelope(species="human"), unplanned_envelope])

    payload = reconciliation_payload(result)

    assert payload["sources"][f"transcriptome|human|{DIGEST}"] == "envelope"
    assert ["mirna_seed", "chicken", DIGEST] in payload["unplanned"]


@pytest.mark.unit
def test_collect_evidence_reads_every_envelope_under_a_directory_and_ignores_the_reconciliation_file(
    tmp_path: Path,
):
    """``collect_evidence`` globs ``*_evidence.json`` only, so ``evidence.json`` itself is not one."""
    write_evidence(
        tmp_path,
        producer=EvidenceProducer.OFFTARGET_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.TRANSCRIPTOME,
            species="human",
            guide_set_digest=DIGEST,
            status=EvidenceStatus.COMPLETE,
        ),
    )
    nested = tmp_path / "nested"
    nested.mkdir()
    write_evidence(
        nested,
        producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.MIRNA_SEED,
            species="mouse",
            guide_set_digest=DIGEST,
            status=EvidenceStatus.COMPLETE,
        ),
    )
    (tmp_path / "evidence.json").write_text("{}")

    envelopes = collect_evidence(tmp_path)

    assert len(envelopes) == 2
    assert [join_key(e.entry) for e in envelopes] == sorted(join_key(e.entry) for e in envelopes)
