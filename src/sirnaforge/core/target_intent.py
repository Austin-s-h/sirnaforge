"""Resolved target intent: what was meant, decided separately from what was found (#101).

Issue #101's first finding is that ``run_complete_workflow`` never asks what the user intended. It
retrieves every same-gene transcript, treats all of them as targets, and then lets two independent
things drift:

1. **The coverage denominator.** ``step1_retrieve_transcripts`` returns only transcripts that are
   protein-coding *and* carry a sequence, and coverage was reported against the enumeration map
   intersected with that survivor set. So an isoform dropped by ORF or sequence filtering silently
   left the denominator, and coverage rose. #101 requires the denominator be frozen *before* that
   filtering, which is why it lives on one frozen :class:`~sirnaforge.models.policy.TargetIntent`
   instance rather than in a mutable attribute nothing prevents a later step from reassigning.
2. **Acceptability.** ``ON_TARGET`` is a statement about *gene* taxonomy. In an isoform-selective
   design a hit on the wrong isoform of the right gene is on-target and unwanted, and the classifier
   is correct in both halves of that sentence. So the four-way taxonomy stays exactly as it is --
   ``core/hit_classification.py`` is not edited by #101 -- and acceptability becomes a second,
   separately published record: :class:`IntentAssessment`, produced here.

This module is the pure half of both jobs. It imports the policy vocabulary, the classifier's result
*types*, and the pairing geometry, and nothing else: no workflow, no configuration, no filesystem, no
network. ``tests/unit/test_target_intent.py`` asserts that import set directly, the way
``tests/unit/test_policy_and_evidence_contracts.py`` already does for the contract modules, because
the value of a pure seam is entirely in its being testable without a workflow instance.

It also takes no variant argument, deliberately. ``VariantWorkflowConfig`` and
``resolve_variants_step`` keep sole ownership of allele targeting; intent describes *which
transcripts*, variant configuration describes *which alleles*, and two objects that cannot see each
other cannot disagree. Intent is resolved before variants are resolved, and neither reads the other.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from enum import Enum

from sirnaforge.core.hit_classification import HitClass, HitClassification
from sirnaforge.core.seed_geometry import normalize_guide_sequence, reverse_complement
from sirnaforge.models.policy import (
    DEFAULT_TARGET_SPECIES,
    CoverageMatch,
    OrthologyAssertion,
    TargetIntent,
    TargetSelectivity,
)

__all__ = [
    "PROTEIN_CODING_BIOTYPE",
    "CoverageObservation",
    "CoverageStatus",
    "IntentAssessment",
    "IntentVerdict",
    "evaluate_intent",
    "observe_coverage",
    "resolve_target_intent",
]

# The one biotype a pan-isoform coverage denominator is restricted to. Named here rather than
# inlined because ``step1_retrieve_transcripts`` already filters on the same literal, and the point
# of freezing a denominator is that the two agree on what they are agreeing about.
PROTEIN_CODING_BIOTYPE = "protein_coding"


class CoverageStatus(str, Enum):
    """Whether a coverage fraction exists, and if not, which absence produced it (#101).

    Three spellings, none of which is a number. #101's rule is that missing annotation yields
    *unknown*, never a fabricated value, and "unknown" has to be distinguishable from "zero" at every
    consumer: an unknown coverage withholds a candidate from a qualified shortlist, whereas a
    coverage of 0.0 rejects it against any positive floor. The two absences are reported separately
    because they have different fixes -- one needs a declared target set, the other needs the
    sequences of transcripts this run never fetched.

    Attributes:
        KNOWN: A denominator was frozen and every member's sequence was available, so the fraction is
            a measurement.
        UNKNOWN_NO_DENOMINATOR: No denominator was frozen, so there is nothing to report coverage
            against. Note this is *not* complete coverage, which is the reading an empty set invites.
        UNKNOWN_MISSING_SEQUENCE: A denominator was frozen but at least one member's sequence is
            absent, so the numerator would be a lower bound. Since the floor is greater-or-equal, a
            lower bound can only wrongly reject, so it is withheld rather than published.
    """

    KNOWN = "known"
    UNKNOWN_NO_DENOMINATOR = "unknown_no_denominator"
    UNKNOWN_MISSING_SEQUENCE = "unknown_missing_sequence"


class IntentVerdict(str, Enum):
    """Why a classified hit is or is not acceptable, beside the class rather than instead of it.

    A separate vocabulary from :class:`~sirnaforge.core.hit_classification.HitClass` because the two
    answer different questions and #101 needs both published: ``on_target`` + ``excluded_isoform`` is
    a readable, non-contradictory pair, whereas folding the second into the first would have required
    a fifth hit class that lies about the gene.

    Attributes:
        INTENDED: On-target and wanted.
        EXCLUDED_ISOFORM: On-target by gene taxonomy, forbidden by an explicit exclusion.
        UNINTENDED_ISOFORM: A known same-gene transcript outside the required set, in an
            isoform-selective run.
        CONSERVED_TARGET: An orthologue in a declared target species.
        EXPLORATORY_ORTHOLOG: An orthologue outside the target species. Descriptive only.
        LIABILITY: Not the target gene at all. Already counted by the existing liability gates;
            intent adds no new rejection here.
        UNKNOWN: The transcript identity or the annotation needed to decide was absent.
    """

    INTENDED = "intended"
    EXCLUDED_ISOFORM = "excluded_isoform"
    UNINTENDED_ISOFORM = "unintended_isoform"
    CONSERVED_TARGET = "conserved_target"
    EXPLORATORY_ORTHOLOG = "exploratory_ortholog"
    LIABILITY = "liability"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class IntentAssessment:
    """One hit's acceptability under one intent, as a record rather than a mutation.

    ``acceptable`` is tri-state, and ``None`` is never coerced to ``False``. The distinction is the
    whole point of the type: an unacceptable hit is evidence *against* a candidate, while an
    undecidable one is a gap in this run's annotation and must reach the gate as ``UNKNOWN``, which
    withholds the candidate from a qualified shortlist without rejecting it. Coercing the second to
    the first would let a thin annotation read as a safety finding -- the same defect as a missing UTR
    annotation masquerading as a clean UTR screen.

    Attributes:
        verdict: Which intent category the hit fell into.
        acceptable: ``True``/``False`` when intent decided; ``None`` when it could not.
        reason: One sentence, never blank, naming what decided it. Same rule as an evidence detail:
            a record whose explanation is empty cannot be audited later.
    """

    verdict: IntentVerdict
    acceptable: bool | None
    reason: str


@dataclass(frozen=True, slots=True)
class CoverageObservation:
    """A coverage measurement, or an explicit refusal to make one.

    Attributes:
        covered: Denominator members this guide was shown to cover.
        denominator: The frozen denominator, unchanged by anything downstream of intent resolution.
        unavailable: Denominator members whose sequence was not available to test against.
        fraction: ``len(covered) / len(denominator)``, or ``None`` whenever that would be a lower
            bound or a division by an empty set.
        match: The criterion ``covered`` was decided by.
    """

    covered: frozenset[str]
    denominator: frozenset[str]
    unavailable: frozenset[str]
    fraction: float | None
    match: CoverageMatch

    @property
    def status(self) -> CoverageStatus:
        """Which of the three spellings this observation is, derived rather than stored."""
        if not self.denominator:
            return CoverageStatus.UNKNOWN_NO_DENOMINATOR
        if self.unavailable:
            return CoverageStatus.UNKNOWN_MISSING_SEQUENCE
        return CoverageStatus.KNOWN


def resolve_target_intent(
    *,
    selectivity: TargetSelectivity,
    annotation_universe: Mapping[str, str | None],
    required_transcript_ids: Iterable[str] = (),
    excluded_transcript_ids: Iterable[str] = (),
    enumeration_inputs: Iterable[str] = (),
    sequence_available: Iterable[str] = (),
    target_species: Iterable[str] = (DEFAULT_TARGET_SPECIES,),
    offtarget_screen_species: Iterable[str] = (),
    annotation_provenance: str | None = None,
    orthology_evidence: Sequence[OrthologyAssertion] = (),
    coverage_match: CoverageMatch = CoverageMatch.EXACT_FULL_SITE,
    unknown_biotype_counts_as_coding: bool = False,
) -> TargetIntent:
    """Freeze one run's intent from what the annotation source offered and what the caller declared.

    Called once, immediately after transcript retrieval and *before* ORF validation or design, and
    stored on a single frozen instance. That ordering is the fix for #101's coverage defect: once the
    denominator is a field on a frozen pydantic model, a later step shrinking it is not a convention
    to be maintained but a raised exception.

    The denominator rule differs by selectivity, and that difference is what lets one unchanged
    ``min_isoform_coverage`` floor mean the right thing in both modes:

    - ``PAN_ISOFORM``: every protein-coding member of the annotation universe. Any isoform of the
      gene is wanted, so coverage is measured against all of them -- including the ones ORF or
      sequence filtering is about to remove, which is precisely the set that used to vanish.
    - ``ISOFORM_SELECTIVE`` with a non-empty required set: the required transcripts. Measuring a
      selective design against every isoform would score a perfectly selective guide as poor
      coverage, and the floor would then reject exactly the design that was asked for.
    - ``ISOFORM_SELECTIVE`` with nothing required: an empty denominator, so coverage is unknown. A
      declared selectivity with no declared targets does not name anything to cover, and inventing
      the pan-isoform denominator here would reinstate the assumption #101 removed.

    Exclusions are subtracted from the denominator in both modes. An exclusion is a prohibition, so
    counting an excluded transcript in the denominator would make covering it improve a candidate's
    score -- the same precedence ("excluded beats required") that :func:`evaluate_intent` applies to
    hits, applied to the denominator.

    Args:
        selectivity: Declared by the caller; there is no default anywhere in #101's path.
        annotation_universe: Every transcript the annotation source offered, mapping
            **version-stripped** transcript ID to its biotype, with ``None`` for "the source
            declared none". Version stripping is the caller's job: three ``_strip_version``
            implementations already exist in this codebase and #101 declines to add a fourth, so this
            builder compares the strings it is given and records any mismatch as unresolved.
        required_transcript_ids: Transcripts a guide is required to cover.
        excluded_transcript_ids: Transcripts a guide must not cover.
        enumeration_inputs: Transcripts guides were actually enumerated on. Recorded as a
            diagnostic, never as the coverage numerator.
        sequence_available: Transcripts whose sequence this run holds. Intersected with the frozen
            denominator, because availability outside the denominator is not coverage's business.
        target_species: Species whose transcripts are meant to be knocked down.
        offtarget_screen_species: Species to screen for liabilities. Independent of the above.
        annotation_provenance: Identity of the annotation source, e.g. ``"ensembl_114"`` or
            ``"input_fasta"``.
        orthology_evidence: Recorded cross-species correspondences.
        coverage_match: Criterion the numerator will be measured with.
        unknown_biotype_counts_as_coding: Whether a universe member with no declared biotype joins a
            pan-isoform denominator. ``False`` by default, so a thin annotation cannot inflate the
            denominator; ``True`` exists for an input-FASTA run, where no biotype was ever available
            and excluding every record would leave the denominator empty and coverage unknown for a
            run whose targets the user supplied by hand.

    Returns:
        The frozen intent. Every set on it is version-stripped in the caller's own spelling.
    """
    required = frozenset(required_transcript_ids)
    excluded = frozenset(excluded_transcript_ids)
    universe = frozenset(annotation_universe)

    # Declared and never offered: recorded, because absence from every other set is indistinguishable
    # from compliance and a typo would otherwise read as a satisfied requirement.
    unresolved = (required | excluded) - universe

    if selectivity is TargetSelectivity.ISOFORM_SELECTIVE and required:
        denominator = required - excluded
    elif selectivity is TargetSelectivity.ISOFORM_SELECTIVE:
        denominator = frozenset()
    else:
        denominator = (
            frozenset(
                transcript_id
                for transcript_id, biotype in annotation_universe.items()
                if biotype == PROTEIN_CODING_BIOTYPE or (biotype is None and unknown_biotype_counts_as_coding)
            )
            - excluded
        )

    return TargetIntent(
        selectivity=selectivity,
        target_species=frozenset(target_species),
        offtarget_screen_species=frozenset(offtarget_screen_species),
        required_transcript_ids=required,
        excluded_transcript_ids=excluded,
        annotation_universe=universe,
        coverage_denominator=denominator,
        enumeration_inputs=frozenset(enumeration_inputs),
        unresolved_target_ids=unresolved,
        coverage_sequence_available=frozenset(sequence_available) & denominator,
        coverage_match=coverage_match,
        annotation_provenance=annotation_provenance,
        orthology_evidence=tuple(orthology_evidence),
    )


def evaluate_intent(
    classification: HitClassification,
    *,
    transcript_id: str | None,
    species: str,
    intent: TargetIntent,
    query_species: str,
) -> IntentAssessment:
    """Decide acceptability from intent, taking the classifier's own verdict as descriptive evidence.

    Pure, and deliberately non-invasive: it takes the classifier's output *object* rather than
    re-deriving a class from the hit (two derivations could disagree and only one would be published),
    never mutates it, and never returns a
    :class:`~sirnaforge.core.hit_classification.HitClass`. ``hit_class`` keeps its meaning on every
    row; the verdict is written beside it.

    Precedence, in order, and pinned by test:

    1. **Excluded beats required.** A transcript named in both sets is ``EXCLUDED_ISOFORM``. An
       exclusion is a safety statement and a requirement is an efficacy statement, and the
       conservative reading of a contradictory declaration is the safety one. It is also the only
       resolution that cannot be reached by accident: silently preferring the requirement would let a
       stale required list defeat a freshly added exclusion.
    2. **Selectivity**, which only decides anything for a same-gene hit whose transcript identity is
       actually known.
    3. **Species**, for orthologue hits.

    Absent identity is ``UNKNOWN``, never unacceptable. ``classify_hit`` reaches ``ON_TARGET`` through
    a gene-ID or symbol match in ``_classify_query_species``, consulting no transcript identity at
    all, so an on-target row can legitimately name a transcript this intent has never heard of. That
    is a gap in annotation, not evidence of an unintended isoform, and reporting it as the latter
    would manufacture rejections out of a thin index.

    Args:
        classification: The classifier's result for this hit, used as evidence and not recomputed.
        transcript_id: Version-stripped transcript ID from the hit row's ``rname``, or ``None`` when
            the row carried none.
        species: Canonical species of the hit.
        intent: The run's frozen intent.
        query_species: Canonical query species. Checked rather than decorative: ``classify_hit``
            returns ``ON_TARGET`` only for the query species, so an on-target row in another species
            is a contradiction, and treating it as intended would let a caller smuggle a
            cross-species on-target claim past intent.

    Returns:
        The assessment. ``acceptable`` is ``None`` only for :attr:`IntentVerdict.UNKNOWN`.
    """
    hit_class = classification.hit_class
    named = transcript_id if transcript_id else None

    if hit_class is HitClass.UNDETERMINED:
        return IntentAssessment(
            verdict=IntentVerdict.UNKNOWN,
            acceptable=None,
            reason=f"hit class is undetermined for {species}, so intent has nothing to decide from",
        )

    if hit_class in (HitClass.REPEAT, HitClass.OFF_TARGET):
        return IntentAssessment(
            verdict=IntentVerdict.LIABILITY,
            acceptable=False,
            reason=f"{hit_class.value} hit in {species} is not the target gene, so it is a liability",
        )

    if hit_class is HitClass.ORTHOLOG:
        if species in intent.target_species:
            return IntentAssessment(
                verdict=IntentVerdict.CONSERVED_TARGET,
                acceptable=True,
                reason=f"orthologue hit in {species}, a declared target species",
            )
        return IntentAssessment(
            verdict=IntentVerdict.EXPLORATORY_ORTHOLOG,
            acceptable=True,
            reason=f"orthologue hit in {species}, which is not a declared target species",
        )

    # ON_TARGET from here: on the target gene by taxonomy, which decides nothing about the isoform.
    return _assess_on_target(named, species=species, intent=intent, query_species=query_species)


def _assess_on_target(
    named: str | None,
    *,
    species: str,
    intent: TargetIntent,
    query_species: str,
) -> IntentAssessment:
    """Acceptability of an on-target hit, where the four-way taxonomy has stopped being informative.

    Split out so the branch structure is readable: an on-target hit has four fates and only one of
    them is "intended". The exclusion check comes first because an exclusion outranks a requirement.
    """
    if species != query_species:
        return IntentAssessment(
            verdict=IntentVerdict.UNKNOWN,
            acceptable=None,
            reason=(
                f"hit is on-target in {species} but the query species is {query_species}; the "
                "classifier cannot produce that pair, so intent declines to decide it"
            ),
        )

    if named is not None and named in intent.excluded_transcript_ids:
        return IntentAssessment(
            verdict=IntentVerdict.EXCLUDED_ISOFORM,
            acceptable=False,
            reason=f"{named} is an excluded transcript, so an on-target hit on it is unacceptable",
        )

    if named is None:
        return _assess_unidentified_on_target(species=species, intent=intent)

    if intent.selectivity is TargetSelectivity.ISOFORM_SELECTIVE and intent.required_transcript_ids:
        return _assess_selective_isoform(named, intent=intent)

    return IntentAssessment(
        verdict=IntentVerdict.INTENDED,
        acceptable=True,
        reason=f"{named} is on the target gene and no declaration excludes it",
    )


def _assess_unidentified_on_target(*, species: str, intent: TargetIntent) -> IntentAssessment:
    """An on-target hit whose transcript this row never named: unknown only if identity could matter."""
    if _requires_transcript_identity(intent):
        return IntentAssessment(
            verdict=IntentVerdict.UNKNOWN,
            acceptable=None,
            reason=(
                f"on-target hit in {species} carries no transcript identity, and this intent names "
                "transcripts, so acceptability is unknown rather than unacceptable"
            ),
        )
    return IntentAssessment(
        verdict=IntentVerdict.INTENDED,
        acceptable=True,
        reason=(
            f"on-target hit in {species} under a pan-isoform intent with no exclusions, "
            "so transcript identity cannot change the verdict"
        ),
    )


def _assess_selective_isoform(named: str, *, intent: TargetIntent) -> IntentAssessment:
    """A named same-gene transcript under an isoform-selective intent that required something.

    A transcript the annotation universe knows about and the required set omits is an unintended
    isoform. A transcript in *neither* is unknown: the intent has never heard of it, and the
    difference between "we did not ask for this isoform" and "we have no annotation for this
    identifier" is exactly the difference between a finding and a gap.
    """
    if named in intent.required_transcript_ids:
        return IntentAssessment(
            verdict=IntentVerdict.INTENDED,
            acceptable=True,
            reason=f"{named} is a required transcript of an isoform-selective intent",
        )
    if named in intent.annotation_universe:
        return IntentAssessment(
            verdict=IntentVerdict.UNINTENDED_ISOFORM,
            acceptable=False,
            reason=f"{named} is a known same-gene transcript outside the required set of an isoform-selective intent",
        )
    return IntentAssessment(
        verdict=IntentVerdict.UNKNOWN,
        acceptable=None,
        reason=(
            f"{named} appears in no declared set and in no annotation universe, so an "
            "isoform-selective intent cannot say whether it was wanted"
        ),
    )


def _requires_transcript_identity(intent: TargetIntent) -> bool:
    """Whether this intent's verdicts can depend on which transcript a hit is on.

    They can when something is excluded (any exclusion needs identity to enforce) or when an
    isoform-selective intent names required transcripts. Otherwise a missing identity is harmless and
    reporting ``UNKNOWN`` for it would fill an ordinary pan-isoform run with undecided rows.
    """
    if intent.excluded_transcript_ids:
        return True
    return intent.selectivity is TargetSelectivity.ISOFORM_SELECTIVE and bool(intent.required_transcript_ids)


def observe_coverage(
    guide_sequence: str,
    intent: TargetIntent,
    sequences: Mapping[str, str],
) -> CoverageObservation:
    """Test one guide's actual complementarity against the frozen denominator.

    This replaces the enumeration map as the coverage numerator. ``_store_guide_to_transcripts``
    records which transcripts a guide was *enumerated on*; intersecting that with the protein-coding
    set answered "was this guide generated from that isoform", which is not the same question as
    "does this guide silence that isoform" and differs on exactly the transcripts enumeration never
    saw. Those transcripts are the ones an isoform-selective design cares most about, so the test
    here runs over the whole denominator, including them.

    The criterion is :attr:`~sirnaforge.models.policy.CoverageMatch.EXACT_FULL_SITE`: the guide's
    full-length reverse complement occurring verbatim in the transcript. Full-length rather than the
    seed, because a seed-level criterion counts off-target silencing as coverage. Reverse complement
    rather than the guide, because a guide is antisense to its target -- searching for the guide
    itself finds passenger-orientation sites at the same rate and would produce a numerator that
    looks plausible and means nothing (see :mod:`sirnaforge.core.seed_geometry`).

    A denominator member is unavailable when the intent recorded no sequence for it *or* this mapping
    does not hold one, and the two are unioned rather than reconciled. Any doubt makes coverage
    unknown, which is the safe direction: ``min_isoform_coverage`` is a greater-or-equal floor, so an
    under-counted numerator can only wrongly reject.

    Args:
        guide_sequence: The guide, RNA or DNA, any case.
        intent: The run's frozen intent, supplying the denominator and the recorded availability.
        sequences: Transcript ID to cDNA sequence, for as many transcripts as this run holds. Keys
            outside the denominator are ignored.

    Returns:
        The observation, whose ``fraction`` is ``None`` unless a denominator was frozen and every
        member of it was available.
    """
    denominator = intent.coverage_denominator
    site = reverse_complement(guide_sequence)

    held: set[str] = set()
    for transcript_id in denominator:
        sequence = sequences.get(transcript_id)
        if not sequence:
            continue
        if intent.coverage_sequence_available and transcript_id not in intent.coverage_sequence_available:
            continue
        held.add(transcript_id)

    unavailable = frozenset(denominator - held)
    covered = frozenset(
        transcript_id for transcript_id in held if site in normalize_guide_sequence(sequences[transcript_id])
    )

    fraction = None if (unavailable or not denominator) else len(covered) / len(denominator)
    return CoverageObservation(
        covered=covered,
        denominator=denominator,
        unavailable=unavailable,
        fraction=fraction,
        match=intent.coverage_match,
    )
