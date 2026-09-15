"""Transcript-seed liability: complementary seed sites in transcript sequence (#101).

The third liability channel, and deliberately a **new module** rather than a new argument to the
miRNA seed scanner. ``off_target._MiRNASeedScanner`` / ``scan_mirna_seed_matches`` answer a different
question -- *does this guide resemble a known miRNA guide?* -- and they answer it by matching the
guide's seed window **forward** against a miRNA database (``off_target._prepare_seed_queries`` slices
``guide[seed_start-1:seed_end]`` and hands that to the backends unreversed). That is correct for
their question and wrong for this one: pointed at a cDNA, a forward window finds guide-*identical*,
passenger-orientation windows. In a random transcriptome those occur at the same expected rate as
real target sites, so the resulting null is indistinguishable from a true one by any count. #101 is
emphatic that the existing scanner must not be repointed and relabelled, so nothing in
``core/off_target.py`` is touched and this module carries its own typed result instead.

What is shared is only arithmetic: :mod:`sirnaforge.core.seed_geometry` owns the antiparallel map
``p(i) = anchor + 2 - i`` and the class table, and every coordinate below is computed through
:func:`~sirnaforge.core.seed_geometry.paired_transcript_position` rather than by open-coding
``anchor - 6``, so a future change to the map cannot pass these tests.

WHAT ONE SITE IS
----------------
A site is one *anchor* -- the transcript position paired with guide position 2 -- not one matched
window. The four seed classes are nested windows sharing that base and nothing else, so a 6mer at
``[102..107]`` and an 8mer at ``[101..108]`` are one site, published at its maximal class. That is
why :attr:`TranscriptSeedSite.anchor_position` exists and why it, not ``site_start``, carries the
dedup identity: keying on ``site_start`` publishes one site as two.

WHAT IS COUNTED, AND SEPARATELY
-------------------------------
Sites, distinct transcripts and distinct genes are three different questions with three different
caps, so they are reported through #100's :class:`~sirnaforge.models.evidence.ObservedCounts` and
never collapsed. A site whose transcript has no resolvable gene is routed to
``unresolved_gene_sites`` rather than dropped, and its existence makes ``distinct_genes`` a declared
lower bound -- that is the whole mechanism by which a gene-level ceiling stops silently discarding
what it could not resolve. Nothing is ever summed across species: each species is its own evidence
unit and its own column suffix.

WHY THE CAP IS NOT OPTIONAL IN PRACTICE
---------------------------------------
A 7mer is expected roughly once per 16 kb, so a full human cDNA scan yields ``O(10^3-10^4)`` sites
per guide -- orders of magnitude more rows than the alignment table. Hence a per-guide cap, applied
deterministically (strongest class first), and a ``CENSORED`` status through the same
:func:`~sirnaforge.core.screening_evidence.censored_counts` helper every other capped channel uses,
so a ceiling read against a truncated count reports UNKNOWN instead of a pass.

WHY A MISSING UTR ANNOTATION IS NOT A CLEAN UTR SCREEN
------------------------------------------------------
``region_scope`` beyond :attr:`SiteRegion.FULL_CDNA` needs per-transcript intervals *in cDNA
coordinates*, and nothing in this repository can supply them: ``models/transcript_annotation.py``
carries genomic ``exons``/``cds`` intervals only, and ``data/annotation_manager.py`` is a
reference-file cache, not a coordinate service. So such a request returns ``FAILED`` with a reason
and zero sites, every count left unobserved. #100's vocabulary has no ``UNKNOWN`` status, and
``FAILED`` + detail is its spelling of "we asked and could not answer"; the consequence at the gate
is then automatic and correct, because the pair never enters ``completed_pairs`` and
``evaluate_gate`` returns UNKNOWN rather than passing on a fabricated zero.
"""

from __future__ import annotations

import re
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from sirnaforge.core.screening_evidence import censored_counts
from sirnaforge.core.seed_geometry import (
    A1_POSITION,
    M8_POSITION,
    SEED_START_POSITION,
    normalize_guide_sequence,
    paired_transcript_position,
    reverse_complement,
    seed_window,
)
from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcript_index import TranscriptGeneIndex
from sirnaforge.models.evidence import (
    EvidenceStatus,
    ObservedCount,
    ObservedCounts,
    ScreeningEvidenceEntry,
)
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

__all__ = [
    "ALL_SEED_CLASSES",
    "COORDINATE_SYSTEM",
    "SEED_CLASS_RANK",
    "SIXMER_END_POSITION",
    "GuideStrand",
    "SeedClass",
    "SeedScanScope",
    "SiteRegion",
    "SiteStrand",
    "TranscriptSeedScanResult",
    "TranscriptSeedSite",
    "scan_transcript_seed_sites",
    "transcript_seed_evidence_entry",
]

#: Guide position the 6mer/7mer-A1 window ends at. ``seed_geometry`` names ``M8_POSITION`` (8) and
#: ``SEED_START_POSITION`` (2); this is the third boundary the class table needs, kept here because
#: it is a property of the class table rather than of the pairing map.
SIXMER_END_POSITION = 7

#: Spelled on every site rather than implied, because the aligner-derived tables in this repository
#: publish 0-based read coordinates and these are 1-based transcript coordinates.
COORDINATE_SYSTEM = "transcript_cdna_1based"

#: Splits an Ensembl-style identifier into its stable part and its version suffix.
_VERSIONED_ID = re.compile(r"^(?P<stable>.+?)\.(?P<version>\d+)$")


class GuideStrand(str, Enum):
    """Which strand of the duplex was submitted to the scan.

    Explicit, and on every row, because 0.7.1 submits guide sequences only:
    ``_prepare_offtarget_input`` writes the guide strand, so passenger-strand liability is simply
    not screened. Recording that per site is the difference between an unscreened strand and a
    strand that screened clean.

    Attributes:
        GUIDE: The antisense/guide strand -- the only value reachable in 0.7.1.
        PASSENGER: The sense/passenger strand, reserved for a screen that submits it.
    """

    GUIDE = "guide"
    PASSENGER = "passenger"


class SiteStrand(str, Enum):
    """Which strand of the reference the site sits on.

    A cDNA reference is written on the transcript's own sense strand, so this is always
    ``TRANSCRIPT_SENSE`` here. A minus-strand spelling would only mean something for a
    genome-coordinate scan, and inventing the member now would let a future reader believe this
    scan considered both.

    Attributes:
        TRANSCRIPT_SENSE: The transcript's sense strand, as written in the reference FASTA.
    """

    TRANSCRIPT_SENSE = "transcript_sense"


class SeedClass(str, Enum):
    """The canonical seed-site classes, as nested windows of the guide's 5' seed.

    Attributes:
        SIXMER: Guide positions 2-7 pair.
        SEVENMER_A1: Guide positions 2-7 pair and the transcript carries an unpaired ``A`` at
            ``p(1)``.
        SEVENMER_M8: Guide positions 2-8 pair.
        EIGHTMER: Guide positions 2-8 pair and the transcript carries an unpaired ``A`` at ``p(1)``.
    """

    SIXMER = "6mer"
    SEVENMER_A1 = "7mer-A1"
    SEVENMER_M8 = "7mer-m8"
    EIGHTMER = "8mer"


#: Ordering used to publish the maximal class at one anchor. ``7mer-m8`` outranks ``7mer-A1``
#: deliberately: the m8 base is a real base pair, whereas the A1 is an unpaired target adenosine, so
#: an m8 site rests on one more base of complementarity than an A1 site of the same width.
SEED_CLASS_RANK: Mapping[SeedClass, int] = {
    SeedClass.SIXMER: 0,
    SeedClass.SEVENMER_A1: 1,
    SeedClass.SEVENMER_M8: 2,
    SeedClass.EIGHTMER: 3,
}

#: Every class, i.e. the default request. Not a wildcard: the scope publishes the frozenset itself,
#: so a run's class membership is readable from its evidence rather than assumed from a default.
ALL_SEED_CLASSES: frozenset[SeedClass] = frozenset(SeedClass)


class SiteRegion(str, Enum):
    """Transcript region a scan was scoped to, or a site was found in.

    Only ``FULL_CDNA`` is answerable in 0.7.1 (see the module docstring): the other members exist so
    a *request* for a region this repository cannot resolve is representable, and can therefore be
    refused with a reason instead of silently answered over the whole cDNA.

    Attributes:
        FULL_CDNA: The entire cDNA as written in the reference FASTA.
        UTR3: The 3' UTR only.
        UTR5: The 5' UTR only.
        CDS: The coding sequence only.
        UNKNOWN: Region membership was not determined.
    """

    FULL_CDNA = "full_cdna"
    UTR3 = "utr3"
    UTR5 = "utr5"
    CDS = "cds"
    UNKNOWN = "unknown"


#: Why each unanswerable region scope was refused, naming the annotation that was absent rather than
#: saying "unsupported scope". A reader has to be able to tell a screen that found nothing from a
#: screen that never ran, and the reason is the only thing carrying that.
_REGION_REFUSALS: Mapping[SiteRegion, str] = {
    SiteRegion.UTR3: "3' UTR scope requested and no UTR annotation was available",
    SiteRegion.UTR5: "5' UTR scope requested and no UTR annotation was available",
    SiteRegion.CDS: "CDS scope requested and no cDNA-coordinate CDS annotation was available",
    SiteRegion.UNKNOWN: "an unknown region scope was requested, which no annotation can satisfy",
}


@dataclass(frozen=True, slots=True)
class TranscriptSeedSite:
    """One deduplicated transcript-seed site, at hit-level identity.

    Attributes:
        guide_id: The screen query name (``SiRNACandidate.screen_query_id``), so a site joins the
            alignment table on the same key rather than on a sequence.
        guide_sequence: The normalised DNA guide, so a row is interpretable without the FASTA.
        queried_strand: Which strand was submitted. See :class:`GuideStrand`.
        species: Canonical species name.
        transcript_id: Version-stripped transcript identifier.
        transcript_version: Version as written in the FASTA header; ``None`` means the header
            carried none, which is different from an unknown version.
        gene_id: Version-stripped gene identifier from the transcript index, or ``None``. A ``None``
            here is counted in ``unresolved_gene_sites``, never dropped.
        site_start: 1-based inclusive start of the *complementary* window.
        site_end: 1-based inclusive end, always equal to ``anchor_position``.
        anchor_position: Transcript position paired with guide position 2. The site identity.
        site_class: The maximal *requested* class realised at this anchor.
        site_strand: Reference strand. See :class:`SiteStrand`.
        region: Region this site was screened in.
        coordinate_system: Always :data:`COORDINATE_SYSTEM`; stated, not implied.
        annotation_provenance: Which reference produced this site, e.g.
            ``ensembl_cdna:<reference_id>``.
    """

    guide_id: str
    guide_sequence: str
    queried_strand: GuideStrand
    species: str
    transcript_id: str
    transcript_version: str | None
    gene_id: str | None
    site_start: int
    site_end: int
    anchor_position: int
    site_class: SeedClass
    site_strand: SiteStrand
    region: SiteRegion
    coordinate_system: str
    annotation_provenance: str

    @property
    def identity(self) -> tuple[str, str, str, str, str | None, int]:
        """The tuple two rows must share to be the same site.

        ``anchor_position`` and not ``site_start``: the classes are nested windows whose starts
        differ by one base while their anchor is shared, so a ``site_start`` identity publishes one
        8mer site as an 8mer plus a 6mer. ``transcript_version`` is included and is ``None``-tolerant
        rather than defaulted, because two releases of one transcript are two sites in two
        references. ``queried_strand`` is included because a guide-strand and a passenger-strand site
        at one anchor are two different liabilities.
        """
        return (
            self.guide_id,
            self.queried_strand.value,
            self.species,
            self.transcript_id,
            self.transcript_version,
            self.anchor_position,
        )


@dataclass(frozen=True, slots=True)
class SeedScanScope:
    """What a scan actually asked, published beside what it found.

    A site count is meaningless without it: the same guide scores an order of magnitude higher
    against a full cDNA set than against a canonical-transcript set, and higher again if 6mers are
    counted. So the classes, the region, the cap and the reference are recorded rather than left to
    a reader's assumption about defaults.

    Attributes:
        species: Canonical species name this scan covers. One species per scope, never aggregated.
        region: Region requested. See :class:`SiteRegion`.
        classes: Classes requested. A class outside this set is not reported even where realised.
        max_sites_per_guide: Retention cap per guide, or ``None`` for uncapped.
        reference_id: Identity of the cDNA reference searched, when the caller knows it.
    """

    species: str
    region: SiteRegion
    classes: frozenset[SeedClass]
    max_sites_per_guide: int | None
    reference_id: str | None


@dataclass(frozen=True, slots=True)
class TranscriptSeedScanResult:
    """What one species' transcript-seed scan found, and whether it could answer at all.

    Attributes:
        sites: Deduplicated sites, in guide order then cap order. Deterministic by construction.
        counts: #100 counts at each separately-reported aggregation unit.
        status: ``COMPLETE``, ``CENSORED`` when the cap truncated, ``FAILED`` when the scan could
            not be performed as scoped.
        detail: Reason, required by :class:`~sirnaforge.models.evidence.ScreeningEvidenceEntry` for
            ``FAILED`` and ``CENSORED``.
        scope: The question that was asked.
        submitted_guides: Guides handed to the scan.
        processed_guides: Guides actually screened. Lower than ``submitted_guides`` when a guide was
            too short to carry the seed windows; carried here rather than inferred, so a partially
            screened guide set is visible in the envelope.
    """

    sites: tuple[TranscriptSeedSite, ...]
    counts: ObservedCounts
    status: EvidenceStatus
    detail: str | None
    scope: SeedScanScope
    submitted_guides: int
    processed_guides: int


@dataclass(frozen=True, slots=True)
class _PreparedGuide:
    """One guide's precomputed search strings, in the target's orientation.

    Both strings are reverse complements, which is the entire point of the channel: what occurs
    verbatim in a cDNA is what the guide pairs with, never the guide's own window. ``None`` means the
    requested class set needs no search of that width.
    """

    guide_id: str
    sequence: str
    sixmer_site: str | None
    m8_site: str | None


def _split_version(identifier: str) -> tuple[str, str | None]:
    """Split ``ENST00000123456.7`` into its stable id and its version, or ``(id, None)``."""
    match = _VERSIONED_ID.match(identifier.strip())
    if match is None:
        return identifier.strip(), None
    return match.group("stable"), match.group("version")


def _iter_fasta_records(fasta_path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` one record at a time, never materialising the file.

    ``FastaUtils.read_fasta`` returns a list, so a human cDNA set would be held whole in memory for
    the duration of a scan. One record at a time is the most a substring search needs, and a cDNA
    reference is exactly the file where the difference matters.
    """
    header: str | None = None
    chunks: list[str] = []
    with Path(fasta_path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            elif header is not None and line:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def _prepare_guides(guides: Mapping[str, str], classes: frozenset[SeedClass]) -> tuple[_PreparedGuide, ...]:
    """Precompute each guide's target-orientation search strings, skipping guides too short to screen.

    A guide shorter than :data:`~sirnaforge.core.seed_geometry.M8_POSITION` cannot carry the m8
    window, and screening it for 6mers alone would publish a class-incomplete zero that reads exactly
    like a clean result. It is therefore dropped from ``processed_guides`` and named in a warning
    instead. Guide lengths in this repository are 19-23 nt, so this is a guard, not a code path.
    """
    wants_sixmer_family = bool(classes & {SeedClass.SIXMER, SeedClass.SEVENMER_A1})
    wants_m8_family = bool(classes & {SeedClass.SEVENMER_M8, SeedClass.EIGHTMER})

    prepared: list[_PreparedGuide] = []
    for guide_id, raw_sequence in guides.items():
        sequence = normalize_guide_sequence(raw_sequence)
        if len(sequence) < M8_POSITION:
            logger.warning(
                f"Guide {guide_id} is {len(sequence)} nt, shorter than the {M8_POSITION}-nt seed window; "
                "it was not screened for transcript-seed sites"
            )
            continue
        prepared.append(
            _PreparedGuide(
                guide_id=guide_id,
                sequence=sequence,
                sixmer_site=(
                    reverse_complement(seed_window(sequence, SEED_START_POSITION, SIXMER_END_POSITION))
                    if wants_sixmer_family
                    else None
                ),
                m8_site=(
                    reverse_complement(seed_window(sequence, SEED_START_POSITION, M8_POSITION))
                    if wants_m8_family
                    else None
                ),
            )
        )
    return tuple(prepared)


def _realised_classes(*, is_m8_window: bool, has_a1: bool) -> frozenset[SeedClass]:
    """Classes realised at one anchor by a match of the given width.

    An m8-window match realises the m8 classes only. It does not also report a 6mer, because the
    6mer window is searched in its own pass and the two candidates are merged by anchor identity --
    that merge is what makes nested classes one site rather than two.
    """
    if is_m8_window:
        return frozenset({SeedClass.EIGHTMER, SeedClass.SEVENMER_M8}) if has_a1 else frozenset({SeedClass.SEVENMER_M8})
    return frozenset({SeedClass.SEVENMER_A1, SeedClass.SIXMER}) if has_a1 else frozenset({SeedClass.SIXMER})


def _site_start_for(site_class: SeedClass, anchor: int) -> int:
    """1-based start of the complementary window a published class rests on.

    Computed through the pairing map rather than as ``anchor - 6``, so the coordinate cannot drift
    away from :mod:`sirnaforge.core.seed_geometry`'s definition of it. The A1 base is *not* part of
    the window: it sits at ``p(1) = anchor + 1``, past the site's end, and is an unpaired adenosine
    rather than a paired base.
    """
    if site_class in {SeedClass.SEVENMER_M8, SeedClass.EIGHTMER}:
        return paired_transcript_position(M8_POSITION, anchor)
    return paired_transcript_position(SIXMER_END_POSITION, anchor)


def _unanswerable_region_detail(region: SiteRegion, species: str) -> str:
    """Why a region-scoped request was refused, naming the annotation that was absent."""
    reason = _REGION_REFUSALS.get(region, f"{region.value} scope requested and no annotation was available")
    return f"{reason} for {species}; no sites were screened"


def _censoring_detail(*, cap: int, discarded: int, pre_cap: int, guides: int) -> str:
    """Why a truncated scan's counts are lower bounds, in the wording the other channels use."""
    return (
        f"per-guide site cap {cap} truncated the scan: {discarded} of {pre_cap} sites were discarded "
        f"across {guides} guide(s), so every reported count is a lower bound"
    )


def _complete_counts(
    sites: tuple[TranscriptSeedSite, ...],
    *,
    resolves_genes: bool,
) -> ObservedCounts:
    """Counts for an untruncated scan, with the gene unit left unobserved when nothing resolved it.

    ``distinct_genes`` is deliberately ``None`` -- unobserved -- when no transcript index was
    supplied. Reporting ``0`` would be a fabricated zero of exactly the kind #100's count model
    exists to forbid: it would pass any gene ceiling on a scan that never attempted gene resolution.
    When an index *was* supplied and some transcripts did not resolve, the count is real but partial,
    so it is emitted as a declared lower bound and the unresolved sites are published beside it.
    """
    unresolved = sum(1 for site in sites if site.gene_id is None)
    resolved_genes = {site.gene_id for site in sites if site.gene_id is not None}
    return ObservedCounts(
        sites=ObservedCount(value=len(sites)),
        distinct_transcripts=ObservedCount(value=len({site.transcript_id for site in sites})),
        distinct_genes=(
            ObservedCount(value=len(resolved_genes), is_lower_bound=unresolved > 0)
            if resolves_genes
            else ObservedCount()
        ),
        unresolved_gene_sites=ObservedCount(value=unresolved),
    )


def _censored_scan_counts(
    sites: tuple[TranscriptSeedSite, ...],
    *,
    pre_cap: int,
    resolves_genes: bool,
) -> ObservedCounts:
    """Counts for a truncated scan, through the shared helper every capped channel uses.

    ``cap=None`` on purpose: the cap in force is *per guide*, and
    :class:`~sirnaforge.models.evidence.ObservedCount` validates ``value <= cap``, so passing the
    per-guide cap for an aggregate over many guides would either raise or require inventing an
    aggregate ceiling (``cap x guides``) that was never configured. The cap that was configured is
    published where it belongs, on :attr:`SeedScanScope.max_sites_per_guide`, and named again in the
    detail.
    """
    unresolved = sum(1 for site in sites if site.gene_id is None)
    resolved_genes = {site.gene_id for site in sites if site.gene_id is not None}
    return censored_counts(
        retained=len(sites),
        cap=None,
        pre_cap=pre_cap,
        distinct_transcripts=len({site.transcript_id for site in sites}),
        distinct_genes=len(resolved_genes) if resolves_genes else None,
        unresolved_gene_sites=unresolved,
    )


def _anchors(sequence: str, search_string: str) -> Iterator[int]:
    """Every 1-based anchor at which ``search_string`` occurs in ``sequence``.

    The anchor is the match's last base in 1-based coordinates, which is ``p(2)`` for both class
    widths: the 6mer window spans ``anchor-5..anchor`` and the m8 window ``anchor-6..anchor``, so
    both end on the anchor. Overlapping occurrences are all yielded -- a step of ``len(match)`` would
    drop a second site in a repeat.
    """
    start = sequence.find(search_string)
    while start != -1:
        yield start + len(search_string)
        start = sequence.find(search_string, start + 1)


def scan_transcript_seed_sites(
    guides: Mapping[str, str],
    cdna_fasta: str | Path,
    *,
    species: str,
    region_scope: SiteRegion = SiteRegion.FULL_CDNA,
    classes: frozenset[SeedClass] = ALL_SEED_CLASSES,
    max_sites_per_guide: int | None = None,
    gene_index: TranscriptGeneIndex | None = None,
    queried_strand: GuideStrand = GuideStrand.GUIDE,
    reference_id: str | None = None,
) -> TranscriptSeedScanResult:
    """Find complementary seed sites for every guide in one species' cDNA reference (#101).

    Streams the FASTA record by record, searching each requested class width for the guide's
    *reverse-complemented* seed window -- the string a guide actually pairs with. Candidates are then
    merged on :attr:`TranscriptSeedSite.identity`, so nested classes at one anchor become one site
    published at its maximal requested class.

    Args:
        guides: Guide id (the screen query name) to guide sequence, RNA or DNA.
        cdna_fasta: The already-materialised cDNA reference. Never downloaded or indexed here: the
            caller shares the reference the transcriptome channel already resolved.
        species: Species name; normalised to canonical form for every site and for the scope.
        region_scope: Region to screen. Anything but :attr:`SiteRegion.FULL_CDNA` is refused with a
            reason, because no cDNA-coordinate UTR annotation exists in this repository.
        classes: Classes to report. A realised class outside this set is not published, and a width
            no requested class needs is not searched.
        max_sites_per_guide: Retention cap per guide. When it truncates, the result is ``CENSORED``.
        gene_index: Optional transcript→gene index. Absent, every site is unresolved and the gene
            unit is reported as unobserved rather than as zero.
        queried_strand: Which strand these guides are. 0.7.1 submits guides only.
        reference_id: Reference identity for provenance, e.g. an Ensembl release tag.

    Returns:
        A :class:`TranscriptSeedScanResult` whose ``status`` distinguishes a clean scan, a truncated
        one and a request that could not be answered at all.
    """
    canonical_species = normalize_species_name(species)
    requested_classes = frozenset(classes)
    scope = SeedScanScope(
        species=canonical_species,
        region=region_scope,
        classes=requested_classes,
        max_sites_per_guide=max_sites_per_guide,
        reference_id=reference_id,
    )
    provenance = f"ensembl_cdna:{reference_id}" if reference_id else f"index_sidecar_fasta:{Path(cdna_fasta)}"

    if region_scope is not SiteRegion.FULL_CDNA:
        return TranscriptSeedScanResult(
            sites=(),
            counts=ObservedCounts(sites=ObservedCount(value=0)),
            status=EvidenceStatus.COMPLETE,
            detail=_unanswerable_region_detail(region_scope, canonical_species),
            scope=scope,
            submitted_guides=len(guides),
            processed_guides=0,
        )

    if not requested_classes:
        return TranscriptSeedScanResult(
            sites=(),
            counts=ObservedCounts(),
            status=EvidenceStatus.FAILED,
            detail=f"no seed classes were requested for {canonical_species}; no sites were screened",
            scope=scope,
            submitted_guides=len(guides),
            processed_guides=0,
        )

    prepared = _prepare_guides(guides, requested_classes)
    if not prepared:
        return TranscriptSeedScanResult(
            sites=(),
            counts=ObservedCounts(),
            status=EvidenceStatus.FAILED,
            detail=(
                f"none of the {len(guides)} submitted guide(s) could be screened for transcript-seed sites "
                f"in {canonical_species}; no sites were screened"
            ),
            scope=scope,
            submitted_guides=len(guides),
            processed_guides=0,
        )

    species_index = gene_index.for_species(canonical_species) if gene_index is not None else None
    by_guide: dict[str, dict[tuple[str, str, str, str, str | None, int], TranscriptSeedSite]] = {
        guide.guide_id: {} for guide in prepared
    }

    for header, raw_sequence in _iter_fasta_records(cdna_fasta):
        sequence = normalize_guide_sequence(raw_sequence)
        transcript_id, transcript_version = _split_version(header.split()[0]) if header.split() else ("", None)
        gene_id = species_index.gene_id_for(transcript_id) if species_index is not None else None
        for guide in prepared:
            for search_string, is_m8_window in ((guide.sixmer_site, False), (guide.m8_site, True)):
                if search_string is None:
                    continue
                for anchor in _anchors(sequence, search_string):
                    a1_position = paired_transcript_position(A1_POSITION, anchor)
                    has_a1 = a1_position <= len(sequence) and sequence[a1_position - 1] == "A"
                    publishable = _realised_classes(is_m8_window=is_m8_window, has_a1=has_a1) & requested_classes
                    if not publishable:
                        continue
                    site_class = max(publishable, key=lambda member: SEED_CLASS_RANK[member])
                    site = TranscriptSeedSite(
                        guide_id=guide.guide_id,
                        guide_sequence=guide.sequence,
                        queried_strand=queried_strand,
                        species=canonical_species,
                        transcript_id=transcript_id,
                        transcript_version=transcript_version,
                        gene_id=gene_id,
                        site_start=_site_start_for(site_class, anchor),
                        site_end=anchor,
                        anchor_position=anchor,
                        site_class=site_class,
                        site_strand=SiteStrand.TRANSCRIPT_SENSE,
                        region=SiteRegion.FULL_CDNA,
                        coordinate_system=COORDINATE_SYSTEM,
                        annotation_provenance=provenance,
                    )
                    existing = by_guide[guide.guide_id].get(site.identity)
                    if existing is None or SEED_CLASS_RANK[site_class] > SEED_CLASS_RANK[existing.site_class]:
                        by_guide[guide.guide_id][site.identity] = site

    retained: list[TranscriptSeedSite] = []
    pre_cap_total = 0
    truncated_guides = 0
    for guide in prepared:
        guide_sites = sorted(
            by_guide[guide.guide_id].values(),
            key=lambda site: (-SEED_CLASS_RANK[site.site_class], site.transcript_id, site.anchor_position),
        )
        pre_cap_total += len(guide_sites)
        if max_sites_per_guide is not None and len(guide_sites) > max_sites_per_guide:
            truncated_guides += 1
            guide_sites = guide_sites[:max_sites_per_guide]
        retained.extend(guide_sites)

    sites = tuple(retained)
    resolves_genes = species_index is not None
    if truncated_guides and max_sites_per_guide is not None:
        return TranscriptSeedScanResult(
            sites=sites,
            counts=_censored_scan_counts(sites, pre_cap=pre_cap_total, resolves_genes=resolves_genes),
            status=EvidenceStatus.CENSORED,
            detail=_censoring_detail(
                cap=max_sites_per_guide,
                discarded=pre_cap_total - len(sites),
                pre_cap=pre_cap_total,
                guides=truncated_guides,
            ),
            scope=scope,
            submitted_guides=len(guides),
            processed_guides=len(prepared),
        )

    return TranscriptSeedScanResult(
        sites=sites,
        counts=_complete_counts(sites, resolves_genes=resolves_genes),
        status=EvidenceStatus.COMPLETE,
        detail=None,
        scope=scope,
        submitted_guides=len(guides),
        processed_guides=len(prepared),
    )


def transcript_seed_evidence_entry(
    result: TranscriptSeedScanResult,
    *,
    guide_set_digest: str,
    reference_id: str | None = None,
) -> ScreeningEvidenceEntry:
    """One #100 evidence entry for a transcript-seed scan, with no new contract.

    The channel needs no change to :func:`~sirnaforge.core.screening_evidence.reconcile` or
    :func:`~sirnaforge.core.screening_evidence.completed_pairs`: the ``(channel, species, digest)``
    join already generalises, and a ``FAILED`` scan keeps its pair out of ``completed_pairs``, which
    is exactly how a transcript-seed ceiling reports UNKNOWN instead of passing on a zero it never
    observed.

    Args:
        result: The scan to report.
        guide_set_digest: Digest of the guide FASTA, from
            :func:`~sirnaforge.core.screening_evidence.guide_set_digest`.
        reference_id: Reference identity to publish; defaults to the scan's own scope.

    Returns:
        The entry, ready for :func:`~sirnaforge.core.screening_evidence.write_evidence`.
    """
    return ScreeningEvidenceEntry(
        channel=ScreeningChannel.TRANSCRIPT_SEED,
        species=result.scope.species,
        reference_id=reference_id if reference_id is not None else result.scope.reference_id,
        guide_set_digest=guide_set_digest,
        status=result.status,
        counts=result.counts,
        submitted_guides=result.submitted_guides,
        processed_guides=result.processed_guides,
        detail=result.detail,
    )
