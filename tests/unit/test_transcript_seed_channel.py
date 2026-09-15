"""The transcript-seed liability channel: orientation, identity and honest counts (#101).

Every assertion about a site is on a **string and a coordinate**, never on a number of matches. That
is not a style preference: the guide's seed window searched *forward* -- which is what
``off_target._prepare_seed_queries`` correctly does against a miRNA database -- finds
guide-identical, passenger-orientation windows in a cDNA at the same expected rate as the real
reverse-complement site. A reversed orientation therefore produces a null that is indistinguishable
by count and wrong by identity, so the composition-matched decoy transcript below exists to catch it.

The rest of the file pins the four things #101 says a count must not be able to hide: nested classes
at one anchor are one site, a gene count over unresolved transcripts is a declared lower bound, a
truncating cap is CENSORED rather than a smaller number, and a UTR-only request with no UTR
annotation is a refusal rather than a clean screen.
"""

import ast
import inspect
from pathlib import Path

import pytest

from sirnaforge.core import transcript_seed as transcript_seed_module
from sirnaforge.core.screening_evidence import (
    EvidenceProducer,
    build_plan,
    collect_evidence,
    completed_pairs,
    guide_set_digest,
    reconcile,
    write_evidence,
)
from sirnaforge.core.seed_geometry import (
    M8_POSITION,
    SEED_START_POSITION,
    complement,
    paired_transcript_position,
    reverse_complement,
    seed_window,
)
from sirnaforge.core.transcript_seed import (
    ALL_SEED_CLASSES,
    COORDINATE_SYSTEM,
    SEED_CLASS_RANK,
    GuideStrand,
    SeedClass,
    SiteRegion,
    SiteStrand,
    scan_transcript_seed_sites,
    transcript_seed_evidence_entry,
)
from sirnaforge.data.transcript_index import TranscriptGeneIndex
from sirnaforge.models.evidence import EvidenceStatus
from sirnaforge.models.policy import ScreeningChannel

# A let-7a-like guide, so the seed windows below are checkable against published let-7 sequences
# rather than only against this code.
GUIDE = "TGAGGTAGTAGGTTGTATAGT"
M8_SITE = "CTACCTC"  # reverse_complement(GUIDE[2..8]) -- what a guide actually pairs with
SIXMER_SITE = "TACCTC"  # reverse_complement(GUIDE[2..7])
FORWARD_M8_WINDOW = "GAGGTAG"  # GUIDE[2..8] itself: a passenger-orientation window, not a site

# 25 bases, chosen so a 7-base core sits at 26..32 and the A1 base at 33, and so the composition
# difference between M8_SITE and FORWARD_M8_WINDOW can be absorbed by the prefix (see DECOY).
PREFIX = "AGGGG" + "T" * 20
ANCHOR = 32

GUIDES = {"guide_1": GUIDE}


def _transcript(core: str, a1_base: str) -> str:
    """A 33-base transcript: 25 bases of prefix, a 7-base core at 26..32, one base at 33."""
    return PREFIX + core + a1_base


# The real site and a decoy of identical length AND identical base composition that carries the
# forward window instead. Composition is matched deliberately: a decoy that merely differed would
# leave "zero sites" attributable to base content rather than to orientation.
M8_TRANSCRIPT = _transcript(M8_SITE, "G")
EIGHTMER_TRANSCRIPT = _transcript(M8_SITE, "A")
SIXMER_TRANSCRIPT = _transcript("A" + SIXMER_SITE, "G")
A1_TRANSCRIPT = _transcript("A" + SIXMER_SITE, "A")
DECOY_TRANSCRIPT = "CCCCT" + "T" * 20 + FORWARD_M8_WINDOW + "G"


def _write_fasta(tmp_path: Path, records: list[tuple[str, str]], name: str = "cdna.fa") -> Path:
    """Write a FASTA whose headers are Ensembl-shaped, so TranscriptGeneIndex can read them too."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / name
    path.write_text("".join(f">{header}\n{sequence}\n" for header, sequence in records))
    return path


def _ensembl_header(transcript_id: str, *, gene_id: str | None, symbol: str = "AAA") -> str:
    """One Ensembl cDNA header, optionally carrying no ``gene:`` field at all."""
    parts = [transcript_id, "cdna", "chromosome:GRCh38:1:1:33:1"]
    if gene_id is not None:
        parts += [f"gene:{gene_id}", "gene_biotype:protein_coding", "transcript_biotype:protein_coding"]
        parts.append(f"gene_symbol:{symbol}")
    return " ".join(parts)


def _single_transcript_fasta(tmp_path: Path, sequence: str) -> Path:
    """A one-record cDNA FASTA carrying ``sequence`` under a versioned transcript id."""
    return _write_fasta(tmp_path, [(_ensembl_header("ENST00000000001.1", gene_id="ENSG00000000001.1"), sequence)])


def _scan(fasta: Path, **kwargs) -> transcript_seed_module.TranscriptSeedScanResult:
    """Scan ``GUIDES`` against ``fasta`` for human, with every class requested unless overridden."""
    kwargs.setdefault("classes", ALL_SEED_CLASSES)
    return scan_transcript_seed_sites(GUIDES, fasta, species="human", **kwargs)


def test_the_search_strings_are_the_reverse_complements_the_docstrings_claim() -> None:
    """The fixtures' own arithmetic, so a wrong constant fails here rather than everywhere."""
    assert seed_window(GUIDE, SEED_START_POSITION, M8_POSITION) == FORWARD_M8_WINDOW
    assert reverse_complement(FORWARD_M8_WINDOW) == M8_SITE
    assert reverse_complement(seed_window(GUIDE, SEED_START_POSITION, 7)) == SIXMER_SITE
    # The decoy is a permutation of the real transcript, so orientation is the only difference.
    assert sorted(DECOY_TRANSCRIPT) == sorted(M8_TRANSCRIPT)
    assert FORWARD_M8_WINDOW in DECOY_TRANSCRIPT
    assert M8_SITE not in DECOY_TRANSCRIPT
    assert SIXMER_SITE not in DECOY_TRANSCRIPT


def test_a_forward_seed_window_is_not_a_target_site(tmp_path: Path) -> None:
    """The orientation test #101 exists for: the decoy carries the guide's own window, not a site.

    Both halves are asserted together on purpose. A scanner pointed at the forward window reports one
    site for the decoy and none for the real transcript, i.e. it inverts this test rather than
    changing a total, which is exactly why a count-based null could not have caught it.
    """
    real = _scan(_single_transcript_fasta(tmp_path / "real", M8_TRANSCRIPT))
    decoy = _scan(_single_transcript_fasta(tmp_path / "decoy", DECOY_TRANSCRIPT))

    assert decoy.sites == ()
    assert decoy.status is EvidenceStatus.COMPLETE
    assert decoy.counts.sites.value == 0
    assert len(real.sites) == 1
    assert real.sites[0].site_class is SeedClass.SEVENMER_M8


def test_a_7mer_m8_site_is_found_at_the_computed_anchor(tmp_path: Path) -> None:
    """Site string and coordinates, not a count: this is the assertion a reversed map fails on."""
    result = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT))

    (site,) = result.sites
    assert site.anchor_position == ANCHOR
    assert site.site_end == ANCHOR
    assert site.site_start == paired_transcript_position(M8_POSITION, ANCHOR)
    assert M8_TRANSCRIPT[site.site_start - 1 : site.site_end] == M8_SITE
    assert site.site_class is SeedClass.SEVENMER_M8
    assert site.coordinate_system == COORDINATE_SYSTEM


def test_the_pairing_map_holds_at_every_position_of_a_found_site(tmp_path: Path) -> None:
    """``transcript[p(i)] == complement(guide[i])`` for i in 2..8, asserted on the site itself.

    A reversed map (``anchor - 2 + i``) leaves ``p(2) == anchor`` intact, so asserting only the anchor
    would pass. Walking the whole window is what makes the direction of the map load-bearing.
    """
    (site,) = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT)).sites

    for guide_position in range(SEED_START_POSITION, M8_POSITION + 1):
        transcript_position = paired_transcript_position(guide_position, site.anchor_position)
        assert M8_TRANSCRIPT[transcript_position - 1] == complement(GUIDE[guide_position - 1])


@pytest.mark.parametrize(
    ("sequence", "expected_class", "expected_window"),
    [
        (SIXMER_TRANSCRIPT, SeedClass.SIXMER, SIXMER_SITE),
        (A1_TRANSCRIPT, SeedClass.SEVENMER_A1, SIXMER_SITE),
        (M8_TRANSCRIPT, SeedClass.SEVENMER_M8, M8_SITE),
        (EIGHTMER_TRANSCRIPT, SeedClass.EIGHTMER, M8_SITE),
    ],
)
def test_each_seed_class_is_pinned_by_one_base(
    tmp_path: Path, sequence: str, expected_class: SeedClass, expected_window: str
) -> None:
    """Four fixtures differing from a neighbour in exactly one base: the m8 base, or the A1 base.

    All four share ``anchor == 32``, which is the property that makes them one identity when they
    co-occur, and the published window follows the class rather than being fixed.
    """
    (site,) = _scan(_single_transcript_fasta(tmp_path, sequence)).sites

    assert site.site_class is expected_class
    assert site.anchor_position == ANCHOR
    assert sequence[site.site_start - 1 : site.site_end] == expected_window


def test_the_a1_base_is_a_literal_adenosine_not_a_complement(tmp_path: Path) -> None:
    """A1 is an unpaired target adenosine, so it is tested as 'A' and never through complement().

    ``complement(GUIDE[1]) == 'A'`` for this guide, so the two rules agree here by coincidence. A
    guide whose first base is not ``T`` would make them disagree, and the class table's wording is
    what decides. Pinned as an explicit statement about the fixture, so a future guide change cannot
    silently rely on the coincidence.
    """
    a1_position = paired_transcript_position(1, ANCHOR)
    assert a1_position == ANCHOR + 1
    assert EIGHTMER_TRANSCRIPT[a1_position - 1] == "A"
    assert M8_TRANSCRIPT[a1_position - 1] == "G"

    (eightmer,) = _scan(_single_transcript_fasta(tmp_path / "eight", EIGHTMER_TRANSCRIPT)).sites
    (m8,) = _scan(_single_transcript_fasta(tmp_path / "m8", M8_TRANSCRIPT)).sites
    assert (eightmer.site_class, m8.site_class) == (SeedClass.EIGHTMER, SeedClass.SEVENMER_M8)


def test_nested_classes_at_one_anchor_are_one_site(tmp_path: Path) -> None:
    """``CTACCTCA`` realises a 6mer, a 7mer-A1, a 7mer-m8 and an 8mer, and publishes one site.

    The 6mer window matches at 27..32 and the m8 window at 26..32 -- two candidates with two
    different starts and one shared anchor. Anchoring the identity on ``site_start`` publishes them
    as two sites, which is why the identity is the anchor and this test counts.
    """
    result = _scan(_single_transcript_fasta(tmp_path, EIGHTMER_TRANSCRIPT))

    assert len(result.sites) == 1
    assert result.counts.sites.value == 1
    (site,) = result.sites
    assert site.site_class is SeedClass.EIGHTMER
    assert SEED_CLASS_RANK[site.site_class] == max(SEED_CLASS_RANK.values())
    assert site.identity == ("guide_1", GuideStrand.GUIDE.value, "human", "ENST00000000001", "1", ANCHOR)


def test_two_records_of_one_transcript_version_are_one_site(tmp_path: Path) -> None:
    """A FASTA carrying the same transcript twice must not double the site count."""
    header = _ensembl_header("ENST00000000001.1", gene_id="ENSG00000000001.1")
    fasta = _write_fasta(tmp_path, [(header, M8_TRANSCRIPT), (header, M8_TRANSCRIPT)])

    result = _scan(fasta)

    assert len(result.sites) == 1
    assert result.counts.distinct_transcripts.value == 1


def test_a_requested_class_set_decides_what_is_published(tmp_path: Path) -> None:
    """Classes are a declared scope, not a default: an 8mer-only request skips a 7mer-m8 anchor.

    Requesting only 6mers over an 8mer anchor still yields a site -- it *is* a 6mer site -- published
    at the maximal *requested* class, so narrowing the scope can never invent a stronger class than
    was asked for.
    """
    m8_fasta = _single_transcript_fasta(tmp_path / "m8", M8_TRANSCRIPT)
    eightmer_fasta = _single_transcript_fasta(tmp_path / "eight", EIGHTMER_TRANSCRIPT)

    assert _scan(m8_fasta, classes=frozenset({SeedClass.EIGHTMER})).sites == ()

    (sixmer_only,) = _scan(eightmer_fasta, classes=frozenset({SeedClass.SIXMER})).sites
    assert sixmer_only.site_class is SeedClass.SIXMER
    assert sixmer_only.anchor_position == ANCHOR
    assert EIGHTMER_TRANSCRIPT[sixmer_only.site_start - 1 : sixmer_only.site_end] == SIXMER_SITE


def test_every_row_says_which_strand_was_queried(tmp_path: Path) -> None:
    """0.7.1 submits guides only, so passenger liability is unscreened and every row says so."""
    (site,) = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT)).sites

    assert site.queried_strand is GuideStrand.GUIDE
    assert site.site_strand is SiteStrand.TRANSCRIPT_SENSE
    assert site.region is SiteRegion.FULL_CDNA


def test_a_gene_cap_over_unresolved_transcripts_is_a_lower_bound(tmp_path: Path) -> None:
    """A site whose transcript has no gene is counted, not dropped, and it bounds the gene count.

    Three transcripts each carry one site: two share a gene, the third's header carries no ``gene:``
    field at all. So the honest report is one resolved gene, *at least*, plus one site whose gene is
    unknown -- which is how a gene-level ceiling stops silently discarding what it could not resolve.
    """
    fasta = _write_fasta(
        tmp_path,
        [
            (_ensembl_header("ENST00000000001.1", gene_id="ENSG00000000001.1"), M8_TRANSCRIPT),
            (_ensembl_header("ENST00000000002.1", gene_id="ENSG00000000001.1"), EIGHTMER_TRANSCRIPT),
            (_ensembl_header("ENST00000000003.1", gene_id=None), SIXMER_TRANSCRIPT),
        ],
    )
    index = TranscriptGeneIndex()
    index.build("human", fasta)

    result = _scan(fasta, gene_index=index)

    assert result.status is EvidenceStatus.COMPLETE
    assert result.counts.sites.value == 3
    assert result.counts.distinct_transcripts.value == 3
    assert result.counts.distinct_genes.value == 1
    assert result.counts.distinct_genes.is_lower_bound is True
    assert result.counts.unresolved_gene_sites.value == 1


def test_a_fully_resolved_gene_count_is_not_a_lower_bound(tmp_path: Path) -> None:
    """The sibling case, so the lower-bound flag means something rather than being always true."""
    fasta = _write_fasta(
        tmp_path,
        [(_ensembl_header("ENST00000000001.1", gene_id="ENSG00000000001.1"), M8_TRANSCRIPT)],
    )
    index = TranscriptGeneIndex()
    index.build("human", fasta)

    result = _scan(fasta, gene_index=index)

    assert result.counts.distinct_genes.value == 1
    assert result.counts.distinct_genes.is_lower_bound is False
    assert result.counts.unresolved_gene_sites.value == 0


def test_no_gene_index_leaves_the_gene_unit_unobserved(tmp_path: Path) -> None:
    """No index means gene resolution was never attempted, and 0 would pass any gene ceiling."""
    result = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT))

    assert result.counts.distinct_genes.value is None
    assert "distinct_genes" not in result.counts.observed_units
    assert result.counts.unresolved_gene_sites.value == 1


def test_a_truncating_cap_is_censored_not_a_smaller_number(tmp_path: Path) -> None:
    """Three sites, a cap of two: the strongest survive, and every count says it is a lower bound."""
    multi = "AGGGG" + "T" * 10 + M8_SITE + "G" + "T" * 10 + M8_SITE + "G" + "T" * 10 + M8_SITE + "A"
    result = _scan(_single_transcript_fasta(tmp_path, multi), max_sites_per_guide=2)

    assert result.status is EvidenceStatus.CENSORED
    assert [site.site_class for site in result.sites] == [SeedClass.EIGHTMER, SeedClass.SEVENMER_M8]
    assert [site.anchor_position for site in result.sites] == [58, 22]
    assert result.counts.sites.value == 2
    assert result.counts.sites.is_lower_bound is True
    assert result.counts.sites.truncated is True
    assert result.detail is not None
    assert "per-guide site cap 2" in result.detail
    assert "1 of 3 sites were discarded" in result.detail
    assert result.scope.max_sites_per_guide == 2


def test_a_cap_that_does_not_bite_stays_complete(tmp_path: Path) -> None:
    """A cap in force is not censorship; only a cap that actually discarded members is."""
    result = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT), max_sites_per_guide=2)

    assert result.status is EvidenceStatus.COMPLETE
    assert result.counts.sites.is_lower_bound is False
    assert result.scope.max_sites_per_guide == 2


def test_a_missing_utr_annotation_is_not_a_clean_utr_screen(tmp_path: Path) -> None:
    """A 3' UTR request with no UTR annotation is a refusal with a reason, never zero sites.

    Returning COMPLETE with zero sites would let a transcript-seed ceiling pass on a screen that
    never happened. FAILED is #100's spelling of "we asked and could not answer": the pair then never
    enters ``completed_pairs``, so the gate reports UNKNOWN.
    """
    result = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT), region_scope=SiteRegion.UTR3)

    assert result.status is EvidenceStatus.FAILED
    assert result.sites == ()
    assert result.counts.observed_units == ()
    assert result.detail is not None
    assert "3' UTR scope requested" in result.detail
    assert "human" in result.detail
    assert result.scope.region is SiteRegion.UTR3


def test_a_guide_too_short_to_carry_the_seed_window_is_not_screened(tmp_path: Path) -> None:
    """A short guide is unprocessed and named, rather than screened for 6mers and reported clean."""
    result = scan_transcript_seed_sites(
        {"short": "TGAGGT"},
        _single_transcript_fasta(tmp_path, M8_TRANSCRIPT),
        species="human",
    )

    assert result.status is EvidenceStatus.FAILED
    assert (result.submitted_guides, result.processed_guides) == (1, 0)
    assert result.counts.observed_units == ()


def test_the_scan_streams_the_reference_rather_than_loading_it(tmp_path: Path) -> None:
    """The FASTA reader is a generator, because a human cDNA set must not be held whole.

    ``FastaUtils.read_fasta`` returns a list and would be the obvious reuse; a cDNA reference is
    exactly the file where that costs gigabytes, so the reader is pinned as lazy and as yielding one
    record at a time.
    """
    fasta = _write_fasta(
        tmp_path,
        [
            (_ensembl_header("ENST00000000001.1", gene_id="ENSG00000000001.1"), M8_TRANSCRIPT),
            (_ensembl_header("ENST00000000002.1", gene_id="ENSG00000000002.1"), DECOY_TRANSCRIPT),
        ],
    )

    assert inspect.isgeneratorfunction(transcript_seed_module._iter_fasta_records)
    records = transcript_seed_module._iter_fasta_records(fasta)
    first = next(iter(records))
    assert first[0].startswith("ENST00000000001.1")
    assert first[1] == M8_TRANSCRIPT


def test_the_channel_does_not_reach_into_the_mirna_scanner_or_the_classifier() -> None:
    """#101 forbids repointing the miRNA seed scanner, so this module cannot even see it.

    An AST-level assertion rather than a monkeypatch: the value of the separation is that the two
    modules' contracts stay independent, and an import is the cheapest way for that to stop being
    true.
    """
    tree = ast.parse(Path(transcript_seed_module.__file__).read_text())
    imported = {node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom) and node.module is not None}

    assert "sirnaforge.core.off_target" not in imported
    assert "sirnaforge.core.hit_classification" not in imported
    assert "sirnaforge.core.seed_geometry" in imported


def test_the_producer_names_the_transcript_seed_scan() -> None:
    """A separate producer, so an envelope's origin is readable from the envelope."""
    assert EvidenceProducer.TRANSCRIPT_SEED_ANALYSIS.value == "transcript_seed_analysis"
    assert EvidenceProducer.TRANSCRIPT_SEED_ANALYSIS is not EvidenceProducer.MIRNA_SEED_ANALYSIS


def test_a_plan_that_does_not_request_the_channel_is_byte_identical() -> None:
    """The opt-in guarantee: defaulting the new argument changes no existing plan's bytes."""
    before = build_plan(guide_set_digest="d" * 16, transcriptome=[("human", "ref")], mirna_species=["human"])
    after = build_plan(
        guide_set_digest="d" * 16,
        transcriptome=[("human", "ref")],
        mirna_species=["human"],
        transcript_seed=(),
    )

    assert before.model_dump_json() == after.model_dump_json()
    assert all(entry.channel is not ScreeningChannel.TRANSCRIPT_SEED for entry in after.entries)


def test_a_transcript_seed_plan_entry_is_appended_after_the_existing_ones() -> None:
    """Appended, never interleaved: entry order for the two older channels is unchanged."""
    plan = build_plan(
        guide_set_digest="d" * 16,
        transcriptome=[("human", "ref")],
        mirna_species=["human"],
        transcript_seed=[("human", "cdna-115")],
    )

    assert [entry.channel for entry in plan.entries] == [
        ScreeningChannel.TRANSCRIPTOME,
        ScreeningChannel.MIRNA_SEED,
        ScreeningChannel.TRANSCRIPT_SEED,
    ]
    assert plan.entries[-1].reference_id == "cdna-115"


def test_a_transcript_seed_envelope_round_trips_through_reconcile_and_completed_pairs(tmp_path: Path) -> None:
    """#100's contract needs no change: the (channel, species, digest) join already generalises."""
    guides_fasta = tmp_path / "guides.fa"
    guides_fasta.write_text(f">guide_1\n{GUIDE}\n")
    digest = guide_set_digest(guides_fasta)
    result = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT))

    entry = transcript_seed_evidence_entry(result, guide_set_digest=digest, reference_id="cdna-115")
    write_evidence(tmp_path / "evidence", producer=EvidenceProducer.TRANSCRIPT_SEED_ANALYSIS, entry=entry)

    plan = build_plan(
        guide_set_digest=digest,
        transcriptome=[],
        mirna_species=[],
        transcript_seed=[("human", "cdna-115")],
    )
    reconciliation = reconcile(plan, collect_evidence(tmp_path / "evidence"))

    assert reconciliation.unplanned == ()
    assert completed_pairs(reconciliation) == frozenset({(ScreeningChannel.TRANSCRIPT_SEED.value, "human")})
    (reconciled,) = reconciliation.evidence.entries
    assert reconciled.counts.sites.value == 1
    assert reconciled.processed_guides == 1


def test_a_refused_scope_keeps_its_pair_out_of_completed_pairs(tmp_path: Path) -> None:
    """The consequence of FAILED, end to end: the gate sees an incomplete pair, not a clean zero."""
    guides_fasta = tmp_path / "guides.fa"
    guides_fasta.write_text(f">guide_1\n{GUIDE}\n")
    digest = guide_set_digest(guides_fasta)
    result = _scan(_single_transcript_fasta(tmp_path, M8_TRANSCRIPT), region_scope=SiteRegion.UTR3)

    entry = transcript_seed_evidence_entry(result, guide_set_digest=digest, reference_id="cdna-115")
    write_evidence(tmp_path / "evidence", producer=EvidenceProducer.TRANSCRIPT_SEED_ANALYSIS, entry=entry)

    plan = build_plan(
        guide_set_digest=digest,
        transcriptome=[],
        mirna_species=[],
        transcript_seed=[("human", "cdna-115")],
    )
    reconciliation = reconcile(plan, collect_evidence(tmp_path / "evidence"))

    assert completed_pairs(reconciliation) == frozenset()


def test_a_censored_scan_carries_the_reason_its_evidence_entry_requires(tmp_path: Path) -> None:
    """``ScreeningEvidenceEntry`` rejects a reasonless CENSORED entry, so the scan derives one."""
    multi = "AGGGG" + "T" * 10 + M8_SITE + "G" + "T" * 10 + M8_SITE + "A"
    result = _scan(_single_transcript_fasta(tmp_path, multi), max_sites_per_guide=1)

    entry = transcript_seed_evidence_entry(result, guide_set_digest="d" * 16)

    assert entry.status is EvidenceStatus.CENSORED
    assert entry.detail is not None and entry.detail.strip()
    assert entry.channel is ScreeningChannel.TRANSCRIPT_SEED
    assert entry.counts.sites.value == 1
