"""Pipeline wiring for #101's transcript-seed channel: the module, the gate and the aggregation.

The channel itself (``core/transcript_seed.py``) and its policy surface are pinned elsewhere. What is
pinned here is everything that decides *whether and how* it runs inside a real off-target run, and
each of these was a way to get the channel wrong without any of its own tests noticing:

1. **Off unless requested.** The scan is uncalibrated and enormous -- a 7mer occurs roughly once per
   16 kb -- so it must be reachable only through an explicit parameter, exactly as the miRNA channel's
   species list gates that channel. An unconditional invocation is the defect.
2. **The reference is the one already on disk.** The whole reason the scan lives beside
   ``OFFTARGET_ANALYSIS`` is that the expensive thing, the cDNA reference, has already been resolved,
   staged and indexed for the alignment channel. A task that resolved or downloaded its own reference
   would double the cost of the run and could silently screen a *different* reference than the
   alignment table it is published beside.
3. **Its own envelope, or a failed unit.** A scan that publishes a sites table but no #100 envelope
   leaves the (transcript_seed, species) pair unaccounted for: nothing then distinguishes "scanned and
   found nothing" from "never ran", and a transcript-seed ceiling would read the zero as a pass.
4. **Counters stay distinct.** ``<species>_transcript_seed_summary.json`` matches the ``*_summary.json``
   glob that feeds ``aggregate_offtarget_results``, so a seed table reaching the transcriptome/miRNA
   staging would have had its site count folded into the miRNA hit totals -- a third channel's number
   reported as more of the first two's.

The ``.nf`` assertions parse the workflow text, as ``tests/unit/test_evidence_plan_threading.py``
already does for the miRNA gate: none of these tests invokes a real ``nextflow`` binary.
"""

import ast
import inspect
import json
import textwrap
from pathlib import Path

import pytest

from sirnaforge.core.screening_evidence import (
    EvidenceProducer,
    collect_evidence,
    read_evidence,
    write_evidence,
)
from sirnaforge.models.evidence import EvidenceStatus, ScreeningEvidenceEntry
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.models.schemas import TranscriptSeedSiteSchema
from sirnaforge.pipeline import nextflow_cli

_WORKFLOWS_DIR = Path(__file__).resolve().parents[2] / "src/sirnaforge/pipeline/nextflow/workflows"
_SUBWORKFLOW = _WORKFLOWS_DIR / "subworkflows/local/sirna_offtarget_analysis.nf"
_SEED_MODULE = _WORKFLOWS_DIR / "modules/local/transcript_seed_analysis.nf"

#: A let-7a-like guide, so the geometry is independently checkable against #101's worked example:
#: ``guide[2..8] == "GAGGTAG"``, whose reverse complement ``"CTACCTC"`` is what a real site contains.
_GUIDE = "TGAGGTAGTAGGTTGTATAGT"

#: 108 bases whose 101..108 are ``CTACCTCA``: one 8mer site at anchor 107 and nothing else.
_TRANSCRIPT_WITH_A_SITE = "G" * 100 + "CTACCTCA"


def _stub_block(text: str) -> str:
    """The ``stub:`` half of one module's source, as a plain string."""
    return text.split("stub:", 1)[1]


def _executable_body(function: object) -> str:
    """One function's code with its docstring and comments removed.

    A prose ban has to be checked against code, not against prose: this function's own docstring says
    in as many words that it never downloads or re-indexes a reference, so a substring scan over the
    raw source would fail on the very sentence that promises the property.
    """
    tree = ast.parse(textwrap.dedent(inspect.getsource(function)))  # type: ignore[arg-type]
    definition = tree.body[0]
    assert isinstance(definition, ast.FunctionDef)
    first = definition.body[0]
    if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
        definition.body.pop(0)
    return ast.unparse(definition)


def _seed_inputs(directory: Path) -> tuple[Path, Path]:
    """One candidates FASTA and one cDNA reference carrying exactly one site for it."""
    candidates = directory / "candidates.fasta"
    candidates.write_text(f">cand_1\n{_GUIDE}\n")
    cdna = directory / "cdna.fasta"
    cdna.write_text(f">ENST00000000001.3 cdna gene:ENSG00000000001.5\n{_TRANSCRIPT_WITH_A_SITE}\n")
    return candidates, cdna


def _write_mirna_batch(directory: Path) -> tuple[Path, Path]:
    """One batch miRNA analysis TSV and summary, named the way MIRNA_SEED_ANALYSIS publishes them."""
    analysis = directory / "batch_mirna_analysis.tsv"
    analysis.write_text(
        "qname\tqseq\tspecies\tdatabase\tmirna_id\tcoord\tstrand\tcigar\tmapq\tas_score\tnm\t"
        "seed_mismatches\tofftarget_score\n"
    )
    summary = directory / "batch_mirna_summary.json"
    summary.write_text(json.dumps({"total_candidates": 0, "total_hits": 0}))
    return analysis, summary


# ---------------------------------------------------------------------------
# 1. Off unless requested
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_seed_channel_is_off_unless_requested():
    """TRANSCRIPT_SEED_ANALYSIS must sit behind params.transcript_seed_enabled and a resolved list.

    The channel ships off because it is uncalibrated and its output is enormous: a 7mer is expected
    roughly once per 16 kb, so an uncapped full-cDNA scan publishes orders of magnitude more rows than
    the alignment table. An unconditional invocation would make every existing run pay for it.

    The gate must also be a *value* test, not a truthiness test: Nextflow hands
    ``--transcript_seed_enabled false`` to Groovy as the non-empty String ``"false"``, which is truthy,
    so ``if (params.transcript_seed_enabled)`` would turn the documented way of asking for the channel
    off into the way of switching it on.
    """
    text = _SUBWORKFLOW.read_text()

    assert "params.transcript_seed_enabled" in text, "the channel must be keyed on its own parameter"
    gate = "if (transcript_seed_enabled && ch_transcript_seed_species_list) {"
    assert gate in text, "the invocation must be gated on the parameter AND a resolved species list"
    assert text.index(gate) < text.index("TRANSCRIPT_SEED_ANALYSIS("), "the gate must precede the call"

    # The literal allow-list, not truthiness. "false" must not be able to enable the channel.
    resolution = text.split("def transcript_seed_enabled", 1)[1].split("\n", 1)[0]
    assert "'true'" in resolution and "in [" in resolution, resolution
    assert "'false'" not in resolution, resolution

    # A request that resolves to zero species runs nothing, the same idiom the miRNA channel uses.
    species_resolution = text.split("def ch_transcript_seed_species_list", 1)[1].split("\n\n", 1)[0]
    assert ".findAll { it }" in species_resolution


@pytest.mark.unit
def test_the_existing_channel_floors_and_the_both_channels_off_abort_are_untouched():
    """#100's three ifEmpty([]) floors and the both-channels-off abort survive the new channel.

    The seed channel mixes into ``ch_all_evidence`` only, so the count of collected channels feeding
    AGGREGATE_RESULTS stays at three and the abort's key -- an empty alignment-input channel with the
    miRNA channel off -- keeps its meaning. Pinned here as well as in
    ``test_evidence_plan_threading.py`` because a fourth ``ch_all_*`` channel is the obvious way to
    wire this in and it would silently change what that file asserts.
    """
    text = _SUBWORKFLOW.read_text()

    assignments = [block for block in text.split("ch_all_")[1:] if "=" in block.split("\n", 1)[0]]
    assert len(assignments) == 3, "the new channel must not add a fourth collected channel"
    for block in assignments:
        assert ".ifEmpty([])" in block.split("AGGREGATE_RESULTS", 1)[0]

    assert "if (!ch_mirna_species_list) {" in text
    assert "Nothing to screen" in text


@pytest.mark.unit
def test_the_seed_tables_never_join_the_transcriptome_or_mirna_aggregation():
    """Only the seed *evidence* mixes into the aggregate; its tables reach their reader separately.

    ``ch_all_analysis``/``ch_all_summary`` feed ``aggregate_offtarget_results``, which sums
    transcriptome and known-miRNA hits. A ``<species>_transcript_seed_summary.json`` staged there is
    matched by that aggregate's own ``*_summary.json`` glob, so mixing the seed tables in would fold a
    seed site count into the miRNA counters (#101).
    """
    text = _SUBWORKFLOW.read_text()

    evidence_assignment = text.split("ch_all_evidence =", 1)[1].split("AGGREGATE_RESULTS", 1)[0]
    assert "ch_transcript_seed_evidence" in evidence_assignment

    for channel in ("ch_all_analysis =", "ch_all_summary ="):
        assignment = text.split(channel, 1)[1].split("\n\n", 1)[0]
        assert "transcript_seed" not in assignment, assignment


# ---------------------------------------------------------------------------
# 2. The already-materialised reference
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_module_takes_the_already_materialised_reference():
    """The scan is handed the staged cDNA FASTA and the staged guide set; it resolves neither.

    Both halves are load-bearing. ``cdna_fasta`` is a REQUIRED argument with no default, so there is
    no signature by which this entry point could be asked to find a reference for itself; and the
    module declares it as a ``path`` input taken from the same ``ch_reference_fastas`` BUILD_BWA_INDEX
    indexes, so the file scanned is the very file the alignment table was produced from. A task that
    downloaded or re-indexed its own copy would double the run's cost and could publish sites from a
    different reference release than the alignments it sits beside.
    """
    signature = inspect.signature(nextflow_cli.transcript_seed_analysis_cli)
    assert "cdna_fasta" in signature.parameters, "the reference must be an argument, not a lookup"
    assert signature.parameters["cdna_fasta"].default is inspect.Parameter.empty, "it must be required"
    assert signature.parameters["candidates_file"].default is inspect.Parameter.empty

    # The executable body only: ``ast.unparse`` drops the docstring and every comment, both of which
    # legitimately talk about the downloading and indexing this function refuses to do.
    code = _executable_body(nextflow_cli.transcript_seed_analysis_cli)
    for forbidden in ("download", "TranscriptomeManager", "build_bwa_index", "urlopen", "requests", "http"):
        assert forbidden not in code, f"the scan must not {forbidden}: the reference is already staged"

    module = _SEED_MODULE.read_text()
    input_block = module.split("input:", 1)[1].split("output:", 1)[0]
    assert "path(cdna_fasta)" in input_block, "declared as path, so Nextflow stages the resolved file"
    assert "path(candidates_fasta)" in input_block

    subworkflow = _SUBWORKFLOW.read_text()
    seed_input = subworkflow.split("ch_transcript_seed_input =", 1)[1].split("TRANSCRIPT_SEED_ANALYSIS(", 1)[0]
    assert "ch_reference_fastas" in seed_input, "the same FASTA channel BUILD_BWA_INDEX consumes"
    assert "candidates_fasta" in seed_input, "the same deduplicated guide set the aligner screens"
    assert "BUILD_BWA_INDEX" not in seed_input, "nothing is re-indexed for this channel"


@pytest.mark.unit
def test_the_module_declares_output_globs_because_the_producer_canonicalises_the_species():
    """``--transcriptome_fastas 'homo_sapiens:...'`` publishes ``human_*`` files, so the globs matter.

    ``transcript_seed_analysis_cli`` names every output after ``normalize_species_name(species)``, the
    same canonical form the counts and the evidence join key use. An output declaration interpolating
    the raw ``${species}`` would therefore fail to match its own task's real output and abort a run
    that had in fact screened correctly -- the same reason the miRNA module declares an evidence glob.
    """
    module = _SEED_MODULE.read_text()
    output_block = module.split("output:", 1)[1].split("when:", 1)[0]
    assert 'path "*_transcript_seed_sites.tsv", emit: sites' in output_block
    assert 'path "*_transcript_seed_summary.json", emit: summary' in output_block
    assert 'path "transcript_seed_*_evidence.json", emit: evidence' in output_block
    assert 'path "versions.yml", emit: versions' in output_block
    assert "${species}_transcript_seed" not in output_block


@pytest.mark.unit
def test_the_producer_names_its_outputs_after_the_canonical_species(tmp_path):
    """The behaviour the globs above exist for, exercised rather than asserted from the text."""
    candidates, cdna = _seed_inputs(tmp_path)
    out = tmp_path / "out"

    nextflow_cli.transcript_seed_analysis_cli(
        species="homo_sapiens",
        cdna_fasta=str(cdna),
        candidates_file=str(candidates),
        output_dir=str(out),
    )

    assert sorted(path.name for path in out.iterdir()) == [
        "human_transcript_seed_sites.tsv",
        "human_transcript_seed_summary.json",
        "transcript_seed_human_evidence.json",
    ]


@pytest.mark.unit
def test_the_default_site_cap_is_the_same_number_the_in_process_scan_uses():
    """Two producers of one artifact must not disagree about the cap that censors it.

    The cap is not cosmetic: when it truncates, the scan reports CENSORED and every count becomes a
    declared lower bound, so a ceiling read against it cannot claim to have been respected. Two
    different defaults would make the same guide set censored by one producer and complete by the
    other. Not imported from ``workflow.py`` because this module runs inside a Nextflow task that has
    no reason to import the design stack; pinned equal here instead.
    """
    from sirnaforge.workflow import _TRANSCRIPT_SEED_SITE_CAP

    assert nextflow_cli.DEFAULT_TRANSCRIPT_SEED_SITE_CAP == _TRANSCRIPT_SEED_SITE_CAP


# ---------------------------------------------------------------------------
# 3. Its own envelope, or a failed unit
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_seed_unit_publishes_its_own_envelope(tmp_path):
    """The scan writes a #100 envelope beside its table, under its own producer and channel.

    Without it the (transcript_seed, species) pair never enters ``completed_pairs``, so nothing
    distinguishes a scan that ran and found nothing from a scan that never ran -- and the three
    transcript-seed ceilings would read the table's zero rows as a clean screen. The producer is
    ``transcript_seed_analysis``, not ``mirna_seed_analysis``: an envelope whose origin is unreadable
    from the envelope is how the two seed questions get conflated.
    """
    candidates, cdna = _seed_inputs(tmp_path)
    out = tmp_path / "out"

    result = nextflow_cli.transcript_seed_analysis_cli(
        species="human",
        cdna_fasta=str(cdna),
        candidates_file=str(candidates),
        output_dir=str(out),
    )

    envelope_path = out / "transcript_seed_human_evidence.json"
    assert envelope_path.exists(), "a scan that publishes no envelope leaves its unit unaccounted for"
    envelope = read_evidence(envelope_path)
    assert envelope is not None
    assert envelope.producer is EvidenceProducer.TRANSCRIPT_SEED_ANALYSIS
    assert envelope.entry.channel is ScreeningChannel.TRANSCRIPT_SEED
    assert envelope.entry.species == "human"
    assert envelope.entry.status is EvidenceStatus.COMPLETE
    assert envelope.entry.counts.sites.value == 1
    assert envelope.entry.submitted_guides == 1
    assert envelope.entry.processed_guides == 1
    # The join key is the guide set, so the envelope can only be attributed to these guides.
    assert envelope.entry.guide_set_digest == nextflow_cli.guide_set_digest(str(candidates))

    assert result["status"] == "complete"
    assert result["sites"] == 1
    assert Path(result["evidence_file"]) == envelope_path


@pytest.mark.unit
def test_the_published_table_has_the_declared_site_shape_and_the_worked_geometry(tmp_path):
    """One shape for two producers, and the orientation is the reverse-complement one.

    The header is asserted against ``TranscriptSeedSiteSchema`` (``strict=True``) rather than against a
    literal list, so this module and ``workflow.py`` cannot drift apart. The row is #101's worked
    example: ``guide[2..8] == "GAGGTAG"``, whose reverse complement ``"CTACCTC"`` occupies 101..107, so
    the anchor is 107 and the 'A' at 108 makes it an 8mer. Getting the geometry backwards -- searching
    for ``"GAGGTAG"`` itself -- would find passenger-orientation windows at the same expected rate, so
    the coordinates, not the count, are what pins it.
    """
    candidates, cdna = _seed_inputs(tmp_path)
    out = tmp_path / "out"

    nextflow_cli.transcript_seed_analysis_cli(
        species="human",
        cdna_fasta=str(cdna),
        candidates_file=str(candidates),
        output_dir=str(out),
    )

    lines = (out / "human_transcript_seed_sites.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    assert tuple(header) == nextflow_cli.TRANSCRIPT_SEED_SITE_COLUMNS
    assert header == list(TranscriptSeedSiteSchema.to_schema().columns)

    row = dict(zip(header, lines[1].split("\t")))
    assert row["site_start"] == "101"
    assert row["site_end"] == "107"
    assert row["anchor_position"] == "107"
    assert row["site_class"] == "8mer"
    assert row["site_strand"] == "transcript_sense"
    assert row["queried_strand"] == "guide"
    assert row["transcript_id"] == "ENST00000000001"
    assert row["transcript_version"] == "3"
    assert row["gene_id"] == "ENSG00000000001"
    assert row["coordinate_system"] == "transcript_cdna_1based"


@pytest.mark.unit
@pytest.mark.parametrize(
    ("scope", "expected_in_detail"),
    [("utr3", "no UTR annotation"), ("not_a_region", "unknown region scope")],
)
def test_a_scope_this_run_cannot_answer_is_a_failed_unit_not_an_empty_table(tmp_path, scope, expected_in_detail):
    """A refused scope publishes a 0-byte table and a FAILED envelope, never a clean zero.

    Nothing in this repository carries UTR intervals in cDNA coordinates, and an unrecognised scope
    string is not silently widened to the whole cDNA. Both refusals must survive the pipeline
    boundary: the module's outputs are non-optional, so *something* has to be written, and a
    header-only table is a scan that ran and found nothing. 0 bytes is the spelling of "did not run" --
    ``workflow.py``'s parser skips a 0-byte table for exactly this reason -- and the FAILED envelope
    keeps the pair out of ``completed_pairs`` so the ceilings report UNKNOWN.
    """
    candidates, cdna = _seed_inputs(tmp_path)
    out = tmp_path / "out"

    result = nextflow_cli.transcript_seed_analysis_cli(
        species="human",
        cdna_fasta=str(cdna),
        candidates_file=str(candidates),
        output_dir=str(out),
        region_scope=scope,
    )

    assert result["status"] == "failed"
    assert expected_in_detail in (result["detail"] or "")
    assert (out / "human_transcript_seed_sites.tsv").stat().st_size == 0
    envelope = read_evidence(out / "transcript_seed_human_evidence.json")
    assert envelope is not None
    assert envelope.entry.status is EvidenceStatus.FAILED
    # Wholly unobserved, not zero: an unobserved count is what makes a ceiling report UNKNOWN.
    assert envelope.entry.counts.sites.value is None


@pytest.mark.unit
def test_an_unreadable_reference_is_a_failed_unit_not_a_screen_that_found_nothing(tmp_path):
    """The analogue of ``offtarget_analysis_cli``'s missing-index branch, for a missing cDNA FASTA.

    A scan over a reference that is not there searches nothing and finds nothing, which is
    byte-identical to a clean screen unless the unit is published as failed.
    """
    candidates, _ = _seed_inputs(tmp_path)
    out = tmp_path / "out"

    result = nextflow_cli.transcript_seed_analysis_cli(
        species="human",
        cdna_fasta=str(tmp_path / "absent.fasta"),
        candidates_file=str(candidates),
        output_dir=str(out),
    )

    assert result["status"] == "failed"
    assert "no readable cDNA reference" in (result["detail"] or "")
    assert (out / "human_transcript_seed_sites.tsv").stat().st_size == 0
    envelope = read_evidence(out / "transcript_seed_human_evidence.json")
    assert envelope is not None and envelope.entry.status is EvidenceStatus.FAILED


@pytest.mark.unit
def test_the_stub_mirrors_every_declared_output_and_never_claims_completion():
    """A ``-stub-run`` must satisfy the non-optional outputs without asserting a screen happened.

    Nothing scanned, so the envelope is ``failed``/``stub``; the table is 0 bytes; and the stub imports
    no ``sirnaforge`` (a stub run has no container and no importable package).
    """
    stub = _stub_block(_SEED_MODULE.read_text())

    assert "_transcript_seed_sites.tsv" in stub
    assert "_transcript_seed_summary.json" in stub
    assert "transcript_seed_stub_evidence.json" in stub
    assert "versions.yml" in stub
    assert '"status": "failed"' in stub
    assert '"producer": "stub"' in stub
    assert '"status": "complete"' not in stub
    assert "python" not in stub.split("versions.yml")[0]


@pytest.mark.unit
def test_reconciliation_reports_a_missing_seed_unit_as_failed(tmp_path, monkeypatch):
    """A requested seed species that published no envelope reconciles FAILED, not silence.

    Two halves. The scan's own envelope has to reach the aggregate's search root, and the plan has to
    name every seed unit the run asked for -- including the fallback plan, which is what a bare
    Nextflow run reconciles against. Without the plan side, the observed envelope arrives as
    ``unplanned``: reported, but with no entry for the species that produced none, so a scan that died
    reads as a channel nobody asked about. Without the envelope side, the species that DID scan
    reconciles failed as well and the distinction the test makes is vacuous.
    """
    monkeypatch.chdir(tmp_path)
    candidates, cdna = _seed_inputs(tmp_path)
    nextflow_cli.transcript_seed_analysis_cli(
        species="human",
        cdna_fasta=str(cdna),
        candidates_file=str(candidates),
        output_dir=str(tmp_path / "seed"),
    )

    result = nextflow_cli.aggregate_results_cli(
        transcriptome_species="",
        output_dir=str(tmp_path / "aggregated"),
        analysis_files=[],
        summary_files=[],
        transcript_seed_species="human,mouse",
    )

    statuses = {(entry["channel"], entry["species"]): entry["status"] for entry in result["evidence"]["entries"]}
    assert statuses[("transcript_seed", "human")] == "complete"
    assert statuses[("transcript_seed", "mouse")] == "failed"


@pytest.mark.unit
def test_mirna_and_transcript_seed_counters_stay_distinct(tmp_path, monkeypatch):
    """A seed table handed to the aggregate is refused by name, not folded into the miRNA counters.

    ``<species>_transcript_seed_summary.json`` matches the ``*_summary.json`` glob that feeds
    ``aggregate_offtarget_results``, and ``species in f.name`` matches it again in the per-species
    staging, so a seed summary reaching this task by any route would have had its site count added to
    the transcriptome and miRNA hit totals. The channel is reported through its own section instead,
    per species and never summed across them.
    """
    monkeypatch.chdir(tmp_path)
    candidates, cdna = _seed_inputs(tmp_path)
    seed = nextflow_cli.transcript_seed_analysis_cli(
        species="human",
        cdna_fasta=str(cdna),
        candidates_file=str(candidates),
        output_dir=str(tmp_path),
    )
    mirna_analysis, mirna_summary = _write_mirna_batch(tmp_path)
    write_evidence(
        tmp_path,
        producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.MIRNA_SEED,
            species="human",
            guide_set_digest=nextflow_cli.guide_set_digest(str(candidates)),
            status=EvidenceStatus.COMPLETE,
        ),
    )

    result = nextflow_cli.aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(tmp_path / "aggregated"),
        mirna_db="toy_db",
        mirna_species="human",
        analysis_files=[str(mirna_analysis), seed["sites_file"]],
        summary_files=[str(mirna_summary), seed["summary_file"]],
        transcript_seed_species="human",
    )

    # The seed files reached neither aggregate's staging directory.
    staged = [path.name for path in (tmp_path / "temp_results").rglob("*") if path.is_file()]
    assert not [name for name in staged if "transcript_seed" in name], staged
    assert result["mirna"]["analysis_files_processed"] == 1
    assert result["analysis_files_processed"] == 1
    assert result["summary_files_processed"] == 1

    # The miRNA aggregate's own numbers are untouched by a channel that found a site.
    combined_mirna = json.loads((tmp_path / "aggregated" / "combined_mirna_summary.json").read_text())
    assert combined_mirna["total_mirna_hits"] == 0
    assert combined_mirna["species_screened"] == ["human"]

    # The seed channel reports separately, per species, with all three units carried through.
    seed_section = result["transcript_seed"]
    assert seed_section["requested_species"] == ["human"]
    human = seed_section["by_species"]["human"]
    assert human["status"] == "complete"
    assert human["counts"]["sites"]["value"] == 1
    assert human["counts"]["distinct_transcripts"]["value"] == 1
    assert human["counts"]["distinct_genes"]["value"] == 1
    assert "mirna" not in json.dumps(seed_section)


@pytest.mark.unit
def test_the_channel_evidence_stager_is_parameterised_and_keeps_the_channels_apart(tmp_path):
    """``_staged_channel_evidence`` selects by channel prefix, so no channel's evidence leaks.

    The miRNA aggregate reads whatever envelopes it finds in the directory it is handed. A stager that
    globbed ``*_evidence.json`` would copy transcript-seed envelopes in there and let a third channel's
    species decide which species the miRNA screen reported.
    """
    for channel in (ScreeningChannel.MIRNA_SEED, ScreeningChannel.TRANSCRIPT_SEED):
        write_evidence(
            tmp_path,
            producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
            entry=ScreeningEvidenceEntry(
                channel=channel,
                species="human",
                guide_set_digest="a" * 16,
                status=EvidenceStatus.COMPLETE,
            ),
        )
    destination = tmp_path / "mirna"
    destination.mkdir()

    mirna = nextflow_cli._staged_mirna_evidence(tmp_path, destination)
    seed = nextflow_cli._staged_channel_evidence(tmp_path, destination, channel=ScreeningChannel.TRANSCRIPT_SEED)

    assert [path.name for path in mirna] == ["mirna_seed_human_evidence.json"]
    assert [path.name for path in seed] == ["transcript_seed_human_evidence.json"]
    assert len(collect_evidence(tmp_path)) == 2
