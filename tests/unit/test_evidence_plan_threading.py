"""Nextflow-side plumbing for #100's expected-plan contract (A7).

Five defects, pinned independently:

1. ``AGGREGATE_RESULTS`` used to learn which species to reconcile from ``ch_reference_indices`` --
   the channel of indices that actually *built* -- so a species whose ``BUILD_BWA_INDEX`` crashed
   silently vanished from the expected list instead of reconciling as failed. The species handed
   to aggregation must come from the subworkflow's own ``expected_species`` ``take:`` argument.
2. A bare ``nextflow run -stub-run`` (no container, no aligner, no importable ``sirnaforge``) must
   still satisfy every non-optional output declared on a real script block, including the new
   #100 evidence outputs, or the run fails on an unmatched output glob.
3. ``NextflowExecutionError`` used to carry only ``str(CalledProcessError)`` in its message,
   leaving ``.stdout``/``.stderr`` at their empty-string defaults even though the subprocess call
   that raised it had captured real output.
4. With BOTH screening channels off, every channel the subworkflow collected was empty. An empty
   ``collect()`` emits nothing at all (unlike ``toList()``), so ``AGGREGATE_RESULTS`` never ran and
   a ``nextflow run`` exited 0 having published no evidence whatsoever -- reached through the
   subworkflow's own documented disable idiom, a value for ``--mirna_species`` that resolves to
   zero species with no transcriptome reference configured.
5. ``aggregate_results_cli`` staged only the analysis and summary TSV/JSON files into the miRNA
   results directory it then hands to ``aggregate_mirna_results``. That function searches the
   directory it is given for per-unit evidence envelopes, so it found none on every real run, its
   no-envelopes fallback fired, and ``species_screened`` reported the full requested list however
   the run had actually gone.

Every fixture here is text/static or a direct function call -- none of these tests invokes a real
``nextflow`` binary.
"""

import asyncio
import json
import subprocess
from pathlib import Path

import pytest

from sirnaforge.core.screening_evidence import EvidenceProducer, write_evidence
from sirnaforge.models.evidence import EvidenceStatus, ScreeningEvidenceEntry
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.pipeline.nextflow.config import NextflowConfig
from sirnaforge.pipeline.nextflow.runner import NextflowExecutionError, NextflowRunner
from sirnaforge.pipeline.nextflow_cli import aggregate_results_cli

_WORKFLOWS_DIR = Path(__file__).resolve().parents[2] / "src/sirnaforge/pipeline/nextflow/workflows"
_SUBWORKFLOW = _WORKFLOWS_DIR / "subworkflows/local/sirna_offtarget_analysis.nf"
_OFFTARGET_MODULE = _WORKFLOWS_DIR / "modules/local/offtarget_analysis.nf"
_MIRNA_MODULE = _WORKFLOWS_DIR / "modules/local/mirna_seed_analysis.nf"
_AGGREGATE_MODULE = _WORKFLOWS_DIR / "modules/local/aggregate_results.nf"
_MAIN_NF = _WORKFLOWS_DIR / "main.nf"
_NEXTFLOW_CONFIG = _WORKFLOWS_DIR / "nextflow.config"


def _stub_block(text: str) -> str:
    """The ``stub:`` half of one process/module's source, as a plain string."""
    return text.split("stub:", 1)[1]


# ---------------------------------------------------------------------------
# 1. expected_species, not ch_reference_indices
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_aggregate_results_is_handed_expected_species_not_the_built_index_channel():
    """The subworkflow's own take: argument, not a derivation from which indices happened to build."""
    text = _SUBWORKFLOW.read_text()

    take_block = text.split("take:", 1)[1].split("main:", 1)[0]
    assert "expected_species" in take_block, "expected_species must be a take: argument"
    assert "evidence_plan" in take_block, "evidence_plan must be a take: argument"

    call_block = text.split("AGGREGATE_RESULTS(", 1)[1].split(")", 1)[0]
    assert "expected_species" in call_block, "AGGREGATE_RESULTS must be called with expected_species"

    # The historical defect: species derived from ch_reference_indices via a map/unique/toList
    # chain named ch_screened_species. That derivation must be gone entirely, not merely unused.
    assert "ch_screened_species" not in text


@pytest.mark.unit
def test_main_nf_passes_transcriptome_species_into_the_subworkflow():
    """main.nf must forward params.transcriptome_species as the subworkflow's expected_species.

    Before #100, params.transcriptome_species had exactly one reader: the log banner. A species
    the caller asked for never reached the subworkflow at all.
    """
    text = _MAIN_NF.read_text()
    call_block = text.split("SIRNA_OFFTARGET_ANALYSIS(", 1)[1].split(")", 1)[0]
    assert "params.transcriptome_species" in call_block
    assert "evidence_plan" in call_block.lower() or "ch_evidence_plan" in call_block


@pytest.mark.unit
def test_evidence_plan_param_is_declared_and_guarded_by_name():
    """--evidence_plan is a real declared param, and a differently-named caller is refused, not silently ignored."""
    config_text = _NEXTFLOW_CONFIG.read_text()
    assert "evidence_plan" in config_text.split("params {", 1)[1].split("\n}", 1)[0]

    main_text = _MAIN_NF.read_text()
    renamed_block = main_text.split("renamed_params = [", 1)[1].split("]", 1)[0]
    assert "evidence_plan" in renamed_block


# ---------------------------------------------------------------------------
# 2. Every stub mirrors every non-optional output of its real script block
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_offtarget_analysis_stub_writes_its_declared_evidence_output():
    """OFFTARGET_ANALYSIS declares a non-optional transcriptome_${species}_evidence.json output."""
    text = _OFFTARGET_MODULE.read_text()
    assert 'path "transcriptome_${species}_evidence.json", emit: evidence' in text

    stub = _stub_block(text)
    assert "transcriptome_${species}_evidence.json" in stub
    assert '"status": "failed"' in stub
    assert '"producer": "stub"' in stub
    assert '"status": "complete"' not in stub


@pytest.mark.unit
def test_mirna_seed_analysis_stub_writes_a_file_matching_its_declared_evidence_glob():
    """MIRNA_SEED_ANALYSIS declares a non-optional mirna_seed_*_evidence.json output."""
    text = _MIRNA_MODULE.read_text()
    assert 'path "mirna_seed_*_evidence.json", emit: evidence' in text

    stub = _stub_block(text)
    assert "mirna_seed_stub_evidence.json" in stub
    assert '"status": "failed"' in stub
    assert '"producer": "stub"' in stub
    assert '"status": "complete"' not in stub


@pytest.mark.unit
def test_aggregate_results_stub_writes_evidence_json():
    """AGGREGATE_RESULTS declares a non-optional evidence.json output.

    The pre-existing pinned assertion (test_offtarget_aggregation_completeness.py's stub test)
    only ever checked the combined_offtargets.tsv header; evidence.json is additive and must not
    disturb that check.
    """
    text = _AGGREGATE_MODULE.read_text()
    assert 'path "evidence.json", emit: evidence' in text

    stub = _stub_block(text)
    assert "evidence.json" in stub
    assert '"status": "complete"' not in stub
    # Untouched: the header line the other test file pins must survive verbatim.
    assert "combined_offtargets.tsv" in stub


@pytest.mark.unit
def test_no_stub_block_imports_sirnaforge():
    """A `-stub-run` has no container, no aligner and no importable sirnaforge package.

    Guards the same property tests/unit/test_offtarget_aggregation_completeness.py pins for
    AGGREGATE_RESULTS specifically, extended here to all three modules this slice touched.
    """
    for module in (_OFFTARGET_MODULE, _MIRNA_MODULE, _AGGREGATE_MODULE):
        stub = _stub_block(module.read_text())
        before_versions = stub.split("versions.yml")[0]
        assert "python" not in before_versions, module


@pytest.mark.unit
def test_aggregate_results_process_declares_evidence_plan_and_evidence_files_inputs():
    """The process must actually accept the plan file and the staged per-unit envelopes.

    Declaring `evidence.json` as an output with no way to see the plan or the envelopes it
    reconciles against would make every real run's reconciliation fall back to "nothing observed".
    """
    text = _AGGREGATE_MODULE.read_text()
    input_block = text.split("input:", 1)[1].split("output:", 1)[0]
    assert "evidence_plan" in input_block
    assert "evidence_files" in input_block
    # Declared as `path`, not `val`: only `path` inputs get staged into this task's own working
    # directory, which is exactly what collect_evidence() globs at aggregation time.
    assert "path evidence_plan" in input_block
    assert "path evidence_files" in input_block


# ---------------------------------------------------------------------------
# 4. Both channels off must not be a silent success
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_every_collected_channel_has_an_empty_list_floor():
    """AGGREGATE_RESULTS must run even when no upstream unit produced a single file (#100).

    ``collect()`` emits nothing for an empty source, which withheld the aggregation step
    altogether -- so no reconciliation was published and the run exited 0 having screened nothing.
    An empty *list* runs the aggregation on nothing, which is what publishes "these units were
    expected and not one of them produced anything". Verified against a standalone stub harness:
    an empty ``collect()`` leaves the consuming process at ``completed=0``, an
    ``ifEmpty([])`` floor runs it.
    """
    text = _SUBWORKFLOW.read_text()

    collected = [
        block
        for block in text.split("ch_all_")[1:]
        # Only the three assignments, not later reads of the same names.
        if "=" in block.split("\n", 1)[0]
    ]
    assert len(collected) == 3, "expected exactly three collected channels feeding AGGREGATE_RESULTS"
    for block in collected:
        assignment = block.split("AGGREGATE_RESULTS", 1)[0]
        assert ".collect()" in assignment
        assert ".ifEmpty([])" in assignment, assignment


@pytest.mark.unit
def test_the_subworkflow_aborts_when_neither_screening_channel_can_run():
    """With the miRNA channel off and no alignment unit, the only honest outcome is a non-zero exit.

    The floor above cannot make this configuration safe on its own: nothing was requested, so the
    expected plan can be empty too, and a reconciliation over an empty plan has no shortfall to
    report -- it is indistinguishable from a clean screen. Reproduced before the fix with
    ``nextflow run main.nf -stub-run --mirna_species ','``: exit 0, ``completed=0``, and not one
    published file under ``aggregated/``. After it, exit 1 with the message below.
    """
    text = _SUBWORKFLOW.read_text()

    guard = text.split("if (!ch_mirna_species_list) {", 1)
    assert len(guard) == 2, "the both-channels-off guard must be keyed on the resolved species list"
    guard_block = guard[1].split("OFFTARGET_ANALYSIS(", 1)[0]
    assert "ifEmpty" in guard_block, "the guard must fire on an empty alignment-input channel"
    assert "error(" in guard_block, "the guard must abort the run, not warn"
    assert "Nothing to screen" in guard_block
    # The guard must sit upstream of the process it protects, or the run reaches aggregation first.
    assert text.index("if (!ch_mirna_species_list) {") < text.index("OFFTARGET_ANALYSIS(\n")


# ---------------------------------------------------------------------------
# 5. The miRNA evidence envelopes reach the directory the miRNA aggregate searches
# ---------------------------------------------------------------------------


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


@pytest.mark.unit
def test_aggregate_results_cli_stages_mirna_evidence_so_screened_species_reflect_it(tmp_path, monkeypatch):
    """A COMPLETE envelope for one requested species must not report the other one screened (#100).

    ``aggregate_mirna_results`` reads envelopes from the search root it is given, and before #100 a
    root with none in it kept the permissive answer: every requested species reported screened.
    Staging the analysis files without their envelopes made that fallback fire on every real run.
    Both halves of the repair are exercised here -- the envelopes are staged, and the search root is
    named -- so the reported species come from what the run actually did.
    """
    monkeypatch.chdir(tmp_path)
    analysis, summary = _write_mirna_batch(tmp_path)
    write_evidence(
        tmp_path,
        producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.MIRNA_SEED,
            species="human",
            guide_set_digest="a" * 16,
            status=EvidenceStatus.COMPLETE,
        ),
    )

    aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(tmp_path / "aggregated"),
        mirna_db="toy_db",
        mirna_species="human,mouse",
        analysis_files=[str(analysis)],
        summary_files=[str(summary)],
    )

    mirna_summary = json.loads((tmp_path / "aggregated" / "combined_mirna_summary.json").read_text())
    assert mirna_summary["species_screened"] == ["human"]
    assert mirna_summary["unscreened_species"] == ["mouse"]


@pytest.mark.unit
def test_a_failed_mirna_species_reaches_the_aggregate_instead_of_being_masked(tmp_path, monkeypatch):
    """The FAILED envelope for an unresolvable database survives staging and is reported (#100).

    This is the end-to-end statement of defects 5/7 together with the backend's unresolvable-database
    repair: ``run_mirna_seed_analysis`` publishes a FAILED envelope naming the species whose database
    could not be resolved, ``aggregate_results_cli`` stages it beside the tables, and the aggregate
    reads it from the search root it was handed. All three have to hold. If any one of them does not,
    the mouse envelope is invisible here and the pre-#100 fallback would have reported *both* species
    screened -- turning the one species that demonstrably failed into a completed screen.

    The COMPLETE human envelope in the same directory is what makes the assertion discriminating
    rather than vacuous: the aggregate is distinguishing the two outcomes, not refusing everything.
    """
    monkeypatch.chdir(tmp_path)
    analysis, summary = _write_mirna_batch(tmp_path)
    for species, status, detail in (
        ("human", EvidenceStatus.COMPLETE, None),
        ("mouse", EvidenceStatus.FAILED, "no miRNA database entry for mouse in toy_db"),
    ):
        write_evidence(
            tmp_path,
            producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
            entry=ScreeningEvidenceEntry(
                channel=ScreeningChannel.MIRNA_SEED,
                species=species,
                guide_set_digest="a" * 16,
                status=status,
                detail=detail,
            ),
        )

    aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(tmp_path / "aggregated"),
        mirna_db="toy_db",
        mirna_species="human,mouse",
        analysis_files=[str(analysis)],
        summary_files=[str(summary)],
    )

    staged = sorted(path.name for path in (tmp_path / "temp_results" / "mirna").glob("*_evidence.json"))
    assert staged == ["mirna_seed_human_evidence.json", "mirna_seed_mouse_evidence.json"]
    mirna_summary = json.loads((tmp_path / "aggregated" / "combined_mirna_summary.json").read_text())
    assert mirna_summary["species_screened"] == ["human"]
    assert mirna_summary["unscreened_species"] == ["mouse"]


@pytest.mark.unit
def test_mirna_evidence_staging_copies_the_envelope_next_to_the_analysis_files(tmp_path, monkeypatch):
    """The envelope lands in the same directory the miRNA aggregate is pointed at, once."""
    monkeypatch.chdir(tmp_path)
    analysis, summary = _write_mirna_batch(tmp_path)
    write_evidence(
        tmp_path,
        producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.MIRNA_SEED,
            species="human",
            guide_set_digest="a" * 16,
            status=EvidenceStatus.COMPLETE,
        ),
    )

    aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(tmp_path / "aggregated"),
        mirna_db="toy_db",
        mirna_species="human",
        analysis_files=[str(analysis)],
        summary_files=[str(summary)],
    )

    staged = sorted((tmp_path / "temp_results" / "mirna").glob("*_evidence.json"))
    assert [path.name for path in staged] == ["mirna_seed_human_evidence.json"]


@pytest.mark.unit
def test_the_mirna_aggregate_is_told_where_the_envelopes_actually_are(tmp_path, monkeypatch):
    """The evidence search root is named explicitly: the task's own directory, not the staging one.

    ``results_dir`` holds copies of the analysis files this call staged; the envelopes are written
    by the upstream tasks and staged into this task's working directory. Naming the root keeps the
    two facts from being conflated (#100).
    """
    monkeypatch.chdir(tmp_path)
    analysis, summary = _write_mirna_batch(tmp_path)
    recorded: dict[str, object] = {}

    def _spy(**kwargs):
        recorded.update(kwargs)
        return Path(kwargs["output_dir"])

    monkeypatch.setattr("sirnaforge.pipeline.nextflow_cli.aggregate_mirna_results", _spy)

    aggregate_results_cli(
        transcriptome_species="human",
        output_dir=str(tmp_path / "aggregated"),
        mirna_db="toy_db",
        mirna_species="human",
        analysis_files=[str(analysis)],
        summary_files=[str(summary)],
    )

    assert Path(str(recorded["evidence_root"])) == Path()
    assert Path(str(recorded["results_dir"])) != Path()


# ---------------------------------------------------------------------------
# config.py needs no change: additional_params already passes any key through generically
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_nextflow_config_passes_evidence_plan_through_generically(tmp_path):
    """--evidence_plan needs no special-cased handling in NextflowConfig.get_nextflow_args.

    additional_params already forwards any key it is given as a bare ``--key value`` pair; this
    pins that ``evidence_plan`` gets that treatment for free, rather than requiring a bespoke
    branch the way ``transcriptome_species`` needed one.
    """
    config = NextflowConfig(profile="local")
    plan_path = tmp_path / "screening_plan.json"
    plan_path.write_text("{}")

    args = config.get_nextflow_args(
        input_file=tmp_path / "input.fasta",
        output_dir=tmp_path / "out",
        screen_species=["human"],
        additional_params={"evidence_plan": str(plan_path)},
    )

    assert "--evidence_plan" in args
    assert args[args.index("--evidence_plan") + 1] == str(plan_path)


# ---------------------------------------------------------------------------
# 3. NextflowExecutionError carries the captured stdout/stderr
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_failed_nextflow_run_carries_its_captured_stdout_and_stderr(tmp_path):
    """A CalledProcessError's own output/stderr must survive onto NextflowExecutionError.

    Before #100 both raise sites in run_offtarget_analysis only interpolated ``str(e)`` into the
    message, leaving ``.stdout``/``.stderr`` at their empty-string defaults -- a caller inspecting
    the exception for diagnostics saw nothing, even though the subprocess itself had captured
    plenty.
    """
    runner = NextflowRunner(NextflowConfig(profile="local", work_dir=tmp_path / "work"))
    input_file = tmp_path / "candidates.fasta"
    input_file.write_text(">cand_1\nACGT\n")

    async def _boom(cmd, env=None):
        raise subprocess.CalledProcessError(1, cmd, output=b"captured stdout", stderr=b"captured stderr")

    runner._run_subprocess = _boom  # type: ignore[method-assign]

    with pytest.raises(NextflowExecutionError) as excinfo:
        asyncio.run(
            runner.run_offtarget_analysis(
                input_file=input_file,
                output_dir=tmp_path / "out",
                screen_species=["human"],
                show_progress=False,
            )
        )

    assert excinfo.value.stdout == "captured stdout"
    assert excinfo.value.stderr == "captured stderr"


@pytest.mark.unit
def test_a_failed_nextflow_run_with_progress_display_also_carries_stdout_and_stderr(tmp_path):
    """The progress-bar raise site (show_progress=True, the default) gets the identical fix."""
    runner = NextflowRunner(NextflowConfig(profile="local", work_dir=tmp_path / "work"))
    input_file = tmp_path / "candidates.fasta"
    input_file.write_text(">cand_1\nACGT\n")

    async def _boom(cmd, env=None):
        raise subprocess.CalledProcessError(1, cmd, output=b"progress stdout", stderr=b"progress stderr")

    runner._run_subprocess = _boom  # type: ignore[method-assign]

    with pytest.raises(NextflowExecutionError) as excinfo:
        asyncio.run(
            runner.run_offtarget_analysis(
                input_file=input_file,
                output_dir=tmp_path / "out",
                screen_species=["human"],
                show_progress=True,
            )
        )

    assert excinfo.value.stdout == "progress stdout"
    assert excinfo.value.stderr == "progress stderr"


@pytest.mark.unit
def test_a_failed_nextflow_run_with_no_captured_output_stays_empty_string(tmp_path):
    """No captured output must decode to "", not crash on a None -> str conversion."""
    runner = NextflowRunner(NextflowConfig(profile="local", work_dir=tmp_path / "work"))
    input_file = tmp_path / "candidates.fasta"
    input_file.write_text(">cand_1\nACGT\n")

    async def _boom(cmd, env=None):
        raise subprocess.CalledProcessError(1, cmd)

    runner._run_subprocess = _boom  # type: ignore[method-assign]

    with pytest.raises(NextflowExecutionError) as excinfo:
        asyncio.run(
            runner.run_offtarget_analysis(
                input_file=input_file,
                output_dir=tmp_path / "out",
                screen_species=["human"],
                show_progress=False,
            )
        )

    assert excinfo.value.stdout == ""
    assert excinfo.value.stderr == ""
