"""Nextflow-side plumbing for #100's expected-plan contract (A7).

Three defects, pinned independently:

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

All fixtures here are text/static -- none of these tests invokes a real ``nextflow`` binary.
"""

import asyncio
import subprocess
from pathlib import Path

import pytest

from sirnaforge.pipeline.nextflow.config import NextflowConfig
from sirnaforge.pipeline.nextflow.runner import NextflowExecutionError, NextflowRunner

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
