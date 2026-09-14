"""Tests for the #100 execution-outcome contract (A5-execution-outcome).

``BwaAnalyzer.analyze_sequences`` and ``scan_mirna_seed_matches`` return only a hit list, so a
hard aligner/backend failure and a clean zero-hit run were byte-identical: both return ``[]``.
``ExecutionOutcome`` plus the ``*_with_outcome`` sibling entry points are what let a caller tell
them apart, along with the submitted/processed guide membership and cap/truncation metadata the
evidence envelope contract requires. ``run_mirna_seed_analysis`` is the batch entry point that
must not let one species' failure discard another species' already-computed results.
"""

import importlib
import json
import subprocess
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.core.off_target import (
    BwaAnalyzer,
    ExecutionOutcome,
    MiRNASeedBackend,
    aggregate_mirna_results,
    parse_fasta_file,
    run_bwa_alignment_analysis,
    run_mirna_seed_analysis,
    scan_mirna_seed_matches_with_outcome,
)
from sirnaforge.core.screening_evidence import EvidenceProducer, read_evidence, write_evidence
from sirnaforge.models.evidence import EvidenceStatus, ScreeningEvidenceEntry
from sirnaforge.models.policy import ScreeningChannel

TEST_DATA_DIR = Path(__file__).parent / "data"
GUIDE = "ATGCGATGCGATGCGATGCGC"  # 21nt, matches test_offtarget_sam_parsing.py's fixture guide.


def _sam_line(*, qname: str, rname: str, pos: int, flag: int = 0, as_score: int = 42) -> str:
    """A single perfect-match SAM record for ``GUIDE`` against one target."""
    fields = [qname, str(flag), rname, str(pos), "60", "21M", "*", "0", "0", GUIDE, "*"]
    fields.append("NM:i:0")
    fields.append(f"AS:i:{as_score}")
    return "\t".join(fields)


class _FakeCompletedProcess:
    """Minimal stand-in for ``subprocess.CompletedProcess``, exposing only ``.stdout``."""

    def __init__(self, stdout: str) -> None:
        self.stdout = stdout


# ---------------------------------------------------------------------------
# ExecutionOutcome itself
# ---------------------------------------------------------------------------


def test_execution_outcome_completed_maps_to_complete_status_and_reports_sites() -> None:
    """A completed, uncapped execution reports COMPLETE and its retained count as an exact site total."""
    outcome = ExecutionOutcome(
        completed=True, submitted=5, processed=5, retained_hits=3, pre_cap_hits=3, cap=None, truncated=False
    )
    assert outcome.status is EvidenceStatus.COMPLETE
    counts = outcome.counts()
    assert counts.sites.value == 3
    assert counts.sites.is_lower_bound is False
    assert counts.sites.truncated is False


def test_execution_outcome_not_completed_maps_to_failed_status_and_observes_nothing() -> None:
    """A failed execution reports FAILED and observed nothing, never a misleading zero."""
    outcome = ExecutionOutcome(
        completed=False,
        submitted=5,
        processed=0,
        retained_hits=0,
        pre_cap_hits=0,
        cap=None,
        truncated=False,
        detail="the aligner exited non-zero",
    )
    assert outcome.status is EvidenceStatus.FAILED
    counts = outcome.counts()
    # A failed execution observed nothing -- not even a zero, which would misreport a clean screen.
    assert counts.sites.value is None
    assert counts.sites.is_observed is False


def test_execution_outcome_truncated_maps_to_censored_status_and_lower_bound_counts() -> None:
    """A capped execution reports CENSORED and a lower-bound, truncated site count."""
    outcome = ExecutionOutcome(
        completed=True, submitted=5, processed=5, retained_hits=2, pre_cap_hits=9, cap=2, truncated=True
    )
    assert outcome.status is EvidenceStatus.CENSORED
    counts = outcome.counts()
    assert counts.sites.value == 2
    assert counts.sites.is_lower_bound is True
    assert counts.sites.truncated is True
    assert counts.sites.cap == 2


# ---------------------------------------------------------------------------
# BwaAnalyzer.analyze_sequences_with_outcome -- red-without-fix #1 and #2
# ---------------------------------------------------------------------------


def test_analyze_sequences_with_outcome_reports_failed_execution_on_aligner_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A hard aligner failure must not be byte-identical to a clean zero-hit run (#100)."""
    monkeypatch.setattr("sirnaforge.core.off_target._get_executable_path", lambda _tool: "/usr/bin/bwa-mem2")

    def _raise(*_args: Any, **_kwargs: Any) -> Any:
        raise subprocess.CalledProcessError(1, ["bwa-mem2"], output="", stderr="fatal: could not open index")

    monkeypatch.setattr("sirnaforge.core.off_target.subprocess.run", _raise)

    analyzer = BwaAnalyzer(index_prefix="/tmp/does-not-matter", mode="transcriptome")
    hits, outcome = analyzer.analyze_sequences_with_outcome({"g1": GUIDE})

    assert hits == []
    assert outcome.completed is False
    assert outcome.status is EvidenceStatus.FAILED
    assert outcome.detail


def test_analyze_sequences_with_outcome_reports_truncation_against_the_return_cap(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Truncation is attributed to the method's own cap, never to ``_filter_and_rank``."""
    monkeypatch.setattr("sirnaforge.core.off_target._get_executable_path", lambda _tool: "/usr/bin/bwa-mem2")
    sam_output = (
        "\n".join(
            [
                _sam_line(qname="cand_1", rname="ENST00000000001", pos=100),
                _sam_line(qname="cand_1", rname="ENST00000000002", pos=200),
                _sam_line(qname="cand_1", rname="ENST00000000003", pos=300),
            ]
        )
        + "\n"
    )
    monkeypatch.setattr(
        "sirnaforge.core.off_target.subprocess.run",
        lambda *_a, **_k: _FakeCompletedProcess(sam_output),
    )

    analyzer = BwaAnalyzer(index_prefix="/tmp/does-not-matter", mode="transcriptome", max_hits=1)
    hits, outcome = analyzer.analyze_sequences_with_outcome({"cand_1": GUIDE})

    assert len(hits) == 1
    assert outcome.truncated is True
    assert outcome.cap == 1
    assert outcome.pre_cap_hits == 3
    assert outcome.retained_hits == 1
    assert outcome.counts().sites.is_lower_bound is True


def test_analyze_sequences_with_outcome_processed_counts_unmapped_records() -> None:
    """An unmapped SAM record still proves the query reached the aligner (#100 ``processed``)."""
    analyzer = BwaAnalyzer(index_prefix="/tmp/does-not-matter", mode="transcriptome")
    mapped = _sam_line(qname="cand_1", rname="ENST00000000001", pos=100)
    unmapped = "\t".join(["cand_2", "4", "*", "0", "0", "*", "*", "0", "0", GUIDE, "*"])
    hits, processed = analyzer._parse_sam_output_with_membership(
        f"{mapped}\n{unmapped}\n", {"cand_1": GUIDE, "cand_2": GUIDE}
    )
    assert len(hits) == 1
    assert processed == frozenset({"cand_1", "cand_2"})


def test_analyze_sequences_unchanged_still_returns_a_bare_list(monkeypatch: pytest.MonkeyPatch) -> None:
    """The legacy method keeps its exact return shape so no existing caller/test is touched."""
    monkeypatch.setattr("sirnaforge.core.off_target._get_executable_path", lambda _tool: "/usr/bin/bwa-mem2")

    def _raise(*_args: Any, **_kwargs: Any) -> Any:
        raise subprocess.CalledProcessError(1, ["bwa-mem2"], output="", stderr="boom")

    monkeypatch.setattr("sirnaforge.core.off_target.subprocess.run", _raise)

    analyzer = BwaAnalyzer(index_prefix="/tmp/does-not-matter", mode="transcriptome")
    assert analyzer.analyze_sequences({"g1": GUIDE}) == []


# ---------------------------------------------------------------------------
# scan_mirna_seed_matches_with_outcome
# ---------------------------------------------------------------------------


def test_scan_mirna_seed_matches_with_outcome_reports_backend_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
    """An unavailable in-process backend is a failed execution, not a hit list of zero."""
    real_import_module = importlib.import_module

    def _fake_import_module(name: str, package: str | None = None) -> Any:
        if name == "ahocorasick":
            raise ImportError("missing test dependency")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", _fake_import_module)

    hits, outcome = scan_mirna_seed_matches_with_outcome(
        {"g1": "ACGTACGTAC"}, {"m1": "ACGTACGTACGTACGTACGTAC"}, backend=MiRNASeedBackend.PYAHOCORASICK
    )

    assert hits == []
    assert outcome.completed is False
    assert outcome.status is EvidenceStatus.FAILED
    assert outcome.detail


def test_scan_mirna_seed_matches_with_outcome_truncates_and_reports_cap() -> None:
    """max_hits truncation on a real multi-hit scan is reported, not silently applied."""
    candidate_sequences = parse_fasta_file(TEST_DATA_DIR / "toy_candidates.fasta")
    mirna_sequences = parse_fasta_file(TEST_DATA_DIR / "toy_mirna_db.fasta")

    uncapped, uncapped_outcome = scan_mirna_seed_matches_with_outcome(
        candidate_sequences, mirna_sequences, backend=MiRNASeedBackend.EXHAUSTIVE_PYTHON
    )
    assert uncapped_outcome.pre_cap_hits > 1, "toy fixtures should produce more than one raw seed hit"
    assert uncapped_outcome.truncated is False

    capped, outcome = scan_mirna_seed_matches_with_outcome(
        candidate_sequences, mirna_sequences, backend=MiRNASeedBackend.EXHAUSTIVE_PYTHON, max_hits=1
    )
    assert len(capped) == 1
    assert outcome.truncated is True
    assert outcome.cap == 1
    assert outcome.pre_cap_hits == uncapped_outcome.pre_cap_hits
    assert outcome.counts().sites.is_lower_bound is True


# ---------------------------------------------------------------------------
# run_mirna_seed_analysis -- red-without-fix #3
# ---------------------------------------------------------------------------


def _patch_mirna_database(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "sirnaforge.data.mirna_manager.MiRNADatabaseManager.get_database",
        lambda *_args, **_kwargs: TEST_DATA_DIR / "toy_mirna_db.fasta",
    )


def _patch_ahocorasick_unavailable_after(monkeypatch: pytest.MonkeyPatch, *, succeed_calls: int) -> None:
    """Make ``import ahocorasick`` fail starting from the ``succeed_calls + 1``th attempt."""
    real_import_module = importlib.import_module
    call_count = {"n": 0}

    def _fake_import_module(name: str, package: str | None = None) -> Any:
        if name == "ahocorasick":
            call_count["n"] += 1
            if call_count["n"] > succeed_calls:
                raise ImportError("missing test dependency")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", _fake_import_module)


def test_run_mirna_seed_analysis_keeps_the_first_species_result_when_the_second_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A backend failure on species 2 must not discard species 1's already-computed results."""
    _patch_mirna_database(monkeypatch)
    _patch_ahocorasick_unavailable_after(monkeypatch, succeed_calls=1)

    output_dir = tmp_path / "out"
    evidence_dir = tmp_path / "evidence"

    result_path = run_mirna_seed_analysis(
        candidates_file=TEST_DATA_DIR / "toy_candidates.fasta",
        candidate_id="toy",
        mirna_db="toy_db",
        mirna_species=["human", "mouse"],
        output_dir=output_dir,
        backend=MiRNASeedBackend.PYAHOCORASICK,
        evidence_dir=evidence_dir,
    )

    assert result_path == output_dir
    # The first species' TSV is on disk, carrying real rows -- not discarded by the second
    # species' failure.
    analysis_tsv = output_dir / "toy_mirna_analysis.tsv"
    assert analysis_tsv.exists()
    assert analysis_tsv.stat().st_size > 0

    human_envelope = read_evidence(evidence_dir / "mirna_seed_human_evidence.json")
    mouse_envelope = read_evidence(evidence_dir / "mirna_seed_mouse_evidence.json")
    assert human_envelope is not None
    assert human_envelope.entry.status is EvidenceStatus.COMPLETE
    assert human_envelope.entry.channel is ScreeningChannel.MIRNA_SEED
    assert mouse_envelope is not None
    assert mouse_envelope.entry.status is EvidenceStatus.FAILED
    assert mouse_envelope.entry.detail


def test_run_mirna_seed_analysis_raises_only_when_every_species_hits_the_same_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The historical all-or-nothing signal survives when literally nothing could be screened."""
    _patch_mirna_database(monkeypatch)
    _patch_ahocorasick_unavailable_after(monkeypatch, succeed_calls=0)

    with pytest.raises(RuntimeError, match="backend 'pyahocorasick' is unavailable"):
        run_mirna_seed_analysis(
            candidates_file=TEST_DATA_DIR / "toy_candidates.fasta",
            candidate_id="toy",
            mirna_db="toy_db",
            mirna_species=["human", "mouse"],
            output_dir=tmp_path / "out",
            backend=MiRNASeedBackend.PYAHOCORASICK,
        )


# ---------------------------------------------------------------------------
# run_bwa_alignment_analysis -- additive evidence_dir wiring
# ---------------------------------------------------------------------------


def _write_candidates_fasta(path: Path) -> Path:
    path.write_text(f">cand_1\n{GUIDE}\n")
    return path


def test_run_bwa_alignment_analysis_writes_complete_evidence_on_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """evidence_dir opts a caller into a real transcriptome evidence envelope."""
    monkeypatch.setattr(
        "sirnaforge.core.off_target.BwaAnalyzer.analyze_sequences_with_outcome",
        lambda _self, _sequences: (
            [],
            ExecutionOutcome(
                completed=True, submitted=1, processed=1, retained_hits=0, pre_cap_hits=0, cap=None, truncated=False
            ),
        ),
    )
    candidates_file = _write_candidates_fasta(tmp_path / "candidate_0001.fasta")
    evidence_dir = tmp_path / "evidence"

    run_bwa_alignment_analysis(
        candidates_file=candidates_file,
        index_prefix=tmp_path / "index" / "human",
        species="human",
        output_dir=tmp_path / "out",
        evidence_dir=evidence_dir,
    )

    envelope = read_evidence(evidence_dir / "transcriptome_human_evidence.json")
    assert envelope is not None
    assert envelope.entry.status is EvidenceStatus.COMPLETE
    assert envelope.entry.channel is ScreeningChannel.TRANSCRIPTOME
    assert envelope.entry.species == "human"


def test_run_bwa_alignment_analysis_writes_failed_evidence_on_aligner_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed execution is published as a FAILED envelope with its detail intact."""
    monkeypatch.setattr(
        "sirnaforge.core.off_target.BwaAnalyzer.analyze_sequences_with_outcome",
        lambda _self, _sequences: (
            [],
            ExecutionOutcome(
                completed=False,
                submitted=1,
                processed=0,
                retained_hits=0,
                pre_cap_hits=0,
                cap=None,
                truncated=False,
                detail="no usable index",
            ),
        ),
    )
    candidates_file = _write_candidates_fasta(tmp_path / "candidate_0001.fasta")
    evidence_dir = tmp_path / "evidence"

    run_bwa_alignment_analysis(
        candidates_file=candidates_file,
        index_prefix=tmp_path / "index" / "human",
        species="human",
        output_dir=tmp_path / "out",
        evidence_dir=evidence_dir,
    )

    envelope = read_evidence(evidence_dir / "transcriptome_human_evidence.json")
    assert envelope is not None
    assert envelope.entry.status is EvidenceStatus.FAILED
    assert envelope.entry.detail == "no usable index"


def test_run_bwa_alignment_analysis_default_path_still_calls_analyze_sequences(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No ``evidence_dir`` means the exact historical call, so existing mocks stay valid."""
    calls: list[str] = []
    monkeypatch.setattr(
        "sirnaforge.core.off_target.BwaAnalyzer.analyze_sequences",
        lambda _self, _sequences: calls.append("analyze_sequences") or [],
    )
    monkeypatch.setattr(
        "sirnaforge.core.off_target.BwaAnalyzer.analyze_sequences_with_outcome",
        lambda _self, _sequences: calls.append("analyze_sequences_with_outcome") or ([], None),
    )
    candidates_file = _write_candidates_fasta(tmp_path / "candidate_0001.fasta")

    run_bwa_alignment_analysis(
        candidates_file=candidates_file,
        index_prefix=tmp_path / "index" / "human",
        species="human",
        output_dir=tmp_path / "out",
    )

    assert calls == ["analyze_sequences"]
    assert not (tmp_path / "out" / "evidence.json").exists()


# ---------------------------------------------------------------------------
# aggregate_mirna_results -- species_screened/unscreened_species bookkeeping
# ---------------------------------------------------------------------------


def test_aggregate_mirna_results_defaults_to_screened_when_no_evidence_exists(tmp_path: Path) -> None:
    """A pre-#100 results directory with no evidence envelope keeps its historical answer."""
    result_dir = aggregate_mirna_results(
        results_dir=tmp_path / "results",
        output_dir=tmp_path / "aggregated",
        mirna_db="toy_db",
        mirna_species="human,mouse",
    )
    summary = json.loads((result_dir / "combined_mirna_summary.json").read_text())
    assert summary["species_screened"] == ["human", "mouse"]
    assert summary["unscreened_species"] == []


def test_aggregate_mirna_results_reports_unscreened_species_from_evidence(tmp_path: Path) -> None:
    """A COMPLETE envelope for one species and none for the other is real positive evidence."""
    results_dir = tmp_path / "results"
    write_evidence(
        results_dir,
        producer=EvidenceProducer.MIRNA_SEED_ANALYSIS,
        entry=ScreeningEvidenceEntry(
            channel=ScreeningChannel.MIRNA_SEED,
            species="human",
            guide_set_digest="a" * 16,
            status=EvidenceStatus.COMPLETE,
        ),
    )

    result_dir = aggregate_mirna_results(
        results_dir=results_dir, output_dir=tmp_path / "aggregated", mirna_db="toy_db", mirna_species="human,mouse"
    )
    summary = json.loads((result_dir / "combined_mirna_summary.json").read_text())
    assert summary["species_screened"] == ["human"]
    assert summary["unscreened_species"] == ["mouse"]
