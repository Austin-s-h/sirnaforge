"""Container integration tests for ZFN design-mode workflow.

These tests validate the complete ZFN pipeline end-to-end via the CLI,
using the known CCR5 ZFN benchmark data (PROGNOS / Fine et al. 2014).

Designed to run inside Docker via ``make docker-build-test``.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pytest
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# CCR5 benchmark constants (from Fine et al. 2014 / PROGNOS)
# ---------------------------------------------------------------------------
# Primary on-target half-sites reported in PROGNOS Table 1
CCR5_LEFT_HALF_SITE = "GTCATCCTCATC"
CCR5_RIGHT_HALF_SITE = "AAACTGCAAAAG"

# Known validated off-target genes (from S10 visible rows)
CCR5_KNOWN_OFFTARGET_GENES = frozenset({"CCR5", "CSNK1G3", "TMOD1"})


def _get_persistent_output_dir(tmp_path: Path, test_name: str) -> Path:
    """Get output directory with persistence support for failure inspection."""
    if Path("/workspace").exists() and os.access("/workspace", os.W_OK):
        output_dir = Path("/workspace") / f"workflow_test_debug_{test_name}"
    else:
        output_dir = tmp_path / test_name

    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _build_ccr5_synthetic_genome(tmp_path: Path) -> Path:
    """Build a synthetic FASTA containing known CCR5 ZFN off-target sites.

    Embeds the on-target site (CCR5) and the novel validated off-target
    near CSNK1G3 into a simple single-chromosome genome.  This is enough
    to exercise the full sliding-window search without needing hg19.
    """
    s10_path = Path(__file__).resolve().parents[1] / "unit" / "data" / "zfn" / "ccr5_s10_visible_rows.csv"
    with s10_path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    ccr5_row = next(r for r in rows if r["Closest gene"] == "CCR5")
    csnk1g3_row = next(r for r in rows if r["Closest gene"] == "CSNK1G3")
    tmod1_row = next(r for r in rows if r["Closest gene"] == "TMOD1")

    def _site_seq(plus: str, minus: str, spacer_len: int) -> str:
        """Construct genomic site: plus + spacer + RC(minus)."""
        return plus.upper() + "A" * spacer_len + str(Seq(minus.upper()).reverse_complement())

    sites = [
        _site_seq(ccr5_row["(+) half-site"], ccr5_row["(−) half-site"], 5),  # on-target
        _site_seq(csnk1g3_row["(+) half-site"], csnk1g3_row["(−) half-site"], 5),  # novel validated
        _site_seq(tmod1_row["(+) half-site"], tmod1_row["(−) half-site"], 5),  # known off-target
    ]

    # Separate sites by random-ish spacers to avoid overlap
    genome_seq = "N" * 50
    for site in sites:
        genome_seq += site + "N" * 100
    genome_seq += "N" * 50

    fasta = tmp_path / "ccr5_synthetic_genome.fa"
    fasta.write_text(f">chr_synthetic\n{genome_seq}\n")
    return fasta


def _build_chr3_scaled_benchmark_genome(tmp_path: Path) -> Path:
    """Build a larger chr3-only synthetic FASTA to stress backend parity in workflow mode."""
    s10_path = Path(__file__).resolve().parents[1] / "unit" / "data" / "zfn" / "ccr5_s10_visible_rows.csv"
    with s10_path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    ccr5_row = next(r for r in rows if r["Closest gene"] == "CCR5")
    csnk1g3_row = next(r for r in rows if r["Closest gene"] == "CSNK1G3")
    tmod1_row = next(r for r in rows if r["Closest gene"] == "TMOD1")

    def _site_seq(plus: str, minus: str, spacer_len: int) -> str:
        return plus.upper() + "A" * spacer_len + str(Seq(minus.upper()).reverse_complement())

    chr3_seq = (
        ("ACGT" * 10_000)
        + _site_seq(ccr5_row["(+) half-site"], ccr5_row["(−) half-site"], 5)
        + ("TGCA" * 8_000)
        + _site_seq(csnk1g3_row["(+) half-site"], csnk1g3_row["(−) half-site"], 5)
        + ("GATC" * 8_000)
        + _site_seq(tmod1_row["(+) half-site"], tmod1_row["(−) half-site"], 5)
        + ("CATG" * 8_000)
        + _site_seq(ccr5_row["(+) half-site"], ccr5_row["(−) half-site"], 6)
        + ("TTAA" * 10_000)
    )

    fasta = tmp_path / "chr3_scaled_benchmark.fa"
    fasta.write_text(f">chr3\n{chr3_seq}\n", encoding="utf-8")
    return fasta


def _normalize_offtarget_csv(offtarget_csv: Path) -> pd.DataFrame:
    """Normalize workflow output CSV for deterministic backend-to-backend dataframe comparisons."""
    with offtarget_csv.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    df = pd.DataFrame(rows)
    keep_columns = [
        "chrom",
        "start_1based",
        "end_1based",
        "strand",
        "orientation",
        "spacer_len",
        "sequence",
        "left_mismatches",
        "right_mismatches",
        "total_mismatches",
        "score",
    ]
    normalized = df[keep_columns].copy()
    numeric_columns = [
        "start_1based",
        "end_1based",
        "spacer_len",
        "left_mismatches",
        "right_mismatches",
        "total_mismatches",
        "score",
    ]
    for column in numeric_columns:
        normalized[column] = pd.to_numeric(normalized[column])
    normalized["score"] = normalized["score"].round(8)
    return normalized.sort_values(
        by=["chrom", "start_1based", "end_1based", "orientation", "total_mismatches", "sequence"]
    ).reset_index(drop=True)


def _run_zfn_backend_workflow(
    *,
    tmp_path: Path,
    backend_label: str,
    genome_fasta: Path,
    search_space_index: Path | None = None,
) -> tuple[pd.DataFrame, Path]:
    """Run one full ZFN workflow backend and return normalized site dataframe + output directory."""
    output_dir = _get_persistent_output_dir(tmp_path, f"zfn_chr3_{backend_label}")
    command = [
        "sirnaforge",
        "workflow",
        "CCR5_ZFN_CHR3",
        "--design-mode",
        "zfn",
        "--zfn-left-half-site",
        CCR5_LEFT_HALF_SITE,
        "--zfn-right-half-site",
        CCR5_RIGHT_HALF_SITE,
        "--zfn-search-space",
        str(genome_fasta),
        "--zfn-search-backend",
        backend_label,
        "--zfn-spacer-lengths",
        "5,6",
        "--zfn-max-mismatches",
        "2",
        "--zfn-algorithm",
        "homology",
        "--output-dir",
        str(output_dir),
    ]
    if search_space_index is not None:
        command.extend(["--zfn-search-space-index", str(search_space_index)])

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )

    if result.returncode != 0:
        print(f"Backend {backend_label} STDOUT:\n{result.stdout}")
        print(f"Backend {backend_label} STDERR:\n{result.stderr}")
        pytest.fail(f"ZFN workflow backend={backend_label} exited with code {result.returncode}")

    offtarget_csv = output_dir / "sirnaforge" / "zfn_offtarget_sites.csv"
    assert offtarget_csv.exists(), f"Missing zfn_offtarget_sites.csv for backend={backend_label}"
    normalized = _normalize_offtarget_csv(offtarget_csv)
    assert not normalized.empty, f"No off-target rows produced for backend={backend_label}"
    return normalized, output_dir


# ---------------------------------------------------------------------------
# Test: ZFN CLI workflow smoke test
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_zfn_workflow_cli_smoke(tmp_path: Path) -> None:
    """Smoke test: broad-mismatch CCR5 workflow runs end-to-end via CLI and produces expected outputs.

    This benchmark-style scenario intentionally uses ``exhaustive_python`` because
    the 4-mismatch search budget is outside the practical pyahocorasick path for
    these real CCR5 half-sites.
    """
    output_dir = _get_persistent_output_dir(tmp_path, "zfn_smoke")
    genome_fasta = _build_ccr5_synthetic_genome(tmp_path)

    result = subprocess.run(
        [
            "sirnaforge",
            "workflow",
            "CCR5_ZFN",
            "--design-mode",
            "zfn",
            "--zfn-left-half-site",
            CCR5_LEFT_HALF_SITE,
            "--zfn-right-half-site",
            CCR5_RIGHT_HALF_SITE,
            "--zfn-search-space",
            str(genome_fasta),
            "--zfn-search-backend",
            "exhaustive_python",
            "--zfn-spacer-lengths",
            "5,6",
            "--zfn-max-mismatches",
            "4",
            "--zfn-algorithm",
            "zfn_v2",
            "--output-dir",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )

    if result.returncode != 0:
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        pytest.fail(f"ZFN workflow exited with code {result.returncode}")

    # ── Verify output directory structure ──────────────────────────
    sirnaforge_dir = output_dir / "sirnaforge"
    logs_dir = output_dir / "logs"

    assert sirnaforge_dir.exists(), "Missing sirnaforge output directory"
    assert logs_dir.exists(), "Missing logs directory"

    # Off-target CSV
    offtarget_csv = sirnaforge_dir / "zfn_offtarget_sites.csv"
    assert offtarget_csv.exists(), "Missing zfn_offtarget_sites.csv"

    with offtarget_csv.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        ot_rows = list(reader)

    assert fieldnames is not None
    expected_cols = {
        "site_id",
        "chrom",
        "start_1based",
        "end_1based",
        "strand",
        "orientation",
        "spacer_len",
        "sequence",
        "left_mismatches",
        "right_mismatches",
        "total_mismatches",
        "score",
    }
    assert expected_cols.issubset(set(fieldnames)), f"Missing columns: {expected_cols - set(fieldnames)}"
    assert len(ot_rows) >= 1, "Expected at least one off-target site"

    # Candidate summary JSON
    candidate_json = sirnaforge_dir / "zfn_candidate_summary.json"
    assert candidate_json.exists(), "Missing zfn_candidate_summary.json"
    cand_data = json.loads(candidate_json.read_text())
    assert "candidates" in cand_data
    assert len(cand_data["candidates"]) >= 1
    cand = cand_data["candidates"][0]
    assert "composite_score" in cand
    assert 0 <= cand["composite_score"] <= 100

    # Workflow summary JSON
    summary_json = logs_dir / "workflow_summary.json"
    assert summary_json.exists(), "Missing workflow_summary.json"
    summary = json.loads(summary_json.read_text())
    assert summary["workflow_mode"] == "zfn"
    assert summary["left_half_site"] == CCR5_LEFT_HALF_SITE
    assert summary["right_half_site"] == CCR5_RIGHT_HALF_SITE
    assert summary["total_workflow_time_s"] >= 0


# ---------------------------------------------------------------------------
# Test: CCR5 known-site recovery
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_zfn_ccr5_known_site_recovery(tmp_path: Path) -> None:
    """Verify that the ZFN workflow recovers known CCR5 off-target sites.

    Uses the validated CCR5 on-target + CSNK1G3 + CCR2 sites embedded
    in a synthetic genome.  The searcher should find a perfect-match
    on-target site and mismatched off-target sites.

    This fixture keeps the original broad mismatch budget from the CCR5
    benchmark context, so it is explicitly scoped to ``exhaustive_python``.
    """
    output_dir = _get_persistent_output_dir(tmp_path, "zfn_ccr5_recovery")
    genome_fasta = _build_ccr5_synthetic_genome(tmp_path)

    result = subprocess.run(
        [
            "sirnaforge",
            "workflow",
            "CCR5_ZFN_RECOVERY",
            "--design-mode",
            "zfn",
            "--zfn-left-half-site",
            CCR5_LEFT_HALF_SITE,
            "--zfn-right-half-site",
            CCR5_RIGHT_HALF_SITE,
            "--zfn-search-space",
            str(genome_fasta),
            "--zfn-search-backend",
            "exhaustive_python",
            "--zfn-spacer-lengths",
            "5,6",
            "--zfn-max-mismatches",
            "4",
            "--zfn-algorithm",
            "homology",
            "--output-dir",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )

    if result.returncode != 0:
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        pytest.fail(f"CCR5 recovery workflow exited with code {result.returncode}")

    # Parse off-target sites
    offtarget_csv = output_dir / "sirnaforge" / "zfn_offtarget_sites.csv"
    assert offtarget_csv.exists()
    with offtarget_csv.open(newline="", encoding="utf-8") as fh:
        sites = list(csv.DictReader(fh))

    assert len(sites) >= 2, f"Expected ≥2 sites (on-target + off-target), got {len(sites)}"

    # The on-target site (CCR5) should have 0 total mismatches
    perfect_sites = [s for s in sites if int(s["total_mismatches"]) == 0]
    assert len(perfect_sites) >= 1, "Expected at least one perfect-match (on-target) site"

    # Off-target sites should have mismatches > 0
    offtarget_sites = [s for s in sites if int(s["total_mismatches"]) > 0]
    assert len(offtarget_sites) >= 1, "Expected at least one off-target site with mismatches"

    # Verify scoring is reasonable: perfect match should score highest
    best_perfect = max(float(s["score"]) for s in perfect_sites)
    worst_offtarget = min(float(s["score"]) for s in offtarget_sites) if offtarget_sites else 0
    assert best_perfect >= worst_offtarget, "On-target site should score at least as well as off-targets"

    # Verify workflow summary reflects results
    summary = json.loads((output_dir / "logs" / "workflow_summary.json").read_text())
    assert summary["predicted_sites_total"] >= 1
    assert summary["composite_score"] > 0


# ---------------------------------------------------------------------------
# Test: All three ZFN scoring algorithms produce results
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_zfn_all_algorithms_produce_results(tmp_path: Path) -> None:
    """Verify all three ZFN scoring algorithms (homology, conserved_g, zfn_v2) work end-to-end."""
    genome_fasta = _build_ccr5_synthetic_genome(tmp_path)

    for algo in ("homology", "conserved_g", "zfn_v2"):
        output_dir = _get_persistent_output_dir(tmp_path, f"zfn_algo_{algo}")

        result = subprocess.run(
            [
                "sirnaforge",
                "workflow",
                "CCR5_ZFN_ALGO",
                "--design-mode",
                "zfn",
                "--zfn-left-half-site",
                CCR5_LEFT_HALF_SITE,
                "--zfn-right-half-site",
                CCR5_RIGHT_HALF_SITE,
                "--zfn-search-space",
                str(genome_fasta),
                "--zfn-spacer-lengths",
                "5",
                "--zfn-max-mismatches",
                "3",
                "--zfn-algorithm",
                algo,
                "--output-dir",
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )

        if result.returncode != 0:
            print(f"Algorithm {algo} STDOUT:\n{result.stdout}")
            print(f"Algorithm {algo} STDERR:\n{result.stderr}")
            pytest.fail(f"ZFN workflow with algorithm={algo} exited with code {result.returncode}")

        # Verify each algorithm produced output
        offtarget_csv = output_dir / "sirnaforge" / "zfn_offtarget_sites.csv"
        assert offtarget_csv.exists(), f"Missing off-target CSV for algorithm={algo}"

        summary = json.loads((output_dir / "logs" / "workflow_summary.json").read_text())
        assert summary["workflow_mode"] == "zfn"
        assert summary["composite_score"] > 0, f"Algorithm {algo} produced zero composite score"


# ---------------------------------------------------------------------------
# Test: ZFN dimer mode (homodimer) produces wider search results
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_zfn_homodimer_mode(tmp_path: Path) -> None:
    """Verify homodimer mode finds additional sites vs heterodimer-only.

    This comparison retains the broad CCR5 benchmark mismatch budget and is
    therefore scoped to ``exhaustive_python`` rather than the default runtime
    backend.
    """
    genome_fasta = _build_ccr5_synthetic_genome(tmp_path)

    results_by_mode: dict[str, int] = {}
    for dimer_mode in ("heterodimer_only", "include_homodimers"):
        output_dir = _get_persistent_output_dir(tmp_path, f"zfn_dimer_{dimer_mode}")

        result = subprocess.run(
            [
                "sirnaforge",
                "workflow",
                "CCR5_ZFN_DIMER",
                "--design-mode",
                "zfn",
                "--zfn-left-half-site",
                CCR5_LEFT_HALF_SITE,
                "--zfn-right-half-site",
                CCR5_RIGHT_HALF_SITE,
                "--zfn-search-space",
                str(genome_fasta),
                "--zfn-search-backend",
                "exhaustive_python",
                "--zfn-spacer-lengths",
                "5,6",
                "--zfn-max-mismatches",
                "4",
                "--zfn-dimer-mode",
                dimer_mode,
                "--output-dir",
                str(output_dir),
            ],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )

        if result.returncode != 0:
            print(f"Dimer {dimer_mode} STDERR:\n{result.stderr}")
            pytest.fail(f"ZFN workflow dimer_mode={dimer_mode} exited with code {result.returncode}")

        summary = json.loads((output_dir / "logs" / "workflow_summary.json").read_text())
        results_by_mode[dimer_mode] = summary["predicted_sites_total"]

    # Homodimer mode should find at least as many sites as heterodimer-only
    assert results_by_mode["include_homodimers"] >= results_by_mode["heterodimer_only"], (
        f"Homodimer mode ({results_by_mode['include_homodimers']}) found fewer sites "
        f"than heterodimer-only ({results_by_mode['heterodimer_only']})"
    )


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_zfn_chr3_backend_dataframe_parity(tmp_path: Path) -> None:
    """Run a larger chr3 workflow across backends and assert dataframe-level parity with visible diagnostics."""
    genome_fasta = _build_chr3_scaled_benchmark_genome(tmp_path)

    fm_bundle_dir = tmp_path / "zfn_chr3_fm_bundle"
    fm_bundle_result = subprocess.run(
        [
            "sirnaforge",
            "_internal",
            "zfn-build-search-index",
            "--genome-fasta",
            str(genome_fasta),
            "--search-backend",
            "fm_index",
            "--output-dir",
            str(fm_bundle_dir),
        ],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    if fm_bundle_result.returncode != 0:
        print(f"fm_index bundle build STDOUT:\n{fm_bundle_result.stdout}")
        print(f"fm_index bundle build STDERR:\n{fm_bundle_result.stderr}")
        pytest.fail(f"Failed to build fm_index bundle (exit code {fm_bundle_result.returncode})")

    backend_frames: dict[str, pd.DataFrame] = {}
    backend_output_dirs: dict[str, Path] = {}
    backend_cases: list[tuple[str, Path | None]] = [
        ("exhaustive_python", None),
        ("pyahocorasick", None),
        ("fm_index", None),
        ("fm_index_persisted", fm_bundle_dir),
    ]

    for backend_label, index_path in backend_cases:
        workflow_backend = "fm_index" if backend_label == "fm_index_persisted" else backend_label
        frame, out_dir = _run_zfn_backend_workflow(
            tmp_path=tmp_path,
            backend_label=workflow_backend,
            genome_fasta=genome_fasta,
            search_space_index=index_path,
        )
        backend_frames[backend_label] = frame
        backend_output_dirs[backend_label] = out_dir
        frame.to_csv(out_dir / "backend_normalized_sites.csv", index=False)

    baseline = backend_frames["exhaustive_python"]
    parity_counts: list[dict[str, str | int]] = []
    for backend_label, frame in backend_frames.items():
        parity_counts.append({"backend": backend_label, "rows": int(len(frame))})
        if backend_label == "exhaustive_python":
            continue

        if not frame.equals(baseline):
            merged = baseline.merge(
                frame,
                how="outer",
                indicator=True,
                on=[
                    "chrom",
                    "start_1based",
                    "end_1based",
                    "strand",
                    "orientation",
                    "spacer_len",
                    "sequence",
                    "left_mismatches",
                    "right_mismatches",
                    "total_mismatches",
                    "score",
                ],
            )
            diff_preview = merged[merged["_merge"] != "both"].head(20)
            diff_path = backend_output_dirs[backend_label] / "backend_parity_diff_preview.csv"
            diff_preview.to_csv(diff_path, index=False)
            preview_text = cast(Any, diff_preview).to_string(index=False)
            pytest.fail(
                "chr3 backend dataframe mismatch for "
                f"{backend_label}. Preview saved to {diff_path}.\n"
                f"Visible preview:\n{preview_text}"
            )

    counts_df = pd.DataFrame(parity_counts)
    counts_path = tmp_path / "zfn_chr3_backend_row_counts.csv"
    counts_df.to_csv(counts_path, index=False)
    counts_text = cast(Any, counts_df).to_string(index=False)
    print(f"ZFN chr3 backend parity row-count summary:\n{counts_text}")
    print(f"Saved row-count summary to: {counts_path}")


# ---------------------------------------------------------------------------
# Test: ZFN error handling - missing half-sites
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_zfn_missing_half_sites_gives_error(tmp_path: Path) -> None:
    """CLI should fail gracefully when --design-mode zfn is used without half-sites."""
    output_dir = tmp_path / "zfn_error"
    output_dir.mkdir()

    result = subprocess.run(
        [
            "sirnaforge",
            "workflow",
            "ZFN_MISSING",
            "--design-mode",
            "zfn",
            "--output-dir",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )

    assert result.returncode != 0, "Should fail when half-sites are not provided"
    combined = result.stdout + result.stderr
    assert "half-site" in combined.lower(), "Error message should mention missing half-sites"
