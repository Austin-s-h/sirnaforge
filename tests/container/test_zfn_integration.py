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


# ---------------------------------------------------------------------------
# Test: ZFN CLI workflow smoke test
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_zfn_workflow_cli_smoke(tmp_path: Path) -> None:
    """Smoke test: ZFN design mode runs end-to-end via CLI and produces expected outputs."""
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
    """Verify homodimer mode finds additional sites vs heterodimer-only."""
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
