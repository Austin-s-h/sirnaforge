"""Regression checks for public download integrity and benchmark score normalization."""

from __future__ import annotations

import hashlib
import io
import runpy
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_download_checksum_failure_preserves_existing_file(tmp_path: Path) -> None:
    """A changed upstream response must not replace a local source file."""
    fetch_file = runpy.run_path(str(ROOT / "scripts/fetch_public_benchmarks.py"))["fetch_file"]
    target = tmp_path / "panel.csv"
    target.write_bytes(b"local copy")
    entry = {
        "path": "panel.csv",
        "url": "https://example.org/panel.csv",
        "sha256": hashlib.sha256(b"expected").hexdigest(),
    }
    with (
        patch("urllib.request.urlopen", return_value=io.BytesIO(b"changed upstream")),
        pytest.raises(ValueError, match="Checksum mismatch"),
    ):
        fetch_file(entry, tmp_path)
    assert target.read_bytes() == b"local copy"
    assert list(tmp_path.iterdir()) == [target]


def test_download_rejects_path_traversal(tmp_path: Path) -> None:
    """Manifest entries cannot write outside the chosen download directory."""
    fetch_file = runpy.run_path(str(ROOT / "scripts/fetch_public_benchmarks.py"))["fetch_file"]
    with pytest.raises(ValueError, match="escapes destination"):
        fetch_file({"path": "../escaped.csv", "url": "https://example.org", "sha256": ""}, tmp_path)


@pytest.mark.parametrize("published_score", [50, None])
def test_linked_components_use_source_profile_weights(tmp_path: Path, published_score: float | None) -> None:
    """Normalize newer profiles without inventing unpublished post-screen scores."""
    load_linked = runpy.run_path(str(ROOT / "scripts/analyze_oligogym_scoring.py"))["load_linked"]
    native = tmp_path / "native"
    candidates = tmp_path / "runs/panel/sirnaforge"
    native.mkdir()
    candidates.mkdir(parents=True)
    pd.DataFrame([{"benchmark_id": "b1", "dataset": "panel", "candidate_id": "c1"}]).to_csv(
        native / "candidate_benchmark_links.csv", index=False
    )
    pd.DataFrame([{"benchmark_id": "b1", "dataset": "panel", "label_processed": 75}]).to_csv(
        native / "observations.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "id": "c1",
                "gc_content": 50,
                "max_poly_run_length": 2,
                "composite_score": published_score,
                "score_off_target": 12.5,
                "score_target_accessibility": 5,
                "score_asymmetry": 15,
                "score_gc_content": 17.5,
            }
        ]
    ).to_csv(candidates / "candidates_all.csv", index=False)
    linked = load_linked(native, tmp_path / "runs", {"panel": "panel"}, np.array([0.25, 0.10, 0.30, 0.35]))
    for term in ("off_target", "target_accessibility", "asymmetry", "gc_content"):
        assert linked.iloc[0][f"component_{term}"] == pytest.approx(50)
    if published_score is None:
        assert pd.isna(linked.iloc[0]["current_score"])
    else:
        assert linked.iloc[0]["current_score"] == pytest.approx(published_score)


def test_representative_reference_does_not_duplicate_longest_mapped_isoform(tmp_path: Path) -> None:
    """Preservation adds alternate isoforms without emitting the longest one twice."""
    prepare = runpy.run_path(str(ROOT / "scripts/prepare_representative_ensembl_transcriptomes.py"))["prepare"]
    source = tmp_path / "input.fa"
    source.write_text(
        ">t1 gene:g1 gene_biotype:protein_coding\nAAAA\n"
        ">t2 gene:g1 gene_biotype:protein_coding\nAA\n"
        ">t3 gene:g2 gene_biotype:protein_coding\nAAA\n"
    )
    mapping = tmp_path / "mapping.csv"
    mapping.write_text("mapping_gene_id\ng1\n")
    output = tmp_path / "output.fa"
    manifest = prepare(source, output, mapping)
    ids = [line.split()[0] for line in output.read_text().splitlines() if line.startswith(">")]
    assert ids == [">t1", ">t2", ">t3"]
    assert manifest["selected_records"] == 3
