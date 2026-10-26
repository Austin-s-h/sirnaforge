"""Real BWA-MEM2 index/alignment smoke used by the isolated Pixi review lane."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from sirnaforge.core.off_target import BwaAnalyzer, build_bwa_index, validate_index_files

pytestmark = [pytest.mark.integration, pytest.mark.requires_tools]


def test_bwa_mem2_index_and_exact_guide_hit(tmp_path: Path) -> None:
    """The locked environment can build a fresh index and preserve an exact guide hit."""
    executable = shutil.which("bwa-mem2")
    if executable is None:
        pytest.skip("bwa-mem2 is not installed in this environment")

    version = subprocess.run([executable, "version"], capture_output=True, text=True, check=False, timeout=10)
    if version.returncode != 0:
        pytest.skip(f"bwa-mem2 is not runnable: {version.stderr.strip()}")

    guide = "ACGTTGCACTGATCGTACG"
    reference = tmp_path / "tiny_transcriptome.fa"
    reference.write_text(f">ENST_PIXI gene_biotype:protein_coding\n{'G' * 50}{guide}{'C' * 50}\n", encoding="utf-8")
    index_prefix = tmp_path / "tiny_transcriptome_index"

    build_bwa_index(reference, index_prefix, aligner="bwa-mem2")

    assert validate_index_files(index_prefix, tool="bwa-mem2")
    hits, outcome = BwaAnalyzer(index_prefix, aligner="bwa-mem2").analyze_sequences_with_outcome({"guide_exact": guide})

    assert outcome.completed
    assert outcome.submitted == outcome.processed == 1
    assert any(
        hit["qname"] == "guide_exact"
        and hit["rname"] == "ENST_PIXI"
        and hit["nm"] == 0
        and hit["cigar"] == f"{len(guide)}M"
        for hit in hits
    )
