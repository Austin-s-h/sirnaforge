"""Container tests for transcriptome filtering and BWA-MEM2 indexing.

These tests are intended to run inside the project Docker image.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from sirnaforge.core.off_target import build_bwa_index, validate_index_files
from sirnaforge.data.transcriptome_filter import TranscriptFilter
from sirnaforge.data.transcriptome_manager import TranscriptomeManager

pytestmark = [pytest.mark.integration, pytest.mark.runs_in_container]


def _skip_if_bwa_mem2_missing() -> None:
    """Skip if BWA-MEM2 is not available in the current environment."""
    try:
        subprocess.run(
            ["bwa-mem2", "version"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except FileNotFoundError:
        pytest.skip("bwa-mem2 not available - run this test in the Docker container")


@pytest.fixture
def tiny_transcriptome(tmp_path: Path) -> Path:
    """Create a tiny transcriptome FASTA for testing."""
    fasta = tmp_path / "tiny_transcriptome.fa"
    fasta.write_text(
        ">ENST001 gene_biotype:protein_coding\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        ">ENST002 gene_biotype:lncRNA\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
        ">ENST003 gene_biotype:protein_coding canonical:1\n"
        "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA\n",
        encoding="utf-8",
    )
    return fasta


@pytest.mark.requires_tools
class TestBWAIndexBuilding:
    """Test BWA-MEM2 index building in container environment."""

    def test_build_index_from_tiny_fasta(self, tiny_transcriptome: Path, tmp_path: Path):
        """Test building BWA-MEM2 index from small FASTA."""
        _skip_if_bwa_mem2_missing()
        index_prefix = tmp_path / "test_index"

        # Build index
        result = build_bwa_index(tiny_transcriptome, index_prefix)

        # Verify index files exist
        assert result == index_prefix
        assert validate_index_files(index_prefix, tool="bwa-mem2")

    def test_index_validation(self, tiny_transcriptome: Path, tmp_path: Path):
        """Test index file validation."""
        _skip_if_bwa_mem2_missing()
        index_prefix = tmp_path / "test_index"

        # Before building, validation should fail
        assert not validate_index_files(index_prefix, tool="bwa-mem2")

        # Build index
        build_bwa_index(tiny_transcriptome, index_prefix)

        # After building, validation should pass
        assert validate_index_files(index_prefix, tool="bwa-mem2")

    def test_index_cache_reuse(self, tiny_transcriptome: Path, tmp_path: Path):
        """Test that cached indices are reused without rebuilding."""
        _skip_if_bwa_mem2_missing()

        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()

        manager = TranscriptomeManager(cache_dir=cache_dir)

        # First call should build index
        result1 = manager.get_custom_transcriptome(tiny_transcriptome, build_index=True)
        assert result1
        assert "index" in result1
        assert result1["index"].exists()

        # Second call should reuse cached index (no rebuild)
        result2 = manager.get_custom_transcriptome(tiny_transcriptome, build_index=True)
        assert result2
        assert "index" in result2
        assert result2["index"].exists()
        assert result1["index"] == result2["index"]


class TestTranscriptomeFiltering:
    """Test transcriptome filtering with real FASTA files."""

    def test_filter_protein_coding(self, tiny_transcriptome: Path, tmp_path: Path):
        """Test filtering for protein-coding transcripts."""
        output = tmp_path / "filtered.fa"
        kept = TranscriptFilter.apply_protein_coding_filter(tiny_transcriptome, output)

        # Should keep 2 protein_coding transcripts
        assert kept == 2

        # Verify output
        content = output.read_text()
        assert "ENST001" in content
        assert "ENST003" in content
        assert "ENST002" not in content

    def test_filter_canonical(self, tiny_transcriptome: Path, tmp_path: Path):
        """Test filtering for canonical transcripts."""
        output = tmp_path / "canonical.fa"
        kept = TranscriptFilter.apply_canonical_filter(tiny_transcriptome, output)

        # Should keep 1 canonical transcript
        assert kept == 1

        # Verify output
        content = output.read_text()
        assert "ENST003" in content

    def test_combined_filters(self, tiny_transcriptome: Path, tmp_path: Path):
        """Test applying multiple filters."""
        output = tmp_path / "filtered.fa"
        filters = ["protein_coding", "canonical_only"]
        kept = TranscriptFilter.apply_combined_filter(tiny_transcriptome, output, filters)

        # Should keep 1 transcript (protein_coding AND canonical)
        assert kept == 1

        # Verify output
        content = output.read_text()
        assert "ENST003" in content

    @pytest.mark.requires_tools
    def test_filtered_transcriptome_indexing(self, tiny_transcriptome: Path, tmp_path: Path):
        """Test that filtered transcriptomes can be indexed."""
        _skip_if_bwa_mem2_missing()

        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()

        # Copy to cache as a "source"
        cached_fasta = cache_dir / "test.fa"
        cached_fasta.write_text(tiny_transcriptome.read_text(encoding="utf-8"), encoding="utf-8")

        # Simulate filtering and indexing
        filtered = cache_dir / "filtered.fa"
        TranscriptFilter.apply_protein_coding_filter(cached_fasta, filtered)

        # Build index on filtered FASTA
        index_prefix = cache_dir / "filtered_index"
        build_bwa_index(filtered, index_prefix)

        # Verify index was built successfully
        assert validate_index_files(index_prefix, tool="bwa-mem2")
