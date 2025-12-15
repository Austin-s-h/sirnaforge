"""Container tests for BWA-MEM2 indexing and transcriptome filtering.

These tests require Docker/container environment with BWA-MEM2 installed.
Run with: make test-release
"""

import pytest

from sirnaforge.core.off_target import build_bwa_index, validate_index_files
from sirnaforge.data.transcriptome_filter import TranscriptFilter
from sirnaforge.data.transcriptome_manager import TranscriptomeManager

pytestmark = [pytest.mark.runs_in_container]


@pytest.fixture
def tiny_transcriptome(tmp_path):
    """Create a tiny transcriptome FASTA for testing."""
    fasta = tmp_path / "tiny_transcriptome.fa"
    fasta.write_text(
        ">ENST001 gene_biotype:protein_coding\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        ">ENST002 gene_biotype:lncRNA\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
        ">ENST003 gene_biotype:protein_coding canonical:1\n"
        "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA\n"
    )
    return fasta


@pytest.mark.integration
@pytest.mark.requires_tools
class TestBWAIndexBuilding:
    """Test BWA-MEM2 index building in container environment."""

    def test_build_index_from_tiny_fasta(self, tiny_transcriptome, tmp_path):
        """Test building BWA-MEM2 index from small FASTA."""
        index_prefix = tmp_path / "test_index"

        # Build index
        result = build_bwa_index(tiny_transcriptome, index_prefix)

        # Verify index files exist
        assert result == index_prefix
        assert validate_index_files(index_prefix, tool="bwa-mem2")

    def test_index_validation(self, tiny_transcriptome, tmp_path):
        """Test index file validation."""
        index_prefix = tmp_path / "test_index"

        # Before building, validation should fail
        assert not validate_index_files(index_prefix, tool="bwa-mem2")

        # Build index
        build_bwa_index(tiny_transcriptome, index_prefix)

        # After building, validation should pass
        assert validate_index_files(index_prefix, tool="bwa-mem2")

    def test_index_cache_reuse(self, tiny_transcriptome, tmp_path, monkeypatch):
        """Test that cached indices are reused without rebuilding."""
        manager = TranscriptomeManager()

        # Use temp directory for cache
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        monkeypatch.setattr(manager, "cache_dir", cache_dir)

        # First call should build index
        result1 = manager.get_custom_transcriptome(tiny_transcriptome, build_index=True)
        assert result1 is not None
        assert "index" in result1

        # Second call should reuse cached index (no rebuild)
        result2 = manager.get_custom_transcriptome(tiny_transcriptome, build_index=True)
        assert result2 is not None
        assert "index" in result2
        assert result1["index"] == result2["index"]


@pytest.mark.integration
class TestTranscriptomeFiltering:
    """Test transcriptome filtering with real FASTA files."""

    def test_filter_protein_coding(self, tiny_transcriptome, tmp_path):
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

    def test_filter_canonical(self, tiny_transcriptome, tmp_path):
        """Test filtering for canonical transcripts."""
        output = tmp_path / "canonical.fa"
        kept = TranscriptFilter.apply_canonical_filter(tiny_transcriptome, output)

        # Should keep 1 canonical transcript
        assert kept == 1

        # Verify output
        content = output.read_text()
        assert "ENST003" in content

    def test_combined_filters(self, tiny_transcriptome, tmp_path):
        """Test applying multiple filters."""
        output = tmp_path / "filtered.fa"
        filters = ["protein_coding", "canonical_only"]
        kept = TranscriptFilter.apply_combined_filter(tiny_transcriptome, output, filters)

        # Should keep 1 transcript (protein_coding AND canonical)
        assert kept == 1

        # Verify output
        content = output.read_text()
        assert "ENST003" in content

    def test_filtered_transcriptome_indexing(self, tiny_transcriptome, tmp_path, monkeypatch):
        """Test that filtered transcriptomes can be indexed."""
        manager = TranscriptomeManager()

        # Use temp directory for cache
        cache_dir = tmp_path / "cache"
        cache_dir.mkdir()
        monkeypatch.setattr(manager, "cache_dir", cache_dir)

        # Copy to cache as a "source"
        cached_fasta = cache_dir / "test.fa"
        cached_fasta.write_text(tiny_transcriptome.read_text())

        # Simulate filtering and indexing
        filtered = cache_dir / "filtered.fa"
        TranscriptFilter.apply_protein_coding_filter(cached_fasta, filtered)

        # Build index on filtered FASTA
        index_prefix = cache_dir / "filtered_index"
        _ = build_bwa_index(filtered, index_prefix)

        # Verify index was built successfully
        assert validate_index_files(index_prefix, tool="bwa-mem2")
