#!/usr/bin/env python3
"""Integration tests for consolidated reference caching system.

Tests the ReferenceManager base class through MiRNADatabaseManager and
TranscriptomeManager to ensure caching, metadata, TTL, and index building
work correctly after refactoring.
"""

import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import pytest

from sirnaforge.data.mirna_manager import MiRNADatabaseManager, MiRNASource
from sirnaforge.data.reference_manager import CacheMetadata, ReferenceManager
from sirnaforge.data.transcriptome_manager import TranscriptomeManager, TranscriptomeSource

# Sample data for testing
SAMPLE_HUMAN_MIRNA = """>hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p
UGAGGUAGUAGGUUGUAUAGUU
>hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p
UAGCUUAUCAGACUGAUGUUGA
>hsa-miR-155-5p MIMAT0000646 Homo sapiens miR-155-5p
UUAAUGCUAAUCGUGAUAGGGGU
"""

SAMPLE_TRANSCRIPTOME = """>ENST00000456328.2 transcript:ENST00000456328.2 gene:ENSG00000290825.1
GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGT
CTCTTAGCCCAGACTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTT
>ENST00000450305.2 transcript:ENST00000450305.2 gene:ENSG00000223972.6
ATGGTCAGCTGGGGTCTAGGGACTGGCTTCAGGGAATGGGGGAGAGGGGGAAGAGGACC
AGGGGCTGGGAGATGGAAGAGGTGGGGGTGGAGTGGGGCAGAGACAGGGCAGGGGCCAG
>ENST00000488147.1 transcript:ENST00000488147.1 gene:ENSG00000227232.5
CTCCCCAGGAAGCCCTGCCTGACCTCCTGGCTCCTGGGATGCTCCTTGGTTGTTGAGCT
TCTCCACTGTGAGGTTTGGTGGGGTTTTGGGTGGAGTGAAGGGGGAGAGAAGGAAGAGG
"""


class TestMiRNACaching:
    """Test caching behavior of MiRNADatabaseManager."""

    @pytest.fixture
    def temp_cache_dir(self):
        """Create temporary cache directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def mirna_manager(self, temp_cache_dir):
        """Create MiRNADatabaseManager with temporary cache."""
        return MiRNADatabaseManager(cache_dir=temp_cache_dir, cache_ttl_days=1)

    @pytest.fixture
    def sample_mirna_file(self, temp_cache_dir):
        """Create a sample miRNA file."""
        mirna_file = temp_cache_dir / "human_mirna.fa"
        mirna_file.write_text(SAMPLE_HUMAN_MIRNA)
        return mirna_file

    def test_cache_metadata_creation(self, mirna_manager):
        """Test that cache metadata is created correctly."""
        # Initially empty
        assert len(mirna_manager.metadata) == 0
        info = mirna_manager.cache_info()
        assert info["total_files"] == 0
        assert info["cached_items"] == []

    def test_cache_metadata_persistence(self, mirna_manager, temp_cache_dir):
        """Test that metadata persists across manager instances."""
        # Create a test cache file with metadata
        test_file = temp_cache_dir / "test_cache.fa"
        test_file.write_text(SAMPLE_HUMAN_MIRNA)

        # Manually add metadata (simulating a cached download)
        source = MiRNASource(
            name="test",
            url="http://test.example.com/test.fa",
            species="human",
            format="fasta",
            description="Test source",
        )

        cache_key = source.cache_key()
        mirna_manager.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=test_file.stat().st_size,
            checksum=mirna_manager._compute_file_checksum(test_file),
            file_path=str(test_file),
        )
        mirna_manager._save_metadata()

        # Create new manager instance
        new_manager = MiRNADatabaseManager(cache_dir=temp_cache_dir)

        # Metadata should be loaded
        assert cache_key in new_manager.metadata
        assert new_manager.metadata[cache_key].source.name == "test"
        assert new_manager.metadata[cache_key].file_path == str(test_file)

    def test_cache_ttl_validation(self, temp_cache_dir):
        """Test that TTL-based cache invalidation works."""
        # Create manager with very short TTL
        manager = MiRNADatabaseManager(cache_dir=temp_cache_dir, cache_ttl_days=0)

        test_file = temp_cache_dir / "test_ttl.fa"
        test_file.write_text(SAMPLE_HUMAN_MIRNA)

        source = MiRNASource(
            name="test_ttl",
            url="http://test.example.com/test.fa",
            species="human",
            format="fasta",
            description="Test TTL source",
        )

        cache_key = source.cache_key()

        # Add metadata with old timestamp
        old_time = datetime.now() - timedelta(days=2)
        manager.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=old_time.isoformat(),
            file_size=test_file.stat().st_size,
            checksum=manager._compute_file_checksum(test_file),
            file_path=str(test_file),
        )

        # Cache should be invalid due to TTL
        assert not manager._is_cache_valid(cache_key)

    def test_cache_checksum_validation(self, mirna_manager, temp_cache_dir):
        """Test that corrupted cache files are detected."""
        test_file = temp_cache_dir / "test_checksum.fa"
        test_file.write_text(SAMPLE_HUMAN_MIRNA)

        source = MiRNASource(
            name="test_checksum",
            url="http://test.example.com/test.fa",
            species="human",
            format="fasta",
            description="Test checksum source",
        )

        cache_key = source.cache_key()
        original_checksum = mirna_manager._compute_file_checksum(test_file)

        mirna_manager.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=test_file.stat().st_size,
            checksum=original_checksum,
            file_path=str(test_file),
        )

        # Cache should be valid initially
        assert mirna_manager._is_cache_valid(cache_key)

        # Modify file (corrupt it)
        test_file.write_text(SAMPLE_HUMAN_MIRNA + "\n>corrupted\nAAAA\n")

        # Cache should now be invalid due to checksum mismatch
        assert not mirna_manager._is_cache_valid(cache_key)

    def test_cache_cleaning(self, temp_cache_dir):
        """Test that old cache files are cleaned properly."""
        manager = MiRNADatabaseManager(cache_dir=temp_cache_dir, cache_ttl_days=1)

        # Create old and new cache files
        old_file = temp_cache_dir / "old_cache.fa"
        new_file = temp_cache_dir / "new_cache.fa"
        old_file.write_text(SAMPLE_HUMAN_MIRNA)
        new_file.write_text(SAMPLE_HUMAN_MIRNA)

        old_source = MiRNASource(
            name="old", url="http://test.com/old.fa", species="human", format="fasta", description="Old source"
        )
        new_source = MiRNASource(
            name="new", url="http://test.com/new.fa", species="human", format="fasta", description="New source"
        )

        old_key = old_source.cache_key()
        new_key = new_source.cache_key()

        # Add metadata
        manager.metadata[old_key] = CacheMetadata(
            source=old_source,
            downloaded_at=(datetime.now() - timedelta(days=5)).isoformat(),
            file_size=old_file.stat().st_size,
            checksum=manager._compute_file_checksum(old_file),
            file_path=str(old_file),
        )

        manager.metadata[new_key] = CacheMetadata(
            source=new_source,
            downloaded_at=datetime.now().isoformat(),
            file_size=new_file.stat().st_size,
            checksum=manager._compute_file_checksum(new_file),
            file_path=str(new_file),
        )

        manager._save_metadata()

        # Clean cache older than 2 days
        manager.clean_cache(older_than_days=2)

        # Old file should be removed, new file should remain
        assert not old_file.exists()
        assert new_file.exists()
        assert old_key not in manager.metadata
        assert new_key in manager.metadata


class TestTranscriptomeCaching:
    """Test caching behavior of TranscriptomeManager."""

    @pytest.fixture
    def temp_cache_dir(self):
        """Create temporary cache directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def transcriptome_manager(self, temp_cache_dir):
        """Create TranscriptomeManager with temporary cache."""
        return TranscriptomeManager(
            cache_dir=temp_cache_dir,
            cache_ttl_days=1,
            auto_build_indices=False,  # Disable for faster tests
        )

    @pytest.fixture
    def sample_transcriptome_file(self, temp_cache_dir):
        """Create a sample transcriptome file."""
        transcriptome_file = temp_cache_dir / "transcriptome.fa"
        transcriptome_file.write_text(SAMPLE_TRANSCRIPTOME)
        return transcriptome_file

    def test_custom_transcriptome_caching(self, transcriptome_manager, sample_transcriptome_file):
        """Test that custom transcriptomes are cached correctly."""
        # Process custom transcriptome
        result = transcriptome_manager.get_custom_transcriptome(sample_transcriptome_file, build_index=False)

        assert result is not None
        assert "fasta" in result
        assert result["fasta"].exists()

        # Check metadata was created
        assert len(transcriptome_manager.metadata) == 1

        # Get cache info
        info = transcriptome_manager.cache_info()
        assert info["total_fasta_files"] >= 1
        assert len(info["cached_transcriptomes"]) == 1

    def test_transcriptome_index_metadata(self, temp_cache_dir):
        """Test that index metadata is stored in extra field."""
        manager = TranscriptomeManager(cache_dir=temp_cache_dir, auto_build_indices=False)

        test_file = temp_cache_dir / "test_transcriptome.fa"
        test_file.write_text(SAMPLE_TRANSCRIPTOME)

        source = TranscriptomeSource(
            name="test_txome",
            url="http://test.com/txome.fa",
            species="human",
            format="fasta",
            description="Test transcriptome",
        )

        cache_key = source.cache_key()

        # Create metadata with index info in extra field
        metadata = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=test_file.stat().st_size,
            checksum=manager._compute_file_checksum(test_file),
            file_path=str(test_file),
            extra={"index_path": str(temp_cache_dir / "test_index"), "index_built_at": datetime.now().isoformat()},
        )

        manager.metadata[cache_key] = metadata
        manager._save_metadata()

        # Reload and verify extra data persists
        new_manager = TranscriptomeManager(cache_dir=temp_cache_dir)
        assert cache_key in new_manager.metadata
        assert new_manager.metadata[cache_key].extra is not None
        assert "index_path" in new_manager.metadata[cache_key].extra

        # Test helper methods
        index_path = new_manager._get_index_path(new_manager.metadata[cache_key])
        assert index_path is not None
        assert "test_index" in str(index_path)

    def test_transcriptome_cache_reuse(self, transcriptome_manager, sample_transcriptome_file):
        """Test that cached transcriptomes are reused on subsequent calls."""
        # First call - should cache
        result1 = transcriptome_manager.get_custom_transcriptome(
            sample_transcriptome_file, build_index=False, cache_name="test_reuse"
        )

        assert result1 is not None
        cached_file1 = result1["fasta"]

        # Second call - should use cache
        result2 = transcriptome_manager.get_custom_transcriptome(
            sample_transcriptome_file, build_index=False, cache_name="test_reuse"
        )

        assert result2 is not None
        cached_file2 = result2["fasta"]

        # Should be the same cached file
        assert cached_file1 == cached_file2

        # Should still only have one cache entry
        assert len(transcriptome_manager.metadata) == 1

    def test_cache_info_structure(self, transcriptome_manager, sample_transcriptome_file):
        """Test that cache_info returns correct structure."""
        # Add a cached transcriptome
        transcriptome_manager.get_custom_transcriptome(sample_transcriptome_file, build_index=False)

        info = transcriptome_manager.cache_info()

        # Verify structure
        assert "cache_directory" in info
        assert "total_fasta_files" in info
        assert "total_size_mb" in info
        assert "index_files" in info
        assert "cache_ttl_days" in info
        assert "auto_build_indices" in info
        assert "cached_transcriptomes" in info

        # Verify values
        assert info["total_fasta_files"] >= 1
        assert info["total_size_mb"] > 0
        assert info["cache_ttl_days"] == 1
        assert info["auto_build_indices"] is False
        assert len(info["cached_transcriptomes"]) >= 1


class TestCrossManagerCompatibility:
    """Test that both managers work correctly with shared base class."""

    @pytest.fixture
    def temp_cache_dir(self):
        """Create temporary cache directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    def test_both_managers_use_same_base(self, temp_cache_dir):
        """Test that both managers inherit from ReferenceManager."""
        mirna_manager = MiRNADatabaseManager(cache_dir=temp_cache_dir / "mirna")
        txome_manager = TranscriptomeManager(cache_dir=temp_cache_dir / "txome")

        # Both should be instances of ReferenceManager
        assert isinstance(mirna_manager, ReferenceManager)
        assert isinstance(txome_manager, ReferenceManager)

        # Both should have base class methods
        assert hasattr(mirna_manager, "_compute_file_checksum")
        assert hasattr(txome_manager, "_compute_file_checksum")
        assert hasattr(mirna_manager, "_is_cache_valid")
        assert hasattr(txome_manager, "_is_cache_valid")
        assert hasattr(mirna_manager, "cache_info")
        assert hasattr(txome_manager, "cache_info")

    def test_separate_cache_directories(self, temp_cache_dir):
        """Test that managers maintain separate caches."""
        mirna_dir = temp_cache_dir / "mirna"
        txome_dir = temp_cache_dir / "txome"

        mirna_manager = MiRNADatabaseManager(cache_dir=mirna_dir)
        txome_manager = TranscriptomeManager(cache_dir=txome_dir)

        assert mirna_manager.cache_dir == mirna_dir
        assert txome_manager.cache_dir == txome_dir
        assert mirna_manager.cache_dir != txome_manager.cache_dir

    def test_metadata_isolation(self, temp_cache_dir):
        """Test that metadata is isolated between managers."""
        mirna_manager = MiRNADatabaseManager(cache_dir=temp_cache_dir / "mirna")
        txome_manager = TranscriptomeManager(cache_dir=temp_cache_dir / "txome")

        # Each should have empty metadata initially
        assert len(mirna_manager.metadata) == 0
        assert len(txome_manager.metadata) == 0

        # Add dummy metadata to one
        source = MiRNASource(
            name="test", url="http://test.com/test.fa", species="human", format="fasta", description="Test"
        )

        mirna_manager.metadata["test_key"] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=100,
            checksum="abc123",
            file_path="/tmp/test.fa",
        )

        # Other manager should remain empty
        assert len(mirna_manager.metadata) == 1
        assert len(txome_manager.metadata) == 0
