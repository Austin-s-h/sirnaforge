"""Tests for unified cache management system.

These tests validate cache operations for both miRNA and transcriptome databases.
They are designed to run with `uv run pytest` without requiring Docker.
"""

import pytest

from sirnaforge.utils.unified_cache import UnifiedCacheManager


@pytest.fixture
def mock_transcriptome_cache(tmp_path, monkeypatch):
    """Create a minimal mock transcriptome cache."""
    cache_dir = tmp_path / "transcriptomes"
    cache_dir.mkdir()

    # Create mock files
    (cache_dir / "test.fa").write_text(">transcript1\nATCG\n")
    (cache_dir / "test_index.amb").write_text("index")

    return cache_dir


@pytest.mark.dev
class TestUnifiedCacheManager:
    """Test suite for UnifiedCacheManager."""

    def test_initialization(self):
        """Test that UnifiedCacheManager initializes correctly."""
        manager = UnifiedCacheManager()
        assert manager.mirna is not None
        assert manager.transcriptome is not None

    def test_get_info_structure(self):
        """Test cache info returns expected structure."""
        manager = UnifiedCacheManager()
        info = manager.get_info()

        for component in ["mirna", "transcriptome"]:
            assert component in info
            assert "cache_directory" in info[component]
            assert "total_files" in info[component]
            assert "total_size_mb" in info[component]

    def test_get_total_stats(self):
        """Test combined statistics."""
        manager = UnifiedCacheManager()
        stats = manager.get_total_stats()

        assert "total_files" in stats
        assert "total_size_mb" in stats
        assert stats["total_files"] >= 0
        assert stats["total_size_mb"] >= 0.0

    def test_clear_transcriptome_dry_run(self, mock_transcriptome_cache, monkeypatch):
        """Test dry run doesn't delete files."""
        manager = UnifiedCacheManager()
        monkeypatch.setattr(manager.transcriptome, "cache_dir", mock_transcriptome_cache)

        initial_files = list(mock_transcriptome_cache.glob("*"))
        results = manager.clear(clear_transcriptome=True, dry_run=True)

        # Files should still exist
        assert len(list(mock_transcriptome_cache.glob("*"))) == len(initial_files)
        assert "dry run" in results["transcriptome"]["status"]

    def test_clear_transcriptome_actual(self, mock_transcriptome_cache, monkeypatch):
        """Test actual clear removes files."""
        manager = UnifiedCacheManager()
        monkeypatch.setattr(manager.transcriptome, "cache_dir", mock_transcriptome_cache)

        results = manager.clear(clear_transcriptome=True, dry_run=False)

        # Files should be deleted
        assert len(list(mock_transcriptome_cache.glob("*"))) == 0
        assert results["transcriptome"]["files_deleted"] == 2
        assert results["transcriptome"]["status"] == "cleared"
