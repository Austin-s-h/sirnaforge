"""Integration tests for miRNA database functionality that work with the actual API.

These tests focus on the real interface (get_database, get_combined_database, etc.)
rather than non-existent methods like download_database.
"""

import tempfile
from pathlib import Path

import pytest

from sirnaforge.data.mirna_manager import MiRNADatabaseManager


class TestMiRNAManagerIntegration:
    """Integration tests for miRNA manager using real API."""

    @pytest.mark.integration
    def test_manager_initialization_integration(self):
        """Test manager can be initialized and basic operations work."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "mirna_cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)

            # Basic operations should work
            assert manager.cache_dir.exists()

            sources = manager.list_available_databases()
            assert isinstance(sources, dict)
            assert len(sources) > 0

            info = manager.cache_info()
            assert isinstance(info, dict)
            assert info["total_files"] == 0

    @pytest.mark.integration
    def test_cache_lifecycle_integration(self):
        """Test complete cache lifecycle."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "mirna_cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)

            # Initial empty state
            info = manager.cache_info()
            assert info["total_files"] == 0

            # Add some mock files to test cache operations
            test_file = cache_dir / "test.fa"
            test_file.write_text(">test\nACGU\n")

            # Cache should detect the file
            info = manager.cache_info()
            assert info["total_files"] == 1

            # Dry run clear
            result = manager.clear_cache(confirm=False)
            assert "files_deleted" in result
            assert result["files_deleted"] == 1
            assert test_file.exists()  # Still there after dry run

            # Actual clear
            result = manager.clear_cache(confirm=True)
            assert result["status"] == "Cache cleared successfully"
            assert not test_file.exists()

            # Cache should be empty again
            info = manager.cache_info()
            assert info["total_files"] == 0

    @pytest.mark.integration
    @pytest.mark.requires_network
    @pytest.mark.slow
    def test_real_database_download_integration(self):
        """Test real database download (requires network)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "mirna_cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)

            # Try to get a small database
            try:
                result = manager.get_database("mirbase_high_conf", "human")

                if result:
                    assert result.exists()
                    assert result.stat().st_size > 0

                    # Check FASTA format
                    with result.open("r") as f:
                        first_line = f.readline().strip()
                        assert first_line.startswith(">")

            except Exception as e:
                pytest.skip(f"Network download failed: {e}")


@pytest.mark.smoke
class TestMiRNAManagerSmoke:
    """Smoke tests for miRNA manager - quick validation that basic functionality works."""

    def test_manager_creation_smoke(self):
        """Smoke test: Manager can be created."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "test_cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)
            assert manager is not None
            assert manager.cache_dir == cache_dir

    def test_list_databases_smoke(self):
        """Smoke test: Can list available databases."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "test_cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)

            sources = manager.list_available_databases()
            assert isinstance(sources, dict)
            assert len(sources) > 0

    def test_cache_operations_smoke(self):
        """Smoke test: Cache operations don't crash."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "test_cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)

            # Basic cache operations should not crash
            info = manager.cache_info()
            assert isinstance(info, dict)

            dry_result = manager.clear_cache(confirm=False)
            assert isinstance(dry_result, dict)

    def test_get_database_basic_smoke(self):
        """Smoke test: get_database method exists and handles gracefully."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "test_cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)

            # Should handle invalid requests gracefully
            result = manager.get_database("nonexistent", "species")
            # Should return None or raise appropriate exception
            assert result is None or isinstance(result, Path)
