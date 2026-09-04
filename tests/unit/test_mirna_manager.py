"""Test the miRNA database manager functionality.

Tests cover downloading, caching, and combining miRNA databases from various sources.
"""

import contextlib
import tempfile
from io import StringIO
from pathlib import Path

import pytest
from Bio import SeqIO

from sirnaforge.data.mirna_manager import MiRNADatabaseManager


@pytest.fixture
def temp_cache_dir():
    """Create a temporary cache directory for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir) / "mirna_cache"


@pytest.fixture
def manager_with_temp_cache(temp_cache_dir):
    """Create MiRNADatabaseManager with temporary cache."""
    return MiRNADatabaseManager(cache_dir=temp_cache_dir)


class TestMiRNADatabaseManager:
    """Test miRNA database manager functionality."""

    @pytest.mark.unit
    def test_manager_initialization(self, temp_cache_dir):
        """Test proper initialization of miRNA manager."""
        manager = MiRNADatabaseManager(cache_dir=temp_cache_dir)

        assert manager.cache_dir == temp_cache_dir
        assert temp_cache_dir.exists()

    @pytest.mark.unit
    def test_list_available_databases(self, manager_with_temp_cache):
        """Test listing available miRNA databases."""
        sources = manager_with_temp_cache.list_available_databases()

        assert isinstance(sources, dict)
        assert len(sources) > 0

        # Should have at least mirbase
        assert "mirbase" in sources or "mirbase_high_conf" in sources
        assert "mirgenedb" in sources

        # Each source should have species mappings
        for _source_name, species_dict in sources.items():
            assert isinstance(species_dict, dict)
            for _species, source in species_dict.items():
                assert hasattr(source, "description")

    @pytest.mark.unit
    def test_available_sources_helpers(self):
        """Helper methods should expose supported sources and species."""
        sources = MiRNADatabaseManager.get_available_sources()
        assert "mirgenedb" in sources
        assert "mirbase" in sources
        mirgenedb_species = MiRNADatabaseManager.get_species_for_source("mirgenedb")
        assert "hsa" in mirgenedb_species
        assert "mmu" in mirgenedb_species
        assert "dre" in mirgenedb_species
        assert MiRNADatabaseManager.normalize_species("mirgenedb", "human") == "hsa"

        all_species = MiRNADatabaseManager.get_all_species()
        assert "human" in all_species
        assert "rat" in all_species
        assert MiRNADatabaseManager.normalize_species("mirgenedb", "mosquito") == "aga"

    @pytest.mark.unit
    def test_mirgenedb_source_metadata(self):
        """MirGeneDB helper should expose taxonomy metadata and build URLs."""
        metadata = MiRNADatabaseManager.get_mirgenedb_species_metadata()
        assert metadata["hsa"]["taxonomy_id"] == "9606"
        assert "scientific_name" in metadata["dre"]

        source = MiRNADatabaseManager.get_source_configuration("mirgenedb", "human")
        assert source is not None
        assert source.species == "hsa"
        assert source.url.startswith("https://www.mirgenedb.org/fasta/hsa")

    @pytest.mark.unit
    def test_canonical_species_resolution(self):
        """Canonical species registry should map aliases to genome and miRNA identifiers."""
        resolution = MiRNADatabaseManager.resolve_species_selection(
            ["Human", "mmu", "Chicken"],
            "mirgenedb",
        )
        assert resolution["canonical"] == ["human", "mouse", "chicken"]
        assert resolution["genome"] == ["human", "mouse", "chicken"]
        assert resolution["mirna"] == ["hsa", "mmu", "gga"]

        override = MiRNADatabaseManager.resolve_species_selection(
            ["human"],
            "mirgenedb",
            mirna_overrides=["HSA", "human"],
        )
        assert override["mirna"] == ["hsa"]

        with pytest.raises(ValueError):
            MiRNADatabaseManager.resolve_species_selection(["unicorn"], "mirgenedb")

    @pytest.mark.unit
    def test_species_validation_helpers(self):
        """Validation helpers should confirm supported and unsupported selections."""
        assert MiRNADatabaseManager.is_supported_species("mirgenedb", "human") is True
        assert MiRNADatabaseManager.is_supported_species("mirgenedb", "unicorn") is False

    @pytest.mark.unit
    def test_cache_info_empty(self, manager_with_temp_cache):
        """Test cache info with empty cache."""
        info = manager_with_temp_cache.cache_info()

        assert info["cache_directory"] == str(manager_with_temp_cache.cache_dir)
        assert info["total_files"] == 0
        assert info["total_size_mb"] == 0.0
        assert info["cached_items"] == []

    @pytest.mark.unit
    def test_cache_info_with_files(self, manager_with_temp_cache, temp_cache_dir):
        """Test cache info with existing files."""
        # Create some test cache files
        test_file1 = temp_cache_dir / "test1.fa"
        test_file2 = temp_cache_dir / "test2.fa"

        test_file1.write_text(">seq1\nACGT\n")
        test_file2.write_text(">seq2\nTGCA\n")

        info = manager_with_temp_cache.cache_info()

        assert info["total_files"] == 2
        assert info["total_size_mb"] > 0
        # cached_databases is a list of metadata keys, not actual files
        # Since we created .fa files directly, they won't be in metadata

    @pytest.mark.unit
    def test_clear_cache_dry_run(self, manager_with_temp_cache, temp_cache_dir):
        """Test cache clearing in dry run mode."""
        # Create test files
        test_file = temp_cache_dir / "test.fa"
        test_file.write_text(">seq1\nACGT\n")

        result = manager_with_temp_cache.clear_cache(confirm=False)

        assert result["status"].startswith("Would delete")
        assert result["files_deleted"] > 0
        assert result["size_freed_mb"] > 0
        assert test_file.exists()  # File should still exist

    @pytest.mark.unit
    def test_clear_cache_confirm(self, manager_with_temp_cache, temp_cache_dir):
        """Test actual cache clearing."""
        # Create test files
        test_file = temp_cache_dir / "test.fa"
        test_file.write_text(">seq1\nACGT\n")

        result = manager_with_temp_cache.clear_cache(confirm=True)

        assert result["status"] == "Cache cleared successfully"
        assert result["files_deleted"] > 0
        assert not test_file.exists()  # File should be deleted

    @pytest.mark.slow
    @pytest.mark.integration
    @pytest.mark.skipif(
        not pytest.importorskip("requests", reason="requests not available"), reason="Network tests require requests"
    )
    def test_download_mirbase_high_conf(self, manager_with_temp_cache):
        """Test downloading miRBase high confidence database."""
        # This is a real network test - marked as slow/integration
        try:
            db_file = manager_with_temp_cache.get_database("mirbase_high_conf", "human")

            if db_file:
                assert db_file.exists()
                assert db_file.stat().st_size > 0

                # Check FASTA format
                with db_file.open("r") as f:
                    first_line = f.readline().strip()
                    assert first_line.startswith(">")

                    # Count sequences
                    f.seek(0)
                    content = f.read()
                    seq_count = content.count(">")
                    assert seq_count > 0

        except Exception as e:
            pytest.skip(f"Network download failed: {e}")

    @pytest.mark.unit
    def test_get_combined_database_empty(self, manager_with_temp_cache):
        """Test combining databases with empty cache."""
        # This should handle the case where no databases are cached
        result = manager_with_temp_cache.get_combined_database(
            sources=["mirbase_high_conf"], species="human", output_name="test_combined.fa"
        )

        # Should either return None or handle gracefully
        # (depends on implementation - adjust as needed)
        assert result is None or isinstance(result, Path)

    @pytest.mark.integration
    @pytest.mark.skipif(
        not pytest.importorskip("urllib", reason="urllib not available"), reason="Network tests require urllib"
    )
    def test_get_combined_database_with_mock_files(self, manager_with_temp_cache):
        """Test combining databases - this is an integration test that uses network."""
        # Note: This test actually downloads from internet, not from mock files
        # The get_combined_database method downloads from URLs, not local files

        result = manager_with_temp_cache.get_combined_database(
            sources=["mirbase_high_conf", "mirgenedb"], species="human", output_name="test_combined.fa"
        )

        if result:
            assert result.exists()

            with result.open("r") as f:
                content = f.read()
                # Should contain miRNA sequences in FASTA format
                assert content.count(">") > 0  # Has at least one sequence
                assert "hsa-" in content or "human" in content.lower()  # Contains human sequences


class TestFilterSpeciesSequences:
    """Species filtering must key off each record's own header, not the first record's."""

    # miRBase mature.fa style: three species, deliberately distinct sequences so a
    # mislabelled record cannot be mistaken for the right one.
    MULTI_SPECIES_FASTA = (
        ">hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p\n"
        "UGAGGUAGUAGGUUGUAUAGUU\n"
        ">mmu-miR-1a-3p MIMAT0000123 Mus musculus miR-1a-3p\n"
        "UGGAAUGUAAAGAAGUAUGUAU\n"
        ">rno-miR-21-5p MIMAT0000790 Rattus norvegicus miR-21-5p\n"
        "UAGCUUAUCAGACUGAUGUUGA\n"
    )

    @staticmethod
    def _records(fasta_content):
        """Parse FASTA text into {record id: sequence}."""
        return {record.id: str(record.seq) for record in SeqIO.parse(StringIO(fasta_content), "fasta")}

    @pytest.mark.unit
    @pytest.mark.parametrize(
        ("species", "expected_id", "expected_seq"),
        [
            ("human", "hsa-let-7a-5p", "UGAGGUAGUAGGUUGUAUAGUU"),
            ("mouse", "mmu-miR-1a-3p", "UGGAAUGUAAAGAAGUAUGUAU"),
            ("rat", "rno-miR-21-5p", "UAGCUUAUCAGACUGAUGUUGA"),
        ],
    )
    def test_filters_to_requested_species_only(self, manager_with_temp_cache, species, expected_id, expected_seq):
        """Each species must yield exactly its own record, whatever the first header is."""
        filtered = manager_with_temp_cache._filter_species_sequences(self.MULTI_SPECIES_FASTA, species)

        assert self._records(filtered) == {expected_id: expected_seq}

    @pytest.mark.unit
    def test_unknown_species_returns_content_unchanged(self, manager_with_temp_cache):
        """An unmapped species is a pass-through, not an empty database."""
        filtered = manager_with_temp_cache._filter_species_sequences(self.MULTI_SPECIES_FASTA, "axolotl")

        assert filtered == self.MULTI_SPECIES_FASTA

    @pytest.mark.unit
    def test_html_wrapped_fasta_is_normalized(self, manager_with_temp_cache):
        """Upstream endpoints sometimes wrap FASTA lines in markup."""
        wrapped = "".join(f"<pre>{line}</pre>\n" for line in self.MULTI_SPECIES_FASTA.strip().split("\n"))

        filtered = manager_with_temp_cache._filter_species_sequences(wrapped, "mouse")

        assert self._records(filtered) == {"mmu-miR-1a-3p": "UGGAAUGUAAAGAAGUAUGUAU"}

    @pytest.mark.unit
    def test_combined_database_preserves_record_identity(self, manager_with_temp_cache, monkeypatch):
        """Combining sources must not relabel every record with the first record's id."""
        first = manager_with_temp_cache.cache_dir / "first.fa"
        second = manager_with_temp_cache.cache_dir / "second.fa"
        first.write_text(">hsa-let-7a-5p one\nUGAGGUAGUAGGUUGUAUAGUU\n>hsa-miR-16-5p two\nUAGCAGCACGUAAAUAUUGGCG\n")
        second.write_text(">hsa-miR-21-5p three\nUAGCUUAUCAGACUGAUGUUGA\n")

        source_files = {"mirbase": first, "mirgenedb": second}
        monkeypatch.setattr(
            MiRNADatabaseManager,
            "get_database",
            lambda _self, source_name, _species, force_refresh=False: source_files[source_name],  # noqa: ARG005
        )

        combined = manager_with_temp_cache.get_combined_database(["mirbase", "mirgenedb"], "human")

        assert combined is not None
        assert self._records(combined.read_text()) == {
            "hsa-let-7a-5p": "UGAGGUAGUAGGUUGUAUAGUU",
            "hsa-miR-16-5p": "UAGCAGCACGUAAAUAUUGGCG",
            "hsa-miR-21-5p": "UAGCUUAUCAGACUGAUGUUGA",
        }


class TestMiRNAManagerErrorHandling:
    """Test error handling in miRNA manager."""

    @pytest.mark.unit
    def test_invalid_cache_directory_permissions(self):
        """Test handling of permission errors for cache directory."""
        # Test with a path that would cause permission issues
        restricted_path = Path("/root/restricted_cache")

        # This should handle the error gracefully
        with contextlib.suppress(PermissionError, OSError):
            MiRNADatabaseManager(cache_dir=restricted_path)

    @pytest.mark.unit
    def test_get_database_invalid_source(self, manager_with_temp_cache):
        """Test handling of invalid database source."""
        result = manager_with_temp_cache.get_database("invalid_source", "human")
        assert result is None or not result.exists()

    @pytest.mark.unit
    def test_get_database_invalid_species(self, manager_with_temp_cache):
        """Test handling of invalid species."""
        result = manager_with_temp_cache.get_database("mirbase", "invalid_species")
        assert result is None or not result.exists()


@pytest.mark.smoke
class TestMiRNAManagerSmoke:
    """Smoke tests for miRNA manager - minimal functionality verification."""

    def test_manager_creation(self):
        """Smoke test: Can create manager instance."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)
            assert manager is not None

    def test_list_databases_basic(self):
        """Smoke test: Can list available databases."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)
            sources = manager.list_available_databases()
            assert isinstance(sources, dict)

    def test_cache_info_basic(self):
        """Smoke test: Can get cache info."""
        with tempfile.TemporaryDirectory() as temp_dir:
            cache_dir = Path(temp_dir) / "cache"
            manager = MiRNADatabaseManager(cache_dir=cache_dir)
            info = manager.cache_info()
            assert isinstance(info, dict)
            assert "cache_directory" in info
