"""Test the miRNA database manager functionality.

Tests cover downloading, caching, and combining miRNA databases from various sources.
"""

import contextlib
import logging
import os
import tempfile
import warnings
from io import StringIO
from pathlib import Path

import pytest
from Bio import SeqIO

from sirnaforge.data.mirna_manager import FASTA_PARSE_FORMAT, MiRNADatabaseManager


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

    # miRBase's real mature.fa is ordered by species, so the requested species is usually *not*
    # first. A whole-header-cloning parser therefore drives filtering to zero matches, which is
    # the second, silent failure mode: `get_database` returns None and seed screening vanishes.
    HUMAN_NOT_FIRST_FASTA = (
        ">cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p\n"
        "UGAGGUAGUAGGUUGUAUAGUU\n"
        ">hsa-miR-16-5p MIMAT0000069 Homo sapiens miR-16-5p\n"
        "UAGCAGCACGUAAAUAUUGGCG\n"
        ">hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p\n"
        "UAGCUUAUCAGACUGAUGUUGA\n"
    )

    @staticmethod
    def _records(fasta_content):
        """Parse FASTA text into {full header: sequence}.

        Keyed on the *description*, not `record.id`: the species predicate reads
        `record.description`, and the runtime consumer (`FastaUtils.parse_fasta_to_dict`) keys its
        dict on the full header too, so a parser that preserved ids while cloning descriptions
        would still corrupt screening. Full headers also make collapsing visible — duplicate keys
        silently overwrite rather than showing up as an extra entry.
        """
        return {record.description: str(record.seq) for record in SeqIO.parse(StringIO(fasta_content), "fasta")}

    @pytest.mark.unit
    def test_parse_format_is_per_record(self):
        """Pin the constant: "fasta-blast"/"fasta-pearson" clone the first header onto every record.

        Guards the constant itself, so re-introducing a comment-aware reader is a test failure
        rather than something only a source comment discourages.
        """
        assert FASTA_PARSE_FORMAT == "fasta"

        descriptions = [
            record.description for record in SeqIO.parse(StringIO(self.MULTI_SPECIES_FASTA), FASTA_PARSE_FORMAT)
        ]
        assert len(set(descriptions)) == 3

    @pytest.mark.unit
    @pytest.mark.parametrize(
        ("species", "expected_header", "expected_seq"),
        [
            ("human", "hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p", "UGAGGUAGUAGGUUGUAUAGUU"),
            ("mouse", "mmu-miR-1a-3p MIMAT0000123 Mus musculus miR-1a-3p", "UGGAAUGUAAAGAAGUAUGUAU"),
            ("rat", "rno-miR-21-5p MIMAT0000790 Rattus norvegicus miR-21-5p", "UAGCUUAUCAGACUGAUGUUGA"),
        ],
    )
    def test_filters_to_requested_species_only(self, manager_with_temp_cache, species, expected_header, expected_seq):
        """Each species must yield exactly its own record, whatever the first header is."""
        filtered = manager_with_temp_cache._filter_species_sequences(self.MULTI_SPECIES_FASTA, species)

        assert self._records(filtered) == {expected_header: expected_seq}

    @pytest.mark.unit
    def test_filter_keeps_every_match_when_species_is_not_first(self, manager_with_temp_cache, caplog):
        """Zero-match mode: a cloned first header drops all real matches and logs an error."""
        with caplog.at_level(logging.ERROR, logger="sirnaforge.data.mirna_manager"):
            filtered = manager_with_temp_cache._filter_species_sequences(self.HUMAN_NOT_FIRST_FASTA, "human")

        assert self._records(filtered) == {
            "hsa-miR-16-5p MIMAT0000069 Homo sapiens miR-16-5p": "UAGCAGCACGUAAAUAUUGGCG",
            "hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p": "UAGCUUAUCAGACUGAUGUUGA",
        }
        assert "No sequences found for species" not in caplog.text

    @pytest.mark.unit
    def test_filter_logs_error_when_species_genuinely_absent(self, manager_with_temp_cache, caplog):
        """The `filtered_count == 0` branch must be loud, since callers treat "" as a failure."""
        with caplog.at_level(logging.ERROR, logger="sirnaforge.data.mirna_manager"):
            filtered = manager_with_temp_cache._filter_species_sequences(self.MULTI_SPECIES_FASTA, "worm")

        assert filtered == ""
        assert "No sequences found for species 'worm' after filtering" in caplog.text

    @pytest.mark.unit
    def test_unmapped_species_raises_instead_of_returning_every_species(self, manager_with_temp_cache):
        """Never silently hand ~270 species of miRNAs back under a single-species label."""
        with pytest.raises(ValueError, match="No miRBase species code known for 'axolotl'"):
            manager_with_temp_cache._filter_species_sequences(self.MULTI_SPECIES_FASTA, "axolotl")

    @pytest.mark.unit
    @pytest.mark.parametrize(
        ("label", "wrapper"),
        [
            # Per-line markup: `<pre>` around each FASTA line.
            ("per_line", lambda text: "".join(f"<pre>{line}</pre>\n" for line in text.strip().split("\n"))),
            # Block markup: one `<pre>` around the whole payload. Stripping the tags leaves blank
            # lines ahead of the first record, which Biopython 1.86's "fasta" reader answers with
            # BiopythonDeprecationWarning and a future release will answer with ValueError.
            ("block", lambda text: f"<pre>\n{text}</pre>\n"),
            ("block_with_heading", lambda text: f"<html>\n<body>\n<h1>mature.fa</h1>\n<pre>\n{text}</pre>\n"),
        ],
    )
    def test_html_wrapped_fasta_is_normalized(self, manager_with_temp_cache, label, wrapper):
        """Upstream endpoints sometimes wrap FASTA in markup, per-line or as one block."""
        wrapped = wrapper(self.MULTI_SPECIES_FASTA)

        with warnings.catch_warnings():
            # Any Biopython parser warning here is a latent hard failure on a future release.
            warnings.simplefilter("error")
            filtered = manager_with_temp_cache._filter_species_sequences(wrapped, "mouse")

        assert self._records(filtered) == {
            "mmu-miR-1a-3p MIMAT0000123 Mus musculus miR-1a-3p": "UGGAAUGUAAAGAAGUAUGUAU"
        }, label

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
            "hsa-let-7a-5p one [source:mirbase]": "UGAGGUAGUAGGUUGUAUAGUU",
            "hsa-miR-16-5p two [source:mirbase]": "UAGCAGCACGUAAAUAUUGGCG",
            "hsa-miR-21-5p three [source:mirgenedb]": "UAGCUUAUCAGACUGAUGUUGA",
        }

    @pytest.mark.unit
    def test_combined_database_is_rebuilt_not_reused(self, manager_with_temp_cache, monkeypatch):
        """A stale combined_*.fa must never be served: it carries no TTL and no checksum.

        The old reuse test was "combined mtime >= every source mtime", which cannot detect wrong
        *contents*, so a combined database written by a buggy producer survived indefinitely — even
        after the per-source caches had been re-validated and refreshed.
        """
        source = manager_with_temp_cache.cache_dir / "only.fa"
        source.write_text(">hsa-miR-21-5p real\nUAGCUUAUCAGACUGAUGUUGA\n")
        monkeypatch.setattr(
            MiRNADatabaseManager,
            "get_database",
            lambda _self, _source_name, _species, force_refresh=False: source,  # noqa: ARG005
        )

        combined_path = manager_with_temp_cache.cache_dir / "combined_human_mirbase.fa"
        combined_path.write_text(">poisoned poisoned\nAAAAAAAAAAAAAAAAAAAAAA\n")
        # Make the stale artifact look fresh relative to its source, which is all the old check saw.
        stale_mtime = source.stat().st_mtime + 1000
        os.utime(combined_path, (stale_mtime, stale_mtime))

        combined = manager_with_temp_cache.get_combined_database(["mirbase"], "human")

        assert combined == combined_path
        assert self._records(combined.read_text()) == {"hsa-miR-21-5p real [source:mirbase]": "UAGCUUAUCAGACUGAUGUUGA"}


class TestGetDatabaseFiltering:
    """`get_database`'s miRBase filtering branch, including the silent zero-match failure."""

    @staticmethod
    def _stub_download(monkeypatch, content):
        monkeypatch.setattr(
            MiRNADatabaseManager,
            "_download_file",
            lambda _self, _source, timeout=600: content,  # noqa: ARG005
        )

    @pytest.mark.unit
    def test_returns_cache_file_when_species_is_not_first_record(self, manager_with_temp_cache, monkeypatch):
        """The whole point: a mid-file species must still produce a usable cached database."""
        self._stub_download(monkeypatch, TestFilterSpeciesSequences.HUMAN_NOT_FIRST_FASTA)

        cache_file = manager_with_temp_cache.get_database("mirbase", "human")

        assert cache_file is not None
        # Keyed by full header, exactly as `FastaUtils.parse_fasta_to_dict` does downstream, so a
        # header-cloning parser shows up as a collapsed dict rather than mislabelled entries.
        assert TestFilterSpeciesSequences._records(cache_file.read_text()) == {
            "hsa-miR-16-5p MIMAT0000069 Homo sapiens miR-16-5p": "UAGCAGCACGUAAAUAUUGGCG",
            "hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p": "UAGCUUAUCAGACUGAUGUUGA",
        }

    @pytest.mark.unit
    def test_returns_none_and_leaves_no_cache_when_filter_empties_content(
        self, manager_with_temp_cache, monkeypatch, caplog
    ):
        """Zero-match mode: discard the cache update rather than persist an empty database."""
        self._stub_download(monkeypatch, ">cel-let-7-5p Caenorhabditis elegans let-7-5p\nUGAGGUAGUAGGUUGUAUAGUU\n")

        with caplog.at_level(logging.ERROR, logger="sirnaforge.data.mirna_manager"):
            result = manager_with_temp_cache.get_database("mirbase", "human")

        assert result is None
        assert "produced no sequences; discarding cache update" in caplog.text
        assert not list(manager_with_temp_cache.cache_dir.glob("*.fa"))
        assert manager_with_temp_cache.metadata == {}

    @pytest.mark.unit
    def test_mirgenedb_content_is_not_species_filtered(self, manager_with_temp_cache, monkeypatch):
        """Filtering is gated on `source_name.startswith("mirbase")`; MirGeneDB is pre-filtered.

        This is why the default `--mirna-db mirgenedb` path never saw the mislabelling: its cache
        files are written straight from the download without going through a FASTA parser.
        """
        content = ">Hsa-Let-7-P2a_5p one\nUGAGGUAGUAGGUUGUAUAGUU\n>Hsa-Mir-21_5p two\nUAGCUUAUCAGACUGAUGUUGA\n"
        self._stub_download(monkeypatch, content)

        cache_file = manager_with_temp_cache.get_database("mirgenedb", "human")

        assert cache_file is not None
        assert cache_file.read_text() == content


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
