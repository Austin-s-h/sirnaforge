"""Tests for chemical modifications metadata system."""

import json
import tempfile
from pathlib import Path

import pytest

from sirnaforge.models.modifications import (
    ChemicalModification,
    ConfirmationStatus,
    Provenance,
    SequenceRecord,
    SourceType,
    StrandMetadata,
    StrandRole,
)
from sirnaforge.modifications import (
    load_metadata,
    merge_metadata_into_fasta,
    parse_chem_mods,
    parse_header,
    parse_provenance,
)


@pytest.mark.unit
@pytest.mark.local_python
class TestChemicalModification:
    """Test ChemicalModification model."""

    def test_create_modification(self):
        """Test creating a chemical modification."""
        mod = ChemicalModification(type="2OMe", positions=[1, 4, 6, 11])
        assert mod.type == "2OMe"
        assert mod.positions == [1, 4, 6, 11]

    def test_modification_to_header_string(self):
        """Test conversion to FASTA header format."""
        mod = ChemicalModification(type="2OMe", positions=[1, 4, 6])
        assert mod.to_header_string() == "2OMe(1,4,6)"

    def test_modification_empty_positions(self):
        """Test modification with no positions."""
        mod = ChemicalModification(type="2F", positions=[])
        assert mod.to_header_string() == "2F()"


@pytest.mark.unit
@pytest.mark.local_python
class TestProvenance:
    """Test Provenance model."""

    def test_create_provenance(self):
        """Test creating provenance information."""
        prov = Provenance(
            source_type=SourceType.PATENT,
            identifier="US10060921B2",
            url="https://patents.google.com/patent/US10060921B2",
        )
        assert prov.source_type == SourceType.PATENT
        assert prov.identifier == "US10060921B2"

    def test_provenance_to_header_string(self):
        """Test conversion to FASTA header format."""
        prov = Provenance(source_type=SourceType.PATENT, identifier="US10060921B2")
        assert prov.to_header_string() == "Patent:US10060921B2"


@pytest.mark.unit
@pytest.mark.local_python
class TestStrandMetadata:
    """Test StrandMetadata model."""

    def test_create_strand_metadata(self):
        """Test creating strand metadata."""
        mods = [ChemicalModification(type="2OMe", positions=[1, 4, 6])]
        prov = Provenance(source_type=SourceType.PATENT, identifier="US10060921B2")

        metadata = StrandMetadata(
            id="patisiran_ttr_guide",
            sequence="AUGGAAUACUCUUGGUUAC",
            overhang="dTdT",
            chem_mods=mods,
            provenance=prov,
            confirmation_status=ConfirmationStatus.PENDING,
        )

        assert metadata.id == "patisiran_ttr_guide"
        assert metadata.sequence == "AUGGAAUACUCUUGGUUAC"
        assert metadata.overhang == "dTdT"
        assert len(metadata.chem_mods) == 1

    def test_strand_metadata_to_fasta_header(self):
        """Test FASTA header generation."""
        mods = [ChemicalModification(type="2OMe", positions=[1, 4, 6])]
        metadata = StrandMetadata(
            id="test_guide",
            sequence="AUCGAUCGAUCGAUCGAUCGA",
            overhang="dTdT",
            chem_mods=mods,
        )

        header = metadata.to_fasta_header(target_gene="TTR", strand_role=StrandRole.GUIDE)
        assert ">test_guide" in header
        assert "Target=TTR" in header
        assert "Role=guide" in header
        assert "Overhang=dTdT" in header
        assert "ChemMods=2OMe(1,4,6)" in header


@pytest.mark.unit
@pytest.mark.local_python
class TestSequenceRecord:
    """Test SequenceRecord model."""

    def test_create_sequence_record(self):
        """Test creating a complete sequence record."""
        metadata = StrandMetadata(
            id="test_guide",
            sequence="AUCGAUCGAUCGAUCGAUCGA",
        )

        record = SequenceRecord(
            target_gene="TTR",
            strand_role=StrandRole.GUIDE,
            metadata=metadata,
        )

        assert record.target_gene == "TTR"
        assert record.strand_role == StrandRole.GUIDE

    def test_sequence_record_to_fasta(self):
        """Test FASTA output from sequence record."""
        metadata = StrandMetadata(
            id="test_guide",
            sequence="AUCGAUCGAUCGAUCGAUCGA",
        )

        record = SequenceRecord(
            target_gene="TTR",
            strand_role=StrandRole.GUIDE,
            metadata=metadata,
        )

        fasta = record.to_fasta()
        assert ">test_guide" in fasta
        assert "AUCGAUCGAUCGAUCGAUCGA" in fasta


@pytest.mark.unit
@pytest.mark.local_python
class TestParsingFunctions:
    """Test parsing helper functions."""

    def test_parse_chem_mods(self):
        """Test parsing ChemMods string."""
        mods = parse_chem_mods("2OMe(1,4,6,11)|2F()")
        assert len(mods) == 2
        assert mods[0].type == "2OMe"
        assert mods[0].positions == [1, 4, 6, 11]
        assert mods[1].type == "2F"
        assert mods[1].positions == []

    def test_parse_provenance(self):
        """Test parsing Provenance string."""
        prov = parse_provenance("Patent:US10060921B2", "https://example.com")
        assert prov is not None
        assert prov.source_type == SourceType.PATENT
        assert prov.identifier == "US10060921B2"
        assert prov.url == "https://example.com"

    def test_parse_header_basic(self):
        """Test parsing basic FASTA header."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        record = SeqRecord(
            Seq("AUCGAUCGAUCGAUCGAUCGA"),
            id="test_001",
            description="test_001 Target=TTR; Role=guide; Confirmed=pending",
        )

        metadata = parse_header(record)
        assert metadata["id"] == "test_001"
        assert metadata["target_gene"] == "TTR"
        assert metadata["strand_role"] == StrandRole.GUIDE
        assert metadata["confirmation_status"] == ConfirmationStatus.PENDING

    def test_parse_header_with_modifications(self):
        """Test parsing header with chemical modifications."""
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        record = SeqRecord(
            Seq("AUCGAUCGAUCGAUCGAUCGA"),
            id="test_002",
            description="test_002 Target=TTR; Role=guide; ChemMods=2OMe(1,4,6)",
        )

        metadata = parse_header(record)
        assert "chem_mods" in metadata
        assert len(metadata["chem_mods"]) == 1
        assert metadata["chem_mods"][0].type == "2OMe"
        assert metadata["chem_mods"][0].positions == [1, 4, 6]


@pytest.mark.unit
@pytest.mark.local_python
class TestMetadataLoading:
    """Test metadata loading and merging."""

    def test_load_metadata_from_json(self):
        """Test loading metadata from JSON file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json_data = {
                "modifications": {
                    "test_guide": {
                        "id": "test_guide",
                        "sequence": "AUCGAUCGAUCGAUCGAUCGA",
                        "overhang": "dTdT",
                        "chem_mods": [{"type": "2OMe", "positions": [1, 4, 6]}],
                    }
                }
            }
            json.dump(json_data, f)
            json_path = f.name

        try:
            metadata = load_metadata(json_path)
            assert "test_guide" in metadata
            assert metadata["test_guide"]["overhang"] == "dTdT"
        finally:
            Path(json_path).unlink()

    def test_merge_metadata_into_fasta(self):
        """Test merging metadata into FASTA file."""
        # Create temp FASTA
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">test_guide\n")
            f.write("AUCGAUCGAUCGAUCGAUCGA\n")
            fasta_path = f.name

        # Create temp JSON
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json_data = {
                "test_guide": {
                    "id": "test_guide",
                    "sequence": "AUCGAUCGAUCGAUCGAUCGA",
                    "target_gene": "TTR",
                    "strand_role": "guide",
                    "overhang": "dTdT",
                    "chem_mods": [{"type": "2OMe", "positions": [1, 4, 6]}],
                    "confirmation_status": "pending",
                }
            }
            json.dump(json_data, f)
            json_path = f.name

        # Create temp output
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            output_path = f.name

        try:
            count = merge_metadata_into_fasta(fasta_path, json_path, output_path)
            assert count == 1

            # Verify output
            from Bio import SeqIO

            records = list(SeqIO.parse(output_path, "fasta"))
            assert len(records) == 1
            assert "Target=TTR" in records[0].description
            assert "ChemMods=2OMe(1,4,6)" in records[0].description
        finally:
            Path(fasta_path).unlink()
            Path(json_path).unlink()
            Path(output_path).unlink()
