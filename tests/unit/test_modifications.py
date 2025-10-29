"""Tests for chemical modifications metadata system."""

import json
import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydantic import ValidationError

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
    save_metadata_json,
)


@pytest.mark.unit
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
class TestParsingFunctions:
    """Test parsing helper functions."""

    def test_parse_chem_mods(self):
        """Test parsing ChemMods string."""
        mods = parse_chem_mods("2OMe(1,4,6,11)+2F()")
        assert len(mods) == 2
        assert mods[0].type == "2OMe"
        assert mods[0].positions == [1, 4, 6, 11]
        assert mods[1].type == "2F"
        assert mods[1].positions == []

    def test_parse_chem_mods_pipe_compatibility(self):
        """Ensure legacy pipe delimiter is still supported."""
        mods = parse_chem_mods("2OMe(1,4,6)|2F()")
        assert len(mods) == 2
        assert {mod.type for mod in mods} == {"2OMe", "2F"}

    def test_parse_provenance(self):
        """Test parsing Provenance string."""
        prov = parse_provenance("Patent:US10060921B2", "https://example.com")
        assert prov is not None
        assert prov.source_type == SourceType.PATENT
        assert prov.identifier == "US10060921B2"
        assert prov.url == "https://example.com"

    def test_parse_header_basic(self):
        """Test parsing basic FASTA header."""
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
            records = list(SeqIO.parse(output_path, "fasta"))
            assert len(records) == 1
            assert "Target=TTR" in records[0].description
            assert "Role=guide" in records[0].description
            assert "ChemMods=2OMe(1,4,6)" in records[0].description
            assert "|" not in records[0].description
        finally:
            Path(fasta_path).unlink()
            Path(json_path).unlink()
            Path(output_path).unlink()


@pytest.mark.unit
class TestPydanticValidation:
    """Test Pydantic validation features."""

    def test_invalid_modification_type(self):
        """Test that empty modification type raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            ChemicalModification(type="", positions=[1, 2, 3])
        assert "Modification type cannot be empty" in str(exc_info.value)

    def test_invalid_positions(self):
        """Test that positions must be >= 1."""
        with pytest.raises(ValidationError) as exc_info:
            ChemicalModification(type="2OMe", positions=[1, 0, -1])
        assert "All positions must be >= 1" in str(exc_info.value)

    def test_duplicate_positions_removed(self):
        """Test that duplicate positions are removed and sorted."""
        mod = ChemicalModification(type="2OMe", positions=[5, 1, 3, 1, 5, 2])
        assert mod.positions == [1, 2, 3, 5]

    def test_invalid_sequence(self):
        """Test that invalid sequence characters raise ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            StrandMetadata(
                id="test",
                sequence="AUCGXYZ",  # Invalid characters
            )
        assert "invalid characters" in str(exc_info.value).lower()

    def test_empty_sequence(self):
        """Test that empty sequence raises ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            StrandMetadata(id="test", sequence="")
        assert "cannot be empty" in str(exc_info.value).lower()

    def test_modification_position_exceeds_sequence_length(self):
        """Test that modification positions beyond sequence length raise error."""
        with pytest.raises(ValidationError) as exc_info:
            StrandMetadata(
                id="test",
                sequence="AUCG",  # Length 4
                chem_mods=[ChemicalModification(type="2OMe", positions=[1, 2, 5])],  # Position 5 > 4
            )
        assert "position" in str(exc_info.value).lower()
        assert "sequence length" in str(exc_info.value).lower()

    def test_sequence_case_normalization(self):
        """Test that sequence is normalized to uppercase."""
        metadata = StrandMetadata(id="test", sequence="aucgaucg")
        assert metadata.sequence == "AUCGAUCG"

    def test_modification_type_whitespace_stripped(self):
        """Test that modification type has whitespace stripped."""
        mod = ChemicalModification(type="  2OMe  ", positions=[1, 2])
        assert mod.type == "2OMe"


@pytest.mark.unit
class TestPydanticJSONSerialization:
    """Test Pydantic JSON serialization and deserialization."""

    def test_load_metadata_with_validation(self):
        """Test that load_metadata validates data using Pydantic."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            # Valid data
            json_data = {
                "test_guide": {
                    "id": "test_guide",
                    "sequence": "AUCGAUCG",
                    "chem_mods": [{"type": "2OMe", "positions": [1, 3, 5]}],
                }
            }
            json.dump(json_data, f)
            json_path = f.name

        try:
            metadata = load_metadata(json_path)
            assert isinstance(metadata["test_guide"], StrandMetadata)
            assert metadata["test_guide"].sequence == "AUCGAUCG"
            assert len(metadata["test_guide"].chem_mods) == 1
        finally:
            Path(json_path).unlink()

    def test_load_metadata_invalid_data(self):
        """Test that load_metadata raises ValidationError for invalid data."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            # Invalid data - modification position exceeds sequence length
            json_data = {
                "test_guide": {
                    "id": "test_guide",
                    "sequence": "AUCG",  # Length 4
                    "chem_mods": [{"type": "2OMe", "positions": [1, 10]}],  # Position 10 > 4
                }
            }
            json.dump(json_data, f)
            json_path = f.name

        try:
            with pytest.raises(ValidationError):
                load_metadata(json_path)
        finally:
            Path(json_path).unlink()

    def test_save_and_load_roundtrip(self):
        """Test that save and load operations preserve data."""
        original_metadata = {
            "test_guide": StrandMetadata(
                id="test_guide",
                sequence="AUCGAUCGAUCG",
                overhang="dTdT",
                chem_mods=[
                    ChemicalModification(type="2OMe", positions=[1, 3, 5]),
                    ChemicalModification(type="2F", positions=[2, 4]),
                ],
                provenance=Provenance(source_type=SourceType.PATENT, identifier="TEST123", url="https://example.com"),
                confirmation_status=ConfirmationStatus.CONFIRMED,
                notes="Test data",
            )
        }

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json_path = Path(f.name)

        try:
            # Save
            save_metadata_json(original_metadata, json_path)

            # Load
            loaded_metadata = load_metadata(json_path)

            # Verify
            assert "test_guide" in loaded_metadata
            loaded = loaded_metadata["test_guide"]
            original = original_metadata["test_guide"]

            assert loaded.id == original.id
            assert loaded.sequence == original.sequence
            assert loaded.overhang == original.overhang
            assert len(loaded.chem_mods) == len(original.chem_mods)
            assert loaded.chem_mods[0].type == original.chem_mods[0].type
            assert loaded.chem_mods[0].positions == original.chem_mods[0].positions
            assert loaded.provenance.source_type == original.provenance.source_type
            assert loaded.confirmation_status == original.confirmation_status
            assert loaded.notes == original.notes
        finally:
            json_path.unlink()

    def test_model_dump_excludes_none(self):
        """Test that model_dump with exclude_none=True works correctly."""
        metadata = StrandMetadata(
            id="test",
            sequence="AUCG",
            # overhang, notes, provenance are None
        )

        dumped = metadata.model_dump(mode="json", exclude_none=True)
        assert "overhang" not in dumped
        assert "notes" not in dumped
        assert "provenance" not in dumped
        assert "id" in dumped
        assert "sequence" in dumped


@pytest.mark.unit
class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_empty_chem_mods_list(self):
        """Test strand with no chemical modifications."""
        metadata = StrandMetadata(id="test", sequence="AUCG", chem_mods=[])
        assert metadata.chem_mods == []
        header = metadata.to_fasta_header()
        assert "ChemMods" not in header

    def test_modification_with_no_positions(self):
        """Test modification annotation without specific positions."""
        mod = ChemicalModification(type="2OMe", positions=[])
        assert mod.to_header_string() == "2OMe()"

    def test_very_long_sequence(self):
        """Test handling of long sequences."""
        long_seq = "AUCG" * 100  # 400 nucleotides
        metadata = StrandMetadata(
            id="long_test",
            sequence=long_seq,
            chem_mods=[ChemicalModification(type="2OMe", positions=[1, 200, 400])],
        )
        assert len(metadata.sequence) == 400
        assert max(metadata.chem_mods[0].positions) == 400

    def test_many_modifications(self):
        """Test strand with multiple modification types."""
        mods = [
            ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9]),
            ChemicalModification(type="2F", positions=[2, 4, 6, 8]),
            ChemicalModification(type="PS", positions=[1, 2, 20, 21]),
            ChemicalModification(type="LNA", positions=[10, 11, 12]),
        ]
        metadata = StrandMetadata(
            id="multi_mod",
            sequence="AUCGAUCGAUCGAUCGAUCGA",  # 21 nt
            chem_mods=mods,
        )
        assert len(metadata.chem_mods) == 4
        header = metadata.to_fasta_header()
        assert "2OMe" in header
        assert "2F" in header
        assert "PS" in header
        assert "LNA" in header

    def test_load_empty_json_file(self):
        """Test loading from empty JSON file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({}, f)
            json_path = f.name

        try:
            metadata = load_metadata(json_path)
            assert metadata == {}
        finally:
            Path(json_path).unlink()

    def test_load_nonexistent_file(self):
        """Test loading from non-existent file."""
        metadata = load_metadata("/tmp/nonexistent_file_12345.json")
        assert metadata == {}

    def test_unicode_in_notes(self):
        """Test handling of Unicode characters in notes."""
        metadata = StrandMetadata(
            id="test",
            sequence="AUCG",
            notes="Special characters: α β γ δ ε — • ™",
        )
        assert "α" in metadata.notes
        # Test serialization
        dumped = metadata.model_dump(mode="json")
        assert "α" in dumped["notes"]
