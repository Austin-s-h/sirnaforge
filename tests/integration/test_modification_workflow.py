"""Integration tests for chemical modification annotation workflow.

These tests validate the end-to-end workflow of:
1. Designing siRNAs
2. Annotating with chemical modifications
3. Exporting with metadata
4. Roundtrip validation
"""

import json
import tempfile
from pathlib import Path

import pytest
from Bio import SeqIO

from sirnaforge.models.modifications import (
    ChemicalModification,
    ConfirmationStatus,
    Provenance,
    SourceType,
    StrandMetadata,
    StrandRole,
)
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.modifications import (
    load_metadata,
    merge_metadata_into_fasta,
    parse_header,
    save_metadata_json,
)


@pytest.mark.integration
class TestModificationWorkflow:
    """Test complete workflow of modification annotation."""

    @pytest.mark.integration
    def test_create_annotate_export_roundtrip(self):
        """Test complete roundtrip: create -> annotate -> export -> load -> validate."""
        # Step 1: Create siRNA candidate (simulating design output)
        candidate = SiRNACandidate(
            id="test_sirna_001",
            transcript_id="ENST00000123456",
            position=100,
            guide_sequence="AUCGAUCGAUCGAUCGAUCGA",
            passenger_sequence="UCGAUCGAUCGAUCGAUCGAU",
            gc_content=52.4,
            length=21,
            asymmetry_score=0.75,
            composite_score=85.2,
        )

        # Step 2: Create modification metadata
        guide_metadata = StrandMetadata(
            id=f"{candidate.id}_guide",
            sequence=candidate.guide_sequence,
            overhang="dTdT",
            chem_mods=[ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9, 11, 13, 15, 17, 19])],
            provenance=Provenance(
                source_type=SourceType.DESIGNED, identifier="sirnaforge_test_v1.0", url="https://example.com/test"
            ),
            confirmation_status=ConfirmationStatus.PENDING,
            notes="Test candidate with standard 2'-O-methyl pattern",
        )

        metadata_dict = {f"{candidate.id}_guide": guide_metadata}

        # Step 3: Save metadata to JSON
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = Path(tmpdir) / "modifications.json"
            save_metadata_json(metadata_dict, json_path)

            # Step 4: Load and validate
            loaded_metadata = load_metadata(json_path)
            assert f"{candidate.id}_guide" in loaded_metadata
            loaded = loaded_metadata[f"{candidate.id}_guide"]

            # Validate all fields preserved
            assert loaded.id == guide_metadata.id
            assert loaded.sequence == guide_metadata.sequence
            assert loaded.overhang == guide_metadata.overhang
            assert len(loaded.chem_mods) == len(guide_metadata.chem_mods)
            assert loaded.chem_mods[0].type == guide_metadata.chem_mods[0].type
            assert loaded.chem_mods[0].positions == guide_metadata.chem_mods[0].positions
            assert loaded.provenance.source_type == guide_metadata.provenance.source_type
            assert loaded.provenance.identifier == guide_metadata.provenance.identifier
            assert loaded.confirmation_status == guide_metadata.confirmation_status
            assert loaded.notes == guide_metadata.notes

    @pytest.mark.integration
    def test_fasta_annotation_workflow(self):
        """Test annotating FASTA file with modifications from JSON."""
        # Create test FASTA
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "candidates.fasta"
            json_path = Path(tmpdir) / "modifications.json"
            output_path = Path(tmpdir) / "annotated.fasta"

            # Write simple FASTA
            with fasta_path.open("w") as f:
                f.write(">sirna_001\n")
                f.write("AUCGAUCGAUCGAUCGAUCGA\n")
                f.write(">sirna_002\n")
                f.write("GCGCGCGCGCGCGCGCGCGCG\n")

            # Create metadata for both sequences
            # Note: StrandMetadata doesn't have target_gene or strand_role fields
            # These are stored in the raw JSON and passed to to_fasta_header
            metadata = {
                "sirna_001": StrandMetadata(
                    id="sirna_001",
                    sequence="AUCGAUCGAUCGAUCGAUCGA",
                    overhang="dTdT",
                    chem_mods=[ChemicalModification(type="2OMe", positions=[1, 3, 5])],
                    confirmation_status=ConfirmationStatus.PENDING,
                ),
                "sirna_002": StrandMetadata(
                    id="sirna_002",
                    sequence="GCGCGCGCGCGCGCGCGCGCG",
                    overhang="UU",
                    chem_mods=[
                        ChemicalModification(type="2OMe", positions=[2, 4, 6]),
                        ChemicalModification(type="PS", positions=[]),
                    ],
                    confirmation_status=ConfirmationStatus.CONFIRMED,
                ),
            }

            # Save metadata with extra fields in raw JSON
            output_data = {
                "sirna_001": {
                    **metadata["sirna_001"].model_dump(mode="json", exclude_none=True),
                    "target_gene": "TP53",
                    "strand_role": "guide",
                },
                "sirna_002": {
                    **metadata["sirna_002"].model_dump(mode="json", exclude_none=True),
                    "target_gene": "BRCA1",
                    "strand_role": "guide",
                },
            }
            with json_path.open("w") as f:
                json.dump(output_data, f)

            # Merge metadata into FASTA
            count = merge_metadata_into_fasta(fasta_path, json_path, output_path)
            assert count == 2

            # Validate output
            records = list(SeqIO.parse(output_path, "fasta"))
            assert len(records) == 2

            # Check first record
            rec1 = records[0]
            assert "Target=TP53" in rec1.description
            assert "Role=guide" in rec1.description
            assert "Confirmed=pending" in rec1.description
            assert "Overhang=dTdT" in rec1.description
            assert "ChemMods=2OMe(1,3,5)" in rec1.description

            # Check second record
            rec2 = records[1]
            assert "Target=BRCA1" in rec2.description
            assert "Role=guide" in rec2.description
            assert "Confirmed=confirmed" in rec2.description
            assert "Overhang=UU" in rec2.description
            assert "2OMe(2,4,6)" in rec2.description
            assert "PS()" in rec2.description

    @pytest.mark.integration
    def test_parse_annotated_fasta(self):
        """Test parsing FASTA with modification annotations."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "annotated.fasta"

            # Create annotated FASTA
            with fasta_path.open("w") as f:
                f.write(
                    ">sirna_001 Target=TTR; Role=guide; Confirmed=confirmed; "
                    "Overhang=dTdT; ChemMods=2OMe(1,4,6,11)+2F(2,5); "
                    "Provenance=Patent:US10060921B2; URL=https://example.com\n"
                )
                f.write("AUCGAUCGAUCGAUCGAUCGA\n")

            # Parse and validate
            records = list(SeqIO.parse(fasta_path, "fasta"))
            assert len(records) == 1

            metadata = parse_header(records[0])
            assert metadata["id"] == "sirna_001"
            assert metadata["target_gene"] == "TTR"
            assert metadata["strand_role"] == StrandRole.GUIDE
            assert metadata["confirmation_status"] == ConfirmationStatus.CONFIRMED
            assert metadata["overhang"] == "dTdT"

            # Check modifications
            assert "chem_mods" in metadata
            mods = metadata["chem_mods"]
            assert len(mods) == 2
            assert mods[0].type == "2OMe"
            assert mods[0].positions == [1, 4, 6, 11]
            assert mods[1].type == "2F"
            assert mods[1].positions == [2, 5]

            # Check provenance
            assert metadata["provenance"] is not None
            prov = metadata["provenance"]
            assert prov.source_type == SourceType.PATENT
            assert prov.identifier == "US10060921B2"
            assert prov.url == "https://example.com"


@pytest.mark.integration
class TestModificationPatterns:
    """Test modification pattern loading and application."""

    def test_load_pattern_file(self):
        """Test loading modification pattern from JSON."""
        pattern_data = {
            "pattern_name": "test_pattern",
            "guide_modifications": {"2OMe": {"positions": [1, 3, 5, 7], "strategy": "custom"}},
            "passenger_modifications": {"2OMe": {"positions": [2, 4, 6, 8], "strategy": "custom"}},
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            pattern_path = Path(tmpdir) / "test_pattern.json"
            with pattern_path.open("w") as f:
                json.dump(pattern_data, f)

            # Load and validate
            with pattern_path.open() as f:
                loaded = json.load(f)
                assert loaded["pattern_name"] == "test_pattern"
                assert "guide_modifications" in loaded
                assert "2OMe" in loaded["guide_modifications"]

    def test_apply_pattern_to_sequence(self):
        """Test applying a modification pattern to a sequence."""
        sequence = "AUCGAUCGAUCGAUCGAUCGA"
        pattern_positions = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]

        # Create metadata with pattern
        metadata = StrandMetadata(
            id="test_seq",
            sequence=sequence,
            chem_mods=[ChemicalModification(type="2OMe", positions=pattern_positions)],
        )

        # Validate positions are within sequence
        assert all(pos <= len(sequence) for pos in pattern_positions)
        assert max(metadata.chem_mods[0].positions) <= len(metadata.sequence)

    def test_multiple_modification_types(self):
        """Test sequence with multiple modification types."""
        sequence = "AUCGAUCGAUCGAUCGAUCGA"
        metadata = StrandMetadata(
            id="multi_mod",
            sequence=sequence,
            chem_mods=[
                ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9]),
                ChemicalModification(type="2F", positions=[2, 4, 6, 8]),
                ChemicalModification(type="PS", positions=[1, 2, 20, 21]),
            ],
        )

        # Validate all modifications
        assert len(metadata.chem_mods) == 3
        assert {mod.type for mod in metadata.chem_mods} == {"2OMe", "2F", "PS"}

        # Check FASTA header includes all modifications
        header = metadata.to_fasta_header(target_gene="TEST", strand_role=StrandRole.GUIDE)
        assert "2OMe" in header
        assert "2F" in header
        assert "PS" in header


@pytest.mark.integration
class TestCandidateMetadataIntegration:
    """Test integration of modifications with SiRNACandidate."""

    def test_candidate_with_metadata(self):
        """Test creating candidate with modification metadata."""
        guide_metadata = StrandMetadata(
            id="test_guide",
            sequence="AUCGAUCGAUCGAUCGAUCGA",
            overhang="dTdT",
            chem_mods=[ChemicalModification(type="2OMe", positions=[1, 3, 5])],
        )

        passenger_metadata = StrandMetadata(
            id="test_passenger",
            sequence="UCGAUCGAUCGAUCGAUCGAU",
            overhang="dTdT",
            chem_mods=[ChemicalModification(type="2OMe", positions=[2, 4, 6])],
        )

        candidate = SiRNACandidate(
            id="test_001",
            transcript_id="ENST00000123456",
            position=100,
            guide_sequence=guide_metadata.sequence,
            passenger_sequence=passenger_metadata.sequence,
            gc_content=52.4,
            length=21,
            asymmetry_score=0.75,
            composite_score=85.2,
            guide_metadata=guide_metadata,
            passenger_metadata=passenger_metadata,
        )

        assert candidate.guide_metadata is not None
        assert candidate.passenger_metadata is not None
        assert len(candidate.guide_metadata.chem_mods) == 1
        assert len(candidate.passenger_metadata.chem_mods) == 1

    def test_candidate_fasta_with_metadata(self):
        """Test FASTA generation with metadata."""
        guide_metadata = StrandMetadata(
            id="test_guide",
            sequence="AUCGAUCGAUCGAUCGAUCGA",
            overhang="dTdT",
            chem_mods=[ChemicalModification(type="2OMe", positions=[1, 3, 5, 7, 9])],
            provenance=Provenance(source_type=SourceType.DESIGNED, identifier="sirnaforge_test"),
        )

        candidate = SiRNACandidate(
            id="test_001",
            transcript_id="ENST00000123456",
            position=100,
            guide_sequence=guide_metadata.sequence,
            passenger_sequence="UCGAUCGAUCGAUCGAUCGAU",
            gc_content=52.4,
            length=21,
            asymmetry_score=0.75,
            composite_score=85.2,
            guide_metadata=guide_metadata,
        )

        # Test without metadata
        fasta_simple = candidate.to_fasta(include_metadata=False)
        assert "ChemMods" not in fasta_simple
        assert "test_001" in fasta_simple

        # Test with metadata
        fasta_annotated = candidate.to_fasta(include_metadata=True)
        assert "ChemMods=2OMe(1,3,5,7,9)" in fasta_annotated
        assert "Overhang=dTdT" in fasta_annotated
        assert "Provenance=Designed:sirnaforge_test" in fasta_annotated


@pytest.mark.integration
class TestRealWorldExamples:
    """Test with real-world modification patterns."""

    def test_standard_2ome_pattern(self):
        """Test standard 2'-O-methyl alternating pattern."""
        sequence = "AUCGAUCGAUCGAUCGAUCGA"  # 21nt
        alternating_positions = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]

        metadata = StrandMetadata(
            id="standard_pattern",
            sequence=sequence,
            overhang="dTdT",
            chem_mods=[ChemicalModification(type="2OMe", positions=alternating_positions)],
            notes="Standard alternating 2'-O-methyl pattern for balanced stability",
        )

        assert len(metadata.chem_mods[0].positions) == 10
        assert all(pos % 2 == 1 for pos in metadata.chem_mods[0].positions)  # All odd positions

    def test_fda_approved_onpattro_pattern(self):
        """Test Patisiran (Onpattro) modification pattern."""
        # From patent US10060921B2
        guide_sequence = "AUGGAAUACUCUUGGUUAC"  # 19nt
        guide_mods = [1, 4, 6, 11, 13, 16, 19]  # Actual Onpattro pattern

        metadata = StrandMetadata(
            id="patisiran_ttr_guide",
            sequence=guide_sequence,
            overhang="dTdT",
            chem_mods=[ChemicalModification(type="2OMe", positions=guide_mods)],
            provenance=Provenance(
                source_type=SourceType.PATENT,
                identifier="US10060921B2",
                url="https://patents.google.com/patent/US10060921B2",
            ),
            confirmation_status=ConfirmationStatus.CONFIRMED,
            notes="FDA-approved Patisiran (Onpattro) for hATTR amyloidosis",
        )

        assert metadata.sequence == guide_sequence
        assert len(metadata.sequence) == 19
        assert max(metadata.chem_mods[0].positions) == 19
        assert metadata.confirmation_status == ConfirmationStatus.CONFIRMED

    def test_minimal_terminal_pattern(self):
        """Test minimal terminal modification pattern."""
        sequence = "AUCGAUCGAUCGAUCGAUCGA"
        terminal_3prime = [19, 20, 21]  # Last 3 positions

        metadata = StrandMetadata(
            id="minimal_pattern",
            sequence=sequence,
            overhang="dTdT",
            chem_mods=[ChemicalModification(type="2OMe", positions=terminal_3prime)],
            notes="Minimal modification for cost-sensitive applications",
        )

        assert metadata.chem_mods[0].positions == [19, 20, 21]
        assert all(pos > 15 for pos in metadata.chem_mods[0].positions)  # All at 3' end


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
