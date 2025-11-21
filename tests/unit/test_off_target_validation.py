"""Unit tests for off-target analysis validation and DataFrame transformations.

Tests the critical data transformation logic that converts BWA-MEM2 output
to validated Pandera/Pydantic models, focusing on column renaming, type coercion,
and schema validation without requiring external tools (BWA-MEM2, databases).
"""

import pandas as pd
import pytest
from pandera.errors import SchemaErrors

from sirnaforge.models.schemas import GenomeAlignmentSchema, MiRNAAlignmentSchema


class TestMiRNADataFrameTransformation:
    """Test miRNA alignment DataFrame transformations and validation."""

    def test_bwa_output_to_mirna_schema_basic(self):
        """Test basic transformation from BWA output to MiRNAAlignmentSchema."""
        # Simulate BWA-MEM2 output with internal fields
        bwa_output = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "rname": ["hsa-miR-21-5p"],  # BWA reference name
                "coord": ["hsa-miR-21-5p:5"],  # Needs parsing
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [0],
                "mismatch_positions": [[]],  # Internal field
                "seed_mismatches": [0],
                "offtarget_score": [0.0],
            }
        )

        # Transform: parse coord, add miRNA fields, remove internal columns
        df = bwa_output.copy()
        df["coord"] = df["coord"].apply(lambda x: int(x.split(":")[-1]) if isinstance(x, str) and ":" in x else int(x))
        df["species"] = "hsa"
        df["database"] = "mirgenedb"
        df["mirna_id"] = df["rname"]

        # Remove internal columns not in schema
        df = df.drop(columns=["rname", "mismatch_positions"])

        # Validate against MiRNAAlignmentSchema
        validated_df = MiRNAAlignmentSchema.validate(df, lazy=True)

        assert len(validated_df) == 1
        assert validated_df["mirna_id"].iloc[0] == "hsa-miR-21-5p"
        assert validated_df["species"].iloc[0] == "hsa"
        assert validated_df["coord"].iloc[0] == 5
        assert "rname" not in validated_df.columns
        assert "mismatch_positions" not in validated_df.columns

    def test_missing_rname_removal_causes_validation_failure(self):
        """Test that failing to remove 'rname' causes schema validation error."""
        # Simulate the BUG: forgetting to remove 'rname'
        bwa_output = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "rname": ["hsa-miR-21-5p"],  # Should be removed!
                "coord": [5],
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [0],
                "seed_mismatches": [0],
                "offtarget_score": [0.0],
                "species": ["hsa"],
                "database": ["mirgenedb"],
                "mirna_id": ["hsa-miR-21-5p"],
            }
        )

        # Should fail because 'rname' is not in schema
        with pytest.raises(SchemaErrors) as exc_info:
            MiRNAAlignmentSchema.validate(bwa_output, lazy=True)

        assert "column 'rname' not in DataFrameSchema" in str(exc_info.value)

    def test_missing_mismatch_positions_removal_causes_validation_failure(self):
        """Test that failing to remove 'mismatch_positions' causes schema validation error."""
        # Simulate the BUG: forgetting to remove 'mismatch_positions'
        bwa_output = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "coord": [5],
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [0],
                "seed_mismatches": [0],
                "offtarget_score": [0.0],
                "species": ["hsa"],
                "database": ["mirgenedb"],
                "mirna_id": ["hsa-miR-21-5p"],
                "mismatch_positions": [[]],  # Should be removed!
            }
        )

        # Should fail because 'mismatch_positions' is not in schema
        with pytest.raises(SchemaErrors) as exc_info:
            MiRNAAlignmentSchema.validate(bwa_output, lazy=True)

        assert "column 'mismatch_positions' not in DataFrameSchema" in str(exc_info.value)

    def test_coord_string_parsing(self):
        """Test coordinate field parsing from 'miRNA_ID:position' format."""
        # BWA outputs coord as "rname:position" string
        bwa_output = pd.DataFrame(
            {
                "qname": ["sirna_001", "sirna_002"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU", "UGAGGUAGUAGGUUGUAUAGU"],
                "rname": ["hsa-miR-21-5p", "hsa-miR-155-5p"],
                "coord": ["hsa-miR-21-5p:5", "hsa-miR-155-5p:10"],  # String format
                "strand": ["+", "+"],
                "cigar": ["7M", "7M"],
                "mapq": [60, 60],
                "as_score": [14, 14],
                "nm": [0, 0],
                "seed_mismatches": [0, 0],
                "offtarget_score": [0.0, 0.0],
            }
        )

        # Transform
        df = bwa_output.copy()
        df["coord"] = df["coord"].apply(lambda x: int(x.split(":")[-1]) if isinstance(x, str) and ":" in x else int(x))
        df["species"] = "hsa"
        df["database"] = "mirgenedb"
        df["mirna_id"] = df["rname"]
        df = df.drop(columns=["rname"])

        validated_df = MiRNAAlignmentSchema.validate(df, lazy=True)

        assert validated_df["coord"].iloc[0] == 5
        assert validated_df["coord"].iloc[1] == 10
        assert validated_df["coord"].dtype == "int64"

    def test_multiple_species_batch_validation(self):
        """Test validation of multi-species miRNA results."""
        # Simulate results from multiple species
        species_dfs = []

        for species in ["hsa", "mmu", "rno"]:
            df = pd.DataFrame(
                {
                    "qname": ["sirna_001"],
                    "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                    "coord": [5],
                    "strand": ["+"],
                    "cigar": ["7M"],
                    "mapq": [60],
                    "as_score": [14],
                    "nm": [0],
                    "seed_mismatches": [0],
                    "offtarget_score": [0.0],
                    "species": [species],
                    "database": ["mirgenedb"],
                    "mirna_id": [f"{species}-miR-21-5p"],
                }
            )
            validated = MiRNAAlignmentSchema.validate(df, lazy=True)
            species_dfs.append(validated)

        # Concatenate all species results
        combined = pd.concat(species_dfs, ignore_index=True)

        assert len(combined) == 3
        assert set(combined["species"]) == {"hsa", "mmu", "rno"}
        assert all(combined["mirna_id"].str.contains("-miR-21-5p"))

    def test_perfect_match_validation(self):
        """Test that perfect matches (nm=0) have correct offtarget_score=0.0 (highest risk)."""
        df = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "coord": [5],
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [0],  # Perfect match
                "seed_mismatches": [0],
                "offtarget_score": [0.0],  # Must be 0.0 for perfect matches (highest risk)
                "species": ["hsa"],
                "database": ["mirgenedb"],
                "mirna_id": ["hsa-miR-21-5p"],
            }
        )

        # Should pass validation
        validated = MiRNAAlignmentSchema.validate(df, lazy=True)
        assert validated["offtarget_score"].iloc[0] == 0.0

    def test_perfect_match_wrong_score_fails_validation(self):
        """Test that perfect match with wrong score fails validation."""
        df = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "coord": [5],
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [0],  # Perfect match
                "seed_mismatches": [0],
                "offtarget_score": [5.0],  # WRONG! Should be 0.0
                "species": ["hsa"],
                "database": ["mirgenedb"],
                "mirna_id": ["hsa-miR-21-5p"],
            }
        )

        # Should fail validation
        with pytest.raises(SchemaErrors) as exc_info:
            MiRNAAlignmentSchema.validate(df, lazy=True)

        assert "validate_perfect_match_score" in str(exc_info.value)


class TestGenomeAlignmentTransformation:
    """Test genome/transcriptome alignment DataFrame transformations."""

    def test_bwa_output_to_genome_schema(self):
        """Test basic transformation from BWA output to GenomeAlignmentSchema."""
        bwa_output = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "rname": ["chr1"],  # Keep for genome analysis
                "coord": ["chr1:123456"],
                "strand": ["+"],
                "cigar": ["21M"],
                "mapq": [60],
                "as_score": [42],
                "nm": [0],
                "mismatch_positions": [[]],  # Internal field
                "seed_mismatches": [0],
                "offtarget_score": [0.0],
            }
        )

        # Transform: parse coord, remove internal columns (but keep rname)
        df = bwa_output.copy()
        df["coord"] = df["coord"].apply(lambda x: int(x.split(":")[-1]) if isinstance(x, str) and ":" in x else int(x))
        df = df.drop(columns=["mismatch_positions"])

        # Validate against GenomeAlignmentSchema
        validated_df = GenomeAlignmentSchema.validate(df, lazy=True)

        assert len(validated_df) == 1
        assert validated_df["rname"].iloc[0] == "chr1"  # rname kept for genome
        assert validated_df["coord"].iloc[0] == 123456
        assert "mismatch_positions" not in validated_df.columns

    def test_genome_schema_requires_rname(self):
        """Test that GenomeAlignmentSchema requires 'rname' field."""
        # Missing rname should fail
        df = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                # "rname": missing!
                "coord": [123456],
                "strand": ["+"],
                "cigar": ["21M"],
                "mapq": [60],
                "as_score": [42],
                "nm": [0],
                "seed_mismatches": [0],
                "offtarget_score": [0.0],
            }
        )

        with pytest.raises(SchemaErrors) as exc_info:
            GenomeAlignmentSchema.validate(df, lazy=True)

        assert "rname" in str(exc_info.value)


class TestEmptyDataFrameHandling:
    """Test handling of empty results (no hits found)."""

    def test_empty_mirna_dataframe_validation(self):
        """Test that empty DataFrames with correct schema pass validation."""
        # Create empty DataFrame with correct schema columns
        empty_df = pd.DataFrame(
            columns=[
                "qname",
                "qseq",
                "species",
                "database",
                "mirna_id",
                "coord",
                "strand",
                "cigar",
                "mapq",
                "as_score",
                "nm",
                "seed_mismatches",
                "offtarget_score",
            ]
        )

        # Should pass validation even though empty
        validated = MiRNAAlignmentSchema.validate(empty_df, lazy=True)
        assert len(validated) == 0
        assert set(validated.columns) == set(empty_df.columns)

    def test_empty_genome_dataframe_validation(self):
        """Test that empty genome DataFrame validates correctly."""
        empty_df = pd.DataFrame(
            columns=[
                "qname",
                "qseq",
                "rname",
                "coord",
                "strand",
                "cigar",
                "mapq",
                "as_score",
                "nm",
                "seed_mismatches",
                "offtarget_score",
            ]
        )

        validated = GenomeAlignmentSchema.validate(empty_df, lazy=True)
        assert len(validated) == 0


class TestSchemaStrictness:
    """Test that schemas correctly reject unexpected columns."""

    def test_extra_columns_rejected_by_mirna_schema(self):
        """Test that strict=True rejects DataFrames with extra columns."""
        df = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "coord": [5],
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [0],
                "seed_mismatches": [0],
                "offtarget_score": [0.0],
                "species": ["hsa"],
                "database": ["mirgenedb"],
                "mirna_id": ["hsa-miR-21-5p"],
                "extra_column": ["should_fail"],  # Extra column!
            }
        )

        # strict=True should reject this
        with pytest.raises(SchemaErrors) as exc_info:
            MiRNAAlignmentSchema.validate(df, lazy=True)

        assert "extra_column" in str(exc_info.value)


class TestSeedMismatchValidation:
    """Test seed_mismatches field validation logic."""

    def test_seed_mismatches_less_than_total_mismatches(self):
        """Test that seed_mismatches <= nm (total mismatches)."""
        df = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "coord": [5],
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [2],  # Total mismatches
                "seed_mismatches": [1],  # Must be <= nm
                "offtarget_score": [15.0],
                "species": ["hsa"],
                "database": ["mirgenedb"],
                "mirna_id": ["hsa-miR-21-5p"],
            }
        )

        # Should pass
        validated = MiRNAAlignmentSchema.validate(df, lazy=True)
        assert validated["seed_mismatches"].iloc[0] <= validated["nm"].iloc[0]

    def test_seed_mismatches_exceeds_total_fails(self):
        """Test that seed_mismatches > nm fails validation."""
        df = pd.DataFrame(
            {
                "qname": ["sirna_001"],
                "qseq": ["UGAGGUAGUAGGUUGUAUAGU"],
                "coord": [5],
                "strand": ["+"],
                "cigar": ["7M"],
                "mapq": [60],
                "as_score": [14],
                "nm": [1],  # Total mismatches
                "seed_mismatches": [2],  # INVALID: > nm!
                "offtarget_score": [15.0],
                "species": ["hsa"],
                "database": ["mirgenedb"],
                "mirna_id": ["hsa-miR-21-5p"],
            }
        )

        # Should fail
        with pytest.raises(SchemaErrors) as exc_info:
            MiRNAAlignmentSchema.validate(df, lazy=True)

        assert "validate_seed_mismatches" in str(exc_info.value)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
