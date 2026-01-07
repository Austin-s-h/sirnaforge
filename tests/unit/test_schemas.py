"""Tests for the modernized pandera schemas."""

import pandas as pd
import pytest

from sirnaforge.models.schemas import OffTargetHitsSchema, ORFValidationSchema, SiRNACandidateSchema


class TestSiRNACandidateSchema:
    """Test the modernized siRNA candidate schema."""

    def test_valid_data_passes_validation(self):
        """Test that valid siRNA candidate data passes schema validation."""
        test_data = {
            "id": ["test_1", "test_2"],
            "transcript_id": ["ENST00000000001", "ENST00000000002"],
            "position": [100, 200],
            "guide_sequence": ["TTTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAAA"],
            "passenger_sequence": ["AAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTT"],
            "gc_content": [0.0, 0.0],
            "asymmetry_score": [0.5, 0.7],
            "paired_fraction": [0.3, 0.4],
            # New thermodynamic/structure fields (nullable)
            "structure": [None, None],
            "mfe": [None, None],
            "duplex_stability_dg": [None, None],
            "duplex_stability_score": [None, None],
            "dg_5p": [None, None],
            "dg_3p": [None, None],
            "delta_dg_end": [None, None],
            "melting_temp_c": [None, None],
            "off_target_count": [2, 1],
            "transcript_hit_count": [1, 1],
            "transcript_hit_fraction": [1.0, 1.0],
            "composite_score": [75.5, 82.3],
            "passes_filters": [True, True],
        }

        df = pd.DataFrame(test_data)
        # Explicitly set dtypes for nullable float columns to avoid object inference
        df = df.astype(
            {
                "mfe": "float64",
                "duplex_stability_dg": "float64",
                "duplex_stability_score": "float64",
                "dg_5p": "float64",
                "dg_3p": "float64",
                "delta_dg_end": "float64",
                "melting_temp_c": "float64",
            }
        )
        result = SiRNACandidateSchema.validate(df)
        extra_columns = {
            "guide_overhang",
            "guide_modifications",
            "passenger_overhang",
            "passenger_modifications",
            "guide_pos1_base",
            "pos1_pairing_state",
            "seed_class",
            "supp_13_16_score",
            "seed_7mer_hits",
            "seed_8mer_hits",
            "seed_hits_weighted",
            "off_target_seed_risk_class",
            "variant_mode",
            "allele_specific",
            "targeted_alleles",
            "overlapped_variants",
        }

        expected_cols = len(test_data.keys()) + len(extra_columns)
        assert result.shape == (2, expected_cols)
        # Order-agnostic column check - now includes modification and miRNA columns
        expected_columns = set(test_data.keys()) | extra_columns
        assert set(result.columns) == expected_columns

    def test_passes_filters_allows_new_reason_strings(self):
        """Ensure passes_filters accepts detailed reason strings from filters."""
        test_data = {
            "id": ["test_1", "test_2"],
            "transcript_id": ["ENST00000000001", "ENST00000000002"],
            "position": [100, 200],
            "guide_sequence": ["TTTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAAAAAAAAAAAAA"],
            "passenger_sequence": ["AAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTT"],
            "gc_content": [0.0, 0.0],
            "asymmetry_score": [0.5, 0.7],
            "paired_fraction": [0.3, 0.4],
            "structure": [None, None],
            "mfe": [None, None],
            "duplex_stability_dg": [None, None],
            "duplex_stability_score": [None, None],
            "dg_5p": [None, None],
            "dg_3p": [None, None],
            "delta_dg_end": [None, None],
            "melting_temp_c": [None, None],
            "off_target_count": [2, 1],
            "transcript_hit_count": [1, 1],
            "transcript_hit_fraction": [1.0, 1.0],
            "composite_score": [75.5, 82.3],
            "passes_filters": [
                "MIRNA_PERFECT_SEED (1 hits)",
                "TRANSCRIPTOME_PERFECT_MATCH (2 hits)",
            ],
        }

        df = pd.DataFrame(test_data)
        df = df.astype(
            {
                "mfe": "float64",
                "duplex_stability_dg": "float64",
                "duplex_stability_score": "float64",
                "dg_5p": "float64",
                "dg_3p": "float64",
                "delta_dg_end": "float64",
                "melting_temp_c": "float64",
            }
        )

        result = SiRNACandidateSchema.validate(df)
        assert result.loc[0, "passes_filters"].startswith("MIRNA_PERFECT_SEED")
        assert result.loc[1, "passes_filters"].startswith("TRANSCRIPTOME_PERFECT_MATCH")

    def test_invalid_sequence_length_fails(self):
        """Test that invalid sequence lengths fail validation."""
        test_data = {
            "id": ["test_1"],
            "transcript_id": ["ENST00000000001"],
            "position": [100],
            "guide_sequence": ["TTTT"],  # Too short
            "passenger_sequence": ["AAAA"],  # Too short
            "gc_content": [0.0],
            "asymmetry_score": [0.5],
            "paired_fraction": [0.3],
            # Include new nullable fields
            "structure": [None],
            "mfe": [None],
            "duplex_stability_dg": [None],
            "duplex_stability_score": [None],
            "dg_5p": [None],
            "dg_3p": [None],
            "delta_dg_end": [None],
            "melting_temp_c": [None],
            "off_target_count": [2],
            "transcript_hit_count": [1],
            "transcript_hit_fraction": [1.0],
            "composite_score": [75.5],
            "passes_filters": [True],
        }

        df = pd.DataFrame(test_data)
        # Explicitly set dtypes for nullable float columns
        df = df.astype(
            {
                "mfe": "float64",
                "duplex_stability_dg": "float64",
                "duplex_stability_score": "float64",
                "dg_5p": "float64",
                "dg_3p": "float64",
                "delta_dg_end": "float64",
                "melting_temp_c": "float64",
            }
        )
        with pytest.raises((Exception, ValueError)):  # Should fail sequence length validation
            SiRNACandidateSchema.validate(df)

    def test_invalid_nucleotides_fails(self):
        """Test that invalid nucleotide characters fail validation."""
        test_data = {
            "id": ["test_1"],
            "transcript_id": ["ENST00000000001"],
            "position": [100],
            "guide_sequence": ["NNNNNNNNNNNNNNNNNNNNX"],  # Contains invalid 'X'
            "passenger_sequence": ["AAAAAAAAAAAAAAAAAAAAA"],
            "gc_content": [0.0],
            "asymmetry_score": [0.5],
            "paired_fraction": [0.3],
            # Include new nullable fields
            "structure": [None],
            "mfe": [None],
            "duplex_stability_dg": [None],
            "duplex_stability_score": [None],
            "dg_5p": [None],
            "dg_3p": [None],
            "delta_dg_end": [None],
            "melting_temp_c": [None],
            "off_target_count": [2],
            "transcript_hit_count": [1],
            "transcript_hit_fraction": [1.0],
            "composite_score": [75.5],
            "passes_filters": [True],
        }

        df = pd.DataFrame(test_data)
        # Explicitly set dtypes for nullable float columns
        df = df.astype(
            {
                "mfe": "float64",
                "duplex_stability_dg": "float64",
                "duplex_stability_score": "float64",
                "dg_5p": "float64",
                "dg_3p": "float64",
                "delta_dg_end": "float64",
                "melting_temp_c": "float64",
            }
        )
        with pytest.raises((Exception, ValueError)):  # Should fail nucleotide validation
            SiRNACandidateSchema.validate(df)


class TestORFValidationSchema:
    """Test the modernized ORF validation schema."""

    def test_valid_orf_data_passes(self):
        """Test that valid ORF analysis data passes validation."""
        test_data = {
            "transcript_id": ["ENST00000000001"],
            "sequence_length": [1000],
            "gc_content": [50.0],
            "orfs_found": [1],
            "has_valid_orf": [True],
            "longest_orf_start": [10],
            "longest_orf_end": [500],
            "longest_orf_length": [490],
            "longest_orf_frame": [0],
            "start_codon": ["ATG"],
            "stop_codon": ["TAA"],
            "orf_gc_content": [45.0],
        }

        df = pd.DataFrame(test_data)
        # Convert to object types for nullable fields as schema expects
        df = df.astype(
            {
                "longest_orf_start": object,
                "longest_orf_end": object,
                "longest_orf_length": object,
                "longest_orf_frame": object,
                "start_codon": object,
                "stop_codon": object,
                "orf_gc_content": object,
            }
        )
        result = ORFValidationSchema.validate(df)

        assert result.shape == (1, 12)
        assert result["transcript_id"].iloc[0] == "ENST00000000001"

    def test_empty_orf_data_passes(self):
        """Test that empty ORF data (no ORFs found) passes validation."""
        test_data = {
            "transcript_id": ["ENST00000000001"],
            "sequence_length": [500],
            "gc_content": [30.0],
            "orfs_found": [0],
            "has_valid_orf": [False],
            "longest_orf_start": [None],
            "longest_orf_end": [None],
            "longest_orf_length": [None],
            "longest_orf_frame": [None],
            "start_codon": [None],
            "stop_codon": [None],
            "orf_gc_content": [None],
        }

        df = pd.DataFrame(test_data)
        result = ORFValidationSchema.validate(df)

        assert result.shape == (1, 12)
        assert not result["has_valid_orf"].iloc[0]


class TestOffTargetHitsSchema:
    """Test the modernized off-target hits schema."""

    def test_valid_offtarget_data_passes(self):
        """Test that valid off-target data passes validation."""
        test_data = {
            "qname": ["guide_1"],
            "target_id": ["target_1"],
            "species": ["homo_sapiens"],
            "chromosome": ["chr1"],
            "position": [1000],
            "strand": ["+"],
            "mismatches": [2],
            "alignment_score": [15.5],
            "offtarget_score": [3.2],
            "target_sequence": ["ATCGATCGATCGATCGATCGA"],
        }

        df = pd.DataFrame(test_data)
        # Convert to object types for nullable fields
        df = df.astype(
            {
                "target_id": object,
                "species": object,
                "chromosome": object,
                "position": object,
                "strand": object,
                "mismatches": object,
                "alignment_score": object,
                "offtarget_score": object,
                "target_sequence": object,
            }
        )
        result = OffTargetHitsSchema.validate(df)

        assert result.shape == (1, 10)
        assert result["qname"].iloc[0] == "guide_1"

    def test_nullable_fields_work(self):
        """Test that nullable fields work correctly."""
        test_data = {
            "qname": ["guide_1"],
            "target_id": [None],
            "species": [None],
            "chromosome": [None],
            "position": [None],
            "strand": [None],
            "mismatches": [None],
            "alignment_score": [None],
            "offtarget_score": [None],
            "target_sequence": [None],
        }

        df = pd.DataFrame(test_data)
        result = OffTargetHitsSchema.validate(df)

        assert result.shape == (1, 10)
        assert pd.isna(result["target_id"].iloc[0])
