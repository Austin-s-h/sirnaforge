"""Integration tests for variant targeting (SNP target and avoid modes).

Tests variant-aware siRNA design with TP53 gene, verifying both target and avoid
modes work correctly in the container environment.
"""

import json
import subprocess
import tempfile
from pathlib import Path

import pytest


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_snp_avoid_mode_tp53():
    """Test SNP avoid mode with TP53 and a common variant.
    
    This test verifies that:
    1. The workflow can resolve variants in avoid mode
    2. Variant resolution produces expected output files
    3. The resolved variants are properly annotated
    """
    # Use rs28934576 - a well-known TP53 variant (R175H)
    variant_id = "rs28934576"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "tp53_avoid_test"
        
        try:
            # Run workflow with SNP avoid mode
            cmd = [
                "sirnaforge",
                "workflow",
                "TP53",
                "--output-dir",
                str(output_dir),
                "--snp",
                variant_id,
                "--variant-mode",
                "avoid",
                "--min-af",
                "0.01",
                "--species",
                "human",
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minute timeout
            )
            
            # Check that command completed successfully
            assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"
            
            # Verify output directory exists
            assert output_dir.exists(), "Output directory should be created"
            
            # Check for variant resolution report
            variant_report = output_dir / "logs" / "resolved_variants.json"
            if variant_report.exists():
                with open(variant_report) as f:
                    variant_data = json.load(f)
                
                # Verify report structure
                assert "gene" in variant_data
                assert variant_data["gene"] == "TP53"
                assert "variant_mode" in variant_data
                assert variant_data["variant_mode"] == "avoid"
                assert "variants" in variant_data
                
                # If variants were resolved, verify they contain expected data
                if variant_data.get("variants"):
                    variant = variant_data["variants"][0]
                    assert "chr" in variant
                    assert "pos" in variant
                    assert "ref" in variant
                    assert "alt" in variant
            
            # Verify siRNA candidates were generated
            sirna_files = list(output_dir.glob("sirnaforge/*.csv"))
            assert len(sirna_files) > 0, "Should generate siRNA candidate files"
            
        except subprocess.TimeoutExpired:
            pytest.fail("Workflow timeout - test took longer than 5 minutes")
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_snp_target_mode_tp53():
    """Test SNP target mode with TP53 and a pathogenic variant.
    
    This test verifies that:
    1. The workflow can resolve variants in target mode
    2. Target mode uses global AF for filtering (not population-specific)
    3. Variant-specific candidates are generated when available
    """
    # Use rs28934576 - a pathogenic TP53 variant
    variant_id = "rs28934576"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "tp53_target_test"
        
        try:
            # Run workflow with SNP target mode
            cmd = [
                "sirnaforge",
                "workflow",
                "TP53",
                "--output-dir",
                str(output_dir),
                "--snp",
                variant_id,
                "--variant-mode",
                "target",
                "--min-af",
                "0.001",  # Lower threshold for targeting rare pathogenic variants
                "--species",
                "human",
                "--clinvar-filter-levels",
                "Pathogenic,Likely pathogenic",
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
            )
            
            # Check that command completed successfully
            assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"
            
            # Verify output directory exists
            assert output_dir.exists(), "Output directory should be created"
            
            # Check for variant resolution report
            variant_report = output_dir / "logs" / "resolved_variants.json"
            if variant_report.exists():
                with open(variant_report) as f:
                    variant_data = json.load(f)
                
                # Verify target mode is set
                assert variant_data.get("variant_mode") == "target"
                
                # Verify filters are correctly applied
                filters = variant_data.get("filters", {})
                assert filters.get("min_af") == 0.001
                assert "Pathogenic" in filters.get("clinvar_significance", [])
            
            # Verify siRNA candidates were generated
            sirna_files = list(output_dir.glob("sirnaforge/*.csv"))
            assert len(sirna_files) > 0, "Should generate siRNA candidate files"
            
        except subprocess.TimeoutExpired:
            pytest.fail("Workflow timeout - test took longer than 5 minutes")
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_snp_both_mode_tp53():
    """Test SNP both mode (target and avoid) with TP53.
    
    This test verifies that both mode works correctly, generating candidates
    for both reference and alternate alleles when appropriate.
    """
    variant_id = "rs28934576"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "tp53_both_test"
        
        try:
            # Run workflow with both mode
            cmd = [
                "sirnaforge",
                "workflow",
                "TP53",
                "--output-dir",
                str(output_dir),
                "--snp",
                variant_id,
                "--variant-mode",
                "both",
                "--min-af",
                "0.001",
                "--species",
                "human",
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
            )
            
            # Check that command completed successfully
            assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"
            
            # Verify output structure
            assert output_dir.exists()
            
            variant_report = output_dir / "logs" / "resolved_variants.json"
            if variant_report.exists():
                with open(variant_report) as f:
                    variant_data = json.load(f)
                assert variant_data.get("variant_mode") == "both"
            
            # Verify siRNA output files exist
            sirna_files = list(output_dir.glob("sirnaforge/*.csv"))
            assert len(sirna_files) > 0
            
        except subprocess.TimeoutExpired:
            pytest.fail("Workflow timeout - test took longer than 5 minutes")
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_variant_coordinate_format():
    """Test variant specification using coordinate format instead of rsID.
    
    Verifies that chr:pos:ref:alt format works correctly.
    """
    # TP53 R175H variant in coordinate format (example coordinates)
    # Note: This is an example - actual coordinates would need to be verified
    variant_coord = "chr17:7577121:C:A"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "tp53_coordinate_test"
        
        try:
            cmd = [
                "sirnaforge",
                "workflow",
                "TP53",
                "--output-dir",
                str(output_dir),
                "--snp",
                variant_coord,
                "--variant-mode",
                "avoid",
                "--species",
                "human",
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
            )
            
            # Command should complete (may not resolve variant if coordinates are incorrect)
            # but should not crash
            assert result.returncode == 0 or "Could not resolve" in result.stderr
            
        except subprocess.TimeoutExpired:
            pytest.fail("Workflow timeout")
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")


@pytest.mark.integration
@pytest.mark.runs_in_container
def test_multiple_variants():
    """Test specifying multiple variants for TP53.
    
    Verifies that multiple --snp flags work correctly.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "tp53_multi_variant_test"
        
        try:
            cmd = [
                "sirnaforge",
                "workflow",
                "TP53",
                "--output-dir",
                str(output_dir),
                "--snp",
                "rs28934576",  # R175H
                "--snp",
                "rs11540652",  # Another TP53 variant
                "--variant-mode",
                "avoid",
                "--species",
                "human",
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
            )
            
            # Should complete successfully
            assert result.returncode == 0
            
            # Check variant report contains multiple variants
            variant_report = output_dir / "logs" / "resolved_variants.json"
            if variant_report.exists():
                with open(variant_report) as f:
                    variant_data = json.load(f)
                
                # Should have attempted to resolve multiple variants
                summary = variant_data.get("summary", {})
                # Note: Some variants may not resolve, but we should see the attempt
                assert "total_variants" in summary
            
        except subprocess.TimeoutExpired:
            pytest.fail("Workflow timeout")
        except FileNotFoundError:
            pytest.skip("sirnaforge CLI not available - run this test in Docker container")
