#!/usr/bin/env python3
"""
Complete siRNAforge Workflow Example

This example demonstrates the integrated workflow that goes from gene query
to off-target analysis using both Python design algorithms and Nextflow pipeline.
"""

import asyncio
import subprocess
import sys
import traceback
from pathlib import Path

from sirnaforge.workflow import run_sirna_workflow

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


async def run_complete_example() -> None:
    """Run complete siRNA design workflow for TP53."""

    print("🧬 Starting Complete siRNAforge Workflow Example")
    print("=" * 60)

    # Configuration
    gene_query = "TP53"
    output_dir = Path(__file__).parent / "workflow_output"

    print(f"Gene Target: {gene_query}")
    print(f"Output Directory: {output_dir}")
    print()

    try:
        # Run the complete workflow
        results = await run_sirna_workflow(
            gene_query=gene_query,
            output_dir=str(output_dir),
            database="ensembl",
            top_n_candidates=20,
            genome_species=["human", "rat", "rhesus"],
            gc_min=30.0,
            gc_max=52.0,
            sirna_length=21,
        )

        print("✅ Workflow completed successfully!")
        print()

        # Display summary
        print("📊 Results Summary:")
        print("-" * 30)

        workflow_config = results.get("workflow_config", {})
        transcript_summary = results.get("transcript_summary", {})
        design_summary = results.get("design_summary", {})
        offtarget_summary = results.get("offtarget_summary", {})

        print(f"Processing time: {workflow_config.get('processing_time', 0):.2f} seconds")
        print(f"Transcripts retrieved: {transcript_summary.get('total_transcripts', 0)}")
        print(f"siRNA candidates generated: {design_summary.get('total_candidates', 0)}")
        print(f"Off-target analysis: {offtarget_summary.get('method', 'basic')}")
        print()

        # Key output files
        print("📁 Key Output Files:")
        print("-" * 20)
        print(f"• Transcripts: {output_dir}/transcripts/{gene_query}_transcripts.fasta")
        print(f"• siRNA (ALL): {output_dir}/sirnaforge/{gene_query}_all.csv")
        print(f"• siRNA (PASS): {output_dir}/sirnaforge/{gene_query}_pass.csv")
        print(f"• Off-target results: {output_dir}/off_target/results/")
        print(f"• Workflow summary: {output_dir}/logs/workflow_summary.json")

        if offtarget_summary.get("method") == "nextflow":
            print(f"• Full off-target report: {output_dir}/off_target/results/offtarget_report.html")

        print()
        print("🎯 Next steps:")
        print("1. Review siRNA candidates in the CSV files (ALL and PASS)")
        print("2. Examine off-target analysis results")
        print("3. Select top candidates for experimental validation")
        print("4. Consider additional specificity testing if needed")

    except Exception as e:
        print(f"❌ Workflow failed: {e}")
        traceback.print_exc()


async def run_nextflow_only_example() -> None:
    """Example of running just the Nextflow off-target analysis."""

    print("\n🔬 Nextflow Off-target Analysis Example")
    print("=" * 45)

    # Create example siRNA sequences
    output_dir = Path(__file__).parent / "nextflow_output"
    output_dir.mkdir(exist_ok=True)

    # Example sequences for TP53
    example_sequences = [
        ("TP53_siRNA_1", "GCAUGAACCGGAGGCCCAUUU"),
        ("TP53_siRNA_2", "GAAUGUGAAUGAACACUGAUU"),
        ("TP53_siRNA_3", "CAUCCCACUACAAGUGUGAUU"),
    ]

    # Write to FASTA
    input_fasta = output_dir / "example_sirnas.fasta"
    with input_fasta.open("w") as f:
        for seq_id, sequence in example_sequences:
            f.write(f">{seq_id}\n{sequence}\n")

    print(f"Created example siRNA file: {input_fasta}")
    print(f"Sequences: {len(example_sequences)}")

    # Run Nextflow command
    nextflow_script = Path(__file__).parent.parent / "nextflow_pipeline" / "main.nf"

    if nextflow_script.exists():
        cmd = [
            "nextflow",
            "run",
            str(nextflow_script),
            "--input",
            str(input_fasta),
            "--outdir",
            str(output_dir / "results"),
            "--genome_species",
            "human",
            "--download_indexes",
            "true",  # This will build indexes if needed
        ]

        print("\n🚀 Running Nextflow command:")
        print(" ".join(cmd))

        try:
            result = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=1800)

            if result.returncode == 0:
                print("✅ Nextflow pipeline completed successfully!")
                print(f"📁 Results available at: {output_dir / 'results'}")
            else:
                print(f"❌ Nextflow failed: {result.stderr}")

        except subprocess.TimeoutExpired:
            print("⏰ Nextflow timed out after 30 minutes")
        except FileNotFoundError:
            print("❌ Nextflow not found. Please install Nextflow to run this example.")
            print("Visit: https://www.nextflow.io/docs/latest/getstarted.html")

    else:
        print(f"❌ Nextflow script not found at: {nextflow_script}")


if __name__ == "__main__":
    print("🧬 siRNAforge Toolkit - Complete Workflow Examples")
    print("=" * 60)

    if len(sys.argv) > 1 and sys.argv[1] == "--nextflow-only":
        # Run only Nextflow example
        asyncio.run(run_nextflow_only_example())
    else:
        # Run complete workflow
        asyncio.run(run_complete_example())
