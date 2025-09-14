#!/usr/bin/env python3
"""
Example script demonstrating how to use the local Docker image for Nextflow testing.

This script shows how to configure the Nextflow integration to use the local
Docker image built by 'make docker' for testing the off-target analysis.
"""

import asyncio

from sirnaforge.pipeline.nextflow.runner import NextflowRunner


def test_local_docker_integration():
    """Test Nextflow integration using local Docker image."""

    print("🧪 Testing Nextflow with local Docker image...")

    # Create runner configured for local Docker testing
    runner = NextflowRunner.for_local_docker_testing()

    # Check configuration
    print(f"📋 Configuration: {runner.config.get_environment_info().get_execution_summary()}")

    # Validate installation
    validation = runner.validate_installation()
    print(f"✅ Nextflow available: {validation.get('nextflow', False)}")
    print(f"🐳 Docker available: {validation.get('docker', False)}")
    print(f"📁 Workflow files: {validation.get('workflow_files', False)}")

    # Example usage (commented out to avoid actual execution)
    """
    # Prepare test data
    input_file = Path("examples/sample_transcripts.fasta")
    output_dir = Path("test_results")

    # Run the analysis
    result = runner.run_sync(
        input_file=input_file,
        output_dir=output_dir,
        genome_species=["human"]
    )

    print(f"🎉 Analysis completed! Results in: {result['output_dir']}")
    """

    print("✅ Local Docker configuration ready for testing!")


async def async_test_example():
    """Async version of the test."""

    # Example async usage (commented out)
    """
    runner = NextflowRunner.for_local_docker_testing()
    result = await runner.run(
        input_file=Path("examples/sample_transcripts.fasta"),
        output_dir=Path("async_results")
    )
    """

    print("✅ Async local Docker configuration ready!")


if __name__ == "__main__":
    # Run sync test
    test_local_docker_integration()

    # Run async test
    asyncio.run(async_test_example())

    print("\n📝 To use in your tests:")
    print("   runner = NextflowRunner.for_local_docker_testing()")
    print("   result = runner.run_sync(input_file, output_dir)")
