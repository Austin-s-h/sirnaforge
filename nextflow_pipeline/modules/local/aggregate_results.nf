process AGGREGATE_RESULTS {
    tag "aggregate"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: params.publish_dir_mode

    conda "${moduleDir}/environment.yml"

    input:
    path analysis_files
    path summary_files
    val genome_species

    output:
    path "combined_*.tsv", emit: combined_analyses, optional: true
    path "combined_summary.json", emit: combined_summary, optional: true
    path "final_summary.txt", emit: final_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run aggregation using CLI function
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import aggregate_results_cli

result = aggregate_results_cli(
    genome_species='${genome_species}',
    output_dir='.'
)

print(f"Aggregation completed: {result['status']}")
print(f"Analysis files processed: {result['analysis_files_processed']}")
print(f"Summary files processed: {result['summary_files_processed']}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo 'not available')
    END_VERSIONS
    """.stripIndent()

    stub:
    """
    touch combined_offtarget_analysis.tsv
    echo '{"status": "completed"}' > combined_summary.json
    echo 'Aggregation completed' > final_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo 'not available')
    END_VERSIONS
    """.stripIndent()
}
