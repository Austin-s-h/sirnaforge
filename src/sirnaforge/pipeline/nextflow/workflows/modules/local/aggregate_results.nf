process AGGREGATE_RESULTS {
    tag "aggregate"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: params.publish_dir_mode

    input:
    path analysis_files
    path summary_files
    val genome_species

    output:
    path "combined_*.tsv", emit: combined_analyses, optional: true
    path "combined_summary.json", emit: combined_summary, optional: true
    path "final_summary.txt", emit: final_summary
    path "analysis_report.html", emit: html_report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import aggregate_results_cli

# Run aggregation
result = aggregate_results_cli(
    genome_species='${genome_species}',
    output_dir='.'
)

print(f'Aggregation status: {result[\"status\"]}')
print(f'Processed {result[\"analysis_files_processed\"]} analysis files')
print(f'Processed {result[\"summary_files_processed\"]} summary files')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch combined_mirna_analysis.tsv
    touch combined_transcriptome_analysis.tsv
    echo '{}' > combined_summary.json
    echo 'Aggregation completed' > final_summary.txt
    touch analysis_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
