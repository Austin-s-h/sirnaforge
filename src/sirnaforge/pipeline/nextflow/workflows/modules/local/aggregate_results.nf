process AGGREGATE_RESULTS {
    tag "aggregate"
    label 'process_low'
    publishDir "${params.outdir}/aggregated", mode: params.publish_dir_mode

    input:
    val analysis_files
    val summary_files
    val genome_species

    output:
    path "combined_*.tsv", emit: combined_analyses, optional: true
    path "combined*.json", emit: combined_summary, optional: true
    path "final*_summary.txt", emit: final_summary
    path "analysis_report.html", emit: html_report, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def analysisFileList = analysis_files instanceof List ? analysis_files.collect { it.toString() } : [analysis_files.toString()]
    def summaryFileList = summary_files instanceof List ? summary_files.collect { it.toString() } : [summary_files.toString()]
    def analysisFilesJson = groovy.json.JsonOutput.toJson(analysisFileList)
    def summaryFilesJson = groovy.json.JsonOutput.toJson(summaryFileList)
    """
    python3 <<'PYEOF'
import sys
import json
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import aggregate_results_cli

mirna_db = '${params.mirna_db ?: 'mirgenedb'}'.strip()
mirna_species = '${params.mirna_species ?: 'chicken,pig,rat,mouse,human,macaque'}'.strip()
analysis_files = json.loads('''${analysisFilesJson}''')
summary_files = json.loads('''${summaryFilesJson}''')

result = aggregate_results_cli(
    genome_species='${genome_species}',
    output_dir='.',
    mirna_db=mirna_db or None,
    mirna_species=mirna_species or None,
    analysis_files=analysis_files,
    summary_files=summary_files,
)

print(f"Aggregation status: {result.get('status', 'unknown')}")
print(f"Processed {result.get('analysis_files_processed', 0)} analysis files")
print(f"Processed {result.get('summary_files_processed', 0)} summary files")
if result.get('mirna'):
    stats = result['mirna']
    print(f"miRNA aggregation: {stats.get('analysis_files_processed', 0)} files using {stats.get('mirna_db')}")
PYEOF

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
