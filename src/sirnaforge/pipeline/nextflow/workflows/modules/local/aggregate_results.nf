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
    def analysisFileList = analysis_files instanceof List ? analysis_files.collect { f -> f.toString() } : [analysis_files.toString()]
    def summaryFileList = summary_files instanceof List ? summary_files.collect { f -> f.toString() } : [summary_files.toString()]
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
    # The stub publishes the real header, including the classification columns, so a stub run and a
    # real run agree on the column set. It publishes no rows: a stub screened nothing.
    # Written literally, in bash: a stub run must exercise the wiring with no container, no aligner
    # and no importable sirnaforge. tests/unit/test_offtarget_aggregation_completeness.py pins this
    # line against OffTargetHit.tsv_header() + CLASSIFICATION_COLUMNS so it cannot drift.
    touch combined_mirna_analysis.tsv
    printf 'qname\\tqseq\\tspecies\\trname\\tcoord\\tstrand\\tcigar\\tmapq\\tas_score\\tnm\\tseed_mismatches\\tofftarget_score\\thit_class\\tmatched_symbol\\tsymbol_lookup_missing\\thit_symbol\\thit_symbol_missing\\tspecies_index_missing\\tortholog_evidence\\n' > combined_offtargets.tsv
    echo '{"status": "stub", "species_screened": [], "unscreened_species": [], "total_results": 0}' > combined_summary.json
    echo 'Aggregation completed' > final_summary.txt
    touch analysis_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
