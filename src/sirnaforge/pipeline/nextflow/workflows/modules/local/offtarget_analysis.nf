process OFFTARGET_ANALYSIS {
    tag "$species"
    label 'process_medium'
    publishDir "${params.outdir}/genome", mode: params.publish_dir_mode

    input:
    tuple val(species), val(index_path), path(candidates_fasta)
    val max_hits
    val bwa_k
    val bwa_T
    val seed_start
    val seed_end

    output:
    path "${species}_analysis.tsv", emit: analysis
    path "${species}_summary.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run off-target analysis for ALL candidates against this genome in one session
    # This is much more efficient: load index once, process all candidates sequentially.
    # The CLI also owns the "this prefix is not a usable index" case: it publishes an EMPTY
    # analysis file plus a failed summary, which the aggregator reports as a per-species
    # rejection instead of a completed screen with zero hits.
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import offtarget_analysis_cli

result = offtarget_analysis_cli(
    species='${species}',
    index_prefix='${index_path}',
    candidates_file='${candidates_fasta}',
    output_dir='.',
    max_hits=${max_hits},
    bwa_k=${bwa_k},
    bwa_T=${bwa_T},
    seed_start=${seed_start},
    seed_end=${seed_end}
)

print(f"Batch analysis for ${species}: {result['status']}")
if result.get('error'):
    print(f"  reason: {result['error']}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """

    stub:
    """
    touch ${species}_analysis.tsv
    echo '{"species": "${species}", "total_candidates": 0, "total_hits": 0}' > ${species}_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """
}
