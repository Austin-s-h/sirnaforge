process OFFTARGET_ANALYSIS {
    tag "$species"
    label 'process_medium'

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
    # This is much more efficient: load index once, process all candidates sequentially
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import run_bwa_alignment_analysis

# Run batch analysis: one BWA session, all candidates
output_path = run_bwa_alignment_analysis(
    candidates_file='${candidates_fasta}',
    index_prefix='${index_path}',
    species='${species}',
    output_dir='.',
    max_hits=${max_hits},
    bwa_k=${bwa_k},
    bwa_T=${bwa_T},
    seed_start=${seed_start},
    seed_end=${seed_end}
)

print(f"Batch analysis completed for ${species}: all candidates processed")
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
