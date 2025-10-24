process OFFTARGET_ANALYSIS {
    tag "${candidate_meta.id}-$species-$index_type"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f':
        'community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f' }"

    input:
    tuple val(candidate_meta), path(candidate_fasta), val(species), val(index_path), val(index_type)
    val max_hits
    val bwa_k
    val bwa_T
    val seed_start
    val seed_end

    output:
    tuple val(candidate_meta), val(species), val(index_type), path("*_analysis.tsv"), path("*_summary.json"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def candidate_id = candidate_meta.id
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import run_bwa_alignment_analysis
import os

print(f'Running ${index_type} analysis for candidate ${candidate_id} against ${species}')
print(f'Using index path: ${index_path}')

# Run analysis
result_dir = run_bwa_alignment_analysis(
    candidates_file='${candidate_fasta}',
    index_prefix='${index_path}',
    species='${species}',
    output_dir='.',
    max_hits=${max_hits},
    bwa_k=${bwa_k},
    bwa_T=${bwa_T},
    seed_start=${seed_start},
    seed_end=${seed_end}
)

print(f'Analysis completed for ${candidate_id}-${species}-${index_type}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """

    stub:
    """
    touch ${candidate_meta.id}_${species}_${index_type}_analysis.tsv
    echo '{"candidate": "${candidate_meta.id}", "species": "${species}", "index_type": "${index_type}"}' > ${candidate_meta.id}_${species}_${index_type}_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """
}
