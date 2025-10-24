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
    tuple val(candidate_meta), val(species), val(index_type), path("*.tsv"), path("*_summary.txt"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def candidate_id = candidate_meta.id
    def output_prefix = "${candidate_id}_${species}_${index_type}"
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import run_comprehensive_offtarget_analysis
import os

print(f'Running ${index_type} analysis for candidate ${candidate_id} against ${species}')
print(f'Using index path: ${index_path}')

# Run comprehensive off-target analysis using the proper Nextflow entrypoint
tsv_path, json_path, summary_path = run_comprehensive_offtarget_analysis(
    species='${species}',
    sequences_file='${candidate_fasta}',
    index_path='${index_path}',
    output_prefix='${output_prefix}',
    bwa_k=${bwa_k},
    bwa_T=${bwa_T},
    max_hits=${max_hits},
    seed_start=${seed_start},
    seed_end=${seed_end}
)

print(f'Analysis completed for ${candidate_id}-${species}-${index_type}')
print(f'Results written to: {tsv_path}, {json_path}, {summary_path}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """

    stub:
    def output_prefix = "${candidate_meta.id}_${species}_${index_type}"
    """
    touch ${output_prefix}.tsv
    echo "Species: ${species}" > ${output_prefix}_summary.txt
    echo "Total sequences analyzed: 0" >> ${output_prefix}_summary.txt
    echo "Total off-target hits: 0" >> ${output_prefix}_summary.txt
    echo "Analysis completed successfully" >> ${output_prefix}_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """
}
