process MIRNA_OFFTARGET_ANALYSIS {
    tag "${candidate_meta.id}-$species"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f':
        'community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f' }"

    input:
    tuple val(candidate_meta), path(candidate_fasta), val(species), val(mirna_index)

    output:
    tuple val(candidate_meta), val(species), path("*_mirna_hits.tsv"), path("*_mirna_hits.json"), path("*_mirna_summary.txt"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def candidate_id = candidate_meta.id
    def output_prefix = "${candidate_id}_${species}"
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import run_mirna_analysis_for_nextflow
import os

print(f'Running miRNA off-target analysis for candidate ${candidate_id} against ${species}')
print(f'Using miRNA index: ${mirna_index}')

# Run miRNA analysis using the dedicated Nextflow entrypoint
tsv_path, json_path, summary_path = run_mirna_analysis_for_nextflow(
    species='${species}',
    sequences_file='${candidate_fasta}',
    mirna_index='${mirna_index}',
    output_prefix='${output_prefix}'
)

print(f'miRNA analysis completed for ${candidate_id}-${species}')
print(f'Results written to: {tsv_path}, {json_path}, {summary_path}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """

    stub:
    def output_prefix = "${candidate_meta.id}_${species}"
    """
    touch ${output_prefix}_mirna_hits.tsv
    echo '[]' > ${output_prefix}_mirna_hits.json
    echo "Species: ${species}" > ${output_prefix}_mirna_summary.txt
    echo "Total miRNA hits: 0" >> ${output_prefix}_mirna_summary.txt
    echo "Analysis completed successfully" >> ${output_prefix}_mirna_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """
}
