process BUILD_BWA_INDEX {
    tag "$species"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0':
        'biocontainers/bwa-mem2:2.2.1--he513fc3_0' }"

    input:
    tuple val(species), path(genome_fasta)

    output:
    tuple val(species), path("${species}_index*"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import build_bwa_index

# Build BWA index
index_prefix = build_bwa_index(
    fasta_file='${genome_fasta}',
    index_prefix='${species}_index'
)

print(f'Built BWA index for ${species}: {index_prefix}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//')
    END_VERSIONS
    """

    stub:
    """
    touch ${species}_index.0123
    touch ${species}_index.amb
    touch ${species}_index.ann
    touch ${species}_index.bwt.2bit.64
    touch ${species}_index.pac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//')
    END_VERSIONS
    """
}
