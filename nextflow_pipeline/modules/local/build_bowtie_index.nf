process BUILD_BOWTIE_INDEX {
    tag "$species"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/bowtie:1.3.1--py39hf95cd2a_0':
        'biocontainers/bowtie:1.3.1--py39hf95cd2a_0' }"

    input:
    tuple val(species), path(genome_fasta)

    output:
    tuple val(species), path("${species}_bowtie_index*"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import build_bowtie_index

# Build Bowtie index
index_prefix = build_bowtie_index(
    fasta_file='${genome_fasta}',
    index_prefix='${species}_bowtie_index'
)

print(f'Built Bowtie index for ${species}: {index_prefix}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bowtie: \$(bowtie --version 2>&1 | head -n1 | sed 's/.*bowtie-align-s version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${species}_bowtie_index.1.ebwt
    touch ${species}_bowtie_index.2.ebwt
    touch ${species}_bowtie_index.3.ebwt
    touch ${species}_bowtie_index.4.ebwt
    touch ${species}_bowtie_index.rev.1.ebwt
    touch ${species}_bowtie_index.rev.2.ebwt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bowtie: \$(bowtie --version 2>&1 | head -n1 | sed 's/.*bowtie-align-s version //')
    END_VERSIONS
    """
}
