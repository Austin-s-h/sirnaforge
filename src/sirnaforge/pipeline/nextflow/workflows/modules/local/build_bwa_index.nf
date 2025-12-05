process BUILD_BWA_INDEX {
    tag "$species"
    label 'process_high'

    input:
    tuple val(species), path(genome_fasta)

    output:
    tuple val(species), path("${species}_index*"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import build_bwa_index_cli

result = build_bwa_index_cli(
    fasta_file='${genome_fasta}',
    species='${species}',
    output_dir='.'
)

print(f"Built BWA index for {result['species']}: {result['index_prefix']}")
PYEOF

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
