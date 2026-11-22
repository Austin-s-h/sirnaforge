process BUILD_BWA_INDEX {
    tag "$species"
    label 'process_high'

    input:
    tuple val(species), path(transcriptome_fasta)

    output:
    tuple val(species), val("${species}${params.aligner == 'bwa-mem3' ? '_bwa_mem3' : ''}_index"), path("${species}${params.aligner == 'bwa-mem3' ? '_bwa_mem3' : ''}_index*"), emit: index
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
    fasta_file='${transcriptome_fasta}',
    species='${species}',
    output_dir='.',
    aligner='${params.aligner ?: 'bwa-mem2'}'
)

print(f"Built BWA index for {result['species']}: {result['index_prefix']}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    ${params.aligner ?: 'bwa-mem2'}: \$(if [ "${params.aligner ?: 'bwa-mem2'}" = "bwa-mem3" ]; then bwa-mem3 version 2>/dev/null | head -n1; else bwa-mem2 version 2>&1 | tail -n1; fi)
END_VERSIONS
    """

    stub:
    """
    touch ${species}${params.aligner == 'bwa-mem3' ? '_bwa_mem3' : ''}_index.amb
    touch ${species}${params.aligner == 'bwa-mem3' ? '_bwa_mem3' : ''}_index.ann
    touch ${species}${params.aligner == 'bwa-mem3' ? '_bwa_mem3' : ''}_index.bwt.2bit.64
    touch ${species}${params.aligner == 'bwa-mem3' ? '_bwa_mem3' : ''}_index.pac

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    ${params.aligner ?: 'bwa-mem2'}: \$(if [ "${params.aligner ?: 'bwa-mem2'}" = "bwa-mem3" ]; then bwa-mem3 version 2>/dev/null | head -n1; else bwa-mem2 version 2>&1 | tail -n1; fi)
END_VERSIONS
    """
}
