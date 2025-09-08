process PREPARE_CANDIDATES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython:e5b315e81e28f4c6':
        'community.wave.seqera.io/library/python_biopython:e5b315e81e28f4c6' }"

    input:
    tuple val(meta), path(candidates_fasta)

    output:
    tuple val(meta), path("validated_candidates.fasta"), path("validation_report.txt"), emit: candidates
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import validate_and_write_sequences

# Validate siRNA candidates
total, valid, errors = validate_and_write_sequences(
    input_file='${candidates_fasta}',
    output_file='validated_candidates.fasta',
    expected_length=21
)

# Write validation report
with open('validation_report.txt', 'w') as f:
    f.write(f'Total candidates: {total}\\n')
    f.write(f'Valid candidates: {valid}\\n')
    f.write(f'Invalid candidates: {total - valid}\\n')
    if errors:
        f.write('\\nErrors:\\n')
        for error in errors:
            f.write(f'  {error}\\n')

print(f'Validated {valid} out of {total} candidates')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch validated_candidates.fasta
    echo "Total candidates: 0" > validation_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
