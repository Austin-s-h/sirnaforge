process PREPARE_CANDIDATES {
    tag "$meta.id"
    label 'process_low'

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
from sirnaforge.pipeline.nextflow_cli import prepare_candidates_cli

# Validate siRNA candidates
result = prepare_candidates_cli(
    input_fasta='${candidates_fasta}',
    output_fasta='validated_candidates.fasta',
    expected_length=21
)

print(f'Validated {result[\"valid\"]} out of {result[\"total\"]} candidates')
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
