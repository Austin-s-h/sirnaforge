process SPLIT_CANDIDATES {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(candidates_fasta)

    output:
    tuple val(meta), path("candidate_*.fasta"), emit: individual_candidates
    path "candidate_manifest.txt", emit: manifest
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import split_candidates_cli

# Split candidates into individual files
result = split_candidates_cli(
    input_fasta='${candidates_fasta}',
    output_dir='.'
)

print(f'Split {result[\"candidate_count\"]} candidates into individual files')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch candidate_0001.fasta
    echo -e "sequence_id\\tfilename\\ntest_seq\\tcandidate_0001.fasta" > candidate_manifest.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
