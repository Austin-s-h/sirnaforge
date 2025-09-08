process SPLIT_CANDIDATES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython:e5b315e81e28f4c6':
        'community.wave.seqera.io/library/python_biopython:e5b315e81e28f4c6' }"

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
from sirnaforge.core.off_target import parse_fasta_sequences
from pathlib import Path

# Parse candidates
sequences = parse_fasta_sequences('${candidates_fasta}')

# Split into individual files
candidate_files = []
for i, (seq_id, sequence) in enumerate(sequences.items()):
    filename = f'candidate_{i+1:04d}.fasta'
    with open(filename, 'w') as f:
        f.write(f'>{seq_id}\\n{sequence}\\n')
    candidate_files.append(f'{seq_id}\\t{filename}')

# Write manifest
with open('candidate_manifest.txt', 'w') as f:
    f.write('sequence_id\\tfilename\\n')
    for entry in candidate_files:
        f.write(f'{entry}\\n')

print(f'Split {len(sequences)} candidates into individual files')
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
