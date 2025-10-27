process MIRNA_SEED_ANALYSIS {
    tag "mirna_batch"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f':
        'community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f' }"

    input:
    path candidates_fasta
    val mirna_db
    val mirna_species

    output:
    path "mirna_analysis.tsv", emit: analysis
    path "mirna_summary.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run miRNA analysis for ALL candidates in one batch - efficient!
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.core.off_target import run_mirna_seed_analysis

print(f'Running batch miRNA seed match analysis')
print(f'miRNA database: ${mirna_db}')
print(f'miRNA species: ${mirna_species}')

# Parse species list
species_list = [s.strip() for s in '${mirna_species}'.split(',') if s.strip()]

# Run batch miRNA analysis - all candidates in one session
output_path = run_mirna_seed_analysis(
    candidates_file='${candidates_fasta}',
    candidate_id='batch',
    mirna_db='${mirna_db}',
    mirna_species=species_list,
    output_dir='.'
)

# Rename output files for consistency
import os
if os.path.exists('batch_mirna_analysis.tsv'):
    os.rename('batch_mirna_analysis.tsv', 'mirna_analysis.tsv')
if os.path.exists('batch_mirna_summary.json'):
    os.rename('batch_mirna_summary.json', 'mirna_summary.json')

print(f'Batch miRNA analysis completed for all candidates')
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch mirna_analysis.tsv
    echo '{"total_candidates": 0, "total_hits": 0}' > mirna_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
