process MIRNA_SEED_ANALYSIS {
    tag "${candidate_meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f':
        'community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f' }"

    input:
    tuple val(candidate_meta), path(candidate_fasta)
    val mirna_db
    val mirna_species

    output:
    tuple val(candidate_meta), path("*_mirna_analysis.tsv"), path("*_mirna_summary.json"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def candidate_id = candidate_meta.id
    def species_str = mirna_species ? mirna_species.join(',') : 'human'
    """
    python3 -c "
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import run_mirna_seed_analysis_cli

print(f'Running miRNA seed match analysis for candidate ${candidate_id}')
print(f'miRNA database: ${mirna_db}')
print(f'miRNA species: ${species_str}')

# Run miRNA seed match analysis
result = run_mirna_seed_analysis_cli(
    candidate_fasta='${candidate_fasta}',
    candidate_id='${candidate_id}',
    mirna_db='${mirna_db}',
    mirna_species='${species_str}',
    output_dir='.'
)

print(f'miRNA analysis completed for ${candidate_id}')
print(f'Found {result[\"total_hits\"]} miRNA seed matches')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${candidate_meta.id}_mirna_analysis.tsv
    echo '{"candidate": "${candidate_meta.id}", "total_hits": 0}' > ${candidate_meta.id}_mirna_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
