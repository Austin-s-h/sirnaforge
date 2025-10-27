process OFFTARGET_ANALYSIS {
    tag "$candidate_id-$species"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f':
        'community.wave.seqera.io/library/python_biopython_pyyaml:a9b2e2e522b05e9f' }"

    input:
    tuple val(candidate_id), val(fasta_content), val(species), val(index_path)
    val max_hits
    val bwa_k
    val bwa_T
    val seed_start
    val seed_end

    output:
    path "${candidate_id}_${species}_analysis.tsv", emit: analysis
    path "${candidate_id}_${species}_summary.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Write candidate FASTA content to file
    echo '${fasta_content}' > candidate.fasta

    # Run off-target analysis using CLI function
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import run_offtarget_analysis_cli

result = run_offtarget_analysis_cli(
    candidate_fasta='candidate.fasta',
    candidate_id='${candidate_id}',
    species='${species}',
    index_prefix='${index_path}',
    output_dir='.',
    max_hits=${max_hits},
    bwa_k=${bwa_k},
    bwa_T=${bwa_T},
    seed_start=${seed_start},
    seed_end=${seed_end}
)

print(f"Analysis completed: {result['candidate_id']}-{result['species']}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """

    stub:
    """
    touch ${candidate_id}_${species}_analysis.tsv
    echo '{"candidate_id": "${candidate_id}", "species": "${species}"}' > ${candidate_id}_${species}_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """
}
