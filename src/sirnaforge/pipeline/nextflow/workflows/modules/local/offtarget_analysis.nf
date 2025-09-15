process OFFTARGET_ANALYSIS {
    tag "${candidate_meta.id}-$species-$index_type"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(candidate_meta), path(candidate_fasta), val(species), val(index_path), val(index_type)
    val max_hits
    val bwa_k
    val bwa_T
    val seed_start
    val seed_end

    output:
    tuple val(candidate_meta), val(species), val(index_type), path("*_analysis.tsv"), path("*_summary.json"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def candidate_id = candidate_meta.id
    """
    python3 -c "
from sirnaforge.core.off_target import run_comprehensive_offtarget_analysis

print(f'Running comprehensive off-target analysis for candidate ${candidate_id} against ${species} (${index_type} index)')
print(f'Using index path: ${index_path}')

# Run comprehensive analysis
tsv_file, json_file, summary_file = run_comprehensive_offtarget_analysis(
    species='${species}',
    sequences_file='${candidate_fasta}',
    index_path='${index_path}',
    output_prefix='${candidate_id}_${species}_${index_type}',
    mode='transcriptome',  # Can be 'transcriptome' or 'mirna_seed'
    bwa_k=${bwa_k},
    bwa_T=${bwa_T},
    max_hits=${max_hits},
    seed_start=${seed_start},
    seed_end=${seed_end}
)

print(f'Analysis completed for ${candidate_id}-${species}-${index_type}')
print(f'Results: TSV={tsv_file}, JSON={json_file}, Summary={summary_file}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
        bowtie: \$(bowtie --version 2>&1 | head -n1 | sed 's/.*bowtie-align-s version //' || echo 'not available')
    END_VERSIONS
    """

    stub:
    """
    touch ${candidate_meta.id}_${species}_${index_type}_analysis.tsv
    echo '{"candidate": "${candidate_meta.id}", "species": "${species}", "index_type": "${index_type}"}' > ${candidate_meta.id}_${species}_${index_type}_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
        bowtie: \$(bowtie --version 2>&1 | head -n1 | sed 's/.*bowtie-align-s version //' || echo 'not available')
    END_VERSIONS
    """
}
