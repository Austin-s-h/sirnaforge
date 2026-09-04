/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ZFN OFF-TARGET SEARCH MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Exhaustive sliding-window ZFN off-target search.
    Wraps the Python-native ZFNOffTargetSearcher for Nextflow orchestration,
    containerization, resource management and -resume semantics.

    The algorithm is NOT BWA-MEM2 alignment — it is exhaustive paired half-site
    scanning with FokI seed-region penalties, spacer constraints, and IUPAC handling.
----------------------------------------------------------------------------------------
*/

process ZFN_OFFTARGET_SEARCH {
    tag "zfn_offtarget:${shard_id}"
    label 'process_high'
    publishDir "${params.outdir}/zfn_offtarget", mode: params.publish_dir_mode

    input:
    tuple val(shard_id), val(shard_chrom), val(core_start_1), val(core_end_1), val(scan_start_1), val(scan_end_1), val(shard_max_mismatches)
    val left_half_site
    val right_half_site
    path genome_fasta
    val algorithm
    val dimer_mode
    val spacer_lengths
    path annotation_file

    output:
    path "offtarget_sites_${shard_id}.csv", emit: sites
    path "candidate_summary_${shard_id}.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def annotation_arg = annotation_file.name != 'NO_ANNOTATION' ? "--annotation-file ${annotation_file} \\\n" : ''
    """
    sirnaforge _internal zfn-search-shard \
      --shard-id '${shard_id}' \
      --shard-chrom '${shard_chrom}' \
            --core-start-1 ${core_start_1} \
            --core-end-1 ${core_end_1} \
      --scan-start-1 ${scan_start_1} \
      --scan-end-1 ${scan_end_1} \
      --shard-max-mismatches ${shard_max_mismatches} \
      --left-half-site '${left_half_site}' \
      --right-half-site '${right_half_site}' \
      --genome-fasta ${genome_fasta} \
      --algorithm '${algorithm}' \
      --dimer-mode '${dimer_mode}' \
      --spacer-lengths '${spacer_lengths}' \
    ${annotation_arg}      --output-sites-csv offtarget_sites_${shard_id}.csv \
      --output-summary-json candidate_summary_${shard_id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sirnaforge: \$(python -c 'import sirnaforge; print(sirnaforge.__version__)' 2>/dev/null || echo 'dev')
    END_VERSIONS
    """

    stub:
    """
    echo "site_id,chrom,start_1based,end_1based,strand,orientation,spacer_len,sequence,left_mismatches,right_mismatches,total_mismatches,score,region,nearest_gene,left_aligned,right_aligned" > offtarget_sites_${shard_id}.csv
    echo '{"candidates": [], "summary": {"candidates": 0, "off_target_sites": 0}}' > candidate_summary_${shard_id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sirnaforge: \$(python -c 'import sirnaforge; print(sirnaforge.__version__)' 2>/dev/null || echo 'dev')
    END_VERSIONS
    """
}
