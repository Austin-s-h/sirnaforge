def zfnShellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

process ZFN_OFFTARGET_SEARCH_FULL {
    tag "zfn_full_search"
    memory { params.max_memory as nextflow.util.MemoryUnit }
    publishDir "${params.outdir}/zfn_offtarget", mode: params.publish_dir_mode

    input:
    path genome_fasta
    path zfn_contract
    path annotation_file
    path annotation_config
    path search_index_file

    output:
    path "offtarget_sites.csv", emit: sites
    path "candidate_summary.json", emit: summary
    path "versions.yml", emit: versions

    script:
    def annotation_arg = annotation_file.name != 'NO_ANNOTATION' ? "--annotation-file ${zfnShellQuote(annotation_file)} \\\n" : ''
    def annotation_config_arg = annotation_config.name != 'NO_ANNOTATION_CONFIG' ? "--annotation-config-json ${zfnShellQuote(annotation_config)} \\\n" : ''
    def index_arg = search_index_file.name != 'NO_SEARCH_INDEX' ? "--search-space-index ${zfnShellQuote(search_index_file)} \\\n" : ''
    """
    sirnaforge _internal zfn-search-full \\
      --genome-fasta ${zfnShellQuote(genome_fasta)} \\
      --design-params-json ${zfnShellQuote(zfn_contract)} \\
    ${annotation_arg}${annotation_config_arg}${index_arg}      --output-sites-csv offtarget_sites.csv \\
      --output-summary-json candidate_summary.json

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    sirnaforge: \$(python -c 'import sirnaforge; print(sirnaforge.__version__)' 2>/dev/null || echo 'dev')
END_VERSIONS
    """

    stub:
    """
    echo "site_id,chrom,start_1based,end_1based,strand,orientation,spacer_len,sequence,left_mismatches,right_mismatches,left_seed_mismatches,right_seed_mismatches,left_mismatch_positions,right_mismatch_positions,total_mismatches,score,dimer_compatible,score_components,region,nearest_gene,left_aligned,right_aligned" > offtarget_sites.csv
    echo '{"schema_version": "zfn_candidate_summary.v1", "search_contract": {}, "candidates": [], "summary": {"off_target_sites": 0}}' > candidate_summary.json

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    sirnaforge: \$(python -c 'import sirnaforge; print(sirnaforge.__version__)' 2>/dev/null || echo 'dev')
END_VERSIONS
    """
}
