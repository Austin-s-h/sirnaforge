/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ZFN AGGREGATE RESULTS MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Merges shard CSVs, deduplicates by genomic coordinates + orientation,
    and globally sorts using deterministic ZFN ranking.
----------------------------------------------------------------------------------------
*/

process ZFN_AGGREGATE_RESULTS {
    tag "zfn_aggregate"
    label 'process_single'
    publishDir "${params.outdir}/zfn_offtarget", mode: params.publish_dir_mode

    input:
    path shard_site_csvs

    output:
    path "offtarget_sites.csv", emit: sites
    path "candidate_summary.json", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    sirnaforge _internal zfn-aggregate-shards \
      --shard-csv-glob 'offtarget_sites_*.csv' \
      --output-sites-csv offtarget_sites.csv \
      --output-summary-json candidate_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo "site_id,chrom,start_1based,end_1based,strand,orientation,spacer_len,sequence,left_mismatches,right_mismatches,total_mismatches,score,region,nearest_gene,left_aligned,right_aligned" > offtarget_sites.csv
    echo '{"candidates": [], "summary": {"off_target_sites": 0, "shards": 0}}' > candidate_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
