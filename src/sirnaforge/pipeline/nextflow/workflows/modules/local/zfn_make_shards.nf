/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ZFN SHARD MANIFEST MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Builds chromosome/chunk shard manifest for scalable ZFN off-target search.
----------------------------------------------------------------------------------------
*/

process ZFN_MAKE_SHARDS {
    tag "zfn_make_shards"
    label 'process_single'
    publishDir "${params.outdir}/zfn_offtarget", mode: params.publish_dir_mode

    input:
    path genome_fasta
    val left_half_site
    val right_half_site
    val spacer_lengths
    val sharding_enabled
    val shard_chunk_mb
    val shard_overlap_bp
    val shard_chromosomes

    output:
    path "zfn_shards.tsv", emit: shards
    path "versions.yml", emit: versions

    script:
    """
    sirnaforge _internal zfn-make-shards \
      --genome-fasta ${genome_fasta} \
      --left-half-site '${left_half_site}' \
      --right-half-site '${right_half_site}' \
      --spacer-lengths '${spacer_lengths}' \
      --max-mismatches ${params.zfn_max_mismatches} \
      --sharding-enabled '${sharding_enabled}' \
      --shard-chunk-mb ${shard_chunk_mb} \
      --shard-overlap-bp ${shard_overlap_bp} \
      --shard-chromosomes '${shard_chromosomes}' \
      --output zfn_shards.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo -e "shard_id\tchrom\tcore_start_1\tcore_end_1\tscan_start_1\tscan_end_1\tmax_mismatches" > zfn_shards.tsv
    echo -e "chr3:1-1000\tchr3\t1\t1000\t1\t1000\t1" >> zfn_shards.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
