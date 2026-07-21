/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ZFN OFF-TARGET ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Orchestrates ZFN pair evaluation and exhaustive off-target search.
    This is scientifically distinct from siRNA off-target (BWA-MEM2 alignment):
    - Input: two IUPAC half-site sequences (9-18 bp each)
    - Algorithm: exhaustive sliding window across whole-genome FASTA
    - Pairing: finds L+R half-site pairs with valid spacer and dimer mode
    - Scoring: FokI seed-region penalties, three algorithm variants
----------------------------------------------------------------------------------------
*/

include { ZFN_MAKE_SHARDS } from '../../modules/local/zfn_make_shards'
include { ZFN_OFFTARGET_SEARCH } from '../../modules/local/zfn_offtarget_search'
include { ZFN_AGGREGATE_RESULTS } from '../../modules/local/zfn_aggregate_results'

workflow ZFN_OFFTARGET_ANALYSIS {
    take:
    left_half_site      // val: left ZFN half-site sequence
    right_half_site     // val: right ZFN half-site sequence
    genome_fasta        // path: genome FASTA for off-target search
    algorithm           // val: scoring algorithm (homology, conserved_g, zfn_v2)
    dimer_mode          // val: dimer mode (heterodimer_only, include_homodimers)
    spacer_lengths      // val: comma-separated allowed spacer lengths
    max_mismatches      // val: max mismatches per half-site
    annotation_file     // path: optional GTF/GFF annotation file (or NO_ANNOTATION sentinel)
    sharding_enabled    // val: enable chromosome/chunk sharding
    shard_chunk_mb      // val: chunk size in MB
    shard_overlap_bp    // val: overlap in bp (auto-raised to safety minimum)
    shard_chromosomes   // val: comma-separated chromosomes (optional)

    main:
    ch_versions = Channel.empty()

    ZFN_MAKE_SHARDS(
        genome_fasta,
        left_half_site,
        right_half_site,
        spacer_lengths,
        sharding_enabled,
        shard_chunk_mb,
        shard_overlap_bp,
        shard_chromosomes
    )
    ch_versions = ch_versions.mix(ZFN_MAKE_SHARDS.out.versions)

    ch_shards = ZFN_MAKE_SHARDS.out.shards
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.shard_id,
                row.chrom,
                row.core_start_1 as int,
                row.core_end_1 as int,
                row.scan_start_1 as int,
                row.scan_end_1 as int,
                row.max_mismatches as int
            )
        }

    ZFN_OFFTARGET_SEARCH(
        ch_shards,
        left_half_site,
        right_half_site,
        genome_fasta,
        algorithm,
        dimer_mode,
        spacer_lengths,
        annotation_file
    )
    ch_versions = ch_versions.mix(ZFN_OFFTARGET_SEARCH.out.versions)

    ZFN_AGGREGATE_RESULTS(
        ZFN_OFFTARGET_SEARCH.out.sites.collect()
    )
    ch_versions = ch_versions.mix(ZFN_AGGREGATE_RESULTS.out.versions)

    emit:
    sites           = ZFN_AGGREGATE_RESULTS.out.sites
    summary         = ZFN_AGGREGATE_RESULTS.out.summary
    versions        = ch_versions
}
