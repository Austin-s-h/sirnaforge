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

include { ZFN_OFFTARGET_SEARCH } from '../../modules/local/zfn_offtarget_search'

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

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run exhaustive ZFN off-target search
    //
    ZFN_OFFTARGET_SEARCH(
        left_half_site,
        right_half_site,
        genome_fasta,
        algorithm,
        dimer_mode,
        spacer_lengths,
        max_mismatches,
        annotation_file
    )
    ch_versions = ch_versions.mix(ZFN_OFFTARGET_SEARCH.out.versions)

    emit:
    sites           = ZFN_OFFTARGET_SEARCH.out.sites
    summary         = ZFN_OFFTARGET_SEARCH.out.summary
    versions        = ch_versions
}
