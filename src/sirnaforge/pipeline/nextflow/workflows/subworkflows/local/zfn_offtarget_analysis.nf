include { ZFN_OFFTARGET_SEARCH_FULL } from '../../modules/local/zfn_offtarget_search_full'

workflow ZFN_OFFTARGET_ANALYSIS {
    take:
    genome_fasta
    zfn_contract
    annotation_file
    annotation_config
    search_index_file

    main:
    ZFN_OFFTARGET_SEARCH_FULL(genome_fasta, zfn_contract, annotation_file, annotation_config, search_index_file)

    emit:
    sites = ZFN_OFFTARGET_SEARCH_FULL.out.sites
    summary = ZFN_OFFTARGET_SEARCH_FULL.out.summary
    versions = ZFN_OFFTARGET_SEARCH_FULL.out.versions
}
