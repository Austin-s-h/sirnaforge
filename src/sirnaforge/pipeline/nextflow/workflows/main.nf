#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sirnaforge/pipeline/nextflow/workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    siRNA Off-Target Analysis Pipeline - Embedded in Python Package
    Github: https://github.com/austin-s-h/sirnaforge
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SIRNA_OFFTARGET_ANALYSIS } from './subworkflows/local/sirna_offtarget_analysis'
include { ZFN_OFFTARGET_ANALYSIS   } from './subworkflows/local/zfn_offtarget_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SIRNAFORGE_OFFTARGET {

    main:
    //
    // Print parameter summary
    //
    log.info """\
        ===============================================
         S I R N A F O R G E   O F F - T A R G E T
        ===============================================
        input                : ${params.input}
        outdir               : ${params.outdir}

        TRANSCRIPTOME SCREEN (OPTIONAL - Resource Intensive)
        transcriptome_fastas : ${params.transcriptome_fastas ?: 'Not provided'}
        transcriptome_indices: ${params.transcriptome_indices ?: 'Not provided'}
        transcriptome_species: ${params.transcriptome_species ?: 'N/A'}

        ANALYSIS PARAMETERS
        max_hits             : ${params.max_hits}
        bwa_k                : ${params.bwa_k}
        bwa_T                : ${params.bwa_T}
        seed_start           : ${params.seed_start}
        seed_end             : ${params.seed_end}

        RESOURCES
        max_memory           : ${params.max_memory}
        """
        .stripIndent()

    //
    // Validate required parameters
    //
    if (!params.input) {
        error "Input FASTA file must be specified with --input"
    }
    if (!file(params.input).exists()) {
        error "Input file does not exist: ${params.input}"
    }

    //
    // Refuse the parameters the #99 rename removed, by name. Nextflow accepts an unknown --param
    // silently, so a stale genome_indices would configure no reference at all and the run would
    // report success having screened nothing.
    //
    // Nextflow stores a hyphenated --genome-indices under the camelCase key genomeIndices, so every
    // supplied key is folded to snake_case before it is compared: checking the snake_case spellings
    // alone let the hyphenated form through, and a stale one screened nothing and reported success.
    def renamed_params = [
        genome_fastas : 'transcriptome_fastas',
        genome_indices: 'transcriptome_indices',
        genome_species: 'transcriptome_species',
    ]
    params.keySet().each { supplied ->
        def folded = supplied.toString().replaceAll('-', '_').replaceAll(/([a-z0-9])([A-Z])/, '$1_$2').toLowerCase()
        if (renamed_params.containsKey(folded)) {
            error "--${supplied} was renamed to --${renamed_params[folded]}: siRNA/miRNA screening references are transcriptomes, and 'genome' now means ZFN genomic DNA."
        }
    }

    //
    // Create input channel - simple file input
    //
    ch_input = channel.fromPath(params.input, checkIfExists: true)

    //
    // Screening references: FASTAs to index here and prefixes already built, in one channel.
    // Both forms are transcriptomes; the tuple's third element says which work each one needs.
    //
    ch_references = channel.empty()

    if (params.transcriptome_fastas) {
        ch_references = ch_references.mix(
            channel.from(params.transcriptome_fastas.split(','))
                .map { entry ->
                    def (species, fasta_path) = entry.split(':')
                    [species.trim(), file(fasta_path.trim(), checkIfExists: true), 'fasta']
                }
        )
    }

    if (params.transcriptome_indices) {
        ch_references = ch_references.mix(
            channel.from(params.transcriptome_indices.split(','))
                .map { entry ->
                    def (species, index_path) = entry.split(':')
                    [species.trim(), index_path.trim(), 'index']
                }
        )
    }

    // Check if ANY transcriptome off-target analysis is enabled
    def has_offtarget_data = params.transcriptome_fastas || params.transcriptome_indices

    // With no reference at all, skip alignment (miRNA-only mode)
    if (!has_offtarget_data) {
        log.info ""
        log.info "=" * 80
        log.info "NOTE: No transcriptome FASTAs or indices provided"
        log.info "Transcriptome off-target analysis: DISABLED"
        log.info "Running lightweight miRNA seed match analysis only (< 1GB RAM)"
        log.info ""
        log.info "To enable transcriptome off-target analysis, provide either:"
        log.info "  --transcriptome_fastas 'species:path,species2:path2' OR"
        log.info "  --transcriptome_indices 'species:index_prefix,species2:index_prefix2'"
        log.info "=" * 80
        log.info ""
    }

    //
    // SUBWORKFLOW: Route by design_mode
    //
    if (params.design_mode == 'zfn') {
        //
        // ZFN mode: exhaustive paired half-site search (not BWA-MEM2)
        //
        def zfn_genome = params.zfn_genome_fasta ? file(params.zfn_genome_fasta, checkIfExists: true) : file('NO_GENOME')
        def zfn_annot  = params.zfn_annotation   ? file(params.zfn_annotation, checkIfExists: true) : file('NO_ANNOTATION')

        ZFN_OFFTARGET_ANALYSIS(
            params.zfn_left_half_site,
            params.zfn_right_half_site,
            zfn_genome,
            params.zfn_algorithm,
            params.zfn_dimer_mode,
            params.zfn_spacer_lengths,
            params.zfn_max_mismatches,
            zfn_annot,
            params.zfn_sharding_enabled,
            params.zfn_shard_chunk_mb,
            params.zfn_shard_overlap_bp,
            params.zfn_shard_chromosomes
        )
    } else {
        //
        // siRNA / miRNA mode: BWA-MEM2 alignment + miRNA seed analysis
        //
        SIRNA_OFFTARGET_ANALYSIS(
            ch_input,
            ch_references,
            params.max_hits,
            params.bwa_k,
            params.bwa_T,
            params.seed_start,
            params.seed_end
        )
    }

    emit:
    // Outputs depend on design_mode; consumers must check which are populated
    combined_analyses    = params.design_mode != 'zfn' ? SIRNA_OFFTARGET_ANALYSIS.out.combined_analyses : channel.empty()
    combined_summary     = params.design_mode != 'zfn' ? SIRNA_OFFTARGET_ANALYSIS.out.combined_summary : channel.empty()
    final_summary        = params.design_mode != 'zfn' ? SIRNA_OFFTARGET_ANALYSIS.out.final_summary : channel.empty()
    zfn_sites            = params.design_mode == 'zfn' ? ZFN_OFFTARGET_ANALYSIS.out.sites : channel.empty()
    zfn_summary          = params.design_mode == 'zfn' ? ZFN_OFFTARGET_ANALYSIS.out.summary : channel.empty()
    versions             = params.design_mode == 'zfn' ? ZFN_OFFTARGET_ANALYSIS.out.versions : SIRNA_OFFTARGET_ANALYSIS.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    SIRNAFORGE_OFFTARGET()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
