#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sirnaforge/nextflow_pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    siRNA Off-Target Analysis Pipeline
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
    log.info """
        ===============================================
         S I R N A F O R G E   O F F - T A R G E T
        ===============================================
        input                : ${params.input}
        outdir               : ${params.outdir}

        TRANSCRIPTOME ALIGNMENT (OPTIONAL - Resource Intensive)
        transcriptome_fastas (legacy --genome-fastas)  : ${params.genome_fastas ?: 'Not provided'}
        transcriptome_indices (legacy --genome-indices): ${params.genome_indices ?: 'Not provided'}
        supplemental_transcriptome_indices             : ${params.transcriptome_indices ?: 'Not provided'}

        MIRNA ANALYSIS
        transcriptome_species (legacy genome_species) : ${params.genome_species}

        max_hits             : ${params.max_hits}
        bwa_k                : ${params.bwa_k}
        bwa_T                : ${params.bwa_T}
        seed_start           : ${params.seed_start}
        seed_end             : ${params.seed_end}
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
    // Create input channel - simple file input
    //
    ch_input = Channel.fromPath(params.input, checkIfExists: true)

    //
    // Genome configurations: combine FASTAs and indices into single channel
    //
    ch_genomes = Channel.empty()

    if (params.genome_fastas) {
        ch_genomes = ch_genomes.mix(
            Channel.fromList(params.genome_fastas.split(','))
                .map { entry ->
                    def (species, fasta_path) = entry.split(':')
                    [species.trim(), file(fasta_path.trim(), checkIfExists: true), 'fasta']
                }
        )
    }

    if (params.genome_indices) {
        ch_genomes = ch_genomes.mix(
            Channel.fromList(params.genome_indices.split(','))
                .map { entry ->
                    def (species, index_path) = entry.split(':')
                    [species.trim(), index_path.trim(), 'index']
                }
        )
    }

    if (params.transcriptome_indices) {
        ch_genomes = ch_genomes.mix(
            Channel.fromList(params.transcriptome_indices.split(','))
                .map { entry ->
                    def (species, index_path) = entry.split(':')
                    [species.trim(), index_path.trim(), 'index']
                }
        )
    }

    def has_offtarget_data = params.genome_fastas || params.genome_indices || params.transcriptome_indices

    // If no references specified, skip transcriptome analysis (miRNA-only mode)
    if (!has_offtarget_data) {
        log.info "No transcriptome FASTAs, transcriptome indices, or supplemental indices provided"
        log.info "Transcriptome off-target analysis: DISABLED"
        log.info "Running lightweight miRNA seed match analysis only (< 1GB RAM)"
        log.info "Provide --genome_fastas (transcriptome), --genome_indices, or --transcriptome_indices to enable alignment"
    }

    //
    // SUBWORKFLOW: Run comprehensive off-target analysis
    //
    SIRNA_OFFTARGET_ANALYSIS(
        ch_input,
        ch_genomes,
        params.max_hits,
        params.bwa_k,
        params.bwa_T,
        params.seed_start,
        params.seed_end
    )

    emit:
    combined_analyses    = SIRNA_OFFTARGET_ANALYSIS.out.combined_analyses
    combined_summary     = SIRNA_OFFTARGET_ANALYSIS.out.combined_summary
    final_summary        = SIRNA_OFFTARGET_ANALYSIS.out.final_summary
    versions            = SIRNA_OFFTARGET_ANALYSIS.out.versions
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
