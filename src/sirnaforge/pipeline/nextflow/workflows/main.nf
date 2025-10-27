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

        GENOME ANALYSIS (OPTIONAL - Resource Intensive)
        genome_fastas        : ${params.genome_fastas ?: 'Not provided (miRNA-only mode)'}
        genome_indices       : ${params.genome_indices ?: 'Not provided'}
        genome_species       : ${params.genome_species}

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
    // Create input channel - simple file input
    //
    ch_input = Channel.fromPath(params.input, checkIfExists: true)

    //
    // Genome configurations: combine FASTAs and indices into single channel
    //
    ch_genomes = Channel.empty()

    if (params.genome_fastas) {
        ch_genomes = ch_genomes.mix(
            Channel.from(params.genome_fastas.split(','))
                .map { entry ->
                    def (species, fasta_path) = entry.split(':')
                    [species.trim(), file(fasta_path.trim(), checkIfExists: true), 'fasta']
                }
        )
    }

    if (params.genome_indices) {
        ch_genomes = ch_genomes.mix(
            Channel.from(params.genome_indices.split(','))
                .map { entry ->
                    def (species, index_path) = entry.split(':')
                    [species.trim(), index_path.trim(), 'index']
                }
        )
    }

    // If no genomes specified, skip genome analysis (miRNA-only mode)
    if (!params.genome_fastas && !params.genome_indices) {
        log.info "No genome FASTAs or indices provided"
        log.info "Genome/transcriptome off-target analysis: DISABLED"
        log.info "Running lightweight miRNA seed match analysis only (< 1GB RAM)"
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
