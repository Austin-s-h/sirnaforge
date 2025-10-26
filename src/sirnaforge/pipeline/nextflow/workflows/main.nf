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
    // Validate required parameters (simplified - Python handles detailed validation)
    //
    if (!params.input) {
        error "Input FASTA file must be specified with --input"
    }
    if (!file(params.input).exists()) {
        error "Input file does not exist: ${params.input}"
    }

    //
    // Create input channels
    //

    // Input siRNA candidates with metadata
    ch_input = Channel.value([
        [id: file(params.input).simpleName, single_end: false],
        file(params.input, checkIfExists: true)
    ])

    // Genome configurations: parse both FASTAs and indices
    // NOTE: Genome/transcriptome analysis is OPTIONAL (resource-intensive: 8-60GB RAM)
    // If not provided, only lightweight miRNA seed match analysis will run
    ch_genomes = Channel.empty()

    // Parse genome FASTAs (for building indices)
    if (params.genome_fastas) {
        ch_fasta_genomes = Channel.from(params.genome_fastas.split(','))
            .map { entry ->
                def (species, fasta_path) = entry.split(':')
                [species.trim(), file(fasta_path.trim(), checkIfExists: true), 'fasta']
            }
        ch_genomes = ch_genomes.mix(ch_fasta_genomes)
        log.info "Genome FASTAs provided - enabling resource-intensive off-target analysis"
    }

    // Parse pre-built indices
    if (params.genome_indices) {
        ch_index_genomes = Channel.from(params.genome_indices.split(','))
            .map { entry ->
                def (species, index_path) = entry.split(':')
                [species.trim(), index_path.trim(), 'index']
            }
        ch_genomes = ch_genomes.mix(ch_index_genomes)
        log.info "Genome indices provided - enabling resource-intensive off-target analysis"
    }

    // If no genomes specified, only lightweight miRNA analysis will run
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
    mirna_results        = SIRNA_OFFTARGET_ANALYSIS.out.mirna_results
    individual_results   = SIRNA_OFFTARGET_ANALYSIS.out.individual_results
    combined_analyses    = SIRNA_OFFTARGET_ANALYSIS.out.combined_analyses
    combined_summary     = SIRNA_OFFTARGET_ANALYSIS.out.combined_summary
    final_summary        = SIRNA_OFFTARGET_ANALYSIS.out.final_summary
    html_report         = SIRNA_OFFTARGET_ANALYSIS.out.html_report
    validation_report    = SIRNA_OFFTARGET_ANALYSIS.out.validation_report
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
