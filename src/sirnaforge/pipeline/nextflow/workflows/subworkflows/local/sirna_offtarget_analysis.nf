/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SIRNA OFF-TARGET ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_BWA_INDEX     } from '../../modules/local/build_bwa_index'
include { MIRNA_SEED_ANALYSIS } from '../../modules/local/mirna_seed_analysis'
include { OFFTARGET_ANALYSIS  } from '../../modules/local/offtarget_analysis'
include { AGGREGATE_RESULTS   } from '../../modules/local/aggregate_results'

workflow SIRNA_OFFTARGET_ANALYSIS {
    take:
    candidates_fasta    // path: input FASTA file
    genomes             // channel: [species, path_or_null, type] where type is 'fasta', 'index'
    max_hits           // val: maximum hits per candidate
    bwa_k              // val: BWA seed length
    bwa_T              // val: BWA minimum score threshold
    seed_start         // val: seed region start
    seed_end           // val: seed region end

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: ALWAYS run miRNA seed match analysis (lightweight, <4GB RAM)
    // Efficient batch mode: one process for all candidates
    //
    MIRNA_SEED_ANALYSIS(
        candidates_fasta,
        params.mirna_db ?: 'mirgenedb',
        params.mirna_species ?: 'chicken,pig,rat,mouse,human,macaque'
    )
    ch_versions = ch_versions.mix(MIRNA_SEED_ANALYSIS.out.versions)

    //
    // CONDITIONAL: Genome/transcriptome off-target analysis
    // Efficient pattern: One alignment session per genome, all candidates processed sequentially
    //
    ch_genome_indices = Channel.empty()

    // Build BWA indices for FASTA files if provided
    genomes
        .filter { species, path, type -> type == 'fasta' }
        .map { species, path, type -> [species, path] }
        .set { ch_genome_fastas }

    if (ch_genome_fastas) {
        BUILD_BWA_INDEX(ch_genome_fastas)
        ch_versions = ch_versions.mix(BUILD_BWA_INDEX.out.versions)

        // Add built indices to channel
        ch_genome_indices = ch_genome_indices.mix(
            BUILD_BWA_INDEX.out.index
                .map { species, index_files ->
                    def index_prefix = index_files[0].toString().replaceAll(/\.[^.]+$/, '')
                    [species, index_prefix]
                }
        )
    }

    // Use existing indices
    ch_genome_indices = ch_genome_indices.mix(
        genomes
            .filter { species, path, type -> type == 'index' }
            .map { species, index_path, type -> [species, index_path] }
    )

    //
    // EFFICIENT PATTERN: One analysis session per genome with all candidates
    // Instead of candidate × genome combinations (e.g., 100 × 3 = 300 processes),
    // we run 3 processes (one per genome), each processing all 100 candidates sequentially
    //
    // Combine genome indices with the candidates FASTA file
    ch_genome_analysis_input = ch_genome_indices
        .combine(candidates_fasta)
        .map { species, index_path, fasta_file ->
            [species, index_path, fasta_file]
        }

    //
    // MODULE: Run off-target analysis once per genome (all candidates in batch)
    //
    OFFTARGET_ANALYSIS(
        ch_genome_analysis_input,
        max_hits,
        bwa_k,
        bwa_T,
        seed_start,
        seed_end
    )
    ch_versions = ch_versions.mix(OFFTARGET_ANALYSIS.out.versions)

    //
    // Collect all analysis results for aggregation
    //
    ch_genome_analysis_files = OFFTARGET_ANALYSIS.out.analysis.toList().ifEmpty { [] }
    ch_genome_summary_files  = OFFTARGET_ANALYSIS.out.summary.toList().ifEmpty { [] }

    ch_mirna_analysis_files = MIRNA_SEED_ANALYSIS.out.analysis.toList()
    ch_mirna_summary_files  = MIRNA_SEED_ANALYSIS.out.summary.toList()

    // Build a SINGLE list of staged inputs so AGGREGATE_RESULTS runs once
    // (Using mix() here would create multiple invocations and can overwrite published outputs.)
    ch_all_analysis = ch_genome_analysis_files
        .combine(ch_mirna_analysis_files)
        .map { genome_files, mirna_files -> genome_files + mirna_files }

    ch_all_summary = ch_genome_summary_files
        .combine(ch_mirna_summary_files)
        .map { genome_files, mirna_files -> genome_files + mirna_files }

    // Extract species list for aggregation
    ch_genome_species = ch_genome_indices
        .map { species, index_path -> species }
        .unique()
        .toList()
        .map { species_list -> species_list.join(',') }
        .ifEmpty { '' }

    //
    // MODULE: Aggregate all results
    //
    AGGREGATE_RESULTS(
        ch_all_analysis,
        ch_all_summary,
        ch_genome_species
    )
    ch_versions = ch_versions.mix(AGGREGATE_RESULTS.out.versions)

    emit:
    combined_analyses    = AGGREGATE_RESULTS.out.combined_analyses
    combined_summary     = AGGREGATE_RESULTS.out.combined_summary
    final_summary        = AGGREGATE_RESULTS.out.final_summary
    versions            = ch_versions
}
