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
    references          // channel: [species, path_or_prefix, type] where type is 'fasta' or 'index'
    max_hits           // val: maximum hits per candidate
    bwa_k              // val: BWA seed length
    bwa_T              // val: BWA minimum score threshold
    seed_start         // val: seed region start
    seed_end           // val: seed region end

    main:
    ch_versions = channel.empty()

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
    // CONDITIONAL: transcriptome off-target analysis
    // Efficient pattern: one alignment session per reference, all candidates processed sequentially
    //
    ch_reference_indices = channel.empty()

    // Build BWA indices for FASTA files if provided
    references
        .filter { _species, _path, type -> type == 'fasta' }
        .map { species, path, _type -> [species, path] }
        .set { ch_reference_fastas }

    if (ch_reference_fastas) {
        BUILD_BWA_INDEX(ch_reference_fastas)
        ch_versions = ch_versions.mix(BUILD_BWA_INDEX.out.versions)

        // Add built indices to channel
        ch_reference_indices = ch_reference_indices.mix(
            BUILD_BWA_INDEX.out.index
                .map { species, index_files ->
                    def index_prefix = index_files[0].toString().replaceAll(/\.[^.]+$/, '')
                    [species, index_prefix]
                }
        )
    }

    // Use existing indices
    ch_reference_indices = ch_reference_indices.mix(
        references
            .filter { _species, _path, type -> type == 'index' }
            .map { species, index_path, _type -> [species, index_path] }
    )

    //
    // EFFICIENT PATTERN: One analysis session per reference with all candidates
    // Instead of candidate × reference combinations (e.g., 100 × 3 = 300 processes),
    // we run 3 processes (one per reference), each processing all 100 candidates sequentially
    //
    // Combine the reference indices with the candidates FASTA file
    ch_reference_analysis_input = ch_reference_indices
        .combine(candidates_fasta)
        .map { species, index_path, fasta_file ->
            [species, index_path, fasta_file]
        }

    //
    // MODULE: Run off-target analysis once per reference (all candidates in batch)
    //
    OFFTARGET_ANALYSIS(
        ch_reference_analysis_input,
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
    // Aggregate once after all upstream analyses complete.
    // mix(...).collect() keeps the shape simple for miRNA-only and transcriptome+miRNA runs alike.
    ch_all_analysis = OFFTARGET_ANALYSIS.out.analysis
        .mix(MIRNA_SEED_ANALYSIS.out.analysis)
        .collect()

    ch_all_summary = OFFTARGET_ANALYSIS.out.summary
        .mix(MIRNA_SEED_ANALYSIS.out.summary)
        .collect()

    // Extract species list for aggregation
    ch_screened_species = ch_reference_indices
        .map { species, _index_path -> species }
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
        ch_screened_species
    )
    ch_versions = ch_versions.mix(AGGREGATE_RESULTS.out.versions)

    emit:
    combined_analyses    = AGGREGATE_RESULTS.out.combined_analyses
    combined_summary     = AGGREGATE_RESULTS.out.combined_summary
    final_summary        = AGGREGATE_RESULTS.out.final_summary
    versions            = ch_versions
}
