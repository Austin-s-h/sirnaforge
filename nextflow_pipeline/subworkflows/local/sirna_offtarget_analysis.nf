/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SIRNA OFF-TARGET ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Simplified off-target analysis with native Nextflow operators
*/

include { BUILD_BWA_INDEX     } from '../../modules/local/build_bwa_index'
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
    // Split FASTA using native Nextflow operator - no separate process needed!
    //
    ch_individual_candidates = candidates_fasta
        .splitFasta(record: [id: true, seqString: true])
        .map { record ->
            def fasta_file = "candidate_${record.id}.fasta"
            def fasta_content = ">${record.id}\n${record.seqString}\n"
            [record.id, fasta_content]
        }

    //
    // Prepare genome indices (build or use existing)
    //
    ch_genome_indices = Channel.empty()

    // Build BWA indices for FASTA files
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
    // Create combinations for parallel processing: each candidate x each genome
    //
    ch_analysis_combinations = ch_individual_candidates
        .combine(ch_genome_indices)
        .map { candidate_id, fasta_content, species, index_path ->
            [candidate_id, fasta_content, species, index_path]
        }

    //
    // MODULE: Run off-target analysis for each candidate-genome combination
    //
    OFFTARGET_ANALYSIS(
        ch_analysis_combinations,
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
    ch_all_analysis = OFFTARGET_ANALYSIS.out.analysis.collect()
    ch_all_summary = OFFTARGET_ANALYSIS.out.summary.collect()

    // Extract species list for aggregation
    ch_genome_species = ch_genome_indices
        .map { species, index_path -> species }
        .unique()
        .collect()
        .map { species_list -> species_list.join(',') }

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
