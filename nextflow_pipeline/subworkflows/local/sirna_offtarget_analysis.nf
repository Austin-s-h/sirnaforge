/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SIRNA OFF-TARGET ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Comprehensive off-target analysis with parallel processing per candidate per genome
*/

include { PREPARE_CANDIDATES  } from '../../modules/local/prepare_candidates'
include { SPLIT_CANDIDATES    } from '../../modules/local/split_candidates'
include { BUILD_BWA_INDEX     } from '../../modules/local/build_bwa_index'
include { OFFTARGET_ANALYSIS  } from '../../modules/local/offtarget_analysis'
include { AGGREGATE_RESULTS   } from '../../modules/local/aggregate_results'

workflow SIRNA_OFFTARGET_ANALYSIS {
    take:
    candidates_fasta    // tuple: [meta, fasta_file]
    genomes             // channel: [species, path_or_null, type] where type is 'fasta', 'index', or 'discover'
    max_hits           // val: maximum hits per candidate
    bwa_k              // val: BWA seed length
    bwa_T              // val: BWA minimum score threshold
    seed_start         // val: seed region start
    seed_end           // val: seed region end

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Validate and prepare siRNA candidates
    //
    PREPARE_CANDIDATES(candidates_fasta)
    ch_versions = ch_versions.mix(PREPARE_CANDIDATES.out.versions)

    //
    // MODULE: Split candidates for parallel processing
    //
    SPLIT_CANDIDATES(PREPARE_CANDIDATES.out.candidates)
    ch_versions = ch_versions.mix(SPLIT_CANDIDATES.out.versions)

    //
    // Prepare genome indices (build or use existing)
    //
    ch_genome_indices = Channel.empty()

    // Build BWA indices for FASTA files
    genomes
        .filter { species, path, type -> type == 'fasta' }
        .map { species, path, type -> [species, path] }
        .set { ch_genome_fastas }

    BUILD_BWA_INDEX(ch_genome_fastas)
    ch_versions = ch_versions.mix(BUILD_BWA_INDEX.out.versions)

    // Add built indices to channel
    BUILD_BWA_INDEX.out.index
        .map { species, index_files ->
            def index_prefix = index_files[0].toString().replaceAll(/\.[^.]+$/, '')
            [species, index_prefix, 'bwa']
        }
        .set { ch_built_bwa_indices }
    ch_genome_indices = ch_genome_indices.mix(ch_built_bwa_indices)

    // Use existing indices
    genomes
        .filter { species, path, type -> type == 'index' }
        .map { species, index_path, type -> [species, index_path, 'bwa'] }
        .set { ch_existing_indices }
    ch_genome_indices = ch_genome_indices.mix(ch_existing_indices)

    //
    // Create combinations for parallel processing: each candidate x each genome
    //
    SPLIT_CANDIDATES.out.individual_candidates
        .flatten()
        .map { file ->
            // Extract candidate metadata from filename
            def candidate_id = file.simpleName.replaceAll(/candidate_/, '').replaceAll(/\.fasta$/, '')
            [[id: candidate_id, file: file.name], file]
        }
        .set { ch_individual_candidates }

    ch_individual_candidates
        .combine(ch_genome_indices)
        .map { candidate_meta, candidate_file, species, index_path, index_type ->
            [candidate_meta, candidate_file, species, index_path, index_type]
        }
        .set { ch_analysis_combinations }

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
    // Collect all analysis results
    //
    OFFTARGET_ANALYSIS.out.results
        .map { candidate_meta, species, analysis_type, analysis_file, summary_file ->
            [analysis_file, summary_file]
        }
        .collect()
        .set { ch_all_results }

    //
    // Extract species list for aggregation
    //
    ch_genome_indices
        .map { species, index_path, index_type -> species }
        .unique()
        .collect()
        .map { species_list -> species_list.join(',') }
        .set { ch_genome_species_list }

    //
    // MODULE: Aggregate all results
    //
    AGGREGATE_RESULTS(
        ch_all_results.map { it[0] }.flatten(),  // analysis files
        ch_all_results.map { it[1] }.flatten(),  // summary files
        ch_genome_species_list
    )
    ch_versions = ch_versions.mix(AGGREGATE_RESULTS.out.versions)

    emit:
    // Individual results for detailed analysis
    individual_results   = OFFTARGET_ANALYSIS.out.results

    // Aggregated results
    combined_analyses    = AGGREGATE_RESULTS.out.combined_analyses
    combined_summary     = AGGREGATE_RESULTS.out.combined_summary
    final_summary        = AGGREGATE_RESULTS.out.final_summary
    html_report         = AGGREGATE_RESULTS.out.html_report

    // Validation and metadata
    validation_report    = PREPARE_CANDIDATES.out.candidates.map { it[2] }
    candidate_manifest   = SPLIT_CANDIDATES.out.manifest

    // Versions
    versions            = ch_versions
}
