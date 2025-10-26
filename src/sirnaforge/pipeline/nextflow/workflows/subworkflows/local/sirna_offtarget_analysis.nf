/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SIRNA OFF-TARGET ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Comprehensive off-target analysis with parallel processing per candidate per genome
*/

include { PREPARE_CANDIDATES  } from '../../modules/local/prepare_candidates'
include { SPLIT_CANDIDATES    } from '../../modules/local/split_candidates'
include { BUILD_BWA_INDEX     } from '../../modules/local/build_bwa_index'
include { MIRNA_SEED_ANALYSIS } from '../../modules/local/mirna_seed_analysis'
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
    // MODULE: ALWAYS run miRNA seed match analysis (lightweight, <1GB RAM)
    //
    SPLIT_CANDIDATES.out.individual_candidates
        .map { meta, candidate_file ->
            def candidate_id = meta?.id ?: candidate_file.simpleName.replaceAll(/candidate_/, '').replaceAll(/\.fasta$/, '')
            [[id: candidate_id, file: candidate_file.name], candidate_file]
        }
        .set { ch_individual_candidates }

    MIRNA_SEED_ANALYSIS(
        ch_individual_candidates,
        params.mirna_db ?: 'mirgenedb',
        params.mirna_species ?: 'human'
    )
    ch_versions = ch_versions.mix(MIRNA_SEED_ANALYSIS.out.versions)

    //
    // CONDITIONAL: Genome/transcriptome off-target analysis
    // Only run if user provides genome FASTAs or indices (resource-intensive: 8-60GB RAM)
    //
    ch_genome_indices = Channel.empty()

    // Build BWA indices for FASTA files
    genomes
        .filter { species, path, type -> type == 'fasta' }
        .map { species, path, type -> [species, path] }
        .set { ch_genome_fastas }

    BUILD_BWA_INDEX(ch_genome_fastas)
    ch_versions = ch_versions.mix(BUILD_BWA_INDEX.out.versions.ifEmpty([]))

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

    // Log genome analysis status
    ch_genome_indices
        .count()
        .subscribe { count ->
            if (count > 0) {
                log.info "Genome/transcriptome off-target analysis: ENABLED (${count} genome(s), resource-intensive: 8-60GB RAM)"
            } else {
                log.info "Genome/transcriptome off-target analysis: SKIPPED (no --genome_fastas or --genome_indices provided)"
                log.info "Running lightweight miRNA seed match analysis only"
            }
        }

    //
    // Create combinations for parallel processing: each candidate x each genome
    //
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
    // Collect all analysis results (gracefully handles empty channel if no genomes)
    //
    OFFTARGET_ANALYSIS.out.results
        .map { candidate_meta, species, analysis_type, analysis_file, summary_file ->
            [analysis_file, summary_file]
        }
        .collect()
        .ifEmpty([[],[]])
        .set { ch_all_results }

    //
    // Extract species list for aggregation (empty if no genomes)
    //
    ch_genome_indices
        .map { species, index_path, index_type -> species }
        .unique()
        .collect()
        .map { species_list -> species_list.join(',') }
        .ifEmpty('')
        .set { ch_genome_species_list }

    //
    // MODULE: Aggregate all results (only runs if genome analysis was performed)
    //
    AGGREGATE_RESULTS(
        ch_all_results.map { it[0] }.flatten().ifEmpty([]),  // analysis files
        ch_all_results.map { it[1] }.flatten().ifEmpty([]),  // summary files
        ch_genome_species_list
    )
    ch_versions = ch_versions.mix(AGGREGATE_RESULTS.out.versions.ifEmpty([]))

    emit:
    // miRNA analysis results (always present)
    mirna_results        = MIRNA_SEED_ANALYSIS.out.results

    // Individual results for detailed analysis (empty if no genome analysis)
    individual_results   = OFFTARGET_ANALYSIS.out.results.ifEmpty([])

    // Aggregated results (empty if no genome analysis)
    combined_analyses    = AGGREGATE_RESULTS.out.combined_analyses.ifEmpty([])
    combined_summary     = AGGREGATE_RESULTS.out.combined_summary.ifEmpty([])
    final_summary        = AGGREGATE_RESULTS.out.final_summary.ifEmpty([])
    html_report         = AGGREGATE_RESULTS.out.html_report.ifEmpty([])

    // Validation and metadata
    validation_report    = PREPARE_CANDIDATES.out.candidates.map { it[2] }
    candidate_manifest   = SPLIT_CANDIDATES.out.manifest

    // Versions
    versions            = ch_versions
}
