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
    expected_species   // val: comma-separated species this run expects to screen (#100). Carried
                        // independently of `references`/the built-index channel below: deriving
                        // "expected" from which indices actually built loses a species whose
                        // BUILD_BWA_INDEX crashed instead of letting it reconcile as failed.
    evidence_plan       // path: serialized #100 ScreeningPlan JSON, or the NO_EVIDENCE_PLAN sentinel

    main:
    ch_versions = channel.empty()

    //
    // #100: gate the miRNA channel on a non-empty *resolved* species list rather than running it
    // unconditionally. Default behaviour is unchanged (the elvis fallback below still fires when
    // params.mirna_species is null or unset -- Nextflow's own CLI parsing collapses a
    // whitespace-only value to the same falsy empty string, so blank does not defeat it either).
    // A caller that wants the channel off must supply a value that is non-empty (so elvis leaves
    // it alone) but resolves to zero species after splitting -- e.g. a bare comma.
    //
    def ch_mirna_species_list = (params.mirna_species ?: 'chicken,pig,rat,mouse,human,macaque')
        .split(',')
        .collect { it.trim() }
        .findAll { it }

    ch_mirna_analysis = channel.empty()
    ch_mirna_summary = channel.empty()
    // MIRNA_SEED_ANALYSIS.out.evidence is a glob ("mirna_seed_*_evidence.json") emitted as one
    // List<Path> per batch invocation; flatten() so it mixes with OFFTARGET_ANALYSIS's one-file-
    // per-species evidence emission as individual paths, not a list nested inside a list.
    ch_mirna_evidence = channel.empty()

    if (ch_mirna_species_list) {
        //
        // MODULE: miRNA seed match analysis (lightweight, <4GB RAM), only when species were
        // actually requested. Efficient batch mode: one process for all candidates.
        //
        MIRNA_SEED_ANALYSIS(
            candidates_fasta,
            params.mirna_db ?: 'mirgenedb',
            ch_mirna_species_list.join(',')
        )
        ch_versions = ch_versions.mix(MIRNA_SEED_ANALYSIS.out.versions)
        ch_mirna_analysis = MIRNA_SEED_ANALYSIS.out.analysis
        ch_mirna_summary = MIRNA_SEED_ANALYSIS.out.summary
        ch_mirna_evidence = MIRNA_SEED_ANALYSIS.out.evidence.flatten()
    }

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
        .mix(ch_mirna_analysis)
        .collect()

    ch_all_summary = OFFTARGET_ANALYSIS.out.summary
        .mix(ch_mirna_summary)
        .collect()

    // Every per-unit #100 evidence envelope, from both channels. Staged as a `path` input (not
    // `val`, unlike analysis/summary above) so AGGREGATE_RESULTS's own reconciliation -- which
    // globs its OWN task working directory -- actually finds them, rather than reconciling
    // against an empty local directory and marking every species failed regardless of outcome.
    ch_all_evidence = OFFTARGET_ANALYSIS.out.evidence
        .mix(ch_mirna_evidence)
        .collect()

    //
    // MODULE: Aggregate all results
    //
    // transcriptome_species is `expected_species` (#100), not something derived from
    // ch_reference_indices: the aggregator must reconcile every species this run intended to
    // screen, including one whose BUILD_BWA_INDEX crashed and therefore never reached
    // ch_reference_indices at all.
    //
    AGGREGATE_RESULTS(
        ch_all_analysis,
        ch_all_summary,
        expected_species,
        evidence_plan,
        ch_all_evidence
    )
    ch_versions = ch_versions.mix(AGGREGATE_RESULTS.out.versions)

    emit:
    combined_analyses    = AGGREGATE_RESULTS.out.combined_analyses
    combined_summary     = AGGREGATE_RESULTS.out.combined_summary
    final_summary        = AGGREGATE_RESULTS.out.final_summary
    versions            = ch_versions
}
