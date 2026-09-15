/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SIRNA OFF-TARGET ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_BWA_INDEX          } from '../../modules/local/build_bwa_index'
include { MIRNA_SEED_ANALYSIS      } from '../../modules/local/mirna_seed_analysis'
include { OFFTARGET_ANALYSIS       } from '../../modules/local/offtarget_analysis'
include { TRANSCRIPT_SEED_ANALYSIS } from '../../modules/local/transcript_seed_analysis'
include { AGGREGATE_RESULTS        } from '../../modules/local/aggregate_results'

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
    // #101: transcript-seed liability, the third screening channel, OPT-IN and off by default.
    //
    // Off by default because it is not calibrated and it is enormous: a 7mer occurs roughly once per
    // 16 kb, so an uncapped full-cDNA scan publishes orders of magnitude more rows than the alignment
    // table. A caller asks for it with --transcript_seed_enabled.
    //
    // The value is compared against a literal allow-list rather than tested for truthiness, because
    // Nextflow hands `--transcript_seed_enabled false` to Groovy as the non-empty String "false" --
    // which is truthy. A truthiness test would therefore turn the documented way of asking for the
    // channel OFF into the way of switching it ON, which is precisely the silent-on failure a
    // default-off channel exists to avoid.
    //
    def transcript_seed_enabled = "${params.transcript_seed_enabled ?: ''}".trim().toLowerCase() in ['true', '1', 'yes', 'on']

    // Which species to scan. Defaults to `expected_species` -- the run's own declared screen list --
    // so asking for the channel does not also require restating the species. Resolved to a list here
    // and gated on it being non-empty, exactly as the miRNA channel is: a request that resolves to
    // zero species must run nothing rather than run everything.
    def ch_transcript_seed_species_list = "${params.transcript_seed_species ?: expected_species ?: ''}"
        .split(',')
        .collect { it.trim() }
        .findAll { it }

    ch_transcript_seed_sites = channel.empty()
    ch_transcript_seed_summary = channel.empty()
    ch_transcript_seed_evidence = channel.empty()

    if (transcript_seed_enabled && ch_transcript_seed_species_list) {
        //
        // MODULE: transcript-seed site scan, one process per species, all candidates in batch.
        //
        // Takes `ch_reference_fastas` -- the SAME already-materialised cDNA FASTAs BUILD_BWA_INDEX
        // indexes -- and `candidates_fasta`, the SAME deduplicated guide set the alignment channel
        // screens. Nothing is downloaded and nothing is re-indexed: the scan is a substring search
        // over a file that is already staged. DSL2 forks a channel for each consumer, so reading
        // ch_reference_fastas here does not take those entries away from BUILD_BWA_INDEX.
        //
        // A species supplied as `transcriptome_indices` (a prebuilt prefix, no FASTA) is absent from
        // this channel and is therefore NOT scanned -- there is no sequence to scan. That species
        // publishes no envelope, so its plan entry reconciles FAILED and the transcript-seed gates
        // report UNKNOWN for it. An index is not a reference the scan can read, and saying so is the
        // whole point of the evidence contract.
        //
        ch_transcript_seed_input = ch_reference_fastas
            .filter { species, _path -> ch_transcript_seed_species_list.contains(species) }
            .combine(candidates_fasta)

        TRANSCRIPT_SEED_ANALYSIS(
            ch_transcript_seed_input,
            params.transcript_seed_scope ?: 'full_cdna',
            // Interpolated straight into Python, so it must always be an integer: `null` would
            // template as a syntax error. 0 or below means uncapped (see the CLI), because a
            // Nextflow `val` cannot carry Python's None.
            params.transcript_seed_max_sites_per_guide ?: 200
        )
        ch_versions = ch_versions.mix(TRANSCRIPT_SEED_ANALYSIS.out.versions)
        ch_transcript_seed_sites = TRANSCRIPT_SEED_ANALYSIS.out.sites
        ch_transcript_seed_summary = TRANSCRIPT_SEED_ANALYSIS.out.summary
        ch_transcript_seed_evidence = TRANSCRIPT_SEED_ANALYSIS.out.evidence
    }

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
    // #100 (defect 11): a run that can screen nothing must not exit 0.
    //
    // With the miRNA channel deliberately off (a value that resolves to zero species, the idiom
    // documented above) AND no alignment unit to run, every channel collected below is empty.
    // `collect()` emits nothing at all for an empty source -- unlike `toList()`, which emits an
    // empty list -- so AGGREGATE_RESULTS never ran, no reconciliation was ever published, and the
    // pipeline reported success having screened nothing: the exact failure this evidence contract
    // exists to prevent, reachable through the subworkflow's own disable idiom.
    //
    // Two changes, and both are needed:
    //
    //  * the three `collect()` calls below get an `ifEmpty([])` floor, so AGGREGATE_RESULTS always
    //    runs and always publishes a reconciliation -- including on the paths where a reference
    //    WAS configured but no analysis unit ran (every entry filtered out by type, every
    //    BUILD_BWA_INDEX crashed), where the aggregate has an expected plan to reconcile against
    //    and can therefore report each unit as failed;
    //
    //  * this guard aborts the run outright when BOTH channels are off, because the floor alone
    //    cannot make that case honest: with nothing requested the expected plan may itself be
    //    empty, and a reconciliation over an empty plan has no shortfall to report -- it is
    //    byte-comparable to a clean screen. A non-zero exit is the only outcome that cannot be
    //    mistaken for success, so that is what this configuration gets.
    //
    // Reassigned rather than forked off `references`, so the check adds no second consumer of an
    // existing channel: this is the one channel whose emptiness means "not one unit will align".
    //
    if (!ch_mirna_species_list) {
        ch_reference_analysis_input = ch_reference_analysis_input.ifEmpty {
            error(
                "Nothing to screen: --mirna_species resolved to zero species and no usable transcriptome " +
                "reference was configured, so neither screening channel can run and no evidence can be " +
                "published. Provide --transcriptome_fastas or --transcriptome_indices, or leave " +
                "--mirna_species at its default."
            )
        }
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
    // ifEmpty([]) is the #100 floor described above: an empty collect() emits nothing, which
    // silently withheld AGGREGATE_RESULTS altogether; an empty list runs it on nothing, which is
    // what publishes "these units were expected and none of them produced anything".
    ch_all_analysis = OFFTARGET_ANALYSIS.out.analysis
        .mix(ch_mirna_analysis)
        .collect()
        .ifEmpty([])

    ch_all_summary = OFFTARGET_ANALYSIS.out.summary
        .mix(ch_mirna_summary)
        .collect()
        .ifEmpty([])

    // Every per-unit #100 evidence envelope, from both channels. Staged as a `path` input (not
    // `val`, unlike analysis/summary above) so AGGREGATE_RESULTS's own reconciliation -- which
    // globs its OWN task working directory -- actually finds them, rather than reconciling
    // against an empty local directory and marking every species failed regardless of outcome.
    //
    // #101: the transcript-seed envelopes mix in HERE and nowhere else. Their evidence has to reach
    // the aggregate's reconciliation, but their sites/summary tables must NOT join ch_all_analysis
    // or ch_all_summary above: those two feed aggregate_offtarget_results, which sums transcriptome
    // and known-miRNA hits. A `<species>_transcript_seed_summary.json` staged there would be picked
    // up by that aggregate's own `*_summary.json` glob and folded into the miRNA counters, turning a
    // separate channel's site count into miRNA hits. The seed tables reach their reader through the
    // module's publishDir instead, where workflow.py's `**/*_transcript_seed_sites.tsv` parser finds
    // them on their own key. `aggregate_results_cli` additionally refuses them by name, so the two
    // channels' counters cannot be conflated from either side.
    //
    ch_all_evidence = OFFTARGET_ANALYSIS.out.evidence
        .mix(ch_mirna_evidence)
        .mix(ch_transcript_seed_evidence)
        .collect()
        .ifEmpty([])

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
    // #101: emitted as their own outputs, never merged into combined_analyses. The transcript-seed
    // table has no cigar, no mapq, no nm and 1-based transcript coordinates, so it is a different
    // artifact rather than more rows of the alignment one. Empty channels when the channel is off.
    transcript_seed_sites   = ch_transcript_seed_sites
    transcript_seed_summary = ch_transcript_seed_summary
    versions            = ch_versions
}
