process OFFTARGET_ANALYSIS {
    tag "$species"
    label 'process_medium'
    publishDir "${params.outdir}/transcriptome", mode: params.publish_dir_mode

    input:
    tuple val(species), val(index_path), path(candidates_fasta)
    val max_hits
    val bwa_k
    val bwa_T
    val seed_start
    val seed_end

    output:
    path "${species}_analysis.tsv", emit: analysis
    path "${species}_summary.json", emit: summary
    // #100: non-optional evidence envelope. offtarget_analysis_cli already writes this file on
    // both its branches (missing-index failure and completed alignment); declaring it as a real
    // output is what makes AGGREGATE_RESULTS's reconciliation actually see it.
    path "transcriptome_${species}_evidence.json", emit: evidence
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run off-target analysis for ALL candidates against this reference in one session
    # This is much more efficient: load index once, process all candidates sequentially.
    # The CLI also owns the "this prefix is not a usable index" case: it publishes an EMPTY
    # analysis file plus a failed summary, which the aggregator reports as a per-species
    # rejection instead of a completed screen with zero hits.
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import offtarget_analysis_cli

result = offtarget_analysis_cli(
    species='${species}',
    index_prefix='${index_path}',
    candidates_file='${candidates_fasta}',
    output_dir='.',
    max_hits=${max_hits},
    bwa_k=${bwa_k},
    bwa_T=${bwa_T},
    seed_start=${seed_start},
    seed_end=${seed_end}
)

print(f"Batch analysis for ${species}: {result['status']}")
if result.get('error'):
    print(f"  reason: {result['error']}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """

    stub:
    """
    # Deliberately empty, not a header: a stub aligned nothing, so the aggregator must report this
    # species as unscreened rather than as screened and clean.
    touch ${species}_analysis.tsv
    echo '{"species": "${species}", "status": "stub", "total_candidates": 0, "total_hits": 0}' > ${species}_summary.json

    # #100: the evidence envelope this species would carry from a real run, mirrored here so a
    # `-stub-run` still satisfies the non-optional output declared above. status "failed" and
    # producer "stub" -- nothing actually aligned, so this must never read as "complete".
    cat <<-'EVIDENCE' > transcriptome_${species}_evidence.json
    {
      "schema_version": "2",
      "producer": "stub",
      "source": "synthesized",
      "entry": {
        "channel": "transcriptome",
        "species": "${species}",
        "reference_id": null,
        "guide_set_digest": "0000000000000000",
        "status": "failed",
        "counts": {
          "sites": {"value": null, "is_lower_bound": false, "cap": null, "truncated": false},
          "distinct_transcripts": {"value": null, "is_lower_bound": false, "cap": null, "truncated": false},
          "distinct_genes": {"value": null, "is_lower_bound": false, "cap": null, "truncated": false},
          "unresolved_gene_sites": {"value": null, "is_lower_bound": false, "cap": null, "truncated": false}
        },
        "submitted_guide_digest": null,
        "submitted_guides": null,
        "processed_guides": null,
        "detail": "stub run: no aligner executed"
      }
    }
    EVIDENCE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1 | sed 's/.*bwa-mem2-//' || echo 'not available')
    END_VERSIONS
    """
}
