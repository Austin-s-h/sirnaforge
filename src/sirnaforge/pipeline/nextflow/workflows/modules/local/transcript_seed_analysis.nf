process TRANSCRIPT_SEED_ANALYSIS {
    tag "$species"
    label 'process_low'
    publishDir "${params.outdir}/transcript_seed", mode: params.publish_dir_mode

    // #101's third liability channel: complementary seed sites in transcript sequence, which is a
    // different question from resemblance to a known miRNA. It runs here, beside OFFTARGET_ANALYSIS,
    // because the cDNA FASTA the alignment channel already materialised is the only input it needs --
    // the same staged `path`, never a second download and never a second index. `label 'process_low'`
    // and no aligner for the same reason: this is a substring search over a file that is already on
    // disk, so it costs one streaming pass, not an alignment session.

    input:
    // The species key, the ALREADY-MATERIALISED cDNA FASTA and the SAME deduplicated candidates
    // FASTA the alignment channel screens. `path` on both files, so Nextflow stages the very files
    // upstream resolved rather than this process resolving anything of its own.
    tuple val(species), path(cdna_fasta), path(candidates_fasta)
    val region_scope
    val max_sites_per_guide

    output:
    // Globs, not "${species}_...", deliberately: `transcript_seed_analysis_cli` names its outputs
    // after the CANONICAL species (normalize_species_name), so a run keyed on 'homo_sapiens' or
    // 'hsa' publishes human_* files. An interpolated output declaration would then fail to match
    // its own task's real output and abort a run that had in fact screened correctly. The miRNA
    // module's evidence glob exists for the same reason.
    path "*_transcript_seed_sites.tsv", emit: sites
    path "*_transcript_seed_summary.json", emit: summary
    // #100 envelope, non-optional: transcript_seed_analysis_cli writes one on every branch it can
    // reach (scanned, censored by the per-guide cap, and refused for want of an answerable scope or
    // a readable reference), so there is no path here that legitimately publishes none.
    path "transcript_seed_*_evidence.json", emit: evidence
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # One streaming pass over the staged cDNA for all candidates at once, the same batching the
    # miRNA and alignment modules use: the reference is the expensive thing, not the guides.
    python3 <<'PYEOF'
import sys
sys.path.insert(0, '${workflow.projectDir}/../src')
from sirnaforge.pipeline.nextflow_cli import transcript_seed_analysis_cli

result = transcript_seed_analysis_cli(
    species='${species}',
    cdna_fasta='${cdna_fasta}',
    candidates_file='${candidates_fasta}',
    output_dir='.',
    region_scope='${region_scope}',
    max_sites_per_guide=${max_sites_per_guide},
)

print(f"Transcript-seed scan for ${species}: {result['status']}")
print(f"  sites: {result['sites']} (cap {result['max_sites_per_guide']} per guide)")
if result.get('detail'):
    print(f"  detail: {result['detail']}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    # Deliberately empty, not a header: a header-only table is a scan that ran and found nothing,
    # and a stub scanned nothing at all. workflow.py's parser skips a 0-byte table for exactly this
    # reason, so the stub cannot be read back as a measured zero.
    touch stub_transcript_seed_sites.tsv
    echo '{"species": "${species}", "status": "stub", "sites": 0, "detail": "stub run: no reference was scanned"}' > stub_transcript_seed_summary.json

    # #100: mirror the non-optional evidence glob. Nothing was scanned, so status is "failed" and
    # producer "stub" -- never "complete" -- which keeps the (transcript_seed, species) pair out of
    # completed_pairs and makes the three transcript-seed ceilings report UNKNOWN rather than pass
    # on a zero a -stub-run invented.
    cat <<-'EVIDENCE' > transcript_seed_stub_evidence.json
    {
      "schema_version": "2",
      "producer": "stub",
      "source": "synthesized",
      "entry": {
        "channel": "transcript_seed",
        "species": "stub",
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
        "detail": "stub run: no reference was scanned"
      }
    }
    EVIDENCE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
