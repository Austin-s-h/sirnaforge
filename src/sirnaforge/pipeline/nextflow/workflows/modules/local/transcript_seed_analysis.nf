process TRANSCRIPT_SEED_ANALYSIS {
    tag "transcript_seed_${species}"
    label 'process_low'
    publishDir "${params.outdir}/transcript_seed", mode: params.publish_dir_mode

    input:
    tuple val(species), path(transcriptome_fasta), path(candidates_fasta)
    val seed_scope
    val site_identity

    output:
    path "transcript_seed_${species}_summary.json", emit: summary
    path "transcript_seed_${species}_evidence.json", emit: evidence
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - <<'PYEOF'
import sys
from pathlib import Path
sys.path.insert(0, '${workflow.projectDir}/../../../../')
from sirnaforge.core.transcript_seed import run_transcript_seed_analysis

result = run_transcript_seed_analysis(
    '${candidates_fasta}',
    '${transcriptome_fasta}',
    species='${species}',
    output_dir='.',
    scope='${seed_scope}',
    site_identity='${site_identity}',
)
print(f"Transcript seed scope: {result['scope']}; sites: {result['site_count']}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    biopython: \$(python -c "import Bio; print(Bio.__version__)")
END_VERSIONS
    """

    stub:
    """
    echo '{"schema_version":"1","channel":"transcript_seed","status":"failed","detail":"stub_run"}' > transcript_seed_${species}_summary.json
    echo '{"schema_version":"2","producer":"stub","source":"synthesized","entry":{"channel":"transcript_seed","species":"${species}","guide_set_digest":"0000000000000000","status":"failed","detail":"stub_run"}}' > transcript_seed_${species}_evidence.json
    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
END_VERSIONS
    """
}
