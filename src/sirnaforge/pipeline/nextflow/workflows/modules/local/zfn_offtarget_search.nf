/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ZFN OFF-TARGET SEARCH MODULE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Exhaustive sliding-window ZFN off-target search.
    Wraps the Python-native ZFNOffTargetSearcher for Nextflow orchestration,
    containerization, resource management and -resume semantics.

    The algorithm is NOT BWA-MEM2 alignment — it is exhaustive paired half-site
    scanning with FokI seed-region penalties, spacer constraints, and IUPAC handling.
----------------------------------------------------------------------------------------
*/

process ZFN_OFFTARGET_SEARCH {
    tag "zfn_offtarget"
    label 'process_high'
    publishDir "${params.outdir}/zfn_offtarget", mode: params.publish_dir_mode

    input:
    val left_half_site
    val right_half_site
    path genome_fasta
    val algorithm
    val dimer_mode
    val spacer_lengths
    val max_mismatches
    path annotation_file

    output:
    path "zfn_offtarget_sites.tsv", emit: sites
    path "zfn_candidate_summary.json", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def annotation_arg = annotation_file.name != 'NO_ANNOTATION' ? "--annotation ${annotation_file}" : ''
    """
    python3 <<'PYEOF'
import sys, json
sys.path.insert(0, '${workflow.projectDir}/../src')

from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher
from sirnaforge.models.zfn import (
    DimerMode,
    GenomicAnnotationConfig,
    ZFNAlgorithm,
    ZFNDesignParameters,
    ZFNHalfSiteConstraints,
    ZFNSpacerConstraints,
)
from sirnaforge.zfn.design import ZFNDesigner

spacer_list = [int(s) for s in '${spacer_lengths}'.split(',')]

params_obj = ZFNDesignParameters(
    left_half_site='${left_half_site}',
    right_half_site='${right_half_site}',
    search_space_fasta='${genome_fasta}',
    search_space_reference=None,
    algorithm=ZFNAlgorithm('${algorithm}'),
    dimer_mode=DimerMode('${dimer_mode}'),
    spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=spacer_list),
    half_site_constraints=ZFNHalfSiteConstraints(max_mismatches=${max_mismatches}),
)

annotation = None
annotation_path = '${annotation_file}'
if annotation_path and annotation_path != 'NO_ANNOTATION':
    annotation = GenomicAnnotationConfig(annotation_path=annotation_path)

designer = ZFNDesigner()
result = designer.evaluate_pair(params=params_obj, annotation=annotation)

# Write off-target sites TSV
result.save_offtargets_csv('zfn_offtarget_sites.tsv')

# Write candidate summary JSON
summary = result.get_summary()
candidates = [c.model_dump(mode='json') for c in result.candidates]
with open('zfn_candidate_summary.json', 'w') as f:
    json.dump({'candidates': candidates, 'summary': summary}, f, indent=2, default=str)

print(f"ZFN off-target search completed: {len(result.off_target_sites)} sites found")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sirnaforge: \$(python -c 'import sirnaforge; print(sirnaforge.__version__)' 2>/dev/null || echo 'dev')
    END_VERSIONS
    """

    stub:
    """
    echo -e "site_id\\tchrom\\tstart_1based\\tend_1based\\tstrand\\torientation\\tspacer_len\\tsequence\\tleft_mismatches\\tright_mismatches\\ttotal_mismatches\\tscore\\tregion\\tnearest_gene\\tleft_aligned\\tright_aligned" > zfn_offtarget_sites.tsv
    echo '{"candidates": [], "summary": {"candidates": 0, "off_target_sites": 0}}' > zfn_candidate_summary.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sirnaforge: \$(python -c 'import sirnaforge; print(sirnaforge.__version__)' 2>/dev/null || echo 'dev')
    END_VERSIONS
    """
}
