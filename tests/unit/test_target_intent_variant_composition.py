"""Variant-targeted candidates compose with the frozen target-intent sequence universe."""

from sirnaforge.core.target_intent import (
    TargetIntentRequest,
    TargetIntentVerdict,
    evaluate_target_intent_guide,
    resolve_target_intent,
)
from sirnaforge.data.base import DatabaseType, TranscriptInfo
from sirnaforge.models.policy import TargetSelectivity
from sirnaforge.models.sirna import SiRNACandidate
from sirnaforge.models.transcript_annotation import Interval, TranscriptAnnotation
from sirnaforge.models.variant import VariantRecord
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig


def test_alternate_variant_guide_passes_intent_on_its_resolved_grch38_allele(tmp_path) -> None:
    """An alt-only guide is checked against the mapped cDNA edit, not reference-only sequence."""
    transcript_id = "ENST00000000001.1"
    sequence = "C" * 25 + "A" + "C" * 24
    variant = VariantRecord(chr="chr1", pos=125, ref="A", alt="G")
    alternate = sequence[:25] + "G" + sequence[26:]
    target = alternate[15:36]
    guide = target.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    transcript = TranscriptInfo(
        transcript_id=transcript_id,
        gene_id="ENSG00000000001",
        sequence=sequence,
        database=DatabaseType.ENSEMBL,
    )
    resolved = resolve_target_intent(
        TargetIntentRequest(
            selectivity=TargetSelectivity.PAN_ISOFORM,
            target_species=frozenset({"human"}),
        ),
        [transcript],
        annotation_transcript_ids=[transcript_id],
        annotation_provenance="GRCh38 synthetic variant fixture",
        enumeration_input_ids=[transcript_id],
        sequence_species="human",
    )
    candidate = SiRNACandidate(
        id="alt-guide",
        transcript_id=transcript_id,
        position=16,
        guide_sequence=guide,
        passenger_sequence=target,
        gc_content=100 * sum(base in "GC" for base in target) / len(target),
        length=21,
        asymmetry_score=0.5,
        allele_specific=True,
        targeted_alleles=["alt"],
        overlapped_variants=[variant.model_dump(mode="json")],
    )
    workflow = SiRNAWorkflow(WorkflowConfig(output_dir=tmp_path / "run", gene_query="fixture"))
    workflow._transcript_annotations[transcript_id] = TranscriptAnnotation(
        transcript_id=transcript_id,
        gene_id="ENSG00000000001",
        seq_region_name="1",
        start=100,
        end=149,
        strand=1,
        exons=[Interval(seq_region_name="1", start=100, end=149, strand=1)],
        provider="fixture",
        reference_choice="GRCh38",
    )

    overrides, unresolved = workflow._target_intent_allele_sequences(candidate, resolved)
    evaluation = evaluate_target_intent_guide(
        guide,
        resolved,
        sequence_overrides=overrides,
        unresolved_allele_transcript_ids=unresolved,
    )

    assert evaluate_target_intent_guide(guide, resolved).verdict is TargetIntentVerdict.FAIL
    assert unresolved == set()
    assert overrides == {"ENST00000000001": alternate}
    assert evaluation.verdict is TargetIntentVerdict.PASS
    assert evaluation.covered_transcript_ids == frozenset({"ENST00000000001"})

    workflow._transcript_annotations.clear()
    overrides, unresolved = workflow._target_intent_allele_sequences(candidate, resolved)
    unknown = evaluate_target_intent_guide(
        guide,
        resolved,
        sequence_overrides=overrides,
        unresolved_allele_transcript_ids=unresolved,
    )
    assert overrides == {}
    assert unresolved == {"ENST00000000001"}
    assert unknown.verdict is TargetIntentVerdict.UNKNOWN
    assert unknown.unresolved_allele_transcript_ids == frozenset({"ENST00000000001"})
