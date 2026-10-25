"""Discover canonical siRNA seed-complementary sites in transcript sequences.

This is a transcript-site channel, separate from known-miRNA seed resemblance scanning.
"""

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio import SeqIO

from sirnaforge.core.screening_evidence import EvidenceProducer, write_evidence
from sirnaforge.models.evidence import EvidenceStatus, ObservedCount, ObservedCounts, ScreeningEvidenceEntry
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.models.transcript_seed import (
    TranscriptSeedClass,
    TranscriptSeedScanResult,
    TranscriptSeedSite,
    TranscriptSeedTranscript,
    TranscriptSequenceScope,
)
from sirnaforge.utils.hashing import file_sha256


@dataclass(frozen=True)
class _SitePattern:
    site_class: TranscriptSeedClass
    sequence: str


_SITE_PRIORITY = {
    TranscriptSeedClass.EIGHTMER: 0,
    TranscriptSeedClass.SEVENMER_M8: 1,
    TranscriptSeedClass.SEVENMER_A1: 2,
    TranscriptSeedClass.SIXMER: 3,
}


def _normalize(sequence: str, *, label: str) -> str:
    """Normalize RNA and DNA to uppercase DNA alphabet and reject ambiguity symbols."""
    normalized = sequence.upper().replace("U", "T")
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        raise ValueError(f"{label} contains unsupported nucleotide symbol(s): {', '.join(invalid)}")
    return normalized


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _patterns(guide: str) -> tuple[_SitePattern, ...]:
    """Build target 5′→3′ patterns for guide positions 2–8 and the A1 geometry."""
    if len(guide) < 8:
        raise ValueError("guide must contain at least 8 nucleotides to define seed positions 2–8")
    seed_2_8 = guide[1:8]
    seed_2_7 = guide[1:7]
    return (
        _SitePattern(TranscriptSeedClass.EIGHTMER, _reverse_complement(seed_2_8) + "A"),
        _SitePattern(TranscriptSeedClass.SEVENMER_M8, _reverse_complement(seed_2_8)),
        _SitePattern(TranscriptSeedClass.SEVENMER_A1, _reverse_complement(seed_2_7) + "A"),
        _SitePattern(TranscriptSeedClass.SIXMER, _reverse_complement(seed_2_7)),
    )


def _find_pattern_hits(target: str, patterns: tuple[_SitePattern, ...]) -> dict[int, _SitePattern]:
    """Return the most specific canonical class at each target start coordinate."""
    hits: dict[int, _SitePattern] = {}
    for pattern in patterns:
        width = len(pattern.sequence)
        for local_start in range(0, len(target) - width + 1):
            if target[local_start : local_start + width] != pattern.sequence:
                continue
            previous = hits.get(local_start)
            if previous is None or _SITE_PRIORITY[pattern.site_class] < _SITE_PRIORITY[previous.site_class]:
                hits[local_start] = pattern
    return hits


def scan_transcript_seed_sites(
    guides: Mapping[str, str],
    transcripts: Sequence[TranscriptSeedTranscript],
    *,
    scope: TranscriptSequenceScope = TranscriptSequenceScope.FULL_CDNA,
) -> TranscriptSeedScanResult:
    """Find canonical 6mer/7mer-A1/7mer-m8/8mer sites for siRNA guide sequences.

    Coordinates are 1-based and always refer to the transcript-oriented cDNA, including for
    3′UTR scans. A missing UTR annotation is reported in ``transcripts_unavailable`` and is
    skipped; it is never interpreted as an annotated empty UTR.
    """
    scope = TranscriptSequenceScope(scope)
    if len({transcript.transcript_id for transcript in transcripts}) != len(transcripts):
        raise ValueError("transcript_id values must be unique within one scan")
    prepared_guides = {
        guide_id: _patterns(_normalize(seq, label=f"guide {guide_id}")) for guide_id, seq in guides.items()
    }
    if not prepared_guides:
        raise ValueError("at least one guide sequence is required")

    sites: dict[tuple[str, str, TranscriptSequenceScope, int, TranscriptSeedClass], TranscriptSeedSite] = {}
    unavailable: list[str] = []
    scanned = 0
    for transcript in transcripts:
        cdna = _normalize(transcript.cdna_sequence, label=f"transcript {transcript.transcript_id} cDNA")
        if scope is TranscriptSequenceScope.THREE_PRIME_UTR:
            if transcript.utr3_sequence is None or transcript.utr3_start is None:
                unavailable.append(transcript.transcript_id)
                continue
            target = _normalize(transcript.utr3_sequence, label=f"transcript {transcript.transcript_id} 3′UTR")
            start_offset = transcript.utr3_start - 1
            if start_offset + len(target) > len(cdna) or cdna[start_offset : start_offset + len(target)] != target:
                raise ValueError(f"transcript {transcript.transcript_id} 3′UTR does not match its cDNA interval")
        else:
            target = cdna
            start_offset = 0
        scanned += 1

        for guide_id, patterns in prepared_guides.items():
            hits_by_start = _find_pattern_hits(target, patterns)
            for local_start, pattern in hits_by_start.items():
                transcript_start = start_offset + local_start + 1
                key = (guide_id, transcript.transcript_id, scope, transcript_start, pattern.site_class)
                sites[key] = TranscriptSeedSite(
                    guide_id=guide_id,
                    transcript_id=transcript.transcript_id,
                    gene_id=transcript.gene_id,
                    scope=scope,
                    site_class=pattern.site_class,
                    start=transcript_start,
                    end=transcript_start + len(pattern.sequence) - 1,
                    target_sequence=pattern.sequence,
                )

    ordered = tuple(
        sorted(
            sites.values(),
            key=lambda site: (site.transcript_id, site.start, site.guide_id, _SITE_PRIORITY[site.site_class]),
        )
    )
    transcripts_with_sites = {site.transcript_id for site in ordered}
    genes_with_sites = {site.gene_id for site in ordered if site.gene_id is not None}
    unresolved_transcripts = {site.transcript_id for site in ordered if site.gene_id is None}
    return TranscriptSeedScanResult(
        scope=scope,
        sites=ordered,
        site_count=len(ordered),
        transcript_count=len(transcripts_with_sites),
        gene_count=len(genes_with_sites),
        unresolved_gene_count=len(unresolved_transcripts),
        transcripts_scanned=scanned,
        transcripts_unavailable=tuple(sorted(unavailable)),
    )


def run_transcript_seed_analysis(
    candidates_file: str | Path,
    transcriptome_fasta: str | Path,
    *,
    species: str,
    output_dir: str | Path,
    scope: TranscriptSequenceScope = TranscriptSequenceScope.FULL_CDNA,
    guide_set_digest: str | None = None,
    reference_id: str | None = None,
    utr3_manifest: str | Path | None = None,
) -> dict[str, object]:
    """Run site discovery over one materialized full transcriptome and write evidence artifacts.

    The FASTA is passed by the workflow after reference resolution; this function never downloads or
    substitutes a target-gene-only sequence set. Plain cDNA FASTA has no explicit 3′UTR feature, so
    a 3′UTR request is published FAILED with unavailable transcript annotations rather than as an
    empty successful scan.
    """
    candidates_path = Path(candidates_file)
    reference_path = Path(transcriptome_fasta)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    guide_records = list(SeqIO.parse(candidates_path, "fasta"))
    transcript_records = list(SeqIO.parse(reference_path, "fasta"))
    guides = {record.id: str(record.seq) for record in guide_records}
    utr_annotations: dict[str, Any] = {}
    utr_checksum = None
    if utr3_manifest is not None:
        utr_path = Path(utr3_manifest)
        manifest_payload = json.loads(utr_path.read_text())
        raw_annotations = manifest_payload.get("transcripts") if isinstance(manifest_payload, dict) else None
        if not isinstance(raw_annotations, dict):
            raise ValueError("3′UTR manifest must contain a 'transcripts' object keyed by transcript ID")
        utr_annotations = raw_annotations
        utr_checksum = file_sha256(utr_path)
    transcripts: list[TranscriptSeedTranscript] = []
    for record in transcript_records:
        annotation = utr_annotations.get(record.id)
        if annotation is not None and not isinstance(annotation, dict):
            raise ValueError(f"3′UTR manifest entry for {record.id} must be an object")
        values = annotation if isinstance(annotation, dict) else {}
        transcripts.append(
            TranscriptSeedTranscript(
                transcript_id=record.id,
                gene_id=values.get("gene_id"),
                cdna_sequence=str(record.seq),
                utr3_sequence=values.get("utr3_sequence"),
                utr3_start=values.get("utr3_start"),
            )
        )
    digest = guide_set_digest or file_sha256(candidates_path)[:16]
    checksum = file_sha256(reference_path)
    scan = scan_transcript_seed_sites(guides, transcripts, scope=scope)
    unavailable = len(scan.transcripts_unavailable)
    detail = None
    status = EvidenceStatus.COMPLETE
    if not guides:
        status = EvidenceStatus.FAILED
        detail = "candidate FASTA contained no guides"
    elif not transcript_records:
        status = EvidenceStatus.FAILED
        detail = "transcriptome FASTA contained no transcripts"
    elif unavailable:
        status = EvidenceStatus.FAILED
        detail = f"{unavailable} transcript(s) lack explicit 3′UTR annotations; no zero-site result was inferred"
    scan_complete = status is EvidenceStatus.COMPLETE

    per_guide: dict[str, dict[str, object]] = {}
    for guide_id in guides:
        guide_sites = [site for site in scan.sites if site.guide_id == guide_id]
        unresolved_gene_site_count = sum(site.gene_id is None for site in guide_sites)
        per_guide[guide_id] = {
            "status": "complete" if scan_complete else "unknown",
            "site_count": len(guide_sites) if scan_complete else None,
            "site_counts_by_class": {
                site_class.value: sum(site.site_class is site_class for site in guide_sites)
                for site_class in TranscriptSeedClass
            },
            "transcript_count": len({site.transcript_id for site in guide_sites}) if scan_complete else None,
            "gene_count": len({site.gene_id for site in guide_sites if site.gene_id is not None})
            if scan_complete
            else None,
            "gene_count_is_lower_bound": unresolved_gene_site_count > 0 if scan_complete else None,
            "unresolved_gene_site_count": unresolved_gene_site_count if scan_complete else None,
        }

    payload: dict[str, object] = {
        "schema_version": "1",
        "channel": ScreeningChannel.TRANSCRIPT_SEED.value,
        "status": status.value,
        "detail": detail,
        "species": species,
        "scope": TranscriptSequenceScope(scope).value,
        "reference_path": str(reference_path),
        "reference_sha256": checksum,
        "utr_annotation_sha256": utr_checksum,
        "guide_set_digest": digest,
        "submitted_guides": len(guides),
        "processed_guides": len(guides),
        "transcripts_scanned": scan.transcripts_scanned,
        "transcripts_unavailable": list(scan.transcripts_unavailable),
        "site_count": scan.site_count if scan_complete else None,
        "transcript_count": scan.transcript_count if scan_complete else None,
        "gene_count": scan.gene_count if scan_complete else None,
        "unresolved_gene_site_count": sum(site.gene_id is None for site in scan.sites) if scan_complete else None,
        "gene_count_is_lower_bound": scan.unresolved_gene_count > 0 if scan_complete else None,
        "sites": [site.model_dump(mode="json") for site in scan.sites],
        "per_guide": per_guide,
    }
    summary_path = out_dir / f"transcript_seed_{species}_summary.json"
    summary_path.write_text(json.dumps(payload, indent=2))
    entry = ScreeningEvidenceEntry(
        channel=ScreeningChannel.TRANSCRIPT_SEED,
        species=species,
        reference_id=reference_id or checksum,
        guide_set_digest=digest,
        status=status,
        counts=ObservedCounts(
            sites=ObservedCount(value=scan.site_count if scan_complete else None),
            distinct_transcripts=ObservedCount(value=scan.transcript_count if scan_complete else None),
            distinct_genes=ObservedCount(
                value=scan.gene_count if scan_complete else None,
                is_lower_bound=scan.unresolved_gene_count > 0 if scan_complete else False,
            ),
            unresolved_gene_sites=ObservedCount(
                value=sum(site.gene_id is None for site in scan.sites) if scan_complete else None
            ),
        ),
        submitted_guide_digest=digest,
        submitted_guides=len(guides),
        processed_guides=len(guides),
        detail=detail,
    )
    evidence_path = write_evidence(
        out_dir,
        producer=EvidenceProducer.TRANSCRIPT_SEED_ANALYSIS,
        entry=entry,
    )
    payload["summary_path"] = str(summary_path)
    payload["evidence_path"] = str(evidence_path)
    return payload
