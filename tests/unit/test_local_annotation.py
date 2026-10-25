"""Local annotation providers preserve source structure, identifiers, and provenance offline."""

from __future__ import annotations

import asyncio
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import AsyncMock

import pytest
from typer.testing import CliRunner

from sirnaforge.cli import app
from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.data.base import (
    DatabaseAccessError,
    DatabaseType,
    GencodeClient,
    GeneInfo,
    SequenceType,
    TranscriptInfo,
)
from sirnaforge.data.gene_search import GeneSearcher
from sirnaforge.data.local_annotation import LocalTranscriptModelClient
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig


@pytest.fixture
def gencode_files(tmp_path: Path) -> tuple[Path, Path]:
    """A negative-strand, two-exon transcript with repeated canonical tags."""
    attributes = (
        'gene_id "ENSG1.3"; transcript_id "ENST1.2"; gene_name "TP53"; '
        'transcript_name "TP53-201"; transcript_type "protein_coding"; '
        'tag "basic"; tag "Ensembl_canonical";'
    )
    lines = [
        "#!genome-build GRCh38",
        "#!annotation-version 49",
        'chr17\tHAVANA\tgene\t100\t400\t.\t-\t.\tgene_id "ENSG1.3"; gene_name "TP53"; gene_type "protein_coding";',
    ]
    for feature, start, end in [
        ("transcript", 100, 400),
        ("exon", 100, 149),
        ("exon", 300, 349),
        ("CDS", 110, 140),
        ("CDS", 310, 340),
    ]:
        lines.append(f"chr17\tHAVANA\t{feature}\t{start}\t{end}\t.\t-\t.\t{attributes}")
    annotation = tmp_path / "release.gtf"
    annotation.write_text("\n".join(lines) + "\n")
    fasta = tmp_path / "release.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">ENST1.2|ENSG1.3|TP53-201|100|\n" + "A" * 100 + "\n")
    return annotation, fasta


@pytest.fixture
def refseq_file(tmp_path: Path) -> Path:
    """A RefSeq GFF3 with accession contigs, encoded attributes, and RNA parents."""
    path = tmp_path / "refseq.gff3"
    path.write_text(
        "##gff-version 3\n#!genome-build GRCh38.p14\n"
        "NC_000017.11\tBestRefSeq\tgene\t100\t400\t.\t+\t.\tID=gene-TP53;Dbxref=GeneID:7157,HGNC:HGNC:11998;Name=TP53;gene_biotype=protein_coding\n"
        "NC_000017.11\tBestRefSeq\tmRNA\t100\t400\t.\t+\t.\tID=rna-NM_000546.6;Parent=gene-TP53;transcript_id=NM_000546.6;product=protein%20p53\n"
        "NC_000017.11\tBestRefSeq\texon\t100\t149\t.\t+\t.\tParent=rna-NM_000546.6\n"
        "NC_000017.11\tBestRefSeq\tCDS\t110\t140\t.\t+\t0\tParent=rna-NM_000546.6\n"
    )
    return path


@pytest.mark.asyncio
async def test_gencode_structure_versions_and_gzip(gencode_files: tuple[Path, Path], tmp_path: Path) -> None:
    """Resolve versionless IDs, preserve exact versions, and retain discontinuous CDS."""
    annotation, _ = gencode_files
    compressed = tmp_path / "release.gtf.gz"
    with gzip.open(compressed, "wb") as handle:
        handle.write(annotation.read_bytes())
    client = LocalTranscriptModelClient(compressed, provider="gencode", assembly="GRCh38", species="human")
    bundle = await client.fetch_by_ids(
        ["ENST1", "ENST1.1", "missing"], species="human", reference=ReferenceChoice.explicit("GRCh38")
    )
    assert bundle.unresolved == ["ENST1.1", "missing"]
    model = bundle.transcripts["ENST1.2"]
    assert model.gene_id == "ENSG1.3"
    assert model.symbol == "TP53"
    assert model.is_canonical
    assert model.transcript_name == "TP53-201"
    assert [interval.start for interval in model.exons] == [300, 100]
    assert [(interval.start, interval.end) for interval in model.cds] == [(310, 340), (110, 140)]
    assert client.annotation_version == "49"
    assert client.source_metadata()["annotation_version"] == "49"
    assert client.source_metadata()["sha256"] == hashlib.sha256(compressed.read_bytes()).hexdigest()


@pytest.mark.asyncio
async def test_refseq_parent_mapping_and_regions(refseq_file: Path) -> None:
    """NCBI GeneID and RNA Parent relationships resolve without chromosome-name guesses."""
    client = LocalTranscriptModelClient(refseq_file, provider="refseq", assembly="GRCh38.p14", species="human")
    reference = ReferenceChoice.explicit("GRCh38.p14")
    for query in ["NM_000546", "NM_000546.6", "7157", "TP53"]:
        bundle = await client.fetch_by_ids([query], species="Homo sapiens", reference=reference)
        model = bundle.transcripts["NM_000546.6"]
        assert not bundle.unresolved
        assert model.gene_id == "7157"
        assert model.symbol == "TP53"
        assert model.biotype == "protein_coding"
        assert model.seq_region_name == "NC_000017.11"
        assert model.chr == "NC_000017.11"
        assert model.exons[0].start == 100
        assert model.provider == "refseq_file"
    bundle = await client.fetch_by_regions(["NC_000017.11:400-401", "17:100-400"], species="human", reference=reference)
    assert bundle.resolved_count == 1
    assert bundle.unresolved == ["17:100-400"]


@pytest.mark.asyncio
@pytest.mark.parametrize("species,assembly", [("mouse", "GRCh38"), ("human", "GRCh37")])
async def test_local_source_rejects_reference_mismatch(
    gencode_files: tuple[Path, Path], species: str, assembly: str
) -> None:
    """A local source cannot silently change species or coordinate assemblies."""
    client = LocalTranscriptModelClient(gencode_files[0], provider="gencode", assembly="GRCh38", species="human")
    with pytest.raises(ValueError):
        await client.fetch_by_ids(["ENST1"], species=species, reference=ReferenceChoice.explicit(assembly))


def test_refseq_header_rejects_wrong_declared_build(refseq_file: Path) -> None:
    """An assembly declared by the annotation itself wins over a caller's mistaken label."""
    with pytest.raises(ValueError, match="File declares assembly"):
        LocalTranscriptModelClient(refseq_file, provider="refseq", assembly="GRCh37", species="human")


@pytest.mark.asyncio
async def test_gencode_gene_and_sequence_lookup(gencode_files: tuple[Path, Path]) -> None:
    """Gene, symbol, and transcript queries use one pinned annotation/FASTA pair."""
    annotation, fasta = gencode_files
    client = GencodeClient(gtf_path=annotation, fasta_path=fasta, assembly="GRCh38")
    for query in ["TP53", "ENSG1", "ENST1"]:
        gene, transcripts = await client.search_gene(query)
        assert gene is not None and gene.gene_id == "ENSG1.3"
        assert len(transcripts) == 1
        assert transcripts[0].sequence == "A" * 100
        assert transcripts[0].is_canonical
    assert await client.get_sequence("ENST1.2") == "A" * 100
    with pytest.raises(ValueError, match="cDNA only"):
        await client.get_sequence("ENST1.2", SequenceType.PROTEIN)


@pytest.mark.asyncio
async def test_gencode_canonical_provenance_is_isoform_specific_and_versioned(
    gencode_files: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A canonical tag belongs to its transcript and its annotation release is retained."""
    annotation, fasta = gencode_files
    source_metadata = LocalTranscriptModelClient(
        annotation, provider="gencode", assembly="GRCh38", species="human"
    ).source_metadata()
    assert source_metadata["annotation_version"] == "49"
    attributes = (
        'gene_id "ENSG1.3"; transcript_id "ENST2.1"; gene_name "TP53"; '
        'transcript_name "TP53-202"; transcript_type "protein_coding"; tag "basic";'
    )
    with annotation.open("a") as handle:
        handle.write(
            "chr17\tHAVANA\ttranscript\t100\t400\t.\t-\t.\t"
            + attributes
            + "\n"
            + "chr17\tHAVANA\texon\t100\t149\t.\t-\t.\t"
            + attributes
            + "\n"
        )
    with gzip.open(fasta, "at") as handle:
        handle.write(">ENST2.1|ENSG1.3|TP53-202|100|\n" + "C" * 50 + "\n")

    async def inline_to_thread(function, /, *args, **kwargs):
        return function(*args, **kwargs)

    monkeypatch.setattr(asyncio, "to_thread", inline_to_thread)
    client = GencodeClient(gtf_path=annotation, fasta_path=fasta, assembly="GRCh38")
    _, transcripts = await client.search_gene("TP53")
    by_id = {transcript.transcript_id: transcript for transcript in transcripts}

    assert by_id["ENST1.2"].is_canonical
    assert by_id["ENST1.2"].canonical_release == "49"
    assert by_id["ENST1.2"].canonical_source == "GENCODE explicit canonical/MANE tags (annotation 49)"
    assert not by_id["ENST2.1"].is_canonical
    assert by_id["ENST2.1"].canonical_source is None
    assert by_id["ENST2.1"].canonical_release is None


@pytest.mark.asyncio
async def test_gencode_rejects_sequence_version_mismatch(gencode_files: tuple[Path, Path], tmp_path: Path) -> None:
    """A matching versionless accession must not conceal mixed reference releases."""
    annotation, _ = gencode_files
    fasta = tmp_path / "wrong.fa"
    fasta.write_text(">ENST1.1\nAAAA\n")
    client = GencodeClient(gtf_path=annotation, fasta_path=fasta, assembly="GRCh38")
    with pytest.raises(DatabaseAccessError, match="exact annotation accessions"):
        await client.search_gene("TP53")


@pytest.mark.asyncio
async def test_local_workflow_annotation_records_provider_and_checksum(
    gencode_files: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Configured local annotation reaches the workflow's public provenance surface."""
    annotation, fasta = gencode_files
    monkeypatch.setenv("SIRNAFORGE_ANNOTATION_PROVIDER", "gencode")
    monkeypatch.setenv("SIRNAFORGE_GENCODE_GTF", str(annotation))
    monkeypatch.setenv("SIRNAFORGE_GENCODE_FASTA", str(fasta))
    monkeypatch.setenv("SIRNAFORGE_GENCODE_ASSEMBLY", "GRCh38")
    workflow = SiRNAWorkflow(
        WorkflowConfig(output_dir=tmp_path / "run", gene_query="TP53", database=DatabaseType.GENCODE)
    )
    await workflow._enrich_transcript_annotations(
        [TranscriptInfo(transcript_id="ENST1.2", gene_id="ENSG1.3", database=DatabaseType.GENCODE)]
    )
    assert workflow._annotation_summary["provider"] == "gencode_file"
    assert workflow._annotation_summary["transcripts_resolved"] == 1
    assert workflow._annotation_provenance()["source"]["sha256"] == hashlib.sha256(annotation.read_bytes()).hexdigest()


@pytest.mark.asyncio
async def test_selected_provider_is_first_in_fallback_order() -> None:
    """Choosing a local provider must not incur an Ensembl request before using it."""
    searcher = GeneSearcher()
    searcher.clients[DatabaseType.GENCODE].search_gene = AsyncMock(
        return_value=(
            GeneInfo(gene_id="ENSG1", database=DatabaseType.GENCODE),
            [TranscriptInfo(transcript_id="ENST1", gene_id="ENSG1", database=DatabaseType.GENCODE)],
        )
    )
    searcher.clients[DatabaseType.ENSEMBL].search_gene = AsyncMock(
        side_effect=DatabaseAccessError("offline", "Ensembl")
    )
    searcher.clients[DatabaseType.REFSEQ].search_gene = AsyncMock(side_effect=DatabaseAccessError("offline", "RefSeq"))
    await searcher.search_gene_with_fallback("TP53", preferred_database=DatabaseType.GENCODE)
    searcher.clients[DatabaseType.GENCODE].search_gene.assert_awaited_once()
    searcher.clients[DatabaseType.ENSEMBL].search_gene.assert_not_awaited()


def test_annotation_cli_writes_structured_json(gencode_files: tuple[Path, Path], tmp_path: Path) -> None:
    """The offline provider can be selected explicitly without environment configuration."""
    output = tmp_path / "annotations.json"
    result = CliRunner().invoke(
        app,
        [
            "annotate",
            "ENST1",
            "--provider",
            "gencode",
            "--annotation-file",
            str(gencode_files[0]),
            "--assembly",
            "GRCh38",
            "--output",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    data = json.loads(output.read_text())
    assert data["source"]["provider"] == "gencode_file"
    assert data["transcripts"]["ENST1.2"]["symbol"] == "TP53"


@pytest.mark.asyncio
async def test_refseq_workflow_enriches_non_ensembl_accessions(
    refseq_file: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """RefSeq RNA accessions reach structural enrichment instead of being filtered out as non-ENST."""
    monkeypatch.setenv("SIRNAFORGE_ANNOTATION_PROVIDER", "refseq")
    monkeypatch.setenv("SIRNAFORGE_REFSEQ_ANNOTATION", str(refseq_file))
    monkeypatch.setenv("SIRNAFORGE_REFSEQ_ASSEMBLY", "GRCh38.p14")
    workflow = SiRNAWorkflow(
        WorkflowConfig(output_dir=tmp_path / "run", gene_query="TP53", database=DatabaseType.REFSEQ)
    )
    await workflow._enrich_transcript_annotations(
        [TranscriptInfo(transcript_id="NM_000546.6", gene_id="7157", database=DatabaseType.REFSEQ)]
    )
    assert workflow._annotation_summary["provider"] == "refseq_file"
    assert workflow._annotation_summary["transcripts_resolved"] == 1
    assert workflow._annotation_provenance()["rest_release"]["state"] == "not_applicable"


@pytest.mark.asyncio
async def test_conflicting_transcript_loci_are_rejected(gencode_files: tuple[Path, Path]) -> None:
    """A duplicate accession on a different locus must not silently overwrite the first record."""
    annotation, _ = gencode_files
    transcript = next(line for line in annotation.read_text().splitlines() if "\ttranscript\t" in line)
    with annotation.open("a") as handle:
        handle.write(transcript.replace("chr17", "chr18") + "\n")
    client = LocalTranscriptModelClient(annotation, provider="gencode", species="human", assembly="GRCh38")
    with pytest.raises(ValueError, match="Conflicting transcript"):
        await client.fetch_by_ids(["ENST1"], species="human", reference=ReferenceChoice.explicit("GRCh38"))


def test_gencode_search_cli_writes_sequences_without_ensembl(
    gencode_files: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A selected, configured GENCODE release works through the public search command."""
    annotation, fasta = gencode_files
    monkeypatch.setenv("SIRNAFORGE_GENCODE_GTF", str(annotation))
    monkeypatch.setenv("SIRNAFORGE_GENCODE_FASTA", str(fasta))
    monkeypatch.setenv("SIRNAFORGE_GENCODE_ASSEMBLY", "GRCh38")
    output = tmp_path / "transcripts.fa"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "sirnaforge.cli",
            "search",
            "TP53",
            "--database",
            "gencode",
            "--no-fallback",
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "ENST1.2" in output.read_text()
    assert "".join(output.read_text().splitlines()[1:]) == "A" * 100
