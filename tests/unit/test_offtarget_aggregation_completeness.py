"""Absent evidence must not publish as a clean screen, and no channel may be counted twice.

Every defect pinned here was reproduced on a real public TP53 screen before it was fixed (#100):
a species whose alignment file the aggregator rejected still reported as screened; a species with
no index vanished before Nextflow or came back as a completed screen with zero hits; and every
miRNA counter was exactly 2x because one file was ingested by a glob and again by a named retry.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

import asyncio
import csv
import json
from pathlib import Path

import pytest

from sirnaforge.core.off_target import aggregate_offtarget_results
from sirnaforge.models.off_target import MiRNAHit, OffTargetHit
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

GUIDE = "ATGCGATGCGATGCGATGCGC"

#: Offline orthologue evidence, so no test here touches Ensembl Compara. See tests/unit/data/README.md.
_ORTHOLOG_MAPPING_FIXTURE = Path(__file__).resolve().parent / "data" / "ortholog_mapping_synthetic.json"

GENOME_COLUMNS = OffTargetHit.tsv_header().split("\t")
MIRNA_COLUMNS = MiRNAHit.tsv_header().split("\t")


def _candidate(candidate_id: str = "cand_1", guide: str = GUIDE) -> SiRNACandidate:
    """A minimal candidate carrying one guide."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=guide,
        passenger_sequence=guide[::-1],
        gc_content=57.1,
        length=len(guide),
        asymmetry_score=0.7,
        composite_score=50.0,
    )


def _workflow(tmp_path: Path, name: str, species: list[str] | None = None) -> SiRNAWorkflow:
    """A TP53 workflow that resolves orthologues from a file rather than over the network."""
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        genome_species=species or ["human"],
        design_params=DesignParameters(),
        ortholog_mapping_file=_ORTHOLOG_MAPPING_FIXTURE,
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_ids = {"ENSG00000000001"}
    workflow._query_gene_symbols = {"TP53"}
    return workflow


def _write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    """Write one alignment table in the producer's column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _genome_row(species: str, rname: str, qname: str = "cand_1", nm: int = 2) -> dict[str, str]:
    """One transcriptome alignment row."""
    return {
        "qname": qname,
        "qseq": GUIDE,
        "species": species,
        "rname": rname,
        "coord": "500",
        "strand": "-",
        "cigar": f"{len(GUIDE)}M",
        "mapq": "60",
        "as_score": "30",
        "nm": str(nm),
        "seed_mismatches": str(nm),
        "offtarget_score": "0.0" if nm == 0 else "1.5",
    }


def _mirna_row(qname: str = "cand_1") -> dict[str, str]:
    """One miRNA seed-match row, scored high enough not to be a high-risk hit."""
    return {
        "qname": qname,
        "qseq": GUIDE,
        "species": "human",
        "database": "mirgenedb",
        "mirna_id": "hsa-miR-21-5p",
        "coord": "3",
        "strand": "+",
        "cigar": "7M",
        "mapq": "30",
        "as_score": "14",
        "nm": "1",
        "seed_mismatches": "1",
        "offtarget_score": "9.0",
    }


# ---------------------------------------------------------------------------
# miRNA double-ingest
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_mirna_analysis_file_is_ingested_exactly_once(tmp_path):
    """``mirna/mirna_analysis.tsv`` matches the glob AND the named retry, so it was read twice.

    Inert under the default thresholds, but ``max_total_offtarget_hits`` gates on
    ``human_transcriptome_hits + mirna_human_total``, so the doubled counter is a real gate input.
    """
    workflow = _workflow(tmp_path, "mirna_once")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    _write_tsv(results_dir / "mirna" / "mirna_analysis.tsv", MIRNA_COLUMNS, [_mirna_row(), _mirna_row()])

    parsed = asyncio.run(workflow._parse_nextflow_results(results_dir))

    assert len(parsed["results"]["cand_1"]["hits"]) == 2, "two rows in the file must count as two hits"


@pytest.mark.unit
def test_mirna_counters_report_the_rows_the_file_holds(tmp_path):
    """The candidate-level and run-level miRNA counters must equal the file's row count, not twice it."""
    workflow = _workflow(tmp_path, "mirna_counters")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    _write_tsv(results_dir / "mirna" / "mirna_analysis.tsv", MIRNA_COLUMNS, [_mirna_row()])

    candidate = _candidate()
    outcome = asyncio.run(workflow._process_nextflow_results([candidate], results_dir, {"status": "completed"}))

    assert candidate.mirna_hits_total == 1
    assert outcome["filtering_stats"]["human_mirna_hits"] == 1
