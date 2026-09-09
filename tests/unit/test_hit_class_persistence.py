"""The four-way class must reach the hit table, and reach it exactly once.

``combined_offtargets.tsv`` used to carry 12 columns and no class, so an on-target isoform
alignment and a genuine liability rendered identically. These tests pin the three added columns
and, more importantly, pin that the per-candidate counters are *derived from* the persisted rows:
two independent computations of the same quantity is the defect this package exists to prevent.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

import asyncio
import csv
import json
from pathlib import Path

import pytest

from sirnaforge.core.hit_annotation import (
    CLASSIFICATION_COLUMNS,
    HIT_CLASS_COLUMN,
    MATCHED_SYMBOL_COLUMN,
    SYMBOL_LOOKUP_MISSING_COLUMN,
    UNKNOWN_SYMBOL,
    accumulate_hit_class,
    annotate_hit_row,
    write_classified_hits,
)
from sirnaforge.core.hit_classification import HitClass, HitClassCounts, HitClassification
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

CLEAN_GUIDE = "ATGCGATGCGATGCGATGCGC"
REPEAT_GUIDE = "TTGCGATGCGATGCGATGCGA"
ORPHAN_GUIDE = "GGGCGATGCGATGCGATGCGT"

TSV_COLUMNS = [
    "qname",
    "qseq",
    "species",
    "rname",
    "coord",
    "strand",
    "cigar",
    "mapq",
    "as_score",
    "nm",
    "seed_mismatches",
    "offtarget_score",
]


def _candidate(candidate_id: str, guide: str, *, repeat_flagged: bool = False) -> SiRNACandidate:
    """A minimal candidate carrying one guide."""
    candidate = SiRNACandidate(
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
    candidate.repeat_flagged = repeat_flagged
    return candidate


def _hit_row(qname: str, qseq: str, species: str, rname: str, nm: int = 0) -> dict[str, str]:
    """One aggregated alignment row, in the aggregator's own column order."""
    return {
        "qname": qname,
        "qseq": qseq,
        "species": species,
        "rname": rname,
        "coord": "500",
        "strand": "-",
        "cigar": f"{len(qseq)}M",
        "mapq": "60",
        "as_score": "30",
        "nm": str(nm),
        "seed_mismatches": str(nm),
        "offtarget_score": "1.5",
    }


def _write_results_dir(results_dir: Path, rows: list[dict[str, str]]) -> Path:
    """Lay out an aggregated Nextflow results directory around one hit table."""
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True, exist_ok=True)
    tsv_path = aggregated / "combined_offtargets.tsv"
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (aggregated / "combined_summary.json").write_text(
        json.dumps(
            {
                "species_analyzed": ["human", "mouse", "rat"],
                "species_file_counts": {"human": 1, "mouse": 1, "rat": 1},
                "missing_species": [],
                "total_results": len(rows),
            }
        )
    )
    return tsv_path


def _workflow(tmp_path: Path, name: str, *, with_symbol_index: bool = True) -> SiRNAWorkflow:
    """A TP53 workflow with a synthetic human/mouse index, or with no index at all."""
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        genome_species=["human", "mouse", "rat"],
        design_params=DesignParameters(),
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_ids = {"ENSG00000000001"}
    workflow._query_gene_symbols = {"TP53"}
    if not with_symbol_index:
        return workflow

    human_reference = tmp_path / f"{name}_human_cdna.fasta"
    human_reference.write_text(
        ">ENST00000000001.2 gene:ENSG00000000001.3 gene_symbol:TP53 transcript_biotype:protein_coding\nACGT\n"
        ">ENST00000000002 gene:ENSG00000000001 gene_symbol:TP53 transcript_biotype:protein_coding\nACGT\n"
        ">ENST00000000009 gene:ENSG00000000009 gene_symbol:OTHER transcript_biotype:protein_coding\nACGT\n"
    )
    mouse_reference = tmp_path / f"{name}_mouse_cdna.fasta"
    mouse_reference.write_text(
        ">ENSMUST00000000002 gene:ENSMUSG00000000002 gene_symbol:Tp53 transcript_biotype:protein_coding\nACGT\n"
        ">ENSMUST00000000099 gene:ENSMUSG00000000099 transcript_biotype:protein_coding\nACGT\n"
    )
    workflow._transcript_index.build("human", human_reference)
    workflow._transcript_index.build("mouse", mouse_reference)
    return workflow


def _mixed_rows() -> list[dict[str, str]]:
    """Rows covering all four classes, a thin annotation, a missing index and an orphan qname."""
    return [
        # Own gene: the enumerated isoform and a sibling only the index recognises.
        _hit_row("cand_clean", CLEAN_GUIDE, "human", "ENST00000000001.2"),
        _hit_row("cand_clean", CLEAN_GUIDE, "human", "ENST00000000002"),
        # An unrelated human transcript: a genuine liability.
        _hit_row("cand_clean", CLEAN_GUIDE, "human", "ENST00000000009"),
        # The mouse ortholog, spelled Tp53.
        _hit_row("cand_clean", CLEAN_GUIDE, "mouse", "ENSMUST00000000002"),
        # A mouse transcript the index carries with no symbol: a thin annotation, not a verdict.
        _hit_row("cand_clean", CLEAN_GUIDE, "mouse", "ENSMUST00000000099"),
        # Rat has no index at all: a reference gap, a different problem from the row above.
        _hit_row("cand_clean", CLEAN_GUIDE, "rat", "ENSRNOT00000000123"),
        # A repeat-flagged guide's unrelated hit is repeat-mediated, not a liability.
        _hit_row("cand_repeat", REPEAT_GUIDE, "human", "ENST00000000009"),
        # A row from a previous run whose candidate is not in this one.
        _hit_row("cand_gone", ORPHAN_GUIDE, "human", "ENST00000000009"),
    ]


def _run(tmp_path: Path, name: str, *, with_symbol_index: bool = True) -> tuple[list[SiRNACandidate], Path, dict]:
    """Screen two candidates against the mixed hit table and return the rewritten TSV."""
    workflow = _workflow(tmp_path, name, with_symbol_index=with_symbol_index)
    results_dir = workflow.config.output_dir / "off_target" / "results"
    tsv_path = _write_results_dir(results_dir, _mixed_rows())

    candidates = [
        _candidate("cand_clean", CLEAN_GUIDE),
        _candidate("cand_repeat", REPEAT_GUIDE, repeat_flagged=True),
    ]
    outcome = asyncio.run(workflow._process_nextflow_results(candidates, results_dir, {"status": "completed"}))
    return candidates, tsv_path, outcome


def _read_tsv(tsv_path: Path) -> list[dict[str, str]]:
    """Read the rewritten hit table back."""
    with tsv_path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.unit
def test_hit_table_gains_the_three_classification_columns(tmp_path):
    """The class, the resolved symbol and the lookup shortfall all reach the TSV."""
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_columns")
    rows = _read_tsv(tsv_path)

    assert len(rows) == len(_mixed_rows())
    for column in (HIT_CLASS_COLUMN, MATCHED_SYMBOL_COLUMN, SYMBOL_LOOKUP_MISSING_COLUMN):
        assert column in rows[0]
    # The class is one of exactly four values on every row, including the orphan.
    assert {row[HIT_CLASS_COLUMN] for row in rows} <= {member.value for member in HitClass}
    assert all(row[HIT_CLASS_COLUMN] for row in rows)


@pytest.mark.unit
def test_original_columns_survive_the_rewrite(tmp_path):
    """Adding columns must not disturb the twelve that were already there."""
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_schema")
    rows = _read_tsv(tsv_path)

    assert list(rows[0].keys())[: len(TSV_COLUMNS)] == TSV_COLUMNS
    for original, rewritten in zip(_mixed_rows(), rows, strict=True):
        assert {column: rewritten[column] for column in TSV_COLUMNS} == original


@pytest.mark.unit
def test_candidate_counters_are_derived_from_the_persisted_rows(tmp_path):
    """The row-level and candidate-level views are one computation, so they cannot disagree.

    This is the assertion the deliverable exists for: every per-candidate hit-class counter is
    reproduced by tallying that candidate's rows in the TSV.
    """
    candidates, tsv_path, _outcome = _run(tmp_path, "out_derived")
    rows = _read_tsv(tsv_path)

    for candidate in candidates:
        own_rows = [row for row in rows if row["qname"] == candidate.id]
        tally = {member.value: 0 for member in HitClass}
        for row in own_rows:
            tally[row[HIT_CLASS_COLUMN]] += 1

        assert candidate.on_target_hits == tally[HitClass.ON_TARGET.value]
        assert candidate.ortholog_hits == tally[HitClass.ORTHOLOG.value]
        assert candidate.repeat_hits == tally[HitClass.REPEAT.value]
        assert candidate.off_target_count == tally[HitClass.OFF_TARGET.value]
        # transcriptome_hits_total counts genuine off-targets at any mismatch level, which is the
        # same population as the off_target rows.
        assert candidate.transcriptome_hits_total == tally[HitClass.OFF_TARGET.value]


@pytest.mark.unit
def test_on_target_and_ortholog_rows_are_labelled_not_counted_as_off_targets(tmp_path):
    """An on-target isoform alignment must be visibly not a liability."""
    candidates, tsv_path, _outcome = _run(tmp_path, "out_taxonomy")
    rows = _read_tsv(tsv_path)
    clean = next(candidate for candidate in candidates if candidate.id == "cand_clean")

    by_rname = {row["rname"]: row for row in rows if row["qname"] == "cand_clean"}
    assert by_rname["ENST00000000001.2"][HIT_CLASS_COLUMN] == HitClass.ON_TARGET.value
    assert by_rname["ENST00000000002"][HIT_CLASS_COLUMN] == HitClass.ON_TARGET.value
    assert by_rname["ENSMUST00000000002"][HIT_CLASS_COLUMN] == HitClass.ORTHOLOG.value
    assert by_rname["ENSMUST00000000002"][MATCHED_SYMBOL_COLUMN] == "TP53"
    assert by_rname["ENST00000000009"][HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value
    assert clean.on_target_hits == 2
    assert clean.ortholog_hits == 1

    repeat_row = next(row for row in rows if row["qname"] == "cand_repeat")
    assert repeat_row[HIT_CLASS_COLUMN] == HitClass.REPEAT.value


@pytest.mark.unit
def test_missing_symbol_is_the_literal_unknown(tmp_path):
    """A symbol we do not have is ``unknown``, never a blank cell that renders as "no gene"."""
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_unknown")
    rows = _read_tsv(tsv_path)

    assert all(row[MATCHED_SYMBOL_COLUMN] for row in rows)
    unresolved = [row for row in rows if row[MATCHED_SYMBOL_COLUMN] == UNKNOWN_SYMBOL]
    assert unresolved, "at least one synthetic row has no resolvable symbol"
    assert all(row[MATCHED_SYMBOL_COLUMN] != "" for row in unresolved)


@pytest.mark.unit
def test_thin_annotation_and_missing_index_are_reported_separately(tmp_path):
    """A transcript with no symbol is a lookup miss; a species with no index is a reference gap."""
    _candidates, tsv_path, outcome = _run(tmp_path, "out_shortfall")
    rows = _read_tsv(tsv_path)

    thin = next(row for row in rows if row["rname"] == "ENSMUST00000000099")
    assert thin[SYMBOL_LOOKUP_MISSING_COLUMN] == "True"
    assert thin[MATCHED_SYMBOL_COLUMN] == UNKNOWN_SYMBOL

    # Rat has no index, so the ortholog check never ran; that is not a symbol lookup miss.
    no_index = next(row for row in rows if row["species"] == "rat")
    assert no_index[SYMBOL_LOOKUP_MISSING_COLUMN] == "False"

    stats = outcome["filtering_stats"]
    assert stats["ortholog_symbol_lookup_misses"] == 1
    assert stats["species_index_misses"] == 1


@pytest.mark.unit
def test_runs_without_a_transcript_symbol_cache(tmp_path):
    """With no index at all, every symbol is ``unknown`` and the run still completes."""
    candidates, tsv_path, outcome = _run(tmp_path, "out_no_cache", with_symbol_index=False)
    rows = _read_tsv(tsv_path)

    assert outcome["status"] == "completed"
    assert all(row[MATCHED_SYMBOL_COLUMN] == UNKNOWN_SYMBOL for row in rows)
    # Nothing can be recognised as the query gene beyond the enumerated transcript ID, and the
    # counters still match the rows they were derived from.
    clean = next(candidate for candidate in candidates if candidate.id == "cand_clean")
    own = [row for row in rows if row["qname"] == "cand_clean"]
    assert clean.on_target_hits == sum(row[HIT_CLASS_COLUMN] == HitClass.ON_TARGET.value for row in own)
    assert clean.off_target_count == sum(row[HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value for row in own)


@pytest.mark.unit
def test_rows_belonging_to_no_candidate_are_still_classified(tmp_path):
    """A stale row must not publish a blank class, and must not enter candidate counters."""
    candidates, tsv_path, outcome = _run(tmp_path, "out_orphan")
    rows = _read_tsv(tsv_path)

    orphan = next(row for row in rows if row["qname"] == "cand_gone")
    assert orphan[HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value
    assert outcome["filtering_stats"]["hits_classified_without_candidate"] == 1
    assert all(candidate.id != "cand_gone" for candidate in candidates)


@pytest.mark.unit
def test_persisting_is_idempotent_across_repeated_reads(tmp_path):
    """Re-reading and re-writing an already-classified table changes nothing.

    The classification columns are appended once, not once per pass, so a re-parsed run directory
    does not accumulate duplicate headers.
    """
    workflow = _workflow(tmp_path, "out_idempotent")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    tsv_path = _write_results_dir(results_dir, _mixed_rows())
    candidates = [
        _candidate("cand_clean", CLEAN_GUIDE),
        _candidate("cand_repeat", REPEAT_GUIDE, repeat_flagged=True),
    ]

    asyncio.run(workflow._process_nextflow_results(candidates, results_dir, {"status": "completed"}))
    first = tsv_path.read_text()
    asyncio.run(workflow._process_nextflow_results(candidates, results_dir, {"status": "completed"}))

    assert tsv_path.read_text() == first


@pytest.mark.unit
def test_counters_follow_the_persisted_row_not_the_classification_object():
    """The single-sourcing invariant, stated directly.

    ``accumulate_hit_class`` reads the row it was handed, so there is no second path by which a
    counter could be incremented for a class the table does not show.
    """
    row: dict[str, object] = {"rname": "ENST00000000009"}
    annotate_hit_row(
        row,
        HitClassification(hit_class=HitClass.ON_TARGET, matched_symbol="TP53", symbol_lookup_missing=False),
    )
    row[HIT_CLASS_COLUMN] = HitClass.OFF_TARGET.value

    counts = HitClassCounts()
    bucket: dict[str, int] = {}
    assert accumulate_hit_class(row, counts, bucket, "human") is HitClass.OFF_TARGET
    assert counts.off_target == 1
    assert counts.on_target == 0
    assert bucket == {"off_target": 1}


@pytest.mark.unit
def test_writer_appends_the_columns_once_when_they_are_already_present(tmp_path):
    """A re-persisted table keeps one header per column, whatever order it arrives in."""
    rows = [
        {
            "qname": "cand_clean",
            HIT_CLASS_COLUMN: HitClass.OFF_TARGET.value,
            MATCHED_SYMBOL_COLUMN: UNKNOWN_SYMBOL,
            SYMBOL_LOOKUP_MISSING_COLUMN: "True",
        }
    ]
    tsv_path = tmp_path / "combined_offtargets.tsv"
    written = write_classified_hits(tsv_path, rows, ["qname", *CLASSIFICATION_COLUMNS])

    assert written == 1
    header = tsv_path.read_text().splitlines()[0].split("\t")
    assert header == ["qname", *CLASSIFICATION_COLUMNS]
