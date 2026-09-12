"""The hit class must reach the hit table, reach every table read, and reach it exactly once.

``combined_offtargets.tsv`` used to carry 12 columns and no class, so an on-target isoform
alignment and a genuine liability rendered identically. These tests pin the added columns and, more
importantly, pin that the per-candidate counters are *derived from* the persisted rows: two
independent computations of the same quantity is the defect this package exists to prevent. The
fallback tests pin the reachable path on which that invariant used to be false — real per-species
rows fed the candidate counters while the published table stayed header-only.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

import asyncio
import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from sirnaforge.core.hit_annotation import (
    CLASSIFICATION_COLUMNS,
    HIT_CLASS_COLUMN,
    HIT_SYMBOL_COLUMN,
    HIT_SYMBOL_MISSING_COLUMN,
    MATCHED_SYMBOL_COLUMN,
    ORTHOLOG_EVIDENCE_CELL_VALUES,
    ORTHOLOG_EVIDENCE_COLUMN,
    ORTHOLOG_EVIDENCE_NOT_APPLICABLE,
    SPECIES_INDEX_MISSING_COLUMN,
    SYMBOL_LOOKUP_MISSING_COLUMN,
    UNKNOWN_SYMBOL,
    HitAnnotator,
    accumulate_hit_class,
    annotate_hit_row,
    write_classified_hits,
)
from sirnaforge.core.hit_classification import (
    ClassificationContext,
    HitClass,
    HitClassCounts,
    HitClassification,
    OrthologEvidence,
    classify_hit,
)
from sirnaforge.data.transcript_index import TranscriptGeneIndex
from sirnaforge.models.schemas import HIT_CLASS_VALUES, AggregatedOffTargetSchema
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate, build_candidate_row
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


#: Offline orthologue evidence, so no test here touches Ensembl Compara. See tests/unit/data/README.md.
_ORTHOLOG_MAPPING_FIXTURE = Path(__file__).resolve().parent / "data" / "ortholog_mapping_synthetic.json"


def _repo_root() -> Path:
    """Repository root, for the published-contract assertions."""
    return Path(__file__).resolve().parents[2]


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
    """One aggregated alignment row, in the aggregator's own column order.

    ``offtarget_score`` follows the aligner's convention -- 0.0 is reserved for a perfect match --
    so these rows satisfy ``GenomeAlignmentSchema`` and its aggregated subclass.
    """
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
        "offtarget_score": "0.0" if nm == 0 else "1.5",
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
    _write_summary(aggregated, len(rows))
    return tsv_path


def _write_summary(aggregated: Path, total: int) -> None:
    """The aggregator's own summary, which is the positive evidence that alignments ran."""
    (aggregated / "combined_summary.json").write_text(
        json.dumps(
            {
                "species_analyzed": ["human", "mouse", "rat"],
                "species_file_counts": {"human": 1, "mouse": 1, "rat": 1},
                "missing_species": [],
                "total_results": total,
            }
        )
    )


def _write_fallback_results_dir(results_dir: Path, rows: list[dict[str, str]]) -> tuple[Path, Path]:
    """Lay out the shape that broke the invariant.

    ``aggregate_offtarget_results`` writes a header-only ``combined_offtargets.tsv`` when
    ``GenomeAlignmentSchema`` rejects every per-species file. The real rows are still published
    under ``transcriptome/``, and the parser's fallback is what reads them.

    Returns:
        (header-only aggregated table, per-species table carrying the rows)
    """
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True, exist_ok=True)
    empty_path = aggregated / "combined_offtargets.tsv"
    with empty_path.open("w", newline="") as handle:
        csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n").writeheader()
    _write_summary(aggregated, len(rows))

    per_species = results_dir / "transcriptome"
    per_species.mkdir(parents=True, exist_ok=True)
    species_path = per_species / "human_analysis.tsv"
    with species_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return empty_path, species_path


def _workflow(tmp_path: Path, name: str, *, with_symbol_index: bool = True) -> SiRNAWorkflow:
    """A TP53 workflow with a synthetic human/mouse index, or with no index at all.

    ``ortholog_mapping_file`` is what keeps these tests off the network. Without it,
    ``_process_nextflow_results`` resolves mouse and rat orthologues over Ensembl Compara, which
    behind a TLS-intercepting proxy cost ~25s of retry backoff per call and made this file 90% of
    the dev tier's runtime. The mapping resolves mouse and says nothing about rat, so both the
    resolved and the unresolved branch still run here (#101 point 4: the fixture uses the file).
    """
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        screen_species=["human", "mouse", "rat"],
        design_params=DesignParameters(),
        ortholog_mapping_file=_ORTHOLOG_MAPPING_FIXTURE,
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
    """Every per-candidate hit-class counter is reproduced by tallying that candidate's rows.

    This is the assertion the deliverable exists for. Note the scope: this fixture always publishes
    what it feeds, so it pins the derivation, not the reconciliation. Per-candidate mis-attribution
    that nets to zero across the table is checked by
    ``test_the_returned_results_map_is_looked_up_by_the_screening_id``, and an unpublished counted hit
    by ``test_a_hit_that_is_never_published_is_reported_as_a_reconciliation_failure``.
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
        assert candidate.undetermined_hits == tally[HitClass.UNDETERMINED.value]
        # off_target_count is the liability population: confirmed off-targets plus the rows whose
        # class could not be decided. transcriptome_hits_total counts the same rows at any
        # mismatch level, so all three views are one tally of the published table.
        liabilities = tally[HitClass.OFF_TARGET.value] + tally[HitClass.UNDETERMINED.value]
        assert candidate.off_target_count == liabilities
        assert candidate.transcriptome_hits_total == liabilities


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
def test_orthology_evidence_comes_from_the_mapping_file_not_the_network(tmp_path):
    """The offline path is the one these tests take, and it publishes gene-ID evidence (#101).

    Every test in this file drives ``_process_nextflow_results`` with three species and cross-species
    hit rows, which used to resolve orthologues over Ensembl Compara -- ~25s of retry backoff each
    behind a TLS-intercepting proxy, and a different verdict depending on whether the network
    answered. The mapping file states the mouse orthologue outright, so the mouse row earns GENE_ID
    evidence instead of the symbol heuristic. Rat is absent from the file, so it stays *unresolved*:
    a species nothing was claimed about, not a checked absence.
    """
    _candidates, tsv_path, outcome = _run(tmp_path, "out_offline_orthology")
    rows = _read_tsv(tsv_path)

    mouse = next(row for row in rows if row["rname"] == "ENSMUST00000000002")
    assert mouse[HIT_CLASS_COLUMN] == HitClass.ORTHOLOG.value
    assert mouse[ORTHOLOG_EVIDENCE_COLUMN] == OrthologEvidence.GENE_ID.value, "the file resolves a gene ID"

    provenance = outcome["filtering_stats"]["orthology"]
    assert provenance["source"] == "ortholog_mapping_file", "an offline run must not claim a REST call"
    assert provenance["gene_ids_by_species"] == {"mouse": ["ENSMUSG00000000002"]}
    assert provenance["unresolved_species"] == ["rat"], "a species the file omits is unchecked, not clean"


@pytest.mark.unit
def test_missing_symbol_is_the_literal_unknown(tmp_path):
    """A symbol we do not have is ``unknown``, never a blank cell that renders as "no gene".

    The row named here is the one the index genuinely cannot resolve, so the assertion cannot be
    satisfied by a row that merely had no symbol *consulted*.
    """
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_unknown")
    rows = _read_tsv(tsv_path)

    assert all(row[MATCHED_SYMBOL_COLUMN] for row in rows)
    unresolvable = next(row for row in rows if row["rname"] == "ENSMUST00000000099")
    assert unresolvable[MATCHED_SYMBOL_COLUMN] == UNKNOWN_SYMBOL
    assert unresolvable[SYMBOL_LOOKUP_MISSING_COLUMN] == "True"


@pytest.mark.unit
def test_matched_symbol_is_class_evidence_not_a_per_hit_gene_name(tmp_path):
    """``matched_symbol`` names the symbol that established the class, and nothing more.

    An ``off_target`` or ``repeat`` row reads ``unknown`` even where the hit species' index does
    resolve a symbol for that transcript, and ``symbol_lookup_missing`` is ``False`` there, so
    nothing on the row says the gene name is absent. The published contract must say so rather than
    promise a resolved gene symbol per row — a consumer that renders ``unknown`` as "unannotated"
    would mislabel fully annotated genes. Widening the column is #101's call on
    ``core/hit_classification.py``.
    """
    workflow = _workflow(tmp_path, "out_symbol_semantics")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    tsv_path = _write_results_dir(results_dir, _mixed_rows())
    candidates = [
        _candidate("cand_clean", CLEAN_GUIDE),
        _candidate("cand_repeat", REPEAT_GUIDE, repeat_flagged=True),
    ]
    asyncio.run(workflow._process_nextflow_results(candidates, results_dir, {"status": "completed"}))
    rows = _read_tsv(tsv_path)

    human_index = workflow._transcript_index.for_species("human")
    assert human_index is not None
    assert human_index.symbol_for("ENST00000000009") == "OTHER"

    resolvable_but_unknown = [
        row
        for row in rows
        if row["species"] == "human"
        and row[MATCHED_SYMBOL_COLUMN] == UNKNOWN_SYMBOL
        and human_index.symbol_for(row["rname"]) is not None
    ]
    assert resolvable_but_unknown, "fixture must exercise a row whose symbol the index resolves"
    assert {row[HIT_CLASS_COLUMN] for row in resolvable_but_unknown} == {
        HitClass.ON_TARGET.value,
        HitClass.REPEAT.value,
        HitClass.OFF_TARGET.value,
    }
    assert all(row[SYMBOL_LOOKUP_MISSING_COLUMN] == "False" for row in resolvable_but_unknown)

    for relative in ("docs/cli_reference.md", "CHANGELOG.md"):
        text = (_repo_root() / relative).read_text()
        assert "matched_symbol` is the resolved gene symbol" not in text, (
            f"{relative} promises a resolved gene symbol per row, which no off_target row carries"
        )
        assert "per-hit gene name" in text, f"{relative} must state that matched_symbol is not one"


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
    assert all(row[HIT_SYMBOL_COLUMN] == UNKNOWN_SYMBOL for row in rows)
    # Nothing can be recognised as the query gene beyond the enumerated transcript ID, and the
    # counters still match the rows they were derived from.
    clean = next(candidate for candidate in candidates if candidate.id == "cand_clean")
    own = [row for row in rows if row["qname"] == "cand_clean"]
    assert clean.on_target_hits == sum(row[HIT_CLASS_COLUMN] == HitClass.ON_TARGET.value for row in own)
    liabilities = sum(row[HIT_CLASS_COLUMN] in {HitClass.OFF_TARGET.value, HitClass.UNDETERMINED.value} for row in own)
    assert clean.off_target_count == liabilities


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
    row: dict[str, object] = {"rname": "ENST00000000009", "species": "human"}
    annotate_hit_row(
        row,
        HitClassification(hit_class=HitClass.ON_TARGET, matched_symbol="TP53", symbol_lookup_missing=False),
        HitAnnotator(index=TranscriptGeneIndex(), query_species="human"),
    )
    row[HIT_CLASS_COLUMN] = HitClass.OFF_TARGET.value

    counts = HitClassCounts()
    bucket: dict[str, int] = {}
    assert accumulate_hit_class(row, counts, bucket, "human") is HitClass.OFF_TARGET
    assert counts.off_target == 1
    assert counts.on_target == 0
    # The shortfall flags are read off the row too; this row was annotated with no index at all.
    assert bucket == {"off_target": 1, SPECIES_INDEX_MISSING_COLUMN: 1}


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


@pytest.mark.unit
def test_zero_hit_table_still_gains_the_classification_columns(tmp_path):
    """A header-only hit table is published with the same schema as a populated one.

    miRNA-only mode and a run whose every species index is missing both publish an aggregated
    table with a header and no rows. If the three columns were conditional on there being a hit,
    a consumer selecting ``hit_class`` would fail on exactly the runs with no liabilities.
    """
    workflow = _workflow(tmp_path, "out_zero_hits")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    tsv_path = _write_results_dir(results_dir, [])
    candidates = [_candidate("cand_clean", CLEAN_GUIDE)]

    asyncio.run(workflow._process_nextflow_results(candidates, results_dir, {"status": "completed"}))

    lines = tsv_path.read_text().splitlines()
    assert lines[0].split("\t") == [*TSV_COLUMNS, *CLASSIFICATION_COLUMNS]
    assert len(lines) == 1, "no data rows may be invented for an empty table"
    assert _read_tsv(tsv_path) == []
    assert candidates[0].off_target_count == 0


# --- P1: the fallback path, on which the headline invariant used to be false ------------------


@pytest.mark.unit
def test_fallback_per_species_rows_reach_the_published_table(tmp_path):
    """Rows the fallback ingests must be published with their class, not counted and discarded.

    When the aggregated table is header-only, the parser globs ``transcriptome/*_analysis.tsv``. Those
    rows used to reach ``off_target_count`` while the published table was rewritten header-only,
    so the hit table authoritatively reported no liabilities beside candidates carrying dozens.
    """
    workflow = _workflow(tmp_path, "out_fallback")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    empty_path, species_path = _write_fallback_results_dir(results_dir, _mixed_rows())
    candidates = [
        _candidate("cand_clean", CLEAN_GUIDE),
        _candidate("cand_repeat", REPEAT_GUIDE, repeat_flagged=True),
    ]

    outcome = asyncio.run(workflow._process_nextflow_results(candidates, results_dir, {"status": "completed"}))

    published = _read_tsv(species_path)
    assert len(published) == len(_mixed_rows())
    assert all(row[HIT_CLASS_COLUMN] in {member.value for member in HitClass} for row in published)

    clean = next(candidate for candidate in candidates if candidate.id == "cand_clean")
    own = [row for row in published if row["qname"] == "cand_clean"]
    liabilities = {HitClass.OFF_TARGET.value, HitClass.UNDETERMINED.value}
    assert clean.off_target_count > 0, "the fallback rows are real evidence and must be counted"
    assert clean.off_target_count == sum(row[HIT_CLASS_COLUMN] in liabilities for row in own)

    # The header-only aggregated table is still republished with the schema, and carries no rows.
    assert _read_tsv(empty_path) == []
    assert empty_path.read_text().splitlines()[0].split("\t") == [*TSV_COLUMNS, *CLASSIFICATION_COLUMNS]
    assert outcome["warnings"] == [] or all("mismatch" not in warning for warning in outcome["warnings"])


@pytest.mark.unit
def test_a_hit_that_is_never_published_is_reported_as_a_reconciliation_failure(tmp_path):
    """The check that would have caught P1, exercised on the pre-fix shape directly.

    A parsed payload whose candidate counters were fed by rows no table holds is exactly what the
    fallback used to produce. It must be loud, not silent.
    """
    parsed = {
        "results": {
            "cand_clean": {
                "off_target_count": 2,
                "off_target_score": 0.0,
                "hits": [_hit_row("cand_clean", CLEAN_GUIDE, "human", "ENST00000000009")] * 2,
            },
            # A miRNA hit has no class of its own and must not be counted against the hit table.
            "cand_repeat": {"hits": [{"qname": "cand_repeat", "mirna_id": "hsa-miR-1"}]},
        },
        "transcriptome_hit_tables": [],
    }

    warnings = SiRNAWorkflow._persist_hit_classifications(parsed)

    assert len(warnings) == 1
    assert "2 transcriptome hit(s) reached the candidate counters" in warnings[0]
    assert "only 0 row(s) were published" in warnings[0]


@pytest.mark.unit
def test_reconciliation_is_silent_when_every_counted_hit_was_published(tmp_path):
    """The aggregated path must not raise a false alarm on a run that is in fact consistent."""
    _candidates, _tsv_path, outcome = _run(tmp_path, "out_reconciled")

    assert all("mismatch" not in warning for warning in outcome["warnings"])


# --- P2: a resolved gene symbol per row, independent of the class ------------------------------


@pytest.mark.unit
def test_hit_symbol_is_resolved_per_row_independent_of_hit_class(tmp_path):
    """``hit_symbol`` names the gene the alignment landed on, whatever the class.

    ``matched_symbol`` is ``unknown`` on 100% of ``off_target`` and ``repeat`` rows by
    construction, so it cannot group anything by gene. ``hit_symbol`` is the column that can.
    """
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_hit_symbol")
    rows = _read_tsv(tsv_path)
    by_key = {(row["qname"], row["rname"]): row for row in rows}

    liability = by_key[("cand_clean", "ENST00000000009")]
    assert liability[HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value
    assert liability[MATCHED_SYMBOL_COLUMN] == UNKNOWN_SYMBOL, "no symbol established this class"
    assert liability[HIT_SYMBOL_COLUMN] == "OTHER", "but the gene it landed on is known"
    assert liability[HIT_SYMBOL_MISSING_COLUMN] == "False"

    repeat = by_key[("cand_repeat", "ENST00000000009")]
    assert repeat[HIT_CLASS_COLUMN] == HitClass.REPEAT.value
    assert repeat[HIT_SYMBOL_COLUMN] == "OTHER"

    on_target = by_key[("cand_clean", "ENST00000000001.2")]
    assert on_target[HIT_CLASS_COLUMN] == HitClass.ON_TARGET.value
    assert on_target[HIT_SYMBOL_COLUMN] == "TP53"

    # The published contract must point a consumer at the column that can group by gene.
    for relative in ("docs/cli_reference.md", "CHANGELOG.md"):
        text = (_repo_root() / relative).read_text()
        assert HIT_SYMBOL_COLUMN in text, f"{relative} does not document {HIT_SYMBOL_COLUMN}"


@pytest.mark.unit
def test_an_unresolvable_hit_symbol_is_unknown_and_flagged(tmp_path):
    """A transcript with no symbol, and a species with no index, both read ``unknown`` and say so.

    14.6% of reference transcripts carry no symbol at all, so ``unknown`` is a common real state
    and must be distinguishable from a resolved absence rather than rendering as "no gene".
    """
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_hit_symbol_unknown")
    rows = _read_tsv(tsv_path)

    thin = next(row for row in rows if row["rname"] == "ENSMUST00000000099")
    assert thin[HIT_SYMBOL_COLUMN] == UNKNOWN_SYMBOL
    assert thin[HIT_SYMBOL_MISSING_COLUMN] == "True"

    no_index = next(row for row in rows if row["species"] == "rat")
    assert no_index[HIT_SYMBOL_COLUMN] == UNKNOWN_SYMBOL
    assert no_index[HIT_SYMBOL_MISSING_COLUMN] == "True"
    assert all(row[HIT_SYMBOL_COLUMN] for row in rows), "never an empty cell"


# --- P4: undetermined, and the persisted reason ------------------------------------------------


@pytest.mark.unit
def test_a_species_with_no_index_is_undetermined_rather_than_an_off_target_verdict(tmp_path):
    """A class we could not tell must not be spelled as a liability verdict."""
    candidates, tsv_path, outcome = _run(tmp_path, "out_undetermined")
    rows = _read_tsv(tsv_path)

    rat = next(row for row in rows if row["species"] == "rat")
    assert rat[HIT_CLASS_COLUMN] == HitClass.UNDETERMINED.value
    assert rat[SPECIES_INDEX_MISSING_COLUMN] == "True"

    human_liability = next(row for row in rows if row["rname"] == "ENST00000000009" and row["qname"] == "cand_clean")
    assert human_liability[HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value
    assert human_liability[SPECIES_INDEX_MISSING_COLUMN] == "False"

    clean = next(candidate for candidate in candidates if candidate.id == "cand_clean")
    assert clean.undetermined_hits == 1
    assert outcome["filtering_stats"]["hit_classes"]["undetermined"] == 1


@pytest.mark.unit
def test_an_undetermined_hit_still_counts_as_a_liability(tmp_path):
    """Absent evidence must not loosen the screen: no gate moves because the class exists.

    Excluding undetermined hits from ``off_target_count`` would mean a run with no transcript
    index passed candidates that the same run fails today, i.e. deleting a reference would
    improve the results.
    """
    candidates, tsv_path, _outcome = _run(tmp_path, "out_undetermined_gated")
    rows = _read_tsv(tsv_path)
    clean = next(candidate for candidate in candidates if candidate.id == "cand_clean")

    own = [row for row in rows if row["qname"] == "cand_clean"]
    confirmed = sum(row[HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value for row in own)
    undecided = sum(row[HIT_CLASS_COLUMN] == HitClass.UNDETERMINED.value for row in own)
    assert undecided == 1, "fixture must exercise an undecidable hit"
    assert clean.off_target_count == confirmed + undecided
    assert clean.transcriptome_hits_total == confirmed + undecided


@pytest.mark.unit
def test_a_table_built_with_no_index_is_distinguishable_from_a_screened_one(tmp_path):
    """The defect P4 names: with no index at all, every row used to read ``off_target``."""
    _screened, screened_path, _outcome = _run(tmp_path, "out_with_index")
    _blind, blind_path, _blind_outcome = _run(tmp_path, "out_no_index", with_symbol_index=False)

    blind_rows = _read_tsv(blind_path)
    assert {row[HIT_CLASS_COLUMN] for row in blind_rows} <= {
        HitClass.UNDETERMINED.value,
        HitClass.ON_TARGET.value,
        HitClass.REPEAT.value,
    }
    assert all(row[SPECIES_INDEX_MISSING_COLUMN] == "True" for row in blind_rows)
    assert not any(row[HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value for row in blind_rows), (
        "an unqualified off_target verdict is exactly what a run with no reference must not publish"
    )
    assert screened_path.read_text() != blind_path.read_text()


@pytest.mark.unit
def test_undetermined_is_never_produced_by_the_pure_classifier():
    """The class is the annotation layer's decision: ``classify_hit`` cannot see the inventory."""
    context = ClassificationContext(
        query_gene_ids=frozenset(),
        query_gene_symbols=frozenset(),
        on_target_transcript_ids=frozenset(),
        query_species="human",
        index=TranscriptGeneIndex(),
        repeat_flagged_guides=frozenset(),
        requested_species=frozenset({"human"}),
    )
    verdict = classify_hit({"rname": "ENST00000000009", "species": "rat"}, CLEAN_GUIDE, context)

    assert verdict.hit_class is HitClass.OFF_TARGET
    annotated: dict[str, object] = {"rname": "ENST00000000009", "species": "rat"}
    assert (
        annotate_hit_row(annotated, verdict, HitAnnotator(index=context.index, query_species="human"))
        is HitClass.UNDETERMINED
    )


# --- P6: a schema for the published table, and its two shapes ---------------------------------


@pytest.mark.unit
def test_the_published_table_validates_against_its_own_schema(tmp_path):
    """The 18-column published shape must validate, which the 12-column strict schema rejects."""
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_schema_valid")
    frame = pd.read_csv(tsv_path, sep="\t", dtype=str)

    validated = AggregatedOffTargetSchema.validate(frame, lazy=True)
    assert len(validated) == len(_mixed_rows())
    assert HIT_CLASS_COLUMN in validated.columns


@pytest.mark.unit
def test_the_twelve_column_producer_shape_also_validates(tmp_path):
    """A table written before the producer owned the classification columns has only 12.

    The producer writes all 19 now, so both entry points agree on the shape (#100). The narrow
    shape stays valid because the per-species ``transcriptome/*_analysis.tsv`` files carry it and so do
    tables published by earlier versions.
    """
    frame = pd.DataFrame([_hit_row("cand_clean", CLEAN_GUIDE, "human", "ENST00000000009")])

    validated = AggregatedOffTargetSchema.validate(frame, lazy=True)
    assert HIT_CLASS_COLUMN not in validated.columns


@pytest.mark.unit
def test_the_schema_hit_class_values_match_the_enum():
    """The schema restates the class vocabulary as literals; a new member must not slip past it."""
    assert set(HIT_CLASS_VALUES) == {member.value for member in HitClass}


# --- P8: the join key #103 calls its top risk --------------------------------------------------


@pytest.mark.unit
def test_screen_query_id_is_the_join_key_from_candidate_rows_to_hit_rows(tmp_path):
    """Two candidates sharing one guide are screened under one query id, and both carry it.

    Joining hit rows on ``id`` attributes a guide's whole hit set to one of its candidate rows
    (median 6, max 34 on the frozen baseline). Joining on the guide sequence works only while the
    spellings are byte-identical, which a U-spelled guide and its T-spelled twin are not.
    """
    workflow = _workflow(tmp_path, "out_join_key")
    t_spelled = _candidate("cand_t", CLEAN_GUIDE)
    u_spelled = _candidate("cand_u", CLEAN_GUIDE.replace("T", "U"))
    distinct = _candidate("cand_repeat", REPEAT_GUIDE)

    asyncio.run(workflow._prepare_offtarget_input([t_spelled, u_spelled, distinct]))

    assert t_spelled.screen_query_id == u_spelled.screen_query_id == "cand_t"
    assert distinct.screen_query_id == "cand_repeat"
    assert t_spelled.guide_sequence != u_spelled.guide_sequence, "the spellings genuinely differ"
    for candidate in (t_spelled, u_spelled, distinct):
        assert build_candidate_row(candidate)["screen_query_id"] == candidate.screen_query_id


@pytest.mark.unit
def test_screen_query_id_matches_the_qname_the_hit_rows_carry(tmp_path):
    """The join is an id match against the published table, not a sequence comparison."""
    _candidates, tsv_path, _outcome = _run(tmp_path, "out_join_qname")
    rows = _read_tsv(tsv_path)
    workflow_candidates = {candidate.id: candidate for candidate in _candidates}

    for candidate in workflow_candidates.values():
        assert candidate.screen_query_id == candidate.id, "no dedup happened in this fixture"
        own = [row for row in rows if row["qname"] == candidate.screen_query_id]
        assert own, f"no hit rows joined to {candidate.screen_query_id}"


@pytest.mark.unit
def test_the_join_key_defaults_to_none_and_serializes_as_none():
    """The field default and its passthrough into the candidate row -- nothing more.

    Named for what it covers: it asserts on a freshly constructed object, so it cannot see whether
    the workflow leaves the key unset for a candidate the aligner never saw. That behaviour is
    pinned by ``test_candidate_never_submitted_is_not_scored_as_clean`` in
    ``test_post_screen_scoring_integrity.py``, which builds the real state.
    """
    assert _candidate("cand_unsubmitted", CLEAN_GUIDE).screen_query_id is None
    assert build_candidate_row(_candidate("cand_unsubmitted", CLEAN_GUIDE))["screen_query_id"] is None


@pytest.mark.unit
def test_the_returned_results_map_is_looked_up_by_the_screening_id(tmp_path):
    """``offtarget_summary.results`` must not report a fabricated zero for a deduplicated candidate.

    This map is serialized into ``logs/workflow_summary.json``. Keyed and looked up by ``id``, the
    non-representative of a deduplicated pair got ``off_target_count: 0`` while
    ``candidates_all.csv`` carried its real count for the same id -- 32,463 of 34,863 ids on the
    frozen baseline. Every other test in this file bypasses ``_prepare_offtarget_input``, so the
    dedup map is empty and this whole class of defect is invisible to them.

    The count is also pinned against the candidate rather than the ingest tally: the ingest counts
    every row it read, on-target rows included, so the two hit rows here must publish
    ``off_target_count: 1``, not 2.

    The map carries counts only. It used to carry the alignment rows too, re-serialised once per
    candidate sharing a guide, which is what made this file 2.9 GiB on a real run; the rows now stay
    in the tables ``detail_files`` names.
    """
    workflow = _workflow(tmp_path, "out_dedup_results_map")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    _write_results_dir(
        results_dir,
        [
            _hit_row("cand_rep", CLEAN_GUIDE, "human", "ENST00000000001.2"),  # own gene: not a liability
            _hit_row("cand_rep", CLEAN_GUIDE, "human", "ENST00000000009"),  # the one liability
        ],
    )

    representative = _candidate("cand_rep", CLEAN_GUIDE)
    duplicate = _candidate("cand_dup", CLEAN_GUIDE)
    asyncio.run(workflow._prepare_offtarget_input([representative, duplicate]))
    assert duplicate.screen_query_id == "cand_rep", "the fixture must genuinely deduplicate"

    outcome = asyncio.run(
        workflow._process_nextflow_results([representative, duplicate], results_dir, {"status": "completed"})
    )

    published = outcome["results"]
    assert set(published) == {"cand_rep", "cand_dup"}, "the map is keyed by candidate id"
    for candidate in (representative, duplicate):
        assert candidate.off_target_count == 1, "one liability, fanned out to both candidates"
        entry = published[candidate.id]
        assert entry["off_target_count"] == candidate.off_target_count
        assert "hits" not in entry, "the summary publishes counts; the rows live in detail_files"

    # The rows are not lost, they are pointed at: the table the parser actually read.
    detail = outcome["detail_files"]
    assert detail["transcriptome"], "the summary must name the table holding the alignments"
    assert all(Path(path).exists() for path in detail["transcriptome"]), (
        "a pointer to a file that does not exist is worse than no pointer"
    )


@pytest.mark.unit
def test_a_row_annotated_in_part_is_repaired_rather_than_republished_with_blank_cells(tmp_path):
    """``is_annotated`` is authoritative over every classification column, not over ``hit_class``.

    A stale or partially annotated table -- a valid ``hit_class`` but no symbols -- counted as
    annotated, so the orphan pass skipped it and the writer republished ``row.get(column, "")`` for
    the five columns it had not forced re-annotation of. The resulting file failed
    ``AggregatedOffTargetSchema`` on ``hit_symbol``, ``hit_symbol_missing`` and
    ``species_index_missing`` -- the schema the same run publishes this table against.

    The row's qname belongs to no candidate in this run, which is what makes the orphan pass the
    only thing that can repair it: a row whose candidate *is* present is re-annotated by the main
    integration loop regardless of what ``is_annotated`` says.
    """
    workflow = _workflow(tmp_path, "out_partial_annotation")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True, exist_ok=True)
    tsv_path = aggregated / "combined_offtargets.tsv"
    header = [*TSV_COLUMNS, *CLASSIFICATION_COLUMNS]
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        # Written by a version that knew only the first three columns: a usable class, blank rest.
        writer.writerow(
            {
                **_hit_row("cand_gone", ORPHAN_GUIDE, "human", "ENST00000000009"),
                HIT_CLASS_COLUMN: HitClass.OFF_TARGET.value,
                MATCHED_SYMBOL_COLUMN: UNKNOWN_SYMBOL,
                SYMBOL_LOOKUP_MISSING_COLUMN: "False",
            }
        )
    _write_summary(aggregated, 1)

    outcome = asyncio.run(
        workflow._process_nextflow_results(
            [_candidate("cand_clean", CLEAN_GUIDE)], results_dir, {"status": "completed"}
        )
    )

    assert outcome["status"] == "completed"
    assert outcome["filtering_stats"]["hits_classified_without_candidate"] == 1, "the orphan pass repaired it"
    republished = _read_tsv(tsv_path)
    assert republished[0][HIT_SYMBOL_COLUMN] == "OTHER", "the row was re-annotated, not republished"
    assert republished[0][HIT_SYMBOL_MISSING_COLUMN] == "False"
    assert republished[0][SPECIES_INDEX_MISSING_COLUMN] == "False"
    assert all(republished[0][column] for column in CLASSIFICATION_COLUMNS), "never an empty cell"
    # The published contract's claim: this file validates against the schema the run declares.
    AggregatedOffTargetSchema.validate(pd.read_csv(tsv_path, sep="\t", dtype=str), lazy=True)


# --- the blank-class case, which used to crash after the table was overwritten ------------------


@pytest.mark.unit
def test_a_blank_hit_class_is_reclassified_rather_than_republished(tmp_path):
    """``is_annotated`` tests the value, not the key.

    A table read back with an empty ``hit_class`` cell used to count as annotated: the blank was
    republished and the next read raised ``ValueError`` from ``HitClass("")`` -- after the file had
    already been overwritten, downgrading a completed screen to ``nextflow_failed``.
    """
    workflow = _workflow(tmp_path, "out_blank_class")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True, exist_ok=True)
    tsv_path = aggregated / "combined_offtargets.tsv"
    header = [*TSV_COLUMNS, *CLASSIFICATION_COLUMNS]
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        # A stale row whose class cell is empty: present as a key, useless as a verdict.
        writer.writerow({**_hit_row("cand_gone", ORPHAN_GUIDE, "human", "ENST00000000009")})
    _write_summary(aggregated, 1)

    outcome = asyncio.run(
        workflow._process_nextflow_results(
            [_candidate("cand_clean", CLEAN_GUIDE)], results_dir, {"status": "completed"}
        )
    )

    assert outcome["status"] == "completed"
    republished = _read_tsv(tsv_path)
    assert republished[0][HIT_CLASS_COLUMN] == HitClass.OFF_TARGET.value
    assert outcome["filtering_stats"]["hits_classified_without_candidate"] == 1
    # Re-reading the repaired table must not raise, and must not change it again.
    before = tsv_path.read_text()
    asyncio.run(
        workflow._process_nextflow_results(
            [_candidate("cand_clean", CLEAN_GUIDE)], results_dir, {"status": "completed"}
        )
    )
    assert tsv_path.read_text() == before


@pytest.mark.unit
def test_an_unclassified_published_row_is_reported_rather_than_silently_skipped(tmp_path):
    """The replacement for the removed "leave the file unchanged" branch.

    That branch could not be reached from its only caller, because the orphan pass annotates every
    row of every table first. Reporting the shortfall is reachable, and says which invariant broke
    instead of leaving a run looking clean.
    """
    tsv_path = tmp_path / "combined_offtargets.tsv"
    row = _hit_row("cand_clean", CLEAN_GUIDE, "human", "ENST00000000009")
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    parsed = {
        "results": {},
        "transcriptome_hit_tables": [{"path": tsv_path, "fieldnames": TSV_COLUMNS, "rows": [row]}],
    }

    warnings = SiRNAWorkflow._persist_hit_classifications(parsed)

    assert len(warnings) == 1
    assert "1 published hit row(s) carry no usable classification" in warnings[0]


@pytest.mark.unit
def test_a_blank_ortholog_evidence_cell_is_reannotated_rather_than_republished(tmp_path):
    """A table written before ``ortholog_evidence`` existed must be repaired, not republished blank.

    ``ortholog_evidence`` joined ``CLASSIFICATION_COLUMNS`` as its seventh member. While
    ``is_annotated`` still checked only the first six, a stale row with every other cell valid counted
    as annotated, so the blank was republished -- and ``AggregatedOffTargetSchema``'s ``isin`` on
    ``ORTHOLOG_EVIDENCE_VALUES`` then rejects the very table the same run publishes.
    """
    workflow = _workflow(tmp_path, "out_blank_evidence")
    results_dir = workflow.config.output_dir / "off_target" / "results"
    aggregated = results_dir / "aggregated"
    aggregated.mkdir(parents=True, exist_ok=True)
    tsv_path = aggregated / "combined_offtargets.tsv"
    header = [*TSV_COLUMNS, *CLASSIFICATION_COLUMNS]
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        # Fully annotated by the pre-ortholog_evidence contract: only that one cell is blank.
        writer.writerow(
            {
                **_hit_row("cand_gone", ORPHAN_GUIDE, "human", "ENST00000000009"),
                HIT_CLASS_COLUMN: HitClass.OFF_TARGET.value,
                MATCHED_SYMBOL_COLUMN: UNKNOWN_SYMBOL,
                SYMBOL_LOOKUP_MISSING_COLUMN: "False",
                HIT_SYMBOL_COLUMN: "OTHER",
                HIT_SYMBOL_MISSING_COLUMN: "False",
                SPECIES_INDEX_MISSING_COLUMN: "False",
                ORTHOLOG_EVIDENCE_COLUMN: "",
            }
        )
    _write_summary(aggregated, 1)

    outcome = asyncio.run(
        workflow._process_nextflow_results(
            [_candidate("cand_clean", CLEAN_GUIDE)], results_dir, {"status": "completed"}
        )
    )

    assert outcome["status"] == "completed"
    republished = _read_tsv(tsv_path)
    assert republished[0][ORTHOLOG_EVIDENCE_COLUMN] == ORTHOLOG_EVIDENCE_NOT_APPLICABLE, "never an empty cell"
    assert all(row[ORTHOLOG_EVIDENCE_COLUMN] in ORTHOLOG_EVIDENCE_CELL_VALUES for row in republished)
