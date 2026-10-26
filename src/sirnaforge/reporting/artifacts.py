"""Versioned, queryable report artifacts staged from projected source tables."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import sqlite3
import tempfile
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS
from sirnaforge.reporting.payload import (
    _HIT_REPORT_COLUMNS,
    _MIRNA_REPORT_COLUMNS,
    REPORT_READ_CHUNK_ROWS,
    ReportInputError,
    _candidate_projection,
    _gate_panel,
    _read_manifest,
)

ReportMode = Literal["auto", "compact_single_file", "paged_single_file", "data_backed"]
ARTIFACT_SCHEMA_VERSION = 1
_MODE_ORDER = ("compact_single_file", "paged_single_file", "data_backed")


@dataclass(frozen=True)
class ReportModeThresholds:
    """Provisional report-size budgets; tune from benchmark artifacts, not biology."""

    compact_max_input_bytes: int = 64 * 1024 * 1024
    compact_max_candidate_rows: int = 50_000
    compact_max_guides: int = 5_000
    paged_max_input_bytes: int = 512 * 1024 * 1024
    paged_max_candidate_rows: int = 500_000
    paged_max_guides: int = 50_000


@dataclass(frozen=True)
class ReportModeFacts:
    """Stable input measurements used for deterministic mode selection."""

    input_bytes: int
    candidate_rows: int
    distinct_guides: int


@dataclass
class ReportArtifactStore:
    """A Parquet-backed report dataset with a versioned manifest and bounded batch reads."""

    root: Path
    manifest: dict[str, Any]

    @classmethod
    def from_run(
        cls,
        run_dir: Path | str,
        artifact_dir: Path | str,
        *,
        mode: ReportMode = "auto",
        thresholds: ReportModeThresholds | None = None,
        policy: Any = None,
    ) -> ReportArtifactStore:
        """Stage report-used candidate/hit fields to Parquet without materializing full source tables."""
        run_root = Path(run_dir)
        destination = Path(artifact_dir)
        if mode not in {"auto", *_MODE_ORDER}:
            raise ValueError(f"unsupported report mode: {mode}")
        resolved_run = run_root.resolve()
        resolved_destination = destination.resolve()
        if resolved_destination == resolved_run or resolved_destination in resolved_run.parents:
            raise ReportInputError(
                f"artifact directory {destination} cannot replace the source run directory {run_root}"
            )
        candidates_path = next(iter(sorted(run_root.glob("**/candidates_all.csv"))), None)
        if candidates_path is None:
            raise ReportInputError(f"no candidates_all.csv under {run_root}")
        caveats: list[str] = []
        manifest_path, run_manifest = _read_manifest(run_root, caveats)
        panel = _gate_panel(policy, run_manifest, caveats)
        candidate_header = _header(candidates_path, ",")
        candidate_columns = _candidate_projection(candidate_header, [f.descriptor for f in panel.filters])
        agg = run_root / "off_target" / "results" / "aggregated"
        hit_path = agg / "combined_offtargets.tsv"
        mirna_path = agg / "combined_mirna_hits.tsv"
        report_input_paths = [candidates_path, *(path for path in (hit_path, mirna_path) if path.exists())]
        if any(
            resolved_destination == source.resolve() or resolved_destination in source.resolve().parents
            for source in report_input_paths
        ):
            raise ReportInputError(f"artifact directory {destination} would replace a report input under {run_root}")

        destination.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix=f".{destination.name}.stage-", dir=destination.parent))
        try:
            tables: dict[str, dict[str, Any]] = {}
            tables["candidates"] = _write_projected_parquet(
                candidates_path, stage / "candidates.parquet", ",", set(candidate_columns)
            )
            if hit_path.exists():
                header = _header(hit_path, "\t")
                missing = [column for column in CLASSIFICATION_COLUMNS if column not in header]
                if missing:
                    raise ReportInputError(
                        f"{hit_path.name} is missing {missing}. The report will not classify hits itself."
                    )
                tables["hits"] = _write_projected_parquet(hit_path, stage / "hits.parquet", "\t", _HIT_REPORT_COLUMNS)
            if mirna_path.exists():
                tables["mirna"] = _write_projected_parquet(
                    mirna_path, stage / "mirna.parquet", "\t", _MIRNA_REPORT_COLUMNS
                )

            candidate_meta = tables["candidates"]
            input_bytes = sum(int(table["source_bytes"]) for table in tables.values())
            facts = ReportModeFacts(
                input_bytes=input_bytes,
                candidate_rows=int(candidate_meta["rows"]),
                distinct_guides=int(candidate_meta["distinct_guides"]),
            )
            selected_mode = resolve_report_mode(mode, facts, thresholds or ReportModeThresholds())
            source_tables = {"candidates": candidates_path, "hits": hit_path, "mirna": mirna_path}
            source_provenance = {
                name: _source_provenance(path) for name, path in source_tables.items() if path.exists()
            }
            if manifest_path is not None:
                source_provenance["manifest"] = _source_provenance(manifest_path)
            report_manifest = {
                "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
                "report_mode": selected_mode,
                "source_run_label": run_root.name,
                "source_manifest": str(manifest_path.relative_to(run_root)) if manifest_path else None,
                "mode_facts": {
                    "input_bytes": facts.input_bytes,
                    "candidate_rows": facts.candidate_rows,
                    "distinct_guides": facts.distinct_guides,
                },
                "mode_thresholds": (thresholds or ReportModeThresholds()).__dict__,
                "tables": tables,
                "inputs": source_provenance,
                "exports": {},
            }
            (stage / "report_manifest.json").write_text(json.dumps(report_manifest, indent=2) + "\n")
            _publish_directory(stage, destination)
            return cls(destination, report_manifest)
        except BaseException:
            shutil.rmtree(stage, ignore_errors=True)
            raise

    def iter_batches(
        self,
        table: str,
        *,
        columns: Sequence[str] | None = None,
        batch_size: int = REPORT_READ_CHUNK_ROWS,
    ) -> Iterator[pa.RecordBatch]:
        """Stream a table's selected columns in bounded Arrow batches."""
        metadata = self.manifest["tables"].get(table)
        if metadata is None:
            raise KeyError(f"report artifact has no {table!r} table")
        parquet = pq.ParquetFile(self.root / metadata["file"])
        yield from parquet.iter_batches(batch_size=batch_size, columns=list(columns) if columns else None)

    def export_table(self, table: str, output_path: Path | str | None = None) -> Path:
        """Write a lossless CSV/TSV view and record its source Parquet checksum in the manifest."""
        metadata = self.manifest["tables"].get(table)
        if metadata is None:
            raise KeyError(f"report artifact has no {table!r} table")
        delimiter = "\t" if table in {"hits", "mirna"} else ","
        suffix = ".tsv" if delimiter == "\t" else ".csv"
        output = Path(output_path) if output_path else self.root / "exports" / f"{table}{suffix}"
        output.parent.mkdir(parents=True, exist_ok=True)
        schema = pq.ParquetFile(self.root / metadata["file"]).schema_arrow
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
            writer.writerow(schema.names)
            for batch in self.iter_batches(table):
                writer.writerows([row.get(column) for column in schema.names] for row in batch.to_pylist())
        export = {
            "file": str(output.relative_to(self.root)) if output.is_relative_to(self.root) else str(output),
            "format": "tsv" if delimiter == "\t" else "csv",
            "table": table,
            "schema_version": metadata["schema_version"],
            "artifact_table_file": metadata["file"],
            "artifact_table_sha256": metadata["sha256"],
            "sha256": _sha256(output),
        }
        exports = self.manifest.setdefault("exports", {})
        exports[table] = export
        (self.root / "report_manifest.json").write_text(json.dumps(self.manifest, indent=2) + "\n")
        return output


def resolve_report_mode(mode: ReportMode, facts: ReportModeFacts, thresholds: ReportModeThresholds) -> str:
    """Resolve auto mode from deterministic input budgets; explicit requests are honored unchanged."""
    if mode != "auto":
        return mode
    compact = (
        facts.input_bytes <= thresholds.compact_max_input_bytes
        and facts.candidate_rows <= thresholds.compact_max_candidate_rows
        and facts.distinct_guides <= thresholds.compact_max_guides
    )
    if compact:
        return _MODE_ORDER[0]
    paged = (
        facts.input_bytes <= thresholds.paged_max_input_bytes
        and facts.candidate_rows <= thresholds.paged_max_candidate_rows
        and facts.distinct_guides <= thresholds.paged_max_guides
    )
    return _MODE_ORDER[1] if paged else _MODE_ORDER[2]


def _header(path: Path, delimiter: str) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle, delimiter=delimiter))


def _write_projected_parquet(path: Path, output: Path, delimiter: str, selected: set[str]) -> dict[str, Any]:
    """Arrow streaming CSV reader -> Parquet writer, with text cells preserving source spellings."""
    header = _header(path, delimiter)
    columns = [column for column in header if column in selected]
    if not columns:
        raise ReportInputError(f"{path.name} has no columns needed by the report")
    read_options = pacsv.ReadOptions(block_size=1 << 20, use_threads=False)
    parse_options = pacsv.ParseOptions(delimiter=delimiter)
    convert_options = pacsv.ConvertOptions(
        include_columns=columns,
        column_types={column: pa.string() for column in columns},
    )
    reader = pacsv.open_csv(
        path, read_options=read_options, parse_options=parse_options, convert_options=convert_options
    )
    schema = reader.schema
    rows = 0
    guide_db_path = output.with_suffix(".guides.sqlite")
    guide_db = sqlite3.connect(guide_db_path) if "guide_sequence" in columns else None
    if guide_db is not None:
        guide_db.execute("CREATE TABLE guides (sequence TEXT PRIMARY KEY) WITHOUT ROWID")
    writer = pq.ParquetWriter(output, schema, compression="zstd")
    try:
        for batch in reader:
            rows += batch.num_rows
            if guide_db is not None:
                values = batch.column(batch.schema.get_field_index("guide_sequence")).to_pylist()
                guide_db.executemany(
                    "INSERT OR IGNORE INTO guides VALUES (?)",
                    ((str(value).strip().upper().replace("T", "U"),) for value in values if value),
                )
            writer.write_batch(batch, row_group_size=REPORT_READ_CHUNK_ROWS)
    finally:
        writer.close()
    distinct_guides = 0
    if guide_db is not None:
        guide_db.commit()
        distinct_guides = int(guide_db.execute("SELECT COUNT(*) FROM guides").fetchone()[0])
        guide_db.close()
        guide_db_path.unlink(missing_ok=True)
    return {
        "file": output.name,
        "format": "parquet",
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "rows": rows,
        "columns": columns,
        "source_file": path.name,
        "source_bytes": path.stat().st_size,
        "source_sha256": _sha256(path),
        "file_bytes": output.stat().st_size,
        "sha256": _sha256(output),
        "distinct_guides": distinct_guides,
    }


def _source_provenance(path: Path) -> dict[str, Any]:
    return {"file": path.name, "bytes": path.stat().st_size, "sha256": _sha256(path)}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _publish_directory(stage: Path, destination: Path) -> None:
    """Publish a complete artifact directory atomically, replacing a prior generated directory."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    backup = destination.with_name(f".{destination.name}.previous")
    if backup.exists():
        raise ReportInputError(f"refusing to overwrite {destination}: prior artifact backup exists at {backup}")
    if destination.exists():
        if destination.is_symlink() or not _is_sirnaforge_artifact_directory(destination):
            raise ReportInputError(
                f"refusing to replace existing directory {destination}: no recognized report manifest found"
            )
        destination.replace(backup)
    try:
        stage.replace(destination)
    except BaseException:
        if backup.exists() and not destination.exists():
            backup.replace(destination)
        raise
    shutil.rmtree(backup, ignore_errors=True)


def _is_sirnaforge_artifact_directory(path: Path) -> bool:
    """Return whether a directory carries a manifest identifying owned report artifacts."""
    manifest_path = path / "report_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return (
        manifest.get("artifact_schema_version") == ARTIFACT_SCHEMA_VERSION
        and manifest.get("report_mode") in _MODE_ORDER
        and isinstance(manifest.get("tables"), dict)
    )
