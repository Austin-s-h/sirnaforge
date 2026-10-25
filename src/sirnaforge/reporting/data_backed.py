"""Static, offline data-backed reader built from normalized report artifacts."""

from __future__ import annotations

import importlib
import json
import math
import shutil
import sqlite3
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from sirnaforge.reporting.artifacts import (
    ARTIFACT_SCHEMA_VERSION,
    ReportArtifactStore,
    ReportMode,
    ReportModeThresholds,
    _sha256,
)
from sirnaforge.reporting.payload import (
    _VERDICT_CODE,
    REPORT_READ_CHUNK_ROWS,
    STATUSES,
    VERDICT_WARN,
    _evaluate,
    _gate_panel,
    _normalise_guide,
    _read_manifest,
    _run_verdict,
    build_payload,
    observed_column,
)
from sirnaforge.reporting.render import write_report

_SHARD_GUIDES = 128


@dataclass(frozen=True)
class ReportWriteResult:
    """Files and selected mode produced by :func:`write_report_from_run`."""

    output: Path
    mode: str
    manifest: dict[str, Any]
    manifest_path: Path
    diagnostics: dict[str, Any]
    diagnostics_path: Path
    artifact_dir: Path | None


def write_report_from_run(
    run_dir: Path | str,
    output: Path | str,
    *,
    mode: ReportMode = "auto",
    artifact_dir: Path | str | None = None,
    thresholds: ReportModeThresholds | None = None,
    policy: Any = None,
) -> ReportWriteResult:
    """Write a report from a run, using local detail shards for the data-backed mode.

    Compact reports continue through the established payload and renderer. The `data_backed` path
    never calls `build_payload`: it stages projected tables, builds one row per guide, and places the
    raw guide detail into on-demand local JavaScript shards beside the HTML.
    """
    output_path = Path(output)
    run_root = Path(run_dir)
    build_started = time.perf_counter()
    phases: dict[str, float] = {}
    if mode == "compact_single_file":
        phase_started = time.perf_counter()
        payload = build_payload(run_root, policy=policy, max_guides=None, diagnostics=True)
        result_manifest = _compact_manifest(run_root, mode, thresholds)
        result_manifest["summary"] = {
            "gene_query": payload.run.get("gene_query", run_root.name),
            "guides": payload.run.get("guides", 0),
            "candidate_rows": payload.run.get("candidate_rows", 0),
            "status_counts": payload.run.get("status_counts", {}),
        }
        phases["payload_build"] = time.perf_counter() - phase_started
        phase_started = time.perf_counter()
        write_report(payload, output_path)
        phases["html_render_and_write"] = time.perf_counter() - phase_started
        return _write_result_metadata(
            output_path, mode, result_manifest, payload.diagnostics or {}, phases, build_started
        )

    phase_started = time.perf_counter()
    artifacts = Path(artifact_dir) if artifact_dir else output_path.with_suffix(".data")
    store = ReportArtifactStore.from_run(run_root, artifacts, mode=mode, thresholds=thresholds, policy=policy)
    phases["projected_input_staging"] = time.perf_counter() - phase_started
    selected = str(store.manifest["report_mode"])
    if selected == "compact_single_file":
        phase_started = time.perf_counter()
        payload = build_payload(run_root, policy=policy, max_guides=None, diagnostics=True)
        store.manifest["summary"] = {
            "gene_query": payload.run.get("gene_query", run_root.name),
            "guides": payload.run.get("guides", 0),
            "candidate_rows": payload.run.get("candidate_rows", 0),
            "status_counts": payload.run.get("status_counts", {}),
        }
        phases["payload_build"] = time.perf_counter() - phase_started
        phase_started = time.perf_counter()
        write_report(payload, output_path)
        phases["html_render_and_write"] = time.perf_counter() - phase_started
        return _write_result_metadata(
            output_path, selected, store.manifest, payload.diagnostics or {}, phases, build_started, store.root
        )
    if selected not in {"data_backed", "paged_single_file"}:
        raise ValueError(f"unsupported report mode: {selected}")

    phase_started = time.perf_counter()
    _derive_normalized_tables(store, run_root, policy=policy)
    phases["normalized_table_build"] = time.perf_counter() - phase_started
    shard_dir = store.root / ".render-shards"
    phase_started = time.perf_counter()
    _write_detail_shards(store, shard_dir)
    phases["detail_shard_write"] = time.perf_counter() - phase_started
    if selected == "paged_single_file":
        detail_directory = None
        detail_loading = "inline_paged_data"
    else:
        detail_directory = output_path.with_name(f"{output_path.stem}.details")
        if detail_directory.exists():
            shutil.rmtree(detail_directory)
        shutil.copytree(shard_dir, detail_directory, dirs_exist_ok=True)
        detail_directory_name = detail_directory.name
        detail_loading = "classic_script_shards_on_demand"
    store.manifest["reader"] = {
        "entrypoint": output_path.name,
        "detail_directory": detail_directory_name if selected == "data_backed" else None,
        "detail_loading": detail_loading,
        "offline_file_origin_supported": True,
    }
    store.manifest["diagnostics"] = {
        "schema_version": 1,
        "phases_seconds": phases,
        "peak_rss_bytes": _peak_rss_bytes(),
        "peak_memory_scope": "process_lifetime",
        "report_mode": "data_backed",
    }
    phase_started = time.perf_counter()
    if selected == "data_backed":
        _write_data_backed_html(store, output_path, detail_directory_name)
    else:
        store.manifest["inline_detail_shards"] = [
            {"index": item["index"], "bytes": item["bytes"], "sha256": item["sha256"]}
            for item in store.manifest.pop("detail_shards", [])
        ]
        _write_paged_html(store, output_path, shard_dir)
    phases["html_render_and_write"] = time.perf_counter() - phase_started
    store.manifest["diagnostics"]["phases_seconds"] = phases
    store.manifest["diagnostics"]["output_html_bytes"] = output_path.stat().st_size
    store.manifest["diagnostics"]["total_seconds"] = time.perf_counter() - build_started
    (store.root / "report_manifest.json").write_text(json.dumps(store.manifest, indent=2) + "\n")
    if shard_dir.exists():
        shutil.rmtree(shard_dir)
    return _write_result_metadata(
        output_path,
        selected,
        store.manifest,
        store.manifest["diagnostics"],
        phases,
        build_started,
        store.root,
    )


def _derive_normalized_tables(store: ReportArtifactStore, run_dir: Path, *, policy: Any = None) -> None:
    """Derive all-guide summaries and normalized isoform/hit summaries with bounded row reads."""
    caveats: list[str] = []
    _manifest_path, manifest = _read_manifest(run_dir, caveats)
    panel = _gate_panel(policy, manifest, caveats)
    identity = manifest.get("target_identity") if isinstance(manifest.get("target_identity"), dict) else {}
    identity_rows_raw = identity.get("transcripts") if isinstance(identity, dict) else []
    identity_rows = identity_rows_raw if isinstance(identity_rows_raw, list) else []
    identity_by_transcript = {
        str(entry["transcript_id"]): entry
        for entry in identity_rows
        if isinstance(entry, dict) and isinstance(entry.get("transcript_id"), str)
    }
    canonical_ids = {
        transcript
        for transcript, entry in identity_by_transcript.items()
        if entry.get("is_canonical") is True
        and entry.get("canonical_source")
        and entry.get("canonical_status") not in {"not_established", "unresolved_input", "unknown"}
    }
    store.manifest["target_identity"] = identity
    populated: set[str] = set()
    temp_db = store.root / ".derive.sqlite"
    db = sqlite3.connect(temp_db)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE candidate_rows (guide TEXT, payload TEXT)")
    db.execute("CREATE INDEX candidate_guide ON candidate_rows(guide)")
    for batch in store.iter_batches("candidates", batch_size=REPORT_READ_CHUNK_ROWS):
        records = batch.to_pylist()
        for record in records:
            for name, value in record.items():
                if value is not None and str(value).strip():
                    populated.add(name)
        db.executemany(
            "INSERT INTO candidate_rows VALUES (?, ?)",
            ((_normalise_guide(row.get("guide_sequence")), json.dumps(row, separators=(",", ":"))) for row in records),
        )
    db.commit()

    descriptors = [gate.descriptor for gate in panel.filters]
    gate_columns = [observed_column(descriptor, populated) for descriptor in descriptors]
    gate_effects: Counter[tuple[int, str]] = Counter()
    status_counts: Counter[str] = Counter()
    guide_records: list[dict[str, Any]] = []
    isoform_sink = _ParquetRecordSink(
        store.root / "isoforms.parquet",
        pa.schema(
            [
                ("guide", pa.string()),
                ("transcript_id", pa.string()),
                ("gene_id", pa.string()),
                ("gene_symbol", pa.string()),
                ("transcript_biotype", pa.string()),
                ("canonical", pa.string()),
                ("canonical_status", pa.string()),
                ("canonical_source", pa.string()),
                ("canonical_release", pa.string()),
                ("candidate_rows", pa.int64()),
            ]
        ),
    )
    candidate_rows_for_detail = db.execute("SELECT guide, payload FROM candidate_rows ORDER BY guide")
    current_guide: str | None = None
    group: list[dict[str, Any]] = []

    def flush_guide(guide: str, rows: list[dict[str, Any]]) -> None:
        frame = pd.DataFrame(rows)
        if frame.empty:
            return
        rank_columns = [column for column in ("composite_score", "design_score") if column in frame]
        ranked = frame.copy()
        for column in rank_columns:
            ranked[column] = pd.to_numeric(ranked[column], errors="coerce")
        best = (
            ranked.sort_values(rank_columns, ascending=False, na_position="last").iloc[0]
            if rank_columns
            else frame.iloc[0]
        )
        gates = [list(_evaluate(desc, best, col)) for desc, col in zip(descriptors, gate_columns, strict=True)]
        fail = _VERDICT_CODE["fail"]
        unknown = _VERDICT_CODE["unknown"]
        n_failed = sum(gate[1] == fail for gate in gates)
        n_unknown = sum(gate[1] == unknown for gate in gates)
        n_warned = sum(gate[1] == VERDICT_WARN for gate in gates)
        run_verdict = _run_verdict(frame)
        undeclared = bool(run_verdict not in (None, "PASS") and not n_failed and not n_unknown)
        status = "fail" if n_failed else "unknown" if n_unknown or undeclared else "warn" if n_warned else "pass"
        status_counts[status] += 1
        failed_gate_indexes = [i for i, gate in enumerate(gates) if gate[1] == fail]
        if len(failed_gate_indexes) == 1:
            gate_effects[(failed_gate_indexes[0], "sole_rejects")] += 1
        for i, gate in enumerate(gates):
            if gate[1] == fail:
                gate_effects[(i, "rejects")] += 1
            if gate[1] == VERDICT_WARN:
                gate_effects[(i, "warns")] += 1
            if gate[1] == unknown:
                gate_effects[(i, "unknowns")] += 1
        iso_seen: set[str] = set()
        for row in rows:
            tx = str(row.get("transcript_id") or "")
            if not tx or tx in iso_seen:
                continue
            iso_seen.add(tx)
            identity_row = identity_by_transcript.get(tx, {})
            canonical_status = identity_row.get("canonical_status") or "unknown"
            canonical_value = (
                "canonical"
                if tx in canonical_ids
                else "other_isoform"
                if identity_row.get("canonical_source")
                and canonical_status not in {"not_established", "unresolved_input", "unknown"}
                else "unknown"
            )
            isoform_sink.append(
                {
                    "guide": guide,
                    "transcript_id": tx,
                    "gene_id": identity_row.get("gene_id") or row.get("gene_id"),
                    "gene_symbol": identity_row.get("gene_symbol") or row.get("gene_symbol") or row.get("gene_name"),
                    "transcript_biotype": identity_row.get("transcript_biotype"),
                    "canonical": canonical_value,
                    "canonical_status": canonical_status,
                    "canonical_source": identity_row.get("canonical_source"),
                    "canonical_release": identity_row.get("canonical_release"),
                    "candidate_rows": sum(1 for item in rows if str(item.get("transcript_id") or "") == tx),
                }
            )
        guide_records.append(
            {
                "guide": guide,
                "candidate_rows": len(rows),
                "design_score": _float_or_none(best.get("design_score")),
                "composite_score": _float_or_none(best.get("composite_score")),
                "run_verdict": run_verdict,
                "canonical_hit": (
                    "hit"
                    if set(iso_seen) & canonical_ids
                    else "not_hit"
                    if canonical_ids and canonical_ids <= iso_seen
                    else "not_established"
                ),
                "status": status,
                "gates": json.dumps(gates, separators=(",", ":")),
                "n_failed_gates": n_failed,
                "n_unknown_gates": n_unknown,
                "n_warned_gates": n_warned,
                "detail_shard": "",
            }
        )

    for raw_guide, raw_payload in candidate_rows_for_detail:
        guide = str(raw_guide)
        if current_guide is not None and guide != current_guide:
            flush_guide(current_guide, group)
            group = []
        current_guide = guide
        group.append(json.loads(raw_payload))
    if current_guide is not None:
        flush_guide(current_guide, group)
    db.close()
    temp_db.unlink(missing_ok=True)

    _write_records(store, "guides", guide_records)
    isoform_sink.close()
    _derive_hit_summary(store)
    store.manifest["summary"] = {
        "guides": len(guide_records),
        "status_counts": {status: status_counts.get(status, 0) for status in STATUSES},
        "filters": [
            {
                "filter_id": desc.filter_id,
                **{key: gate_effects[(index, key)] for key in ("rejects", "sole_rejects", "warns", "unknowns")},
            }
            for index, desc in enumerate(descriptors)
        ],
        "candidate_rows": sum(row["candidate_rows"] for row in guide_records),
    }
    store.manifest["tables"].update(_table_metadata(store, "guides"))
    store.manifest["tables"].update(_table_metadata(store, "isoforms"))
    (store.root / "report_manifest.json").write_text(json.dumps(store.manifest, indent=2) + "\n")


def _derive_hit_summary(store: ReportArtifactStore) -> None:
    """Aggregate hit rows in SQLite so memory is bounded by the read batch, not the hit table."""
    db_path = store.root / ".hits.sqlite"
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE hits (guide TEXT,species TEXT,hit_class TEXT,nm_band TEXT,symbol TEXT)")
    if "hits" in store.manifest["tables"]:
        available = set(store.manifest["tables"]["hits"]["columns"])
        selected = [
            name for name in ("qseq", "species", "hit_class", "nm", "hit_symbol", "matched_symbol") if name in available
        ]
        for batch in store.iter_batches("hits", columns=selected):
            rows = []
            for row in batch.to_pylist():
                nm = _float_or_none(row.get("nm"))
                band = "unknown" if nm is None else "0" if nm == 0 else "1" if nm == 1 else "2" if nm == 2 else ">=3"
                rows.append(
                    (
                        _normalise_guide(row.get("qseq")),
                        str(row.get("species") or "query"),
                        str(row.get("hit_class") or "undetermined"),
                        band,
                        str(row.get("hit_symbol") or row.get("matched_symbol") or "unknown"),
                    )
                )
            db.executemany("INSERT INTO hits VALUES (?,?,?,?,?)", rows)
    db.commit()
    sink = _ParquetRecordSink(
        store.root / "hit_summary.parquet",
        pa.schema(
            [
                ("guide", pa.string()),
                ("species", pa.string()),
                ("hit_class", pa.string()),
                ("nm_band", pa.string()),
                ("symbol", pa.string()),
                ("count", pa.int64()),
            ]
        ),
    )
    cursor = db.execute(
        "SELECT guide,species,hit_class,nm_band,symbol,COUNT(*) FROM hits GROUP BY 1,2,3,4,5 ORDER BY 1,2,3,4,5"
    )
    hit_batch: list[tuple[Any, ...]]
    while hit_batch := cursor.fetchmany(10_000):
        for guide, species, hit_class, band, symbol, count in hit_batch:
            sink.append(
                {
                    "guide": guide,
                    "species": species,
                    "hit_class": hit_class,
                    "nm_band": band,
                    "symbol": symbol,
                    "count": count,
                }
            )
    sink.close()
    db.close()
    db_path.unlink(missing_ok=True)
    store.manifest["tables"].update(_table_metadata(store, "hit_summary"))


def _write_records(store: ReportArtifactStore, name: str, records: list[dict[str, Any]]) -> None:
    path = store.root / f"{name}.parquet"
    schema = pa.Table.from_pylist(records[:1]).schema if records else pa.schema([("_empty", pa.string())])
    sink = _ParquetRecordSink(path, schema)
    for record in records:
        sink.append(record)
    sink.close()


class _ParquetRecordSink:
    """Write Python records incrementally using a fixed Arrow schema."""

    def __init__(self, path: Path, schema: pa.Schema, *, batch_size: int = 10_000) -> None:
        self.path = path
        self.schema = schema
        self.batch_size = batch_size
        self.buffer: list[dict[str, Any]] = []
        self.writer = pq.ParquetWriter(path, schema, compression="zstd")

    def append(self, record: dict[str, Any]) -> None:
        self.buffer.append(record)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if self.buffer:
            self.writer.write_table(pa.Table.from_pylist(self.buffer, schema=self.schema))
            self.buffer.clear()

    def close(self) -> None:
        self.flush()
        self.writer.close()


def _table_metadata(store: ReportArtifactStore, name: str) -> dict[str, dict[str, Any]]:
    path = store.root / f"{name}.parquet"
    parquet = pq.ParquetFile(path)
    return {
        name: {
            "file": path.name,
            "format": "parquet",
            "schema_version": ARTIFACT_SCHEMA_VERSION,
            "rows": parquet.metadata.num_rows,
            "columns": parquet.schema_arrow.names,
            "file_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
    }


def _write_detail_shards(store: ReportArtifactStore, shard_dir: Path) -> None:
    """Write local JS data shards whose script files are only loaded when a guide is opened."""
    if shard_dir.exists():
        shutil.rmtree(shard_dir)
    shard_dir.mkdir(parents=True)
    guides = (
        pq.read_table(store.root / store.manifest["tables"]["guides"]["file"], columns=["guide"]).column(0).to_pylist()
    )
    guide_shard = {str(guide): i // _SHARD_GUIDES for i, guide in enumerate(guides)}
    db_path = store.root / ".detail-shards.sqlite"
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE detail (shard INTEGER, guide TEXT, kind TEXT, payload TEXT)")
    detail_sources = (
        ("candidates", "guide_sequence", "candidates"),
        ("hits", "qseq", "hits"),
        ("mirna", "qseq", "mirna"),
    )
    for table, guide_column, kind in detail_sources:
        if table not in store.manifest["tables"]:
            continue
        for batch in store.iter_batches(table, batch_size=REPORT_READ_CHUNK_ROWS):
            entries = []
            for row in batch.to_pylist():
                guide = _normalise_guide(row.get(guide_column))
                shard_no = guide_shard.get(guide)
                if shard_no is not None:
                    entries.append((shard_no, guide, kind, json.dumps(row, separators=(",", ":"))))
            db.executemany("INSERT INTO detail VALUES (?,?,?,?)", entries)
    db.execute("CREATE INDEX detail_shard ON detail(shard,guide)")
    db.commit()
    shard_count = (len(guides) + _SHARD_GUIDES - 1) // _SHARD_GUIDES
    for shard_no in range(shard_count):
        data: dict[str, dict[str, list[Any]]] = {}
        for guide, kind, payload in db.execute(
            "SELECT guide,kind,payload FROM detail WHERE shard=? ORDER BY guide,kind", (shard_no,)
        ):
            data.setdefault(str(guide), {"candidates": [], "hits": [], "mirna": []})[str(kind)].append(
                json.loads(payload)
            )
        name = f"guide-{shard_no:05d}.js"
        target = shard_dir / name
        text = json.dumps(data, separators=(",", ":")).replace("<", "\\u003c").replace("&", "\\u0026")
        target.write_text(
            f"window.SF_DETAIL_SHARDS=window.SF_DETAIL_SHARDS||{{}};window.SF_DETAIL_SHARDS[{shard_no}]={text};\n"
        )
    db.close()
    db_path.unlink(missing_ok=True)
    store.manifest["detail_shards"] = [
        {"index": index, "file": path.name, "bytes": path.stat().st_size, "sha256": _sha256(path)}
        for index, path in enumerate(sorted(shard_dir.glob("*.js")))
    ]
    guide_rows = pq.read_table(store.root / store.manifest["tables"]["guides"]["file"]).to_pylist()
    for index, row in enumerate(guide_rows):
        row["detail_shard"] = index // _SHARD_GUIDES
    _write_records(store, "guides", guide_rows)
    store.manifest["tables"].update(_table_metadata(store, "guides"))


def _write_data_backed_html(store: ReportArtifactStore, output: Path, shard_dir: str) -> None:
    guides = pq.read_table(store.root / store.manifest["tables"]["guides"]["file"]).to_pylist()
    summary = store.manifest["summary"]
    for i, row in enumerate(guides):
        row["shard"] = i // _SHARD_GUIDES
    manifest = json.dumps(store.manifest, separators=(",", ":")).replace("<", "\\u003c").replace("&", "\\u0026")
    data = json.dumps(guides, separators=(",", ":")).replace("<", "\\u003c").replace("&", "\\u0026")
    html = _DATA_BACKED_TEMPLATE.replace("__GUIDES__", data).replace("__MANIFEST__", manifest)
    html = html.replace("__MODE__", str(store.manifest["report_mode"]))
    html = html.replace("__SHARD_DIR__", shard_dir).replace("__STATUS_COUNTS__", json.dumps(summary["status_counts"]))
    html = html.replace("<!--__INLINE_SHARDS__-->", "")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)


def _write_paged_html(store: ReportArtifactStore, output: Path, shard_dir: Path) -> None:
    """Inline every detail shard, while bounding the visible guide list to one page."""
    guides = pq.read_table(store.root / store.manifest["tables"]["guides"]["file"]).to_pylist()
    for index, row in enumerate(guides):
        row["shard"] = index // _SHARD_GUIDES
    data = json.dumps(guides, separators=(",", ":")).replace("<", "\\u003c").replace("&", "\\u0026")
    manifest = json.dumps(store.manifest, separators=(",", ":")).replace("<", "\\u003c").replace("&", "\\u0026")
    scripts = "\n".join(path.read_text() for path in sorted(shard_dir.glob("*.js")))
    html = _DATA_BACKED_TEMPLATE.replace("__GUIDES__", data).replace("__MANIFEST__", manifest)
    html = html.replace("__MODE__", str(store.manifest["report_mode"]))
    html = html.replace("__SHARD_DIR__", "").replace(
        "__STATUS_COUNTS__", json.dumps(store.manifest["summary"]["status_counts"])
    )
    html = html.replace("<!--__INLINE_SHARDS__-->", f"<script>\n{scripts}\n</script>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)


def _compact_manifest(run_dir: Path, mode: str, thresholds: ReportModeThresholds | None) -> dict[str, Any]:
    return {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "report_mode": mode,
        "source_run_label": run_dir.name,
        "mode_thresholds": (thresholds or ReportModeThresholds()).__dict__,
        "mode_facts": None,
        "tables": {},
        "inputs": {},
        "exports": {},
    }


def _write_result_metadata(
    output: Path,
    mode: str,
    manifest: dict[str, Any],
    builder_diagnostics: dict[str, Any],
    phases: dict[str, float],
    build_started: float,
    artifact_dir: Path | None = None,
) -> ReportWriteResult:
    diagnostics = dict(builder_diagnostics)
    diagnostics["schema_version"] = 1
    diagnostics["report_mode"] = mode
    diagnostics["mode_selection"] = {
        "mode": mode,
        "thresholds": manifest.get("mode_thresholds"),
        "facts": manifest.get("mode_facts"),
    }
    diagnostics["phases_seconds"] = {**diagnostics.get("phases_seconds", {}), **phases}
    diagnostics["peak_memory_scope"] = diagnostics.get("peak_memory_scope", "process_lifetime")
    diagnostics["peak_rss_bytes"] = diagnostics.get("peak_rss_bytes") or _peak_rss_bytes()
    diagnostics["output_html_bytes"] = output.stat().st_size
    diagnostics["total_seconds"] = time.perf_counter() - build_started
    manifest["report_mode"] = mode
    manifest["diagnostics"] = diagnostics
    manifest_path = (artifact_dir / "report_manifest.json") if artifact_dir else output.with_suffix(".manifest.json")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    diagnostics_path = output.with_suffix(".diagnostics.json")
    diagnostics_path.write_text(json.dumps(diagnostics, indent=2) + "\n")
    return ReportWriteResult(output, mode, manifest, manifest_path, diagnostics, diagnostics_path, artifact_dir)


def _float_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _peak_rss_bytes() -> int | None:
    try:
        resource_module = importlib.import_module("resource")
    except ImportError:  # pragma: no cover - resource unavailable on Windows
        return None
    value = int(resource_module.getrusage(resource_module.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


_DATA_BACKED_TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>siRNAforge data-backed report</title><style>body{font:14px system-ui;margin:0;color:#1a1d21;background:#fbfbfc}header{padding:18px 24px;background:white;border-bottom:1px solid #ddd}main{display:grid;grid-template-columns:minmax(420px,45%) 1fr;height:calc(100vh - 112px)}section{overflow:auto;padding:12px}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border-bottom:1px solid #ddd;padding:5px;text-align:left}.pill{padding:2px 6px;border-radius:8px;background:#eee}.fail{color:#b91c1c}.pass{color:#15803d}.warn{color:#92400e}.unknown{color:#6b7280}button{font:inherit;cursor:pointer}#details{white-space:normal}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style></head>
<body><header><h1>siRNAforge report · __MODE__</h1><div id="counts"></div><div id="gatecounts"></div><label>Search guides <input id="search" type="search"></label><label>Status <select id="status"><option value="">All</option><option>pass</option><option>warn</option><option>unknown</option><option>fail</option></select></label><label>Gate <select id="gate"><option value="">Any gate</option></select></label><label>Gate result <select id="gate-result"><option value="">Any result</option><option value="0">pass</option><option value="1">fail</option><option value="2">unknown</option><option value="3">not evaluated</option><option value="4">warn</option></select></label><span id="visible-count"></span><button id="previous">Previous page</button><span id="page-info"></span><button id="next">Next page</button></header>
<main><section><table><thead><tr><th>Guide</th><th>Status</th><th>Rows</th><th>Design</th><th>Composite</th></tr></thead><tbody id="rows"></tbody></table></section><section id="details" aria-live="polite"><p>Select any guide to load its detail shard.</p></section></main>
<!--__INLINE_SHARDS__-->
<script>
const PAGE_SIZE=50;let currentPage=0;const G=__GUIDES__, M=__MANIFEST__, SHARD_DIR='__SHARD_DIR__', STATUS_COUNTS=__STATUS_COUNTS__;
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function renderRows(){const q=document.querySelector('#search').value.trim().toUpperCase(),st=document.querySelector('#status').value;
 const gate=document.querySelector('#gate').selectedIndex-1,result=document.querySelector('#gate-result').value;
 const shown=G.filter(g=>(!q||g.guide.includes(q))&&(!st||g.status===st)&&
   (gate<0||!result||JSON.parse(g.gates)[gate]?.[1]===Number(result)));
 const pageCount=Math.max(1,Math.ceil(shown.length/PAGE_SIZE));currentPage=Math.min(currentPage,pageCount-1);
 const pageRows=shown.slice(currentPage*PAGE_SIZE,(currentPage+1)*PAGE_SIZE);
 document.querySelector('#page-info').textContent=`Page ${currentPage+1} of ${pageCount}`;
 document.querySelector('#previous').disabled=currentPage===0;document.querySelector('#next').disabled=currentPage>=pageCount-1;
 document.querySelector('#visible-count').textContent=`${shown.length} guides match current filters`;
 document.querySelector('#rows').innerHTML=pageRows.map(g=>`<tr><td><button data-guide="${esc(g.guide)}">${esc(g.guide)}</button></td><td class="${g.status}">${g.status}</td><td>${g.candidate_rows}</td><td>${g.design_score??'—'}</td><td>${g.composite_score??'—'}</td></tr>`).join('');
 document.querySelectorAll('[data-guide]').forEach(b=>b.addEventListener('click',()=>showGuide(b.dataset.guide)));
}
function loadShard(i){return new Promise((resolve,reject)=>{if(window.SF_DETAIL_SHARDS?.[i])return resolve();const s=document.createElement('script');s.src=`${SHARD_DIR}/guide-${String(i).padStart(5,'0')}.js`;s.onload=resolve;s.onerror=()=>reject(Error('Could not load local detail shard'));document.head.appendChild(s);});}
async function showGuide(guide){const row=G.find(g=>g.guide===guide);if(!row)return;const panel=document.querySelector('#details');panel.textContent='Loading guide detail…';try{await loadShard(row.shard);const d=window.SF_DETAIL_SHARDS[row.shard][guide];panel.innerHTML=`<h2>${esc(guide)} · ${esc(row.status)}</h2><p>${row.candidate_rows} candidate rows · ${d.hits.length} classified hit rows · ${d.mirna.length} miRNA rows</p><h3>Candidate rows</h3><pre>${esc(JSON.stringify(d.candidates,null,2))}</pre><h3>Off-target rows</h3><pre>${esc(JSON.stringify(d.hits,null,2))}</pre><h3>miRNA rows</h3><pre>${esc(JSON.stringify(d.mirna,null,2))}</pre>`;}catch(e){panel.textContent=e.message;}}
document.querySelector('#counts').textContent=Object.entries(STATUS_COUNTS).map(([k,v])=>`${k}: ${v}`).join(' · ')+` · total ${G.length}`;
document.querySelector('#gatecounts').textContent=(M.summary.filters||[]).map(f=>`${f.filter_id}: ${f.rejects} reject, ${f.warns} warn, ${f.unknowns} unknown`).join(' · ');
for(const f of M.summary.filters){const o=document.createElement('option');o.value=f.filter_id;o.textContent=f.filter_id;document.querySelector('#gate').appendChild(o);}
document.querySelector('#previous').addEventListener('click',()=>{currentPage=Math.max(0,currentPage-1);renderRows()});document.querySelector('#next').addEventListener('click',()=>{currentPage++;renderRows()});
document.querySelector('#search').addEventListener('input',()=>{currentPage=0;renderRows()});document.querySelector('#status').addEventListener('change',()=>{currentPage=0;renderRows()});document.querySelector('#gate').addEventListener('change',()=>{currentPage=0;renderRows()});document.querySelector('#gate-result').addEventListener('change',()=>{currentPage=0;renderRows()});renderRows();
</script></body></html>"""
