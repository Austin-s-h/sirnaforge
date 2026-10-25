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
from html import escape
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
from sirnaforge.reporting.evidence import (
    EVIDENCE_STATES,
    FILTER_UNITS,
    GATE_REASON_IDS,
    WARNING_REGISTRY,
    conservation_evidence,
    guide_mirna_evidence,
    mirna_evidence,
    report_evidence,
)
from sirnaforge.reporting.links import canonical_identity_link, source_link
from sirnaforge.reporting.payload import (
    _VERDICT_CODE,
    REPORT_READ_CHUNK_ROWS,
    STATUSES,
    VERDICT_WARN,
    _evaluate,
    _flag,
    _gate_panel,
    _normalise_guide,
    _read_manifest,
    _recorded_verdict,
    _run_verdict,
    build_payload,
    observed_column,
)
from sirnaforge.reporting.render import write_report

_SHARD_GUIDES = 128
_GUIDE_SCHEMA = pa.schema(
    [
        ("guide", pa.string()),
        ("mirna_evidence_state", pa.string()),
        ("mirna_warning_id", pa.string()),
        ("unresolved_evidence", pa.bool_()),
        ("report_evidence", pa.string()),
        ("candidate_rows", pa.int64()),
        ("design_score", pa.float64()),
        ("composite_score", pa.float64()),
        ("run_verdict", pa.string()),
        ("canonical_hit", pa.string()),
        ("target_genes", pa.string()),
        ("ortholog", pa.string()),
        ("off_target_screened", pa.bool_()),
        ("status", pa.string()),
        ("gates", pa.string()),
        ("n_failed_gates", pa.int64()),
        ("n_unknown_gates", pa.int64()),
        ("n_warned_gates", pa.int64()),
        ("detail_shard", pa.int64()),
    ]
)


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
    guide_index_file: str | None = None
    if selected == "paged_single_file":
        detail_directory = None
        detail_loading = "inline_paged_data"
    else:
        phase_started = time.perf_counter()
        guide_index = _write_guide_index(store, shard_dir)
        phases["guide_index_write"] = time.perf_counter() - phase_started
        guide_index_file = str(guide_index["file"])
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
        "guide_index_file": guide_index_file,
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
    reference_evidence = mirna_evidence(manifest)
    store.manifest["reference_evidence"] = {
        key: value for key, value in reference_evidence.items() if key != "guide_sequences"
    }
    store.manifest["conservation_evidence"] = conservation_evidence(manifest)
    store.manifest["evidence_warning_registry"] = WARNING_REGISTRY
    store.manifest["evidence_states"] = EVIDENCE_STATES
    store.manifest["gate_reason_ids"] = GATE_REASON_IDS
    store.manifest["gate_panel"] = {
        "source": panel.source,
        "profile": panel.profile_name,
        "run_mode": panel.run_mode,
        "filters": [
            {**gate.model_dump(mode="json"), "units": FILTER_UNITS.get(gate.setting_key, "not declared")}
            for gate in panel.filters
        ],
    }
    identity_value = manifest.get("target_identity")
    identity: dict[str, Any] = identity_value if isinstance(identity_value, dict) else {}
    identity_resolvers = identity.get("resolvers")
    identity_resolvers = identity_resolvers if isinstance(identity_resolvers, dict) else {}
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
        and entry.get("canonical_release")
        and entry.get("canonical_status") not in {"not_established", "unresolved_input", "unknown"}
    }
    store.manifest["target_identity"] = identity
    populated: set[str] = set()
    transcript_register: set[str] = set()
    temp_db = store.root / ".derive.sqlite"
    db = sqlite3.connect(temp_db)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE candidate_rows (guide TEXT, payload TEXT)")
    db.execute("CREATE INDEX candidate_guide ON candidate_rows(guide)")
    for batch in store.iter_batches("candidates", batch_size=REPORT_READ_CHUNK_ROWS):
        records = batch.to_pylist()
        for record in records:
            if record.get("transcript_id"):
                transcript_register.add(str(record["transcript_id"]))
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
    for gate, column in zip(store.manifest["gate_panel"]["filters"], gate_columns, strict=True):
        gate["read_column"] = column
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
                ("transcript_link", pa.string()),
                ("gene_link", pa.string()),
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
        guide_genes: dict[str, dict[str, str]] = {}
        for row in rows:
            tx = str(row.get("transcript_id") or "")
            identity_row = identity_by_transcript.get(tx, {})
            gene_id = str(identity_row.get("gene_id") or row.get("gene_id") or "").strip()
            gene_symbol = str(
                identity_row.get("gene_symbol") or row.get("gene_symbol") or row.get("gene_name") or ""
            ).strip()
            gene_entry = _target_gene_entry(gene_id, gene_symbol)
            guide_genes[gene_entry["key"]] = gene_entry
            if not tx or tx in iso_seen:
                continue
            iso_seen.add(tx)
            canonical_status = identity_row.get("canonical_status") or "unknown"
            canonical_value = (
                "canonical"
                if tx in canonical_ids
                else "other_isoform"
                if identity_row.get("canonical_source")
                and identity_row.get("canonical_release")
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
                    "transcript_link": json.dumps(
                        canonical_identity_link(
                            tx,
                            identity_row.get("canonical_source"),
                            identity_row.get("canonical_release"),
                            resolver=identity_resolvers.get("ensembl"),
                        )
                        if identity_row.get("canonical_release")
                        else {
                            "identifier": tx,
                            "source": "ensembl",
                            "release": None,
                            "url": None,
                            "reason": "canonical_release_unavailable",
                        },
                        separators=(",", ":"),
                    ),
                    "gene_link": json.dumps(
                        canonical_identity_link(
                            str(identity_row.get("gene_id") or row.get("gene_id") or ""),
                            identity_row.get("canonical_source"),
                            identity_row.get("canonical_release"),
                            resolver=identity_resolvers.get("ensembl"),
                        )
                        if identity_row.get("canonical_release")
                        else {
                            "identifier": str(identity_row.get("gene_id") or row.get("gene_id") or ""),
                            "source": "ensembl",
                            "release": None,
                            "url": None,
                            "reason": "canonical_release_unavailable",
                        },
                        separators=(",", ":"),
                    ),
                    "candidate_rows": sum(1 for item in rows if str(item.get("transcript_id") or "") == tx),
                }
            )
        diagnostic = guide_mirna_evidence(reference_evidence, guide)
        report_diagnostic = report_evidence(
            gates,
            diagnostic,
            undeclared_rejection=undeclared,
            recorded_states=[_recorded_verdict(best, desc.filter_id) for desc in descriptors],
        )
        guide_records.append(
            {
                "guide": guide,
                "mirna_evidence_state": diagnostic["state"],
                "mirna_warning_id": diagnostic["warning_id"],
                "unresolved_evidence": report_diagnostic["unresolved"],
                "report_evidence": json.dumps(report_diagnostic, separators=(",", ":")),
                "candidate_rows": len(rows),
                "design_score": _float_or_none(best.get("design_score")),
                "composite_score": _float_or_none(best.get("composite_score")),
                "run_verdict": run_verdict,
                "canonical_hit": (
                    "hit"
                    if set(iso_seen) & canonical_ids
                    else "not_hit"
                    if canonical_ids and canonical_ids <= transcript_register
                    else "not_established"
                ),
                "target_genes": json.dumps([guide_genes[key] for key in sorted(guide_genes)], separators=(",", ":")),
                "ortholog": "{}",
                "off_target_screened": _flag(best.get("off_target_screened")),
                "status": status,
                "gates": json.dumps(gates, separators=(",", ":")),
                "n_failed_gates": n_failed,
                "n_unknown_gates": n_unknown,
                "n_warned_gates": n_warned,
                "detail_shard": None,
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

    conservation_by_guide = _derive_hit_summary(store)
    for record in guide_records:
        record["ortholog"] = json.dumps(conservation_by_guide.get(record["guide"], {}), separators=(",", ":"))
    _write_records(store, "guides", guide_records)
    isoform_sink.close()
    target_gene_entries = {entry["key"]: entry for row in guide_records for entry in json.loads(row["target_genes"])}
    for identity_row in identity_by_transcript.values():
        entry = _target_gene_entry(identity_row.get("gene_id"), identity_row.get("gene_symbol"))
        target_gene_entries.setdefault(entry["key"], entry)
    known_target_genes = [entry for key, entry in target_gene_entries.items() if key != "__unresolved__"]
    store.manifest["summary"] = {
        "guides": len(guide_records),
        "unresolved_guides": sum(row["unresolved_evidence"] for row in guide_records),
        "passes_without_unresolved_evidence": sum(
            row["status"] == "pass" and not row["unresolved_evidence"] for row in guide_records
        ),
        "mirna_unresolved_guides": sum(
            row["mirna_evidence_state"] not in {"complete", "not_requested"} for row in guide_records
        ),
        "evidence_warning_counts": dict(
            Counter(
                warning_id for row in guide_records for warning_id in json.loads(row["report_evidence"])["warning_ids"]
            )
        ),
        "status_counts": {status: status_counts.get(status, 0) for status in STATUSES},
        "filters": [
            {
                "filter_id": desc.filter_id,
                **{key: gate_effects[(index, key)] for key in ("rejects", "sole_rejects", "warns", "unknowns")},
            }
            for index, desc in enumerate(descriptors)
        ],
        "candidate_rows": sum(row["candidate_rows"] for row in guide_records),
        "target_gene_options": sorted(known_target_genes, key=lambda entry: (entry["label"].casefold(), entry["key"])),
        "gene_filter_available": len(known_target_genes) > 1 or "__unresolved__" in target_gene_entries,
        "unresolved_target_identity": "__unresolved__" in target_gene_entries,
    }
    store.manifest["tables"].update(_table_metadata(store, "guides"))
    store.manifest["tables"].update(_table_metadata(store, "isoforms"))
    (store.root / "report_manifest.json").write_text(json.dumps(store.manifest, indent=2) + "\n")


def _derive_hit_summary(store: ReportArtifactStore) -> dict[str, dict[str, dict[str, Any]]]:
    """Aggregate hit rows in SQLite so memory is bounded by the read batch, not the hit table."""
    db_path = store.root / ".hits.sqlite"
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE hits (guide TEXT,species TEXT,hit_class TEXT,nm_band TEXT,symbol TEXT)")
    best_orthologs: dict[str, dict[str, tuple[tuple[float, float], dict[str, Any]]]] = {}
    if "hits" in store.manifest["tables"]:
        available = set(store.manifest["tables"]["hits"]["columns"])
        selected = [
            name
            for name in (
                "qseq",
                "species",
                "hit_class",
                "nm",
                "seed_mismatches",
                "rname",
                "hit_symbol",
                "matched_symbol",
            )
            if name in available
        ]
        for batch in store.iter_batches("hits", columns=selected):
            rows = []
            for row in batch.to_pylist():
                nm = _float_or_none(row.get("nm"))
                seed = _float_or_none(row.get("seed_mismatches"))
                guide = _normalise_guide(row.get("qseq"))
                species = str(row.get("species") or "query").strip() or "query"
                if str(row.get("hit_class") or "").strip().casefold() == "ortholog" and nm is not None:
                    rank = (nm, seed if seed is not None else math.inf)
                    current = best_orthologs.setdefault(guide, {}).get(species)
                    if current is None or rank < current[0]:
                        best_orthologs[guide][species] = (
                            rank,
                            {
                                "nm": int(nm) if nm.is_integer() else None,
                                "seed_mismatches": int(seed) if seed is not None and seed.is_integer() else None,
                                "transcript_id": str(row.get("rname") or "") or None,
                                "gene_symbol": str(row.get("hit_symbol") or row.get("matched_symbol") or "") or None,
                            },
                        )
                band = "unknown" if nm is None else "0" if nm == 0 else "1" if nm == 1 else "2" if nm == 2 else ">=3"
                rows.append(
                    (
                        guide,
                        species,
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
    return {
        guide: {species: value for species, (_rank, value) in by_species.items()}
        for guide, by_species in best_orthologs.items()
    }


def _write_records(store: ReportArtifactStore, name: str, records: list[dict[str, Any]]) -> None:
    path = store.root / f"{name}.parquet"
    schema = (
        _GUIDE_SCHEMA
        if name == "guides"
        else pa.Table.from_pylist(records[:1]).schema
        if records
        else pa.schema([("_empty", pa.string())])
    )
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
    target_identity = store.manifest.get("target_identity")
    resolver_value = (
        target_identity.get("resolvers")
        if isinstance(target_identity, dict) and isinstance(target_identity.get("resolvers"), dict)
        else None
    )
    resolvers: dict[str, Any] = resolver_value if isinstance(resolver_value, dict) else {}
    db_path = store.root / ".detail-shards.sqlite"
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA journal_mode=OFF")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE detail (shard INTEGER, guide TEXT, kind TEXT, payload TEXT)")
    detail_sources = (
        ("candidates", "guide_sequence", "candidates"),
        ("hits", "qseq", "hits"),
        ("mirna", "qseq", "mirna"),
        ("isoforms", "guide", "isoforms"),
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
                    if kind == "mirna":
                        identifier = next(
                            (str(row[key]) for key in ("mirna_id", "rname", "mirna", "name") if row.get(key)),
                            "",
                        )
                        row["source_link"] = source_link(
                            identifier,
                            "mirgenedb",
                            resolver=resolvers.get("mirgenedb"),
                        )
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
            data.setdefault(str(guide), {"candidates": [], "hits": [], "mirna": [], "isoforms": []})[str(kind)].append(
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


def _write_guide_index(store: ReportArtifactStore, shard_dir: Path) -> dict[str, Any]:
    """Publish the summary rows as a checksummed script beside data-backed HTML.

    Keeping the guide index out of the HTML preserves file-origin loading while avoiding a
    guide-count-proportional JSON tail in the entrypoint. Detail shards remain separate and are
    still loaded only after a reader selects a guide.
    """
    guides = pq.read_table(store.root / store.manifest["tables"]["guides"]["file"]).to_pylist()
    for index, row in enumerate(guides):
        row["shard"] = index // _SHARD_GUIDES
    target = shard_dir / "guide-index.js"
    payload = json.dumps(guides, separators=(",", ":")).replace("<", "\\u003c").replace("&", "\\u0026")
    target.write_text(f"window.SF_GUIDE_INDEX={payload};\n")
    metadata = {
        "file": target.name,
        "format": "javascript",
        "rows": len(guides),
        "bytes": target.stat().st_size,
        "sha256": _sha256(target),
    }
    store.manifest["guide_index"] = metadata
    return metadata


def _write_data_backed_html(store: ReportArtifactStore, output: Path, shard_dir: str) -> None:
    summary = store.manifest["summary"]
    manifest = json.dumps(store.manifest, separators=(",", ":")).replace("<", "\\u003c").replace("&", "\\u0026")
    guide_index_file = store.manifest["reader"]["guide_index_file"]
    guide_index_script = f'<script src="{shard_dir}/{guide_index_file}"></script>'
    html = _DATA_BACKED_TEMPLATE.replace("__GUIDES__", "window.SF_GUIDE_INDEX||[]")
    html = html.replace("__GUIDE_INDEX_SCRIPT__", guide_index_script).replace("__MANIFEST__", manifest)
    html = html.replace("__MODE__", str(store.manifest["report_mode"]))
    html = html.replace("__EVIDENCE_PANEL__", _evidence_panel(store.manifest))
    unresolved_count = store.manifest["summary"]["unresolved_guides"]
    html = html.replace("__UNRESOLVED_CONTROL__", _unresolved_control(unresolved_count))
    html = html.replace("__GENE_FILTER_CONTROL__", _gene_filter_control(store.manifest["summary"]))
    html = html.replace(
        "__CONSERVATION_FILTER_CONTROL__", _conservation_filter_control(store.manifest["conservation_evidence"])
    )
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
    html = html.replace("__GUIDE_INDEX_SCRIPT__", "")
    html = html.replace("__MODE__", str(store.manifest["report_mode"]))
    html = html.replace("__EVIDENCE_PANEL__", _evidence_panel(store.manifest))
    unresolved_count = store.manifest["summary"]["unresolved_guides"]
    html = html.replace("__UNRESOLVED_CONTROL__", _unresolved_control(unresolved_count))
    html = html.replace("__GENE_FILTER_CONTROL__", _gene_filter_control(store.manifest["summary"]))
    html = html.replace(
        "__CONSERVATION_FILTER_CONTROL__", _conservation_filter_control(store.manifest["conservation_evidence"])
    )
    html = html.replace("__SHARD_DIR__", "").replace(
        "__STATUS_COUNTS__", json.dumps(store.manifest["summary"]["status_counts"])
    )
    html = html.replace("<!--__INLINE_SHARDS__-->", f"<script>\n{scripts}\n</script>")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)


def _unresolved_control(count: int) -> str:
    hidden = " hidden" if count == 0 else ""
    return f'<label{hidden}><input id="unresolved-only" type="checkbox">Show unresolved evidence ({count})</label>'


def _gene_filter_control(summary: dict[str, Any]) -> str:
    """Render a gene filter only when the staged artifact confirms multiple or unresolved targets."""
    if not summary.get("gene_filter_available"):
        return ""
    options = ['<option value="">All target genes</option>']
    options.extend(
        f'<option value="{escape(str(item["key"]), quote=True)}">{escape(str(item["label"]))}</option>'
        for item in summary.get("target_gene_options", [])
    )
    if summary.get("unresolved_target_identity"):
        options.append('<option value="__unresolved__">Unknown target gene</option>')
    return (
        '<label><input type="checkbox" id="geneenable" aria-controls="targetgene" '
        'aria-expanded="false">Enable target gene filtering</label>'
        '<label for="targetgene">Target gene</label>'
        '<select id="targetgene" aria-label="Target gene" disabled>'
        + "".join(options)
        + '</select><span id="gene-state" role="status" aria-live="polite"></span>'
    )


def _conservation_filter_control(evidence: dict[str, Any]) -> str:
    """Render accessible per-species conservation filters backed by published orthology evidence."""
    species = evidence.get("denominator_species") or list((evidence.get("per_species") or {}).keys())
    species = list(dict.fromkeys(str(item) for item in species if item))
    if not species:
        return '<p id="conservation-unavailable">No conservation species were recorded for this run.</p>'
    controls = "".join(
        f'<label><input type="checkbox" data-conservation value="{escape(item, quote=True)}"> {escape(item)}</label>'
        for item in species
    )
    return (
        '<fieldset id="conservation-filters" aria-describedby="cons-state consdiagnostic">'
        "<legend>Conserved in</legend>"
        + controls
        + '<label><input type="checkbox" id="cons-seed" checked> Allow one mismatch only outside guide positions 2–8</label>'
        '<div id="cons-state" role="status" aria-live="polite"></div>'
        '<div id="consdiagnostic" role="status" aria-live="polite"></div></fieldset>'
    )


def _evidence_panel(manifest: dict[str, Any]) -> str:
    """Render shared evidence copy and an independently scrollable policy inspector."""
    summary = manifest["summary"]
    registry = manifest["evidence_warning_registry"]
    notices = []
    for warning_id, count in summary["evidence_warning_counts"].items():
        warning = registry[warning_id]
        notices.append(
            f'<li data-warning-id="{escape(warning_id)}"><b>{escape(warning["title"])}</b> '
            f"({count} guides)<p>{escape(warning.get('summary_body', warning['body']))}</p>"
            f"<p>{escape(warning['remediation'])}</p></li>"
        )
    evidence = manifest["reference_evidence"]
    sources = []
    for source in evidence["sources"]:
        references = "; ".join(
            ", ".join(f"{key}: {value}" for key, value in reference.items()) for reference in source["references"]
        )
        reasons = "; ".join(
            str(source[phase]["reason"]) for phase in ("acquisition", "screening") if source[phase].get("reason")
        )
        sources.append(
            f"<li>{escape(str(source['source']))} / {escape(str(source['species']))}: "
            f"{escape(source['state'])}<p>{escape(reasons)}</p><p>{escape(references)}</p></li>"
        )
    notice_html = ""
    if notices:
        notice_html = (
            '<details id="evidence-notices"><summary>Evidence notices</summary><ul>'
            + "".join(notices)
            + "</ul></details>"
        )
    panel = manifest["gate_panel"]
    filters = []
    for gate in panel["filters"]:
        descriptor = gate["descriptor"]
        scope = json.dumps(descriptor["scope"], sort_keys=True)
        filters.append(
            f"<details><summary>{escape(descriptor['filter_id'].replace('_', ' '))}</summary><dl>"
            f"<dt>Setting</dt><dd>{escape(gate['setting_key'])}</dd>"
            f"<dt>Threshold</dt><dd>{escape(descriptor['comparator'])} {escape(str(descriptor['threshold']))}</dd>"
            f"<dt>Units</dt><dd>{escape(gate['units'])}</dd>"
            f"<dt>Action</dt><dd>{escape(descriptor['action'])}</dd>"
            f"<dt>Stage</dt><dd>{escape(descriptor['stage'])}</dd>"
            f"<dt>Scope</dt><dd>{escape(scope)}</dd>"
            f"<dt>Definition</dt><dd>{escape(gate['definition'])}</dd></dl></details>"
        )
    return (
        '<aside aria-label="Evidence and thresholds">'
        "<h2>Evidence</h2><p>Gate results use the published policy and observations. View filters do not change the run. "
        "Unresolved evidence is listed separately; absent matches do not establish a clean screen.</p>"
        + notice_html
        + (
            "<details><summary>miRNA reference provenance</summary><ul>" + "".join(sources) + "</ul></details>"
            if sources
            else ""
        )
        + "<h2>Thresholds</h2><p>Source: "
        + escape(panel["source"])
        + " · profile: "
        + escape(panel["profile"])
        + "</p>"
        + "<p>Pass: criterion met. Warn: finding retained. Fail: criterion exceeded. "
        "Unknown: active gate lacks evidence. Not evaluated: gate was off or unconfigured.</p>"
        + "".join(filters)
        + "</aside>"
    )


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


def _target_gene_entry(gene_id: Any, gene_symbol: Any) -> dict[str, str]:
    """Build one stable gene-filter key from available identity, preserving unresolved targets."""
    identifier = str(gene_id or "").strip()
    symbol = str(gene_symbol or "").strip()
    if identifier:
        return {
            "key": f"id:{identifier}",
            "label": f"{symbol} · {identifier}" if symbol else identifier,
        }
    if symbol:
        return {"key": f"symbol:{symbol.casefold()}", "label": symbol}
    return {"key": "__unresolved__", "label": "Unknown target gene"}


def _peak_rss_bytes() -> int | None:
    try:
        resource_module = importlib.import_module("resource")
    except ImportError:  # pragma: no cover - resource unavailable on Windows
        return None
    value = int(resource_module.getrusage(resource_module.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


_DATA_BACKED_TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>siRNAforge data-backed report</title><style>html,body{height:100%}body{display:flex;flex-direction:column;font:14px system-ui;margin:0;color:#1a1d21;background:#fbfbfc}header{flex:0 0 auto;padding:18px 24px;background:white;border-bottom:1px solid #ddd}main{display:grid;grid-template-columns:minmax(220px,25%) minmax(350px,35%) 1fr;grid-template-rows:minmax(0,1fr);flex:1 1 auto;min-height:0}section,aside{min-height:0;overflow:auto;padding:12px}aside{border-right:1px solid #ddd;font-size:12px;overflow-wrap:anywhere}dd{margin-left:12px}details{margin:8px 0}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border-bottom:1px solid #ddd;padding:5px;text-align:left}.pill{padding:2px 6px;border-radius:8px;background:#eee}.fail{color:#b91c1c}.pass{color:#15803d}.warn{color:#92400e}.unknown{color:#6b7280}button{font:inherit;cursor:pointer}#details{white-space:normal}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style></head>
<body><header><h1>siRNAforge report · __MODE__</h1><div id="counts"></div><div id="gatecounts"></div><label>Search guides <input id="search" type="search"></label><label>Status <select id="status"><option value="">All</option><option>pass</option><option>warn</option><option>unknown</option><option>fail</option></select></label><label>Gate <select id="gate"><option value="">Any gate</option></select></label><label>Gate result <select id="gate-result"><option value="">Any result</option><option value="0">pass</option><option value="1">fail</option><option value="2">unknown</option><option value="3">not evaluated</option><option value="4">warn</option></select></label>__GENE_FILTER_CONTROL____CONSERVATION_FILTER_CONTROL____UNRESOLVED_CONTROL__<span id="visible-count"></span><button id="previous">Previous page</button><span id="page-info"></span><button id="next">Next page</button></header>
<main>__EVIDENCE_PANEL__<section><table><thead><tr><th>Guide</th><th>Status</th><th>Rows</th><th>Design</th><th>Composite</th></tr></thead><tbody id="rows"></tbody></table></section><section id="details" aria-live="polite"><p>Select any guide to load its detail shard.</p></section></main>
<!--__INLINE_SHARDS__-->
__GUIDE_INDEX_SCRIPT__
<script>
const PAGE_SIZE=50;let currentPage=0;const G=__GUIDES__, M=__MANIFEST__, SHARD_DIR='__SHARD_DIR__', STATUS_COUNTS=__STATUS_COUNTS__;
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function renderRows(){const q=document.querySelector('#search').value.trim().toUpperCase(),st=document.querySelector('#status').value;
 const gate=document.querySelector('#gate').selectedIndex-1,result=document.querySelector('#gate-result').value;
 const unresolved=document.querySelector('#unresolved-only')?.checked||false;
 const geneEnabled=document.querySelector('#geneenable')?.checked||false,gene=document.querySelector('#targetgene')?.value||'';
 const selectedSpecies=[...F.conserved];
 const shown=G.filter(g=>(!q||g.guide.includes(q))&&(!st||g.status===st)&&(!unresolved||g.unresolved_evidence)&&
   (!geneEnabled||!gene||(JSON.parse(g.target_genes||'[]').some(item=>item.key===gene)))&&
   selectedSpecies.every(sp=>{const d=(JSON.parse(g.ortholog||'{}'))[sp];return d&&d.nm!==null&&d.nm<=1&&(d.nm!==1||!F.conservedSeed||d.seed_mismatches===0);})&&
   (gate<0||!result||JSON.parse(g.gates)[gate]?.[1]===Number(result)));
 const pageCount=Math.max(1,Math.ceil(shown.length/PAGE_SIZE));currentPage=Math.min(currentPage,pageCount-1);
 const pageRows=shown.slice(currentPage*PAGE_SIZE,(currentPage+1)*PAGE_SIZE);
 document.querySelector('#page-info').textContent=`Page ${currentPage+1} of ${pageCount}`;
 document.querySelector('#previous').disabled=currentPage===0;document.querySelector('#next').disabled=currentPage>=pageCount-1;
 document.querySelector('#visible-count').textContent=`${shown.length} guides match current filters`;
 const geneState=document.querySelector('#gene-state');if(geneState)geneState.textContent=geneEnabled&&gene?`Target gene filter active: ${gene}`:'';
 const conservationDiagnostic=document.querySelector('#consdiagnostic');if(conservationDiagnostic)conservationDiagnostic.textContent=conservationFilterDiagnostic();
 const conservationState=document.querySelector('#cons-state');if(conservationState)conservationState.textContent=F.conserved.size?`Conservation filter active: ${[...F.conserved].sort().join(', ')}${F.conservedSeed?'; one mismatch permitted only outside positions 2–8':'; up to one mismatch including the seed'}`:'';
 document.querySelector('#rows').innerHTML=pageRows.map(g=>`<tr><td><button data-guide="${esc(g.guide)}">${esc(g.guide)}</button></td><td class="${g.status}">${g.status}</td><td>${g.candidate_rows}</td><td>${g.design_score??'—'}</td><td>${g.composite_score??'—'}</td></tr>`).join('');
 document.querySelectorAll('[data-guide]').forEach(b=>b.addEventListener('click',()=>showGuide(b.dataset.guide)));
}
const F={conserved:new Set(),conservedSeed:true};
const CONSERVATION_SPECIES=(M.conservation_evidence?.denominator_species?.length?M.conservation_evidence.denominator_species:Object.keys(M.conservation_evidence?.per_species||{})).filter((sp,i,all)=>all.indexOf(sp)===i);
function conservationDiagnostic(g,sp){
 const d=(JSON.parse(g.ortholog||'{}'))[sp];if(d)return {unavailable:false,label:`ortholog alignment recorded (${d.nm} mismatches; seed ${d.seed_mismatches??'unknown'})`};
 const evidence=M.conservation_evidence||{},info=(evidence.per_species||{})[sp];
 if(evidence.attempted===false)return {unavailable:true,label:'orthology mapping intentionally skipped'};
 if(!info)return {unavailable:true,label:'orthology evidence not established'};
 if(info.state==='unavailable')return {unavailable:true,label:({no_orthologue:'no orthologue found',unrecognized_identifier:'target identifier not recognized',mapping_not_supplied:'orthology mapping not supplied',unsupported_target_species:'target species is unsupported by the mapping source',unsupported_query_species:'query species is unsupported by the mapping source',budget_exhausted:'orthology lookup budget was exhausted',transport_error:'orthology service transport failed',error_response:'orthology service returned an error response',http_error:'orthology service returned an HTTP error',request_error:'orthology request failed',invalid_response:'orthology service returned an invalid response',orthologues_found:'orthologues were found, but no completed alignment was published'}[info.reason]||'orthology unavailable')};
 if(!info.screening_complete||!g.off_target_screened)return {unavailable:true,label:'screening completion not established for this guide'};
 if(info.state==='resolved_empty')return {unavailable:false,label:'no ortholog mapping found'};
 if(info.state==='resolved')return {unavailable:false,label:'no conserved site found in completed screen'};
 return {unavailable:true,label:'orthology evidence not established'};
}
function conservationFilterDiagnostic(){
 if(!F.conserved.size||G.some(g=>[...F.conserved].every(sp=>{const d=(JSON.parse(g.ortholog||'{}'))[sp];return d&&d.nm!==null&&d.nm<=1&&(d.nm!==1||!F.conservedSeed||d.seed_mismatches===0); })))return '';
 const unavailable=[...F.conserved].map(sp=>{const n=G.filter(g=>conservationDiagnostic(g,sp).unavailable).length;return n?`${sp}: evidence unavailable or not established for ${n} guides`:'';}).filter(Boolean);
 return unavailable.length?'No guides match the selected conservation filters. '+unavailable.join('; ')+'. This does not establish biological absence.':'No guides match the selected conservation filters in the published evidence.';
}
function conservationCard(row){const evidence=M.conservation_evidence||{},ortholog=JSON.parse(row.ortholog||'{}');
 const entries=CONSERVATION_SPECIES.map(sp=>{const d=ortholog[sp];if(d)return `<li><b>${esc(sp)}</b>: ${d.nm} mismatches; seed ${esc(d.seed_mismatches??'unknown')}; transcript ${esc(d.transcript_id||'unknown')}; gene ${esc(d.gene_symbol||'unknown')}</li>`;
  return `<li><b>${esc(sp)}</b>: ${esc(conservationDiagnostic(row,sp).label)}</li>`;}).join('');
 return `<h3>Cross-species conservation</h3><p>Source: ${esc(evidence.source||'not established')} · release ${esc(evidence.release||'unknown')}.</p><ul>${entries||'<li>No conservation species were recorded.</li>'}</ul>`;
}
function restoreFilterState(){const toggle=document.querySelector('#geneenable'),select=document.querySelector('#targetgene'),params=new URLSearchParams((window.location?.hash||'').replace(/^#/,''));
 if(toggle&&select){toggle.checked=params.get('gene-enabled')==='1';select.value=params.get('gene')||'';select.disabled=!toggle.checked;toggle.setAttribute('aria-expanded',toggle.checked?'true':'false');}
 const tokens=(params.get('c')||'').split(',').filter(Boolean);F.conserved=new Set(tokens.filter(sp=>CONSERVATION_SPECIES.includes(sp)));F.conservedSeed=!tokens.includes('seed-relaxed');
 document.querySelectorAll('[data-conservation]').forEach(box=>{box.checked=F.conserved.has(box.value);});const seed=document.querySelector('#cons-seed');if(seed)seed.checked=F.conservedSeed;
}
function syncFilterState(){if(!window.location)return;const toggle=document.querySelector('#geneenable'),select=document.querySelector('#targetgene'),params=new URLSearchParams();
 if(toggle?.checked){params.set('gene-enabled','1');if(select?.value)params.set('gene',select.value);}
 const tokens=[...F.conserved].sort();if(!F.conservedSeed)tokens.push('seed-relaxed');if(tokens.length)params.set('c',tokens.join(','));window.location.hash=params.toString();}
function loadShard(i){return new Promise((resolve,reject)=>{if(window.SF_DETAIL_SHARDS?.[i])return resolve();const s=document.createElement('script');s.src=`${SHARD_DIR}/guide-${String(i).padStart(5,'0')}.js`;s.onload=resolve;s.onerror=()=>reject(Error('Could not load local detail shard'));document.head.appendChild(s);});}
function sourceAnchor(link,label){
 if(!link||!link.url)return esc(label);
 try{const url=new URL(link.url),host=url.hostname;
  if(url.protocol!=='https:'||url.username||url.password||url.port
   ||!(/^(?:(?:www|e[0-9]+|[a-z]{3}[0-9]{4}\.archive)\.)?ensembl\.org$/.test(host)
       ||/^(?:www\.)?mirgenedb\.org$/.test(host)))return esc(label);
 }catch(_e){return esc(label);}
 return `<a href="${esc(link.url)}" target="_blank" rel="noopener noreferrer" title="Open source record · release ${esc(link.release||'unknown')}">${esc(label)}</a>`;
}
function renderMirnaRows(rows){
 return rows.map(m=>{const name=m.mirna_id||m.rname||m.mirna||m.name||'unnamed';
  return `<tr><td>${sourceAnchor(m.source_link,name)}</td><td>${esc(m.database||'source unknown')}</td><td>${esc(m.coord??'—')}</td><td>${esc(m.seed_mismatches??'—')}</td></tr>`;
 }).join('');
}
function gateEvidenceCopy(row){
 const diagnostic=JSON.parse(row.report_evidence||'{}'), triples=JSON.parse(row.gates);
 const notices=(diagnostic.gates||[]).filter(e=>e.state==='unavailable'||e.state==='not_established').map(e=>{
  const gate=M.gate_panel.filters[e.filter_index], f=gate.descriptor,w=M.evidence_warning_registry[e.reason_id];
  if(!w)return '';const values={filter_id:f.filter_id,column:f.column,read_column:gate.read_column||f.column};
  const raw=(triples[e.filter_index][0]!=null&&w.body_if_measured)||w.body;
  const body=raw.replace(/\{(filter_id|column|read_column)\}/g,(_m,k)=>values[k]);
  return `<li data-warning-id="${esc(e.reason_id)}"><b>${esc(f.filter_id)}: ${esc(w.title)}</b><p>${esc(body)}</p><p>${esc(w.remediation)}</p></li>`;
 });
 if((diagnostic.warning_ids||[]).includes('UNDECLARED_RUN_REJECTION')){
  const w=M.evidence_warning_registry.UNDECLARED_RUN_REJECTION;
  notices.push(`<li data-warning-id="UNDECLARED_RUN_REJECTION"><b>${esc(w.title)}</b><p>${esc(w.body)}</p><p>${esc(w.remediation)}</p></li>`);
 }
 return notices.length?`<h3>Unresolved gate evidence</h3><ul>${notices.join('')}</ul>`:'';
}
function evidenceCopy(row){const warning=M.evidence_warning_registry[row.mirna_warning_id];return warning?`<div data-warning-id="${esc(row.mirna_warning_id)}"><h3>${esc(warning.title)}</h3><p>${esc(warning.body)}</p><p>${esc(warning.remediation)}</p></div>`:'<p>miRNA screening completed for this submitted guide; the retained rows describe matches under its declared reference and settings.</p>';}
function candidateMap(candidates,transcript,info={}){
 const rows=(candidates||[]).filter(c=>String(c.transcript_id||'')===transcript&&Number.isFinite(Number(c.position)));
 if(!rows.length)return '<p>No candidate positions are recorded for this transcript.</p>';
 const canonical=info.canonical==='canonical'?'Canonical':info.canonical==='other_isoform'?'Other isoform':'Unknown canonical status';
 const identity=[info.gene_symbol||info.gene_id||'Unknown gene',transcript,canonical].join(' · ');
 const W=640,H=120,L=32,R=12,T=12,B=28,positions=rows.map(c=>Number(c.position));
 const lo=positions.reduce((a,b)=>Math.min(a,b),Infinity),hi=positions.reduce((a,b)=>Math.max(a,b),-Infinity),span=hi-lo||1;
 const scores=rows.map(c=>Number(c.composite_score??c.design_score)).filter(Number.isFinite);
 const minScore=scores.length?scores.reduce((a,b)=>Math.min(a,b),Infinity):0,maxScore=scores.length?scores.reduce((a,b)=>Math.max(a,b),-Infinity):1,scoreSpan=maxScore-minScore||1;
 const dots=rows.map(c=>{const raw=Number(c.composite_score??c.design_score),x=L+(W-L-R)*(Number(c.position)-lo)/span,y=scores.length?H-B-(H-T-B)*(Number.isFinite(raw)?(raw-minScore)/scoreSpan:0.5):(T+H-B)/2;
  return `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="3"><title>Position ${esc(c.position)}; score ${Number.isFinite(raw)?esc(raw):'unavailable'}</title></circle>`;}).join('');
 return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Candidate positions: ${esc(identity)}"><title>Candidate positions: ${esc(identity)}</title><line x1="${L}" y1="${H-B}" x2="${W-R}" y2="${H-B}" stroke="#64748b"/><text x="${L}" y="${H-6}" font-size="10">${lo}</text><text x="${W-R}" y="${H-6}" text-anchor="end" font-size="10">${hi}</text><text x="${W/2}" y="${H-6}" text-anchor="middle" font-size="10">transcript position</text><g fill="#2563eb">${dots}</g></svg>`;
}
function candidateContext(candidates,isoforms){
 const transcripts=[...new Set((candidates||[]).map(c=>String(c.transcript_id||'')).filter(Boolean))];
 if(!transcripts.length)return '<h3>Candidate positions</h3><p>Candidate transcript identity is unavailable.</p>';
 const metadata=new Map((isoforms||[]).map(i=>[String(i.transcript_id||''),i]));
 const options=transcripts.map(tx=>{const i=metadata.get(tx)||{},canonical=i.canonical==='canonical'?'Canonical':i.canonical==='other_isoform'?'Other isoform':'Unknown';return `<option value="${esc(tx)}">${esc(i.gene_symbol||i.gene_id||'Unknown gene')} · ${esc(tx)} · ${canonical}</option>`;}).join('');
 const first=transcripts[0],info=metadata.get(first)||{};
 return `<section><h3>Candidate position map</h3><label for="candidate-transcript">Transcript</label><select id="candidate-transcript">${options}</select><p id="candidate-map-context">${esc([info.gene_symbol||info.gene_id||'Unknown gene',first,info.canonical==='canonical'?'Canonical':info.canonical==='other_isoform'?'Other isoform':'Unknown canonical status'].join(' · '))}</p><div id="candidate-map">${candidateMap(candidates,first,info)}</div></section>`;
}
async function showGuide(guide){const row=G.find(g=>g.guide===guide);if(!row)return;const panel=document.querySelector('#details');panel.textContent='Loading guide detail…';try{await loadShard(row.shard);const d=window.SF_DETAIL_SHARDS[row.shard][guide];
 const isoforms=(d.isoforms||[]).map(i=>{const identity=i.canonical==='canonical'?'Canonical':i.canonical==='other_isoform'?'Other isoform':'Unknown';
  const context=[i.gene_symbol||i.gene_id?esc(i.gene_symbol||'symbol unknown')+' · '+sourceAnchor(JSON.parse(i.gene_link||'{}'),i.gene_id||'gene ID unknown'):'Target gene unknown',
   sourceAnchor(JSON.parse(i.transcript_link||'{}'),i.transcript_id),esc(i.transcript_biotype||'biotype unknown')].join('<br>');
  return `<tr><td>${context}</td><td>${identity}<br>${esc(i.canonical_source||'authority unknown')} · release ${esc(i.canonical_release||'unknown')}<br>status: ${esc(i.canonical_status||'unknown')}</td><td>${i.candidate_rows}</td></tr>`;
 }).join('');
 const canonical=row.canonical_hit==='hit'?'hit':row.canonical_hit==='not_hit'?'not hit':'not established';
 panel.innerHTML=`<h2>${esc(guide)} · ${esc(row.status)}</h2>${gateEvidenceCopy(row)}${evidenceCopy(row)}<h3>Canonical transcript</h3><p>Guide is ${canonical} on the canonical transcript.</p><h3>Target transcript identity</h3>${isoforms?`<table><thead><tr><th>Gene · transcript</th><th>Canonical status and provenance</th><th>Candidate rows</th></tr></thead><tbody>${isoforms}</tbody></table>`:'<p>No transcript identity is established for this guide.</p>'}${candidateContext(d.candidates,d.isoforms)}<p>${row.candidate_rows} candidate rows · ${d.hits.length} classified hit rows · ${d.mirna.length} miRNA rows</p><h3>Named miRNA matches</h3>${d.mirna.length?`<table><thead><tr><th>miRNA</th><th>Database</th><th>Coordinate</th><th>Seed mismatches</th></tr></thead><tbody>${renderMirnaRows(d.mirna)}</tbody></table>`:'<p>No named miRNA matches were recorded.</p>'}<h3>Candidate rows</h3><pre>${esc(JSON.stringify(d.candidates,null,2))}</pre><h3>Off-target rows</h3><pre>${esc(JSON.stringify(d.hits,null,2))}</pre><h3>miRNA rows</h3><pre>${esc(JSON.stringify(d.mirna,null,2))}</pre>`;panel.innerHTML=panel.innerHTML.replace('<h3>Named miRNA matches</h3>',conservationCard(row)+'<h3>Named miRNA matches</h3>');const txSelect=document.querySelector('#candidate-transcript');if(txSelect)txSelect.addEventListener('change',()=>{const tx=txSelect.value,info=(d.isoforms||[]).find(i=>String(i.transcript_id||'')===tx)||{};document.querySelector('#candidate-map').innerHTML=candidateMap(d.candidates,tx,info);document.querySelector('#candidate-map-context').textContent=[info.gene_symbol||info.gene_id||'Unknown gene',tx,info.canonical==='canonical'?'Canonical':info.canonical==='other_isoform'?'Other isoform':'Unknown canonical status'].join(' · ');});}catch(e){panel.textContent=e.message;}}
document.querySelector('#counts').textContent='Gate results: '+Object.entries(STATUS_COUNTS).map(([k,v])=>`${k}: ${v}`).join(' · ')+` · total ${G.length} · ${M.summary.passes_without_unresolved_evidence} passes without unresolved evidence`;
document.querySelector('#gatecounts').textContent=(M.summary.filters||[]).map(f=>`${f.filter_id}: ${f.rejects} reject, ${f.warns} warn, ${f.unknowns} unknown`).join(' · ');
for(const f of M.summary.filters){const o=document.createElement('option');o.value=f.filter_id;o.textContent=f.filter_id;document.querySelector('#gate').appendChild(o);}
document.querySelector('#previous').addEventListener('click',()=>{currentPage=Math.max(0,currentPage-1);renderRows()});document.querySelector('#next').addEventListener('click',()=>{currentPage++;renderRows()});
document.querySelector('#unresolved-only').addEventListener('change',()=>{currentPage=0;renderRows()});
document.querySelector('#search').addEventListener('input',()=>{currentPage=0;renderRows()});document.querySelector('#status').addEventListener('change',()=>{currentPage=0;renderRows()});document.querySelector('#gate').addEventListener('change',()=>{currentPage=0;renderRows()});document.querySelector('#gate-result').addEventListener('change',()=>{currentPage=0;renderRows()});
const geneToggle=document.querySelector('#geneenable'),geneSelect=document.querySelector('#targetgene');restoreFilterState();
if(geneToggle&&geneSelect){geneToggle.addEventListener('change',()=>{geneSelect.disabled=!geneToggle.checked;geneToggle.setAttribute('aria-expanded',geneToggle.checked?'true':'false');syncFilterState();currentPage=0;renderRows();});geneSelect.addEventListener('change',()=>{syncFilterState();currentPage=0;renderRows();});}
document.querySelectorAll('[data-conservation]').forEach(box=>box.addEventListener('change',()=>{box.checked?F.conserved.add(box.value):F.conserved.delete(box.value);syncFilterState();currentPage=0;renderRows();}));
const conservationSeed=document.querySelector('#cons-seed');if(conservationSeed)conservationSeed.addEventListener('change',()=>{F.conservedSeed=conservationSeed.checked;syncFilterState();currentPage=0;renderRows();});
if(window.addEventListener)window.addEventListener('hashchange',()=>{restoreFilterState();currentPage=0;renderRows();});
renderRows();
</script></body></html>"""
