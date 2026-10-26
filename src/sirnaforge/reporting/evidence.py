"""Offline evidence diagnostics shared by report payloads and rendering.

Completion is a producer claim about submitted guides, never an inference from hit counts.
The warning identifiers are stable; display text can evolve without changing exported reasons.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

WARNING_REGISTRY: dict[str, dict[str, Any]] = {
    "MIRNA_ACQUISITION_FAILED": {
        "severity": "warning",
        "title": "miRNA reference unavailable",
        "body": "Reference acquisition failed. Missing matches are not evidence of zero miRNA seed matches.",
        "remediation": "Restore the reference or supply a verified local reference, then repeat screening.",
        "default_visible": True,
    },
    "MIRNA_SCREENING_FAILED": {
        "severity": "warning",
        "title": "miRNA screening failed",
        "body": "Screening did not complete. Any published matches are partial evidence.",
        "remediation": "Inspect the recorded failure and repeat miRNA screening.",
        "default_visible": True,
    },
    "MIRNA_SCREENING_SKIPPED": {
        "severity": "info",
        "title": "miRNA screening skipped",
        "body": "This run did not perform the requested miRNA screening.",
        "remediation": "Run miRNA screening before interpreting missing matches as a negative result.",
        "default_visible": True,
    },
    "MIRNA_SCREENING_NOT_REQUESTED": {
        "severity": "info",
        "title": "miRNA screening not requested",
        "body": "miRNA screening was not requested for this run.",
        "remediation": "Enable miRNA screening if this evidence is needed.",
        "default_visible": True,
    },
    "MIRNA_EVIDENCE_UNKNOWN": {
        "severity": "warning",
        "title": "miRNA screening not established",
        "body": "This run does not establish completed miRNA screening for this guide. Missing matches are not a zero-hit result.",
        "remediation": "Repeat screening with provenance enabled to establish reference and guide coverage.",
        "default_visible": True,
    },
}


EVIDENCE_STATES = ("pass", "warn", "fail", "not_evaluated", "unavailable", "not_established")
# Report labels only; comparisons continue to use the run's descriptors and observations.
FILTER_UNITS = {
    "gc_min": "% GC",
    "gc_max": "% GC",
    "max_poly_runs": "nt",
    "max_repeat_transcript_fraction": "fraction (0–1)",
    "max_paired_fraction": "fraction (0–1)",
    "min_asymmetry_score": "score",
    "min_empirical_score": "score",
    "min_isoform_coverage": "fraction (0–1)",
    "max_off_target_count": "hits",
    "max_transcriptome_hits_0mm": "hits",
    "max_transcriptome_hits_1mm": "hits",
    "max_transcriptome_hits_2mm": "hits",
    "max_transcriptome_seed_perfect": "hits",
    "max_transcript_seed_sites": "sites",
    "max_mirna_perfect_seed": "hits",
    "max_mirna_1mm_seed": "hits",
    "fail_on_high_risk_mirna": "hits",
    "max_total_offtarget_hits": "hits",
}
GATE_REASON_IDS = {
    1: "GATE_COLUMN_UNAVAILABLE",
    2: "GATE_VALUE_UNAVAILABLE",
    3: "GATE_FILTER_OFF",
    4: "GATE_THRESHOLD_UNDECLARED",
    5: "GATE_RUN_UNKNOWN",
    6: "GATE_RUN_NOT_EVALUATED",
    7: "GATE_RUN_UNAVAILABLE",
    8: "GATE_RUN_NOT_ESTABLISHED",
    9: "GATE_RECORDED_CONFLICT",
}
_GATE_WARNINGS = {
    "GATE_COLUMN_UNAVAILABLE": (
        "Gate observation unavailable",
        "the run exports neither {filter_id}_observed nor {column}",
    ),
    "GATE_VALUE_UNAVAILABLE": ("Gate value unavailable", "{read_column} is empty"),
    "GATE_FILTER_OFF": ("Gate not evaluated", "filter is off, and the run exports no value for {column}"),
    "GATE_THRESHOLD_UNDECLARED": ("Gate not configured", "no threshold declared, and no value exported"),
    "GATE_RUN_UNKNOWN": (
        "Gate evidence not established",
        "the run recorded unknown for this gate; no value to re-compare",
    ),
    "GATE_RUN_NOT_EVALUATED": ("Active gate not evaluated", "the run did not evaluate this gate"),
    "GATE_RUN_UNAVAILABLE": (
        "Gate evidence unavailable",
        "the run recorded unavailable evidence for this gate; no value to re-compare",
    ),
    "GATE_OBSERVATIONS_UNDECIDED": (
        "Gate observations undecided",
        "{read_column} is exported, but the run evaluated this gate for none of the {n_guides} guides in this run, so no threshold can decide it",
    ),
    "GATE_RUN_NOT_ESTABLISHED": (
        "Gate evidence not established",
        "the run did not establish this gate's evidence; no value to re-compare",
    ),
    "GATE_RECORDED_CONFLICT": (
        "Gate outcome conflicts with observations",
        "The recorded gate outcome disagrees with its exported observation and policy. The report cannot reconstruct a consistent decision.",
    ),
    "UNDECLARED_RUN_REJECTION": (
        "Run rejection cannot be reconstructed",
        "The run rejected this guide for a reason outside the report's declared gate observations.",
    ),
}
for _warning_id, (_title, _body) in _GATE_WARNINGS.items():
    _inactive = _warning_id in {"GATE_FILTER_OFF", "GATE_THRESHOLD_UNDECLARED"}
    WARNING_REGISTRY[_warning_id] = {
        "severity": "info" if _inactive else "warning",
        "title": _title,
        "body": _body,
        "remediation": "Inspect the run policy and evidence export; repeat the required analysis before treating this gate as satisfied."
        if not _inactive
        else "Inspect the recorded observation and policy if this gate should be enabled.",
        "default_visible": not _inactive,
    }
WARNING_REGISTRY["GATE_COLUMN_UNAVAILABLE"]["summary_body"] = (
    "One or more active gates have no exported observation column."
)
WARNING_REGISTRY["GATE_VALUE_UNAVAILABLE"]["summary_body"] = (
    "One or more active gates have no recorded observation value."
)
WARNING_REGISTRY["GATE_FILTER_OFF"]["body_if_measured"] = "filter is off; value measured, nothing acts on it"
WARNING_REGISTRY["GATE_THRESHOLD_UNDECLARED"]["body_if_measured"] = (
    "no threshold declared; value measured, nothing acts on it"
)
WARNING_REGISTRY["GATE_RUN_NOT_EVALUATED"]["body_if_measured"] = (
    "the run did not evaluate this gate; value measured, no verdict applied"
)


def normalize_evidence_state(value: object) -> str:
    """Normalize report-only aliases; an unrecognized producer value cannot become a pass."""
    key = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    if key == "unknown":
        return "not_established"
    return key if key in EVIDENCE_STATES else "not_established"


def report_evidence(
    gates: Sequence[Sequence[Any]],
    mirna: Mapping[str, Any],
    *,
    undeclared_rejection: bool = False,
    recorded_states: Sequence[str | None] = (),
) -> dict[str, Any]:
    """Describe missing evidence separately from recorded gate outcomes and threshold findings."""
    entries = []
    warning_ids: set[str] = set()
    for index, (_value, verdict, reason) in enumerate(gates):
        reason_id = GATE_REASON_IDS.get(reason)
        state = {0: "pass", 1: "fail", 3: "not_evaluated", 4: "warn"}.get(verdict)
        if state is None:
            state = "unavailable" if reason in {1, 2, 7} else "not_established"
            if reason_id:
                warning_ids.add(reason_id)
        entry = {"filter_index": index, "state": state, "reason_id": reason_id}
        if index < len(recorded_states) and recorded_states[index] is not None:
            entry["recorded_state"] = normalize_evidence_state(recorded_states[index])
            entry["recorded_state_raw"] = recorded_states[index]
        entries.append(entry)
    if undeclared_rejection:
        warning_ids.add("UNDECLARED_RUN_REJECTION")
    mirna_state = mirna.get("state", "unknown")
    mirna_warning = mirna.get("warning_id")
    if isinstance(mirna_warning, str):
        warning_ids.add(mirna_warning)
    unresolved = (
        undeclared_rejection
        or any(entry["state"] in {"unavailable", "not_established"} for entry in entries)
        or mirna_state not in {"complete", "not_requested"}
    )
    return {"schema_version": 1, "gates": entries, "warning_ids": sorted(warning_ids), "unresolved": bool(unresolved)}


def _mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _phase_status(value: object) -> str:
    status = _mapping(value).get("status")
    if isinstance(status, str):
        status = status.strip().lower().replace("-", "_").replace(" ", "_")
    status = (
        {"unavailable": "failed", "fail": "failed", "not_established": "unknown", "not_evaluated": "skipped"}.get(
            status, status
        )
        if isinstance(status, str)
        else status
    )
    return (
        status
        if isinstance(status, str) and status in {"complete", "failed", "skipped", "not_requested", "unknown"}
        else "unknown"
    )


def _phase_metadata(value: object) -> dict[str, Any]:
    """Keep stable diagnostics and corpus identity; free-form producer details stay in the manifest."""
    data = _mapping(value)
    return {k: data[k] for k in ("status", "reason", "guide_set_digest", "submitted_guides") if k in data}


def _source_state(source: Mapping[str, Any]) -> str:
    acquisition = _phase_status(source.get("acquisition"))
    screening = _phase_status(source.get("screening"))
    if acquisition == "failed":
        return "acquisition_failed"
    if screening == "failed":
        return "screening_failed"
    if screening == "not_requested":
        return "not_requested"
    if screening == "skipped" or acquisition == "skipped":
        return "skipped"
    return "complete" if screening == "complete" and acquisition == "complete" else "unknown"


def _reference_provenance(value: object) -> list[dict[str, Any]]:
    """Flatten combined provenance, keeping identities/checksums rather than external URLs."""
    entries = value if isinstance(value, list) else [value]
    result = []
    for item in entries:
        data = _mapping(item)
        if not data:
            continue
        fields = {
            k: data[k]
            for k in (
                "source",
                "release",
                "species",
                "acquisition_mode",
                "raw_sha256",
                "filtered_sha256",
                "acquired_at",
            )
            if isinstance(data.get(k), str)
        }
        if fields:
            result.append(fields)
        if isinstance(data.get("sources"), list):
            result.extend(_reference_provenance(data["sources"]))
    return result


def mirna_evidence(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Read the versioned producer contract without trusting legacy empty hit files."""
    contract = _mapping(manifest.get("reference_evidence"))
    version = contract.get("schema_version")
    data = _mapping(contract.get("mirna")) if type(version) is int and version == 1 else {}
    sources = data.get("sources")
    sources = [dict(s) for s in sources if isinstance(s, Mapping)] if isinstance(sources, list) else []
    state = _source_state(data)
    if sources:
        states = [state, *(_source_state(s) for s in sources)]
        # An incomplete requested source prevents a complete aggregate even if another source ran.
        state = next(
            (
                s
                for s in ("acquisition_failed", "screening_failed", "skipped", "not_requested", "unknown")
                if s in states
            ),
            state,
        )
    screening = _mapping(data.get("screening"))
    submitted = screening.get("guide_sequences")
    return {
        "schema_version": 1,
        "state": state,
        "source": data.get("source"),
        "release": data.get("release"),
        "requested_species": data.get("requested_species") if isinstance(data.get("requested_species"), list) else [],
        "acquisition": _phase_metadata(data.get("acquisition")),
        "screening": _phase_metadata(screening),
        "guide_sequences": frozenset(str(g).strip().upper().replace("T", "U") for g in submitted)
        if isinstance(submitted, list)
        else frozenset(),
        "sources": [
            {
                "source": s.get("source"),
                "species": s.get("species"),
                "release": s.get("release"),
                "acquisition": _phase_metadata(s.get("acquisition")),
                "screening": _phase_metadata(s.get("screening")),
                "references": _reference_provenance(s.get("references")),
                "state": _source_state(s),
            }
            for s in sources
        ],
    }


def guide_mirna_evidence(evidence: Mapping[str, Any], guide: str) -> dict[str, Any]:
    """Return a guide's diagnostic, checking explicit producer submission coverage."""
    state = str(evidence["state"])
    if state == "complete" and guide not in evidence["guide_sequences"]:
        state = "unknown"
    warning_id = {
        "acquisition_failed": "MIRNA_ACQUISITION_FAILED",
        "screening_failed": "MIRNA_SCREENING_FAILED",
        "skipped": "MIRNA_SCREENING_SKIPPED",
        "not_requested": "MIRNA_SCREENING_NOT_REQUESTED",
        "unknown": "MIRNA_EVIDENCE_UNKNOWN",
    }.get(state)
    return {"state": state, "warning_id": warning_id}


def conservation_evidence(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Read mapping and screening coverage separately from the existing workflow provenance."""
    provenance = _mapping(manifest.get("provenance"))
    databases = _mapping(provenance.get("databases"))
    orthology = _mapping(databases.get("orthology"))
    scope = _mapping(orthology.get("conservation_scope"))
    unscreened = _mapping(scope.get("unscreened_in_denominator")).get("value")
    coverage_known = isinstance(unscreened, list)
    missing_species = unscreened if isinstance(unscreened, list) else []
    records = orthology.get("per_species")
    records = [s for s in records if isinstance(s, Mapping)] if isinstance(records, list) else []
    per_species = {}
    for record in records:
        species = record.get("species")
        if not isinstance(species, str):
            continue
        status = record.get("status")
        state = {"resolved": "resolved", "resolved_empty": "resolved_empty", "unresolved": "unavailable"}.get(
            status if isinstance(status, str) else "", "unknown"
        )
        per_species[species] = {
            "state": state,
            "reason": record.get("reason") if isinstance(record.get("reason"), str) else None,
            "screening_complete": coverage_known
            and species in (scope.get("denominator_species") or [])
            and species not in missing_species,
            "gene_ids": record.get("gene_ids") if isinstance(record.get("gene_ids"), list) else [],
        }
    return {
        "source": orthology.get("source"),
        "release": _mapping(orthology.get("compara_release")).get("value"),
        "attempted": orthology.get("attempted"),
        "reason": orthology.get("reason"),
        "per_species": per_species,
        "denominator_species": scope.get("denominator_species") or [],
    }
