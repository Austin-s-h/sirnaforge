"""Offline evidence diagnostics shared by report payloads and rendering.

Completion is a producer claim about submitted guides, never an inference from hit counts.
The warning identifiers are stable; display text can evolve without changing exported reasons.
"""

from __future__ import annotations

from collections.abc import Mapping
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


def _mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _phase_status(value: object) -> str:
    status = _mapping(value).get("status")
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
    return "complete" if screening == "complete" else "unknown"


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
    data = _mapping(contract.get("mirna")) if contract.get("schema_version") == 1 else {}
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
