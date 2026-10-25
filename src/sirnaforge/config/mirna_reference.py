"""Stage explicit immutable miRNA references and consume task-produced evidence."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

from sirnaforge.data.mirna_manager import load_pinned_references


def stage_reference_manifest(manifest: str | Path, destination: Path) -> Path:
    """Resolve local paths before Nextflow stages them in an isolated task directory."""
    references = load_pinned_references(manifest)
    destination.mkdir(parents=True, exist_ok=True)
    entries = []
    for source, species_pins in references.items():
        for species, pin in species_pins.items():
            target = destination / f"{source}_{species}_{pin.sha256}.fa"
            shutil.copyfile(pin.path, target)
            entries.append(
                {
                    "source": source,
                    "species": species,
                    "release": pin.release,
                    "sha256": pin.sha256,
                    "original_url": pin.original_url,
                    "path": str(target.resolve()),
                }
            )
    target_manifest = destination / "references.json"
    target_manifest.write_text(json.dumps({"schema_version": 1, "references": entries}, indent=2))
    return target_manifest.resolve()


def read_reference_evidence(
    run_dir: Path,
    *,
    source: str,
    species: list[str],
    skipped: bool,
    allow_artifact: bool = True,
    expected_digest: str | None = None,
) -> dict[str, Any]:
    """Read only the published task artifact, never infer reference use from a host cache."""
    artifact = run_dir / "off_target" / "results" / "mirna" / "mirna_reference_evidence.json"
    if not skipped and allow_artifact and artifact.is_file():
        try:
            payload = json.loads(artifact.read_text())
            if not isinstance(payload, dict) or not isinstance(payload.get("mirna"), dict):
                raise ValueError("Malformed reference evidence")
            mirna = payload["mirna"]
            if (
                type(payload.get("schema_version")) is int
                and payload.get("schema_version") == 1
                and mirna.get("source") == source
                and set(mirna.get("requested_species", [])) == set(species)
                and (expected_digest is None or mirna.get("screening", {}).get("guide_set_digest") == expected_digest)
            ):
                return payload
        except (OSError, ValueError, KeyError, TypeError):
            pass
    status = "skipped" if skipped else "unknown"
    reason = "screening_intentionally_skipped" if skipped else "producer_reference_evidence_unavailable"
    return {
        "schema_version": 1,
        "mirna": {
            "source": source,
            "requested_species": species,
            "sources": [],
            "acquisition": {"status": status, "reason": reason},
            "screening": {"status": status, "reason": reason, "guide_sequences": []},
        },
    }
