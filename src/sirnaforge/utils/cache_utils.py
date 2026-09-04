"""Shared cache helpers for siRNAforge.

Centralizes cache path resolution so every subsystem honors the same
SIRNAFORGE_CACHE_DIR / XDG cache layout, and the producer-version stamping that
lets us discard artifacts a known-bad writer left on disk.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile
from collections.abc import Mapping, Sequence
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

_CACHE_ROOT_SENTINEL = object()
_RESOLVED_CACHE_SUBDIRS: dict[str, Path] = {}
_RESOLUTION_MODES: dict[str, str] = {}

logger = logging.getLogger(__name__)


def resolve_cache_subdir(subdir: str, *, override: str | os.PathLike[str] | None = None) -> Path:
    """Resolve a writable cache directory for the requested subdir.

    The lookup order matches ReferenceManager: explicit override, env vars,
    XDG, $HOME/.cache, workspace-local fallback, then temp dir.
    """
    if not subdir:
        raise ValueError("subdir must be provided")

    resolution_mode = "override" if override is not None else "auto"

    candidates: list[Path] = []
    if override is not None:
        candidates.append(Path(override))
    else:
        env_override = os.getenv("SIRNAFORGE_CACHE_DIR")
        if env_override:
            candidates.append(Path(env_override) / subdir)

        xdg_cache = os.getenv("XDG_CACHE_HOME")
        if xdg_cache:
            candidates.append(Path(xdg_cache) / "sirnaforge" / subdir)

        candidates.append(Path.home() / ".cache" / "sirnaforge" / subdir)
        candidates.append(Path.cwd() / ".sirnaforge_cache" / subdir)
        candidates.append(Path(tempfile.gettempdir()) / "sirnaforge" / subdir)

    first_error: Exception | None = None
    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            prior = _RESOLVED_CACHE_SUBDIRS.get(subdir)
            prior_mode = _RESOLUTION_MODES.get(subdir)
            should_warn = (
                prior is not None and prior != candidate and prior_mode == "auto" and resolution_mode == "auto"
            )
            if should_warn:
                logger.warning(
                    "Multiple cache roots detected for subdir '%s': '%s' then '%s'. "
                    "This can cause repeated downloads across runs.",
                    subdir,
                    prior,
                    candidate,
                )
            _RESOLVED_CACHE_SUBDIRS[subdir] = candidate
            _RESOLUTION_MODES[subdir] = resolution_mode
            return candidate
        except PermissionError as exc:
            if first_error is None:
                first_error = exc
            continue

    if override is not None and first_error is not None:
        raise first_error
    raise RuntimeError(f"Cannot create cache directory for '{subdir}'")


def stable_cache_key(payload: Mapping[str, Any]) -> str:
    """Create a deterministic sha256 digest for arbitrary JSON-serializable data."""
    normalized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest()[:16]


# ---------------------------------------------------------------------------
# Producer versions
# ---------------------------------------------------------------------------
#
# A checksum recorded when a file is written can only prove the bytes did not rot
# afterwards; it can never prove the bytes were *right*, because a buggy writer
# stamps its own bad output. So each artifact class also records the version of
# the code that produced it, and that version is compared on load: bump the entry
# below in the same change that fixes a producer, and every artifact the old
# producer wrote becomes a cache miss and is regenerated exactly once.
#
# Classes are deliberately fine-grained. Invalidation is scoped to the artifact
# class whose producer actually moved, so a miRNA parsing fix (a few MB to
# re-download) never forces a user to re-fetch a multi-GB transcriptome, and an
# index rebuild never forces a re-download of the FASTA it was built from.

LEGACY_PRODUCER_VERSION = "1.0"
"""Version attributed to artifacts written before producer stamping existed.

`CacheMetadata.version` has always defaulted to this string, so pre-existing
metadata reads back as "produced by the pre-versioning writer" without needing a
migration step.
"""

PRODUCER_VERSION_FIELD = "producer_version"
"""Field/column name for row-oriented caches (e.g. variants.parquet).

File-backed caches record the producer version in `CacheMetadata.version`;
row-backed caches have no such document and should carry this column instead so
both kinds validate through `is_producer_version_current`.
"""

# Derived-artifact classes have no cache subdirectory to name them, so they get
# constants; manager-backed classes are named by their cache subdir.
ARTIFACT_MIRNA_COMBINED = "mirna_combined"
ARTIFACT_TRANSCRIPTOME_FILTERED = "transcriptome_filtered"
ARTIFACT_TRANSCRIPTOME_INDEX = "transcriptome_index"
ARTIFACT_VARIANT_ROWS = "variants"

PRODUCER_VERSIONS: dict[str, str] = {
    # Per-species filtered miRNA FASTAs. Bumped because the pre-2.0 writer parsed
    # miRBase with Biopython's "fasta-blast" reader, which copies the first
    # record's id/description onto every record, so cached "species-filtered"
    # FASTAs can hold the wrong sequences entirely.
    "mirna": "2.0",
    # combined_*.fa, derived from the above and therefore contaminated by the
    # same parsing defect.
    ARTIFACT_MIRNA_COMBINED: "2.0",
    # Raw transcriptome downloads are upstream bytes: no producer of ours shapes
    # their content, so they stay at the legacy version and are never discarded
    # by a version bump. Re-downloading tens of GB must be a deliberate act.
    "transcriptomes": LEGACY_PRODUCER_VERSION,
    # Genome/annotation downloads are raw upstream bytes for the same reason.
    "genomes": LEGACY_PRODUCER_VERSION,
    "annotations": LEGACY_PRODUCER_VERSION,
    # Filtered transcriptome FASTAs are derived, but the filter output itself has
    # not changed; legacy entries are repaired by metadata bookkeeping rather
    # than by re-filtering.
    ARTIFACT_TRANSCRIPTOME_FILTERED: LEGACY_PRODUCER_VERSION,
    # BWA-MEM2 indices. These are bound to the checksum of the FASTA they were
    # built from via an artifact stamp, so a version bump here is only needed if
    # the way we *build* indices changes; rebuilding costs CPU, not bandwidth.
    ARTIFACT_TRANSCRIPTOME_INDEX: LEGACY_PRODUCER_VERSION,
    # Rows in variants.parquet. Bumped because rows written by the pre-2.0 writer
    # could lose population allele frequencies and could be overwritten by an
    # unrelated key, so their contents cannot be trusted.
    ARTIFACT_VARIANT_ROWS: "2.0",
}


def current_producer_version(artifact_class: str) -> str:
    """Return the producer version to record for `artifact_class`.

    Unregistered classes report the legacy version: registering a class is what
    opts it in, so a new cache namespace can never invalidate itself by accident.
    """
    return PRODUCER_VERSIONS.get(artifact_class, LEGACY_PRODUCER_VERSION)


def normalize_producer_version(recorded: object) -> str:
    """Coerce a recorded producer version into a comparable string.

    Missing values (None, NaN from a parquet column that did not exist yet, an
    empty string) all mean the same thing: the artifact predates stamping.
    """
    if isinstance(recorded, str) and recorded:
        return recorded
    return LEGACY_PRODUCER_VERSION


def is_producer_version_current(artifact_class: str, recorded: object) -> bool:
    """Whether an artifact's recorded producer version is the current one."""
    return normalize_producer_version(recorded) == current_producer_version(artifact_class)


def log_discarded_artifact(artifact_class: str, artifact: Path | str, reason: str) -> None:
    """Warn that a stale cache artifact is being dropped.

    Deliberately WARNING, not DEBUG: this is the one line that explains why a run
    suddenly spends time regenerating or re-downloading a reference.
    """
    logger.warning(
        "♻️  Discarding stale cached %s artifact '%s' (%s). It will be regenerated automatically "
        "once, which may mean re-downloading the reference; run `sirnaforge cache --clear` if you "
        "would rather reclaim the space yourself.",
        artifact_class,
        artifact,
        reason,
    )


# ---------------------------------------------------------------------------
# Artifact stamps for derived files
# ---------------------------------------------------------------------------
#
# Derived artifacts (a combined FASTA, a BWA-MEM2 index) have no entry in a
# manager's metadata document, so they carry a sidecar JSON stamp instead. The
# stamp records the producer version *and* a fingerprint of the inputs, which is
# what makes a derived artifact unable to outlive the content it was built from.

ARTIFACT_STAMP_SUFFIX = ".sirnaforge-cache.json"


def artifact_stamp_path(artifact: Path) -> Path:
    """Return the sidecar stamp path for a derived artifact or index prefix."""
    return artifact.with_name(artifact.name + ARTIFACT_STAMP_SUFFIX)


def fingerprint_inputs(inputs: Sequence[Path]) -> dict[str, str]:
    """MD5 each input file, keyed by position and basename.

    Position keeps duplicate basenames distinct; basenames rather than full paths
    keep the fingerprint stable when a cache directory is re-homed (containers
    mount it elsewhere). Callers that already hold a digest for a large file
    should pass it directly rather than paying for a second full read.
    """
    digests: dict[str, str] = {}
    for position, path in enumerate(inputs):
        hash_md5 = hashlib.md5()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                hash_md5.update(chunk)
        digests[f"{position}:{Path(path).name}"] = hash_md5.hexdigest()
    return digests


def write_artifact_stamp(
    artifact_class: str,
    artifact: Path,
    *,
    inputs: Mapping[str, str] | None = None,
    extra: Mapping[str, Any] | None = None,
) -> Path:
    """Record the producer version and input fingerprint for a derived artifact."""
    stamp = {
        "artifact_class": artifact_class,
        PRODUCER_VERSION_FIELD: current_producer_version(artifact_class),
        # Basename only: the stamp lives beside the artifact, so recording an
        # absolute path would just invalidate the stamp whenever the cache moves.
        "artifact": artifact.name,
        "inputs": dict(inputs or {}),
        "stamped_at": datetime.now().isoformat(),
        "extra": dict(extra or {}),
    }
    stamp_path = artifact_stamp_path(artifact)
    stamp_path.parent.mkdir(parents=True, exist_ok=True)
    with stamp_path.open("w", encoding="utf-8") as handle:
        json.dump(stamp, handle, indent=2)
    return stamp_path


def read_artifact_stamp(artifact: Path) -> dict[str, Any] | None:
    """Read a derived artifact's stamp, or None when it is absent/unreadable."""
    stamp_path = artifact_stamp_path(artifact)
    if not stamp_path.exists():
        return None
    try:
        with stamp_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (json.JSONDecodeError, OSError) as exc:
        logger.warning("Unreadable cache stamp %s: %s", stamp_path, exc)
        return None
    return data if isinstance(data, dict) else None


def discard_artifact_stamp(artifact: Path) -> None:
    """Remove a derived artifact's stamp so it cannot vouch for a rebuilt file."""
    artifact_stamp_path(artifact).unlink(missing_ok=True)


def is_artifact_stamp_current(  # noqa: PLR0911
    artifact_class: str,
    artifact: Path,
    *,
    inputs: Mapping[str, str] | None = None,
    max_age_days: int | None = None,
) -> bool:
    """Whether a derived artifact may be reused.

    Callers must separately confirm the artifact itself is present and complete
    (a BWA index is several files); this answers only "was it produced by the
    current code from exactly these inputs, recently enough".
    """
    stamp = read_artifact_stamp(artifact)
    if stamp is None:
        # No stamp and no artifact is an ordinary first run, not a discard.
        if artifact.exists():
            log_discarded_artifact(artifact_class, artifact, "written before cache stamping existed")
        return False

    if stamp.get("artifact_class") != artifact_class:
        log_discarded_artifact(
            artifact_class, artifact, f"stamped as '{stamp.get('artifact_class')}' rather than '{artifact_class}'"
        )
        return False

    recorded_version = stamp.get(PRODUCER_VERSION_FIELD)
    if not is_producer_version_current(artifact_class, recorded_version):
        log_discarded_artifact(
            artifact_class,
            artifact,
            f"producer version {normalize_producer_version(recorded_version)} != "
            f"{current_producer_version(artifact_class)}",
        )
        return False

    if inputs is not None and dict(stamp.get("inputs") or {}) != dict(inputs):
        log_discarded_artifact(artifact_class, artifact, "the inputs it was built from have changed")
        return False

    if max_age_days is not None:
        stamped_at = stamp.get("stamped_at")
        try:
            age = datetime.now() - datetime.fromisoformat(str(stamped_at))
        except (TypeError, ValueError):
            log_discarded_artifact(artifact_class, artifact, f"unparseable stamp timestamp {stamped_at!r}")
            return False
        if age > timedelta(days=max_age_days):
            logger.info(
                "Cached %s artifact '%s' is %d days old (TTL %d); regenerating",
                artifact_class,
                artifact,
                age.days,
                max_age_days,
            )
            return False

    return True
