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
"""Key under which an artifact stamp records the version that produced it.

Manager-backed caches keep the same value in `CacheMetadata.version`; sidecar
stamps have no such document, so they carry this field. Naming it once keeps the
writer (`write_artifact_stamp`) and the reader (`is_artifact_stamp_current`) from
drifting apart. Only artifact classes that are actually validated somewhere are
registered below - a constant nobody checks would read as a guarantee we do not
give.
"""

# Derived-artifact classes have no cache subdirectory to name them, so they get
# constants; manager-backed classes are named by their cache subdir.
ARTIFACT_MIRNA_COMBINED = "mirna_combined"
ARTIFACT_TRANSCRIPTOME_INDEX = "transcriptome_index"

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
    # BWA-MEM2 indices. These are bound to the checksum of the FASTA they were
    # built from via an artifact stamp, so a version bump here is only needed if
    # the way we *build* indices changes; rebuilding costs CPU, not bandwidth.
    ARTIFACT_TRANSCRIPTOME_INDEX: LEGACY_PRODUCER_VERSION,
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
# stamp records three independent things, because each answers a question the
# others cannot:
#
#   producer version - "was this written by code we have since fixed?"
#   input fingerprint - "was this built from the content we hold now?"
#   output fingerprint - "are the artifact's own bytes still the bytes we wrote?"
#
# The output fingerprint is what stops a stamp vouching for an artifact that was
# truncated or corrupted *after* stamping (an interrupted copy, a full disk, a
# half-written rebuild): input-only fingerprints match happily in that case, and a
# stamped-but-broken artifact is then served for its whole TTL.

ARTIFACT_STAMP_SUFFIX = ".sirnaforge-cache.json"

_DIGEST_CHUNK_BYTES = 1024 * 1024

STAMP_FULL_DIGEST_MAX_BYTES = 64 * 1024 * 1024
"""Largest output we are willing to read end-to-end on every cache hit.

A combined miRNA FASTA is a few MB, so it is digested whole for a few
milliseconds. A BWA-MEM2 index or a transcriptome is gigabytes: reading it in full
on every hit would cost more than the cache saves, so those fall back to
`DIGEST_MODE_SAMPLED`.
"""

STAMP_SAMPLE_BYTES = 4 * 1024 * 1024
"""Head/tail window digested for outputs above the full-digest budget."""

DIGEST_MODE_FULL = "full"
DIGEST_MODE_SAMPLED = "sampled"


def artifact_stamp_path(artifact: Path) -> Path:
    """Return the sidecar stamp path for a derived artifact or index prefix."""
    return artifact.with_name(artifact.name + ARTIFACT_STAMP_SUFFIX)


def _digest_mode_for(size: int, full_digest_max_bytes: int) -> str:
    """Pick the digest strength that `size` bytes can afford on every cache hit."""
    return DIGEST_MODE_FULL if size <= full_digest_max_bytes else DIGEST_MODE_SAMPLED


def _digest_output(path: Path, size: int, mode: str) -> str:
    """Digest `path` using `mode`, which the stamp records so reads can repeat it."""
    hasher = hashlib.md5()
    with path.open("rb") as handle:
        if mode == DIGEST_MODE_FULL:
            for chunk in iter(lambda: handle.read(_DIGEST_CHUNK_BYTES), b""):
                hasher.update(chunk)
            return hasher.hexdigest()

        # Sampled: bind the digest to the length as well as the head and tail bytes,
        # so a same-length swap of head for tail cannot collide.
        hasher.update(str(size).encode("ascii"))
        hasher.update(handle.read(STAMP_SAMPLE_BYTES))
        handle.seek(max(size - STAMP_SAMPLE_BYTES, 0))
        hasher.update(handle.read(STAMP_SAMPLE_BYTES))
    return hasher.hexdigest()


def fingerprint_outputs(
    outputs: Sequence[Path], *, full_digest_max_bytes: int | None = None
) -> dict[str, dict[str, Any]]:
    """Fingerprint the bytes an artifact is actually made of.

    Keyed by basename so the fingerprint survives a re-homed cache directory, and
    so a multi-file artifact (a BWA index prefix has four members) records one
    entry per member and therefore proves its own completeness on read.

    Every member records its exact size; the digest is full below
    `STAMP_FULL_DIGEST_MAX_BYTES` and sampled above it. Size alone already catches
    every truncation, which is the failure mode that actually happens.
    """
    budget = STAMP_FULL_DIGEST_MAX_BYTES if full_digest_max_bytes is None else full_digest_max_bytes
    fingerprints: dict[str, dict[str, Any]] = {}
    for path in outputs:
        resolved = Path(path)
        try:
            size = resolved.stat().st_size
            mode = _digest_mode_for(size, budget)
            digest = _digest_output(resolved, size, mode)
        except OSError as exc:
            # Never abort a build over bookkeeping; an output we could not read is
            # simply left out, and a stamp with no outputs is not reusable.
            logger.warning("Could not fingerprint cache output %s: %s", resolved, exc)
            continue
        fingerprints[resolved.name] = {"size": size, "digest": digest, "digest_mode": mode}
    return fingerprints


def _stale_member_reason(member: Path, name: str, recorded: Any) -> str | None:
    """Why one recorded member no longer matches what is on disk, or None if it does.

    Ordered cheapest first: presence, then size (one `stat`, and the check that
    catches every truncation), then the digest in the mode the stamp recorded.
    """
    if not isinstance(recorded, Mapping):
        return f"its stamp holds an unreadable fingerprint for '{name}'"
    if not member.exists():
        return f"'{name}' is missing"

    size = member.stat().st_size
    if size != recorded.get("size"):
        return f"'{name}' is {size} bytes on disk, not the {recorded.get('size')} recorded (truncated or rewritten)"

    mode = recorded.get("digest_mode")
    if mode not in (DIGEST_MODE_FULL, DIGEST_MODE_SAMPLED):
        return f"its stamp records an unknown digest mode {mode!r} for '{name}'"
    if _digest_output(member, size, str(mode)) != recorded.get("digest"):
        return f"'{name}' no longer matches its recorded {mode} digest (corrupted after it was stamped)"

    return None


def _stale_output_reason(artifact: Path, outputs: Mapping[str, Any]) -> str | None:
    """Why an artifact's bytes no longer match its stamp, or None when they do."""
    if not outputs:
        return "its stamp records no fingerprint of the artifact's own bytes"

    for name, recorded in outputs.items():
        reason = _stale_member_reason(artifact.parent / str(name), str(name), recorded)
        if reason is not None:
            return reason

    return None


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
    outputs: Sequence[Path] | None = None,
    extra: Mapping[str, Any] | None = None,
) -> Path:
    """Record the producer version, input fingerprint and output bytes of an artifact.

    `outputs` lists the files the artifact consists of, and defaults to the artifact
    itself. Multi-file artifacts must pass their members explicitly: a BWA index
    prefix is not a file, and the four files behind it are what has to be verified.
    """
    stamp = {
        "artifact_class": artifact_class,
        PRODUCER_VERSION_FIELD: current_producer_version(artifact_class),
        # Basename only: the stamp lives beside the artifact, so recording an
        # absolute path would just invalidate the stamp whenever the cache moves.
        "artifact": artifact.name,
        "inputs": dict(inputs or {}),
        "outputs": fingerprint_outputs([artifact] if outputs is None else outputs),
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

    Answers "was it produced by the current code, from exactly these inputs, recently
    enough - and are its own bytes still the bytes we stamped". The last clause is
    what makes a stamp unable to vouch for an artifact that was truncated or
    corrupted after it was written, so callers may reuse on a True answer alone.
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

    # Everything above describes how the artifact was made; only this checks that
    # what is on disk now is still what was made. Without it a truncated artifact
    # keeps a "current" stamp and is served until the TTL expires.
    stale_output = _stale_output_reason(artifact, stamp.get("outputs") or {})
    if stale_output is not None:
        log_discarded_artifact(artifact_class, artifact, stale_output)
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
