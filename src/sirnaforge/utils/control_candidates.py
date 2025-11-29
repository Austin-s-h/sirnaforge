"""Helpers for injecting control candidates into design results."""

from __future__ import annotations

from sirnaforge.models.sirna import DesignResult, SiRNACandidate

DIRTY_CONTROL_LABEL = "DIRTY_CONTROL"
DIRTY_CONTROL_SUFFIX = "__DIRTY_CONTROL"


def inject_dirty_controls(design_result: DesignResult, count: int = 2) -> list[SiRNACandidate]:
    """Append clearly labelled "dirty" control candidates to a design result.

    "Dirty" controls are *not* synthetic sequences. We intentionally reuse
    the lowest-scoring rejected candidates from the same run so every control
    is a plausible design that simply failed QC (GC range, off-target score,
    etc.). Carrying these borderline guides forward guarantees that the
    downstream off-target pipeline receives at least a couple of sequences
    that should light up in the reports, making it easier for users to
    confirm the workflow is wired correctly.

    The clones retain their original guide/passenger sequences and metadata
    but are marked with :attr:`~sirnaforge.models.sirna.SiRNACandidate.FilterStatus.DIRTY_CONTROL`.
    They are appended to both ``candidates`` and ``top_candidates`` so that
    subsequent stages (CSV export, Nextflow off-target analysis, etc.) see
    the controls exactly where they expect normal guides. IDs are suffixed
    with ``__DIRTY_CONTROL_<n>`` for clarity.

    Args:
        design_result: Aggregated design result to augment in-place.
        count: Number of controls to add (default: 2).

    Returns:
        List of control candidates that were appended. Returns an empty list
        if the design result already contains controls or there are no
        candidates to copy.
    """
    if count <= 0:
        return []

    # Avoid duplicating controls when the workflow is re-run on the same object.
    if any(DIRTY_CONTROL_SUFFIX in candidate.id for candidate in design_result.candidates):
        return []

    rejected_pool = [c for c in getattr(design_result, "rejected_candidates", []) if DIRTY_CONTROL_SUFFIX not in c.id]
    if not rejected_pool:
        return []

    worst_candidates = sorted(rejected_pool, key=lambda c: getattr(c, "composite_score", 0.0))[:count]
    if not worst_candidates:
        return []
    dirty_controls: list[SiRNACandidate] = []

    for idx, base_candidate in enumerate(worst_candidates, start=1):
        clone = base_candidate.model_copy(deep=True)
        clone.id = f"{base_candidate.id}{DIRTY_CONTROL_SUFFIX}_{idx}"
        clone.passes_filters = SiRNACandidate.FilterStatus.DIRTY_CONTROL
        issues = list(getattr(base_candidate, "quality_issues", []) or [])
        if DIRTY_CONTROL_LABEL not in issues:
            issues.append(DIRTY_CONTROL_LABEL)
        clone.quality_issues = issues
        dirty_controls.append(clone)

    if not dirty_controls:
        return []

    design_result.candidates.extend(dirty_controls)
    design_result.top_candidates.extend(dirty_controls)
    design_result.total_candidates += len(dirty_controls)

    return dirty_controls
