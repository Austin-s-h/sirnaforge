"""Experimental-status notice for the ZFN arm.

This lives in its own module rather than in ``sirnaforge.zfn.__init__`` because the
package ``__init__`` imports ``design``, and ``design`` has to emit the notice too --
reaching back into the package from one of its own submodules would be circular.

The notice is emitted at most once per process. Every ZFN entry point calls
:func:`emit_zfn_experimental_warning`, so whichever one a caller happens to reach first
carries the text, and the ones layered beneath it (CLI -> workflow -> designer) do not
repeat it. It is a statement about the module's status, not about a result, so repeating
it per candidate or per run would only train users to skip it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from sirnaforge.utils.logging_utils import get_logger

if TYPE_CHECKING:  # pragma: no cover - import only needed for the annotation
    from rich.console import Console

logger = get_logger(__name__)

ZFN_EXPERIMENTAL_ISSUE_URL = "https://github.com/Austin-s-h/sirnaforge/issues/82"

# The right half-site that actually matches the published CCR5 on-target locus, i.e. the
# reverse complement of the plus-strand text every ZFN doc and notebook still passes.
# Named so docs, warnings and tests cannot drift from each other.
ZFN_CCR5_WORKING_RIGHT_HALF_SITE = "CTTTTGCAGTTT"

# Deliberately enumerates the defects a user cannot notice from the output alone: a wrong
# number that looks plausible is more dangerous than a crash, and the two "silent" entries
# below (field inversion, orientation) both produce artifacts that read as successful runs.
ZFN_EXPERIMENTAL_WARNING = (
    "The ZFN module is EXPERIMENTAL and has known unfixed defects, tracked in "
    f"{ZFN_EXPERIMENTAL_ISSUE_URL}. Do not use ZFN results for decisions without "
    "independent validation.\n"
    "Wrong without saying so:\n"
    "  * Half-site orientation. With the default require_opposite_strands=True the right "
    "half-site must be supplied as the reverse complement of its published plus-strand "
    "text. The CCR5 pair used throughout the docs (GTCATCCTCATC / AAACTGCAAAAG) therefore "
    "matches nothing at all, including its own on-target site; pass "
    f"--zfn-right-half-site {ZFN_CCR5_WORKING_RIGHT_HALF_SITE} instead.\n"
    "  * FokI seed-region and polarity weighting are applied at the wrong end of the right "
    "half-site, so the strict seed budget guards the bases farthest from FokI rather than "
    "the nearest.\n"
    "  * Off-target region classification can report a site inside a large containing gene "
    "as 'intergenic', so the exonic/promoter tallies behind the pass/fail filters "
    "undercount and a pair can PASS filters it should fail.\n"
    "  * worst_site_score and best_offtarget_score are inverted in every exported "
    "artifact: worst_site_score holds the minimum and best_offtarget_score the maximum "
    "site score, but among off-targets the highest-scoring site is the most dangerous one.\n"
    "Fails loudly:\n"
    "  * The default pyahocorasick backend rejects --zfn-max-mismatches 4 on a 12 bp "
    "half-site (5,498,165 candidate patterns against a 1,000,000 safety limit) and raises "
    "ValueError. Use --zfn-search-backend exhaustive_python for mismatch budgets above 3."
)

_notice_emitted = False


def emit_zfn_experimental_warning(console: Console | None = None) -> bool:
    """Announce the ZFN arm's experimental status, at most once per process.

    Args:
        console: Rich console to render a highlighted notice on. Omit it for library
            callers -- the warning still reaches the log, and an application that has
            configured no handlers gets it on stderr via logging's last-resort handler.

    Returns:
        True if this call emitted the notice, False if an earlier call already did.
    """
    global _notice_emitted  # noqa: PLW0603 - process-wide once-per-run latch
    if _notice_emitted:
        return False
    _notice_emitted = True

    logger.warning(ZFN_EXPERIMENTAL_WARNING)
    if console is not None:
        # Imported lazily so importing sirnaforge.zfn does not pull in rich.
        from rich.panel import Panel  # noqa: PLC0415

        console.print(
            Panel.fit(
                f"⚠️  [bold yellow]EXPERIMENTAL: ZFN module[/bold yellow]\n{ZFN_EXPERIMENTAL_WARNING}",
                title="Experimental Feature",
                border_style="yellow",
            )
        )
    return True


def reset_zfn_experimental_warning() -> None:
    """Re-arm the once-per-process latch. For tests that exercise several ZFN entry points."""
    global _notice_emitted  # noqa: PLW0603 - see emit_zfn_experimental_warning
    _notice_emitted = False
