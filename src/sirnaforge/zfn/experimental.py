"""Experimental-status notice for the ZFN arm.

This lives in its own module rather than in ``sirnaforge.zfn.__init__`` because the
package ``__init__`` imports ``design``, and ``design`` has to emit the notice too --
reaching back into the package from one of its own submodules would be circular.

The notice is emitted at most once per process. Every ZFN entry point calls
:func:`emit_zfn_experimental_warning`, so whichever one a caller happens to reach first
carries the text, and the ones layered beneath it (CLI -> workflow -> designer) do not
repeat it. It is a statement about the module's status, not about a result, so repeating
it per candidate or per run would only train users to skip it.

"Once" means once *as the user sees it*, not once per log record. When a rich console is
supplied, the log record and the panel would otherwise both land on the same terminal
stream and print the whole notice twice -- ``sirnaforge.utils.logging_utils`` attaches a
``StreamHandler(sys.stdout)`` to the root logger on first ``get_logger`` call, and rich
writes to ``sys.stdout`` too. :func:`emit_zfn_experimental_warning` therefore keeps the
log record (so log files, ``caplog`` and library callers still get it) but mutes, for that
one record only, the handlers writing to the console's own stream.

Where the notice lands with no console: wherever the host application's logging sends
``WARNING``. Under sirnaforge's own :func:`~sirnaforge.utils.logging_utils.get_logger`
that is **stdout**, not stderr -- because a handler is always installed, logging's
last-resort stderr handler never fires.
"""

from __future__ import annotations

import logging
from contextlib import contextmanager
from typing import TYPE_CHECKING

from sirnaforge.utils.logging_utils import get_logger

if TYPE_CHECKING:  # pragma: no cover - imports only needed for annotations
    from collections.abc import Iterator

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


class _DropEverything(logging.Filter):
    """Filter that rejects every record, used to mute one handler for one call."""

    def filter(self, record: logging.LogRecord) -> bool:  # noqa: ARG002, D102 - signature is logging's
        return False


_DROP_EVERYTHING = _DropEverything()


def _effective_handlers(log: logging.Logger) -> Iterator[logging.Handler]:
    """Yield the handlers a record logged on ``log`` would actually reach."""
    current: logging.Logger | None = log
    while current is not None:
        yield from current.handlers
        current = current.parent if current.propagate else None


@contextmanager
def _console_log_output_muted(console: Console) -> Iterator[None]:
    """Silence log handlers writing to ``console``'s stream for the duration of the block.

    The rich panel already shows the notice on that stream; without this the log record
    prints the same text again just above the panel. Handlers pointing anywhere else -- a
    log file, ``caplog``'s capture handler -- are untouched, so the record is not lost.
    """
    try:
        target = console.file
    except Exception:  # pragma: no cover - defensive: a stub console without .file
        target = None

    muted: list[logging.Handler] = []
    if target is not None:
        for handler in _effective_handlers(logger):
            if getattr(handler, "stream", None) is target:
                handler.addFilter(_DROP_EVERYTHING)
                muted.append(handler)
    try:
        yield
    finally:
        for handler in muted:
            handler.removeFilter(_DROP_EVERYTHING)


def emit_zfn_experimental_warning(console: Console | None = None) -> bool:
    """Announce the ZFN arm's experimental status, at most once per process.

    The notice is logged at ``WARNING`` on every path, so it reaches log files and any
    handler the host application configured. It is rendered to the terminal exactly once:
    with a ``console`` the user sees the rich panel, without one they see whatever the log
    handlers show (stdout, under sirnaforge's own logging setup).

    Args:
        console: Rich console to render a highlighted notice on. Omit it for library
            callers -- the warning still reaches the log.

    Returns:
        True if this call emitted the notice, False if an earlier call already did.
    """
    global _notice_emitted  # noqa: PLW0603 - process-wide once-per-run latch
    if _notice_emitted:
        return False
    _notice_emitted = True

    if console is None:
        logger.warning(ZFN_EXPERIMENTAL_WARNING)
        return True

    # Imported lazily so importing sirnaforge.zfn does not pull in rich.
    from rich.panel import Panel  # noqa: PLC0415

    with _console_log_output_muted(console):
        logger.warning(ZFN_EXPERIMENTAL_WARNING)
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
