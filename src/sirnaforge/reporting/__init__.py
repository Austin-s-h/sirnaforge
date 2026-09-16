"""Single-file HTML reporting over a finished run directory (issue #103).

Three modules with one responsibility each: :mod:`payload` reads a run directory and returns typed
data, :mod:`render` turns that data into one self-contained HTML file, and :mod:`quilt` registers that
file -- and the evidence behind it -- in ``quilt_summarize.json``. The split is what lets the CLI and
the pipeline emit the same artifacts from the same entry point.

Neither `payload` nor `render` classifies a hit, applies a threshold of its own, or reaches the network.
Gate descriptors come from the resolved run policy and ``hit_class`` comes from the published table, so
the report cannot disagree with the run it describes.
"""

from sirnaforge.reporting.payload import (
    PAYLOAD_SCHEMA_VERSION,
    GuideEntry,
    ReportInputError,
    ReportPayload,
    build_payload,
)
from sirnaforge.reporting.quilt import quilt_summarize_entries, write_quilt_summarize
from sirnaforge.reporting.render import render_html, write_report

__all__ = [
    "PAYLOAD_SCHEMA_VERSION",
    "GuideEntry",
    "ReportInputError",
    "ReportPayload",
    "build_payload",
    "quilt_summarize_entries",
    "render_html",
    "write_quilt_summarize",
    "write_report",
]
