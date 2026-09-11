"""Single-file HTML reporting over a finished run directory (issue #103).

Two modules with one responsibility each: :mod:`payload` reads a run directory and returns typed
data, :mod:`render` turns that data into one self-contained HTML file. The split is what lets the CLI
and the pipeline emit the same artifact from the same entry point.

Neither module classifies a hit, applies a threshold of its own, or reaches the network. Gate
descriptors come from the resolved run policy and ``hit_class`` comes from the published table, so the
report cannot disagree with the run it describes.
"""

from sirnaforge.reporting.payload import (
    PAYLOAD_SCHEMA_VERSION,
    GuideEntry,
    ReportInputError,
    ReportPayload,
    build_payload,
)
from sirnaforge.reporting.render import render_html, write_report

__all__ = [
    "PAYLOAD_SCHEMA_VERSION",
    "GuideEntry",
    "ReportInputError",
    "ReportPayload",
    "build_payload",
    "render_html",
    "write_report",
]
