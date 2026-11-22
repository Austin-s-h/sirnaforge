"""Prepare synthetic report fixtures for actual-browser evidence acceptance checks."""

import argparse
import json
import tempfile
from pathlib import Path

import pandas as pd

from sirnaforge.reporting import write_report_from_run
from tests.unit.test_reporting_data_backed import _base4
from tests.unit.test_reporting_payload import _write_fully_evidenced_run


def main() -> None:
    """Write resolved and mixed-unavailable reports in each supported reader mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    guides = tuple(_base4(index) for index in range(300))
    for state in ("resolved", "unresolved"):
        with tempfile.TemporaryDirectory(prefix="report-browser-", dir=args.output) as temporary:
            run = _write_fully_evidenced_run(Path(temporary), guides=guides)
            path = run / "sirnaforge/manifest.json"
            manifest = json.loads(path.read_text())
            manifest["reference_evidence"] = {
                "schema_version": 1,
                "mirna": {
                    "source": "mirbase",
                    "acquisition": {"status": "complete"},
                    "screening": {"status": "complete", "guide_sequences": guides},
                },
            }
            path.write_text(json.dumps(manifest))
            if state == "unresolved":
                path = run / "sirnaforge/candidates_all.csv"
                frame = pd.read_csv(path)
                frame["max_off_target_count_verdict"] = ["unavailable" if i % 2 else "pass" for i in range(len(frame))]
                frame.to_csv(path, index=False)
            for mode in ("compact_single_file", "paged_single_file", "data_backed"):
                write_report_from_run(run, args.output / f"{mode}-{state}.html", mode=mode)


if __name__ == "__main__":
    main()
