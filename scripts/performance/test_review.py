"""Meaningful harness checks, runnable without bioinformatics tools or pytest."""

# ruff: noqa: D101, D102

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from review import differences, exact_oracle, run, sam_headers, sam_records


class PerformanceHarnessTests(unittest.TestCase):
    def test_order_ignored_but_missing_hits_and_tags_detected(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            left, right = Path(folder) / "left.sam", Path(folder) / "right.sam"
            header = "@SQ\tSN:ref\tLN:100\n"
            first = "q\t0\tref\t2\t0\t19M\t*\t0\t0\tACGT\t*\tNM:i:0\tMD:Z:19"
            second = "q\t256\tref\t30\t0\t19M\t*\t0\t0\tACGT\t*\tNM:i:1\tMD:Z:8A10"
            left.write_text(header + first + "\n" + second + "\n")
            right.write_text(header + second + "\n" + first + "\n")
            self.assertTrue(differences(sam_records(left), sam_records(right))["equal"])
            self.assertEqual(sam_headers(left), sam_headers(right))
            right.write_text(header + first + "\n")
            self.assertEqual(differences(sam_records(left), sam_records(right))["left_only_count"], 1)
            right.write_text(header + first.replace("NM:i:0", "NM:i:2") + "\n" + second + "\n")
            self.assertFalse(differences(sam_records(left), sam_records(right))["equal"])

    def test_nonzero_and_missing_command_are_not_measurement_success(self) -> None:
        self.assertEqual(run([sys.executable, "-c", "raise SystemExit(2)"], 5)["status"], "failed")
        # GNU time can report missing executable as nonzero rather than OSError.
        self.assertIn(run(["/definitely/missing/tool"], 5)["status"], {"unavailable", "failed"})

    def test_timeout_is_reported(self) -> None:
        self.assertEqual(run([sys.executable, "-c", "import time; time.sleep(2)"], 0.05)["status"], "timeout")

    def test_malformed_sam_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.sam"
            path.write_text("truncated\trow\n")
            with self.assertRaises(ValueError):
                sam_records(path)

    def test_independent_exact_oracle_detects_shared_omission(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            reference, queries, sam = [Path(folder) / name for name in ("reference.fa", "query.fa", "hits.sam")]
            reference.write_text(">ref\nACGTTTACGT\n")
            queries.write_text(">q\nACGT\n")
            sam.write_text("q\t0\tref\t1\t0\t4M\t*\t0\t0\tACGT\t*\tNM:i:0\n")
            result = exact_oracle(reference, queries, sam)
            # Palindromic query has both orientations at both exact positions.
            self.assertEqual(result["expected_count"], 4)
            self.assertEqual(len(result["missing"]), 3)
            self.assertFalse(result["equal"])


if __name__ == "__main__":
    unittest.main()
