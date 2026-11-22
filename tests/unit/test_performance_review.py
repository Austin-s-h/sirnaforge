"""Regression tests for the BWA performance comparison result policy."""

from typing import Any

from scripts.performance.review import bwa_comparison_succeeded


def test_reviewed_sam_differences_do_not_mask_execution_failures() -> None:
    """Allow reported SAM differences only after complete backend execution."""
    result: dict[str, Any] = {
        "backends": {"mem2": {"status": "ok"}, "mem3": {"status": "ok"}},
        "parity": {"equal": False, "reference_headers_equal": True},
    }

    assert not bwa_comparison_succeeded(result)
    assert bwa_comparison_succeeded(result, allow_sam_differences=True)

    result["backends"]["mem3"]["status"] = "alignment_failed"
    assert not bwa_comparison_succeeded(result, allow_sam_differences=True)
