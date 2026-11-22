"""Published notes must describe the selected version, not static marketing claims."""

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / ".github" / "scripts" / "generate_release_notes.py"
spec = importlib.util.spec_from_file_location("generate_release_notes", SCRIPT)
assert spec and spec.loader
notes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(notes)


def test_release_notes_use_exact_changelog_and_valid_branch_links(tmp_path, monkeypatch):
    """Release text comes from the selected notes and points at master."""
    (tmp_path / "CHANGELOG.md").write_text("## [0.8.0] - 2026-11-22\n\n- Actual fix\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("VERSION", "0.8.0")
    monkeypatch.setenv("REPOSITORY_LOWER", "austin-s-h/sirnaforge")
    monkeypatch.setenv("OWNER_LOWER", "austin-s-h")
    monkeypatch.setenv("NAME_LOWER", "sirnaforge")
    monkeypatch.setenv("BASELINE_TAG", "v0.7.1")
    monkeypatch.delenv("GITHUB_OUTPUT", raising=False)
    notes.main()
    output = (tmp_path / "release_notes.md").read_text(encoding="utf-8")
    assert "- Actual fix" in output
    assert "/blob/main/" not in output
    assert "90% duplex binding weight" not in output
    assert "sirnaforge==0.8.0" in output
