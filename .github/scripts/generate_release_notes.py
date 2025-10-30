#!/usr/bin/env python3
"""Generate release notes content for the siRNAforge release workflow."""

from __future__ import annotations

import os
import re
from pathlib import Path


def _extract_changelog_section(changelog_text: str, version: str) -> list[str]:
    pattern = rf"^## \[{re.escape(version)}\].*?(?=^## \[|\Z)"
    match = re.search(pattern, changelog_text, flags=re.MULTILINE | re.DOTALL)
    if not match:
        return ["- Initial release or changelog entry missing", ""]

    lines = match.group(0).splitlines()
    # Drop the repeated header and any horizontal rules used as separators.
    body = [line for line in lines[1:] if line.strip() != "---"]
    if body and body[-1].strip():
        body.append("")
    return body


def _collect_reference_links(changelog_text: str) -> list[str]:
    references = [line for line in changelog_text.splitlines() if re.match(r"^\[[^]]+\]:", line)]
    return references + [""] if references else []


def main() -> None:
    version = os.environ["VERSION"]
    py_files = os.environ.get("PY_FILES", "0")
    test_files = os.environ.get("TEST_FILES", "0")
    total_lines = os.environ.get("TOTAL_LINES", "0")
    cli_commands = os.environ.get("CLI_COMMANDS", "0")
    docker_size = os.environ.get("DOCKER_SIZE", "Unknown")

    changelog_path = Path("CHANGELOG.md")
    changelog_text = changelog_path.read_text(encoding="utf-8") if changelog_path.exists() else ""

    lines: list[str] = [f"## 📋 What's New in v{version}", ""]

    if changelog_text:
        lines.extend(_extract_changelog_section(changelog_text, version))
    else:
        lines.extend(["- No changelog available for this version", ""])

    lines.extend(
        [
            "## 🔧 Technical Details & Project Stats",
            "",
            "**Python Support:** 3.9, 3.10, 3.11, 3.12",
            "**Package Manager:** uv (ultra-fast Python package management)",
            "**Architecture:** Modern async/await with Pydantic models",
            "**Container:** Multi-stage Docker build with conda bioinformatics stack",
            "**Pipeline:** Nextflow integration for scalable execution",
            "",
            "**📊 Project Statistics:**",
            f"- 🐍 **Python Files:** {py_files} source files (~{total_lines} lines of code)",
            f"- 🧪 **Test Files:** {test_files} test files with comprehensive coverage",
            f"- 🎯 **CLI Commands:** {cli_commands} commands available (workflow, search, design, validate, config, cache, version)",
            f"- 🐳 **Docker Image:** {docker_size}",
            "- ⚡ **Dependencies:** Managed by uv for 10-100x faster installs than pip",
            "",
        ]
    )

    lines.extend(
        [
            "## 🧪 Quality Assurance",
            "",
            "This release has been validated through our comprehensive testing pipeline:",
            "- ✅ **Unit Tests** - Core algorithm validation",
            "- ✅ **Integration Tests** - End-to-end workflow testing",
            "- ✅ **Docker Tests** - Container functionality verification",
            "- ✅ **Smoke Tests** - Critical path validation (must pass for release)",
            "- ✅ **Code Quality** - Ruff, MyPy, and Black formatting",
            "",
        ]
    )

    if changelog_text:
        lines.extend(_collect_reference_links(changelog_text))

    output = "\n".join(lines).rstrip() + "\n"
    Path("release_notes.md").write_text(output, encoding="utf-8")

    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with Path(github_output).open("a", encoding="utf-8") as handle:
            handle.write("content<<EOF\n")
            handle.write(output)
            handle.write("EOF\n")

    print(f"Generated release notes for v{version} (saved to release_notes.md)")


if __name__ == "__main__":
    main()
