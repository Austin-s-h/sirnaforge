#!/usr/bin/env python3
"""Generate complete release notes for the siRNAforge release workflow."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path


def _extract_changelog_section(changelog_text: str, version: str) -> list[str]:
    pattern = rf"^## \[{re.escape(version)}\].*?(?=^## \[|\Z)"
    match = re.search(pattern, changelog_text, flags=re.MULTILINE | re.DOTALL)
    if not match:
        return ["- Initial release or changelog entry missing", ""]

    lines = match.group(0).splitlines()
    body = [line for line in lines if line.strip() != "---"]
    if body and body[-1].strip():
        body.append("")
    return body


def _extract_unreleased_section(changelog_text: str) -> list[str]:
    pattern = r"^## \[Unreleased\].*?(?=^## \[|\Z)"
    match = re.search(pattern, changelog_text, flags=re.MULTILINE | re.DOTALL)
    if not match:
        return []

    lines = match.group(0).splitlines()
    body = [line for line in lines if line.strip() != "---"]
    if body and body[-1].strip():
        body.append("")
    return body


def _list_changelog_versions(changelog_text: str) -> list[str]:
    versions: list[str] = []
    for match in re.finditer(r"^## \[([^\]]+)\]", changelog_text, flags=re.MULTILINE):
        version = match.group(1).strip()
        if re.match(r"^[0-9]+\.[0-9]+\.[0-9]+", version):
            versions.append(version)
    return versions


def _version_key(version: str) -> tuple[int, int, int]:
    match = re.match(r"^(\d+)\.(\d+)\.(\d+)", version)
    if not match:
        return (0, 0, 0)
    return (int(match.group(1)), int(match.group(2)), int(match.group(3)))


def _get_latest_git_tag() -> str | None:
    try:
        result = subprocess.run(
            ["git", "tag", "--list", "v*", "--sort=-v:refname"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None

    for line in result.stdout.splitlines():
        tag = line.strip()
        if tag:
            return tag
    return None


def _select_versions_since_tag(changelog_text: str, current_version: str) -> list[str]:
    versions = _list_changelog_versions(changelog_text)
    if not versions:
        return [current_version]

    baseline_tag = os.environ.get("BASELINE_TAG") or _get_latest_git_tag()
    baseline_version = baseline_tag.lstrip("v") if baseline_tag else None

    current_key = _version_key(current_version)
    baseline_key = _version_key(baseline_version) if baseline_version else None

    selected: list[str] = []
    for v in versions:
        vk = _version_key(v)
        if vk <= current_key and (baseline_key is None or vk > baseline_key):
            selected.append(v)

    if not selected:
        return [current_version]

    selected.sort(key=_version_key, reverse=True)
    return selected


def _collect_reference_links(changelog_text: str) -> list[str]:
    references = [line for line in changelog_text.splitlines() if re.match(r"^\[[^]]+\]:", line)]
    return references + [""] if references else []


def _generate_installation_section(version: str, registry: str, repo_lower: str) -> list[str]:
    return [
        "## Installation",
        "",
        "### Docker (recommended)",
        "```bash",
        "# Pull the latest release",
        f"docker pull {registry}/{repo_lower}:{version}",
        "",
        "# Quick test - should complete in ~2 seconds",
        f"docker run --rm {registry}/{repo_lower}:{version} sirnaforge --help",
        "```",
        "",
        "### Python package",
        "```bash",
        "# Via pip",
        f"pip install sirnaforge=={version}",
        "",
        "# Via uv (recommended for speed)",
        f"uv add sirnaforge=={version}",
        "",
        "# Verify installation",
        "sirnaforge --help",
        "```",
        "",
        "### Development",
        "```bash",
        "# Clone repository",
        f"git clone https://github.com/{repo_lower}.git",
        "cd sirnaforge",
        "",
        "# Setup using the Makefile (uses uv under the hood)",
        "make dev",
        "",
        "# Verify with quick tests",
        "make test-dev  # ~15 seconds",
        "```",
        "",
    ]


def _generate_usage_examples(version: str, registry: str, repo_lower: str) -> list[str]:
    return [
        "## Usage",
        "",
        "### Basic workflow (Docker)",
        "```bash",
        "# Complete gene-to-siRNA workflow",
        "docker run --rm -v $(pwd):/workspace -w /workspace \\",
        f"  {registry}/{repo_lower}:{version} \\",
        "  sirnaforge workflow TP53 --output-dir results --genome-species human",
        "",
        "# Custom transcript file",
        "docker run --rm -v $(pwd):/workspace -w /workspace \\",
        f"  {registry}/{repo_lower}:{version} \\",
        "  sirnaforge design transcripts.fasta --output results.csv",
        "```",
        "",
        "### Python API",
        "```python",
        "from sirnaforge import SiRNADesigner, GeneSearcher",
        "",
        "# Search for gene transcripts",
        'searcher = GeneSearcher(species="human")',
        'transcripts = searcher.search_gene("TP53")',
        "",
        "# Design siRNAs",
        "designer = SiRNADesigner()",
        "candidates = designer.design_from_transcripts(transcripts)",
        "```",
        "",
        "### Commands",
        "```bash",
        "sirnaforge workflow   # Complete pipeline: gene → siRNAs",
        "sirnaforge search     # Gene/transcript search",
        "sirnaforge design     # siRNA candidate generation",
        "sirnaforge validate   # Input file validation",
        "sirnaforge config     # Show configuration",
        "sirnaforge cache      # Manage miRNA databases",
        "sirnaforge version    # Version information",
        "```",
        "",
    ]


def _generate_key_features(version: str) -> list[str]:
    return [
        f"## Key features in v{version}",
        "",
        "- **Smart siRNA design** - Thermodynamic scoring with 90% duplex binding weight",
        "- **Off-target analysis** - BWA-MEM2 alignment across human/rat/rhesus genomes",
        "- **ViennaRNA integration** - Secondary structure prediction for accuracy",
        "- **Pandera data schemas** - Type-safe output validation and formatting",
        "- **uv package manager** - Fast dependency resolution",
        "- **Production Docker** - Pre-built images with all bioinformatics tools",
        "- **Nextflow pipeline** - Scalable execution with automatic parallelization",
        "",
    ]


def _generate_testing_section() -> list[str]:
    return [
        "## Testing & quality",
        "",
        "**This release passed comprehensive validation:**",
        "- **Unit tests** - Core algorithm validation",
        "- **Local Python tests** - Fastest development iteration",
        "- **Docker smoke tests** - Critical functionality verification",
        "- **Integration tests** - End-to-end workflow validation",
        "- **Code quality** - Ruff formatting and MyPy type checking",
        "",
    ]


def _generate_resources_section(
    registry: str, repo_lower: str, owner_lower: str, name_lower: str, version: str
) -> list[str]:
    return [
        "## Resources",
        "",
        "### Documentation",
        f"- [**Full documentation**](https://{owner_lower}.github.io/{name_lower}) - Complete user guide",
        f"- [**Quick start**](https://github.com/{repo_lower}#-quick-start) - Get running in minutes",
        f"- [**API reference**](https://{owner_lower}.github.io/{name_lower}/api) - Python API documentation",
        f"- [**Development guide**](https://github.com/{repo_lower}/blob/main/CONTRIBUTING.md) - Contributing instructions",
        "",
        "### Container images",
        f"- **Versioned:** `{registry}/{repo_lower}:{version}`",
        f"- **Latest:** `{registry}/{repo_lower}:latest`",
        f"- **Registry:** [{registry}/{repo_lower}]({registry}/{repo_lower})",
        "",
        "### Support",
        f"- [**Source code**](https://github.com/{repo_lower})",
        f"- [**Issues**](https://github.com/{repo_lower}/issues) - Bug reports & feature requests",
        f"- [**Discussions**](https://github.com/{repo_lower}/discussions) - Community support",
        f"- [**Changelog**](https://github.com/{repo_lower}/blob/main/CHANGELOG.md) - Version history",
        f"- [**All releases**](https://github.com/{repo_lower}/releases) - Previous versions",
        "",
    ]


def _generate_verification_section(version: str, registry: str, repo_lower: str) -> list[str]:
    return [
        "---",
        "",
        "**Quick verification:**",
        "```bash",
        "# Test Docker image (should complete in ~2 seconds)",
        f"docker run --rm {registry}/{repo_lower}:{version} sirnaforge version",
        "",
        f"# Expected output: siRNAforge v{version}",
        "```",
    ]


def main() -> None:
    """Generate complete release notes content."""
    version = os.environ["VERSION"]

    # Get repository info from environment
    registry = os.environ.get("REGISTRY", "ghcr.io")
    repo_lower = os.environ.get("REPOSITORY_LOWER", "")
    owner_lower = os.environ.get("OWNER_LOWER", "")
    name_lower = os.environ.get("NAME_LOWER", "sirnaforge")

    changelog_path = Path("CHANGELOG.md")
    changelog_text = changelog_path.read_text(encoding="utf-8") if changelog_path.exists() else ""

    # Start with title
    lines: list[str] = [
        f"# 🧬 siRNAforge v{version}",
        "",
        "**Comprehensive siRNA design toolkit with multi-species off-target analysis**",
        "",
    ]

    # Add changelog section(s)
    if changelog_text:
        selected_versions = _select_versions_since_tag(changelog_text, version)
        for v in selected_versions:
            section = _extract_changelog_section(changelog_text, v)
            if section == ["- Initial release or changelog entry missing", ""]:
                # Common for prereleases: CHANGELOG.md often only has an [Unreleased]
                # section until the final version is cut.
                unreleased = _extract_unreleased_section(changelog_text)
                if unreleased:
                    lines.extend(unreleased)
                else:
                    lines.extend(section)
            else:
                lines.extend(section)
    else:
        lines.extend([f"## What's new in v{version}", "", "- Initial release", ""])

    # Add all major sections
    lines.extend(_generate_installation_section(version, registry, repo_lower))
    lines.extend(_generate_usage_examples(version, registry, repo_lower))
    lines.extend(_generate_key_features(version))
    lines.extend(_generate_testing_section())
    lines.extend(_generate_resources_section(registry, repo_lower, owner_lower, name_lower, version))
    lines.extend(_generate_verification_section(version, registry, repo_lower))

    # Add reference links if any
    if changelog_text:
        refs = _collect_reference_links(changelog_text)
        if refs:
            lines.extend(["", ""])
            lines.extend(refs)

    output = "\n".join(lines).rstrip() + "\n"
    Path("release_notes.md").write_text(output, encoding="utf-8")

    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with Path(github_output).open("a", encoding="utf-8") as handle:
            # Preferred output name used by release workflow.
            handle.write("content<<EOF\n")
            handle.write(output)
            handle.write("EOF\n")
            # Backwards-compat output name (older workflows may consume this).
            handle.write("changelog<<EOF\n")
            handle.write(output)
            handle.write("EOF\n")

    print(f"Generated complete release notes for v{version} (saved to release_notes.md)")


if __name__ == "__main__":
    main()
