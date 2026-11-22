"""Release identity must be consistent before publication."""

import hashlib
import importlib.util
import zipfile
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github" / "scripts" / "release_contract.py"
spec = importlib.util.spec_from_file_location("release_contract", SCRIPT)
assert spec and spec.loader
release_contract = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release_contract)


@pytest.mark.parametrize("version", ["0.8.0", "0.8.0rc1", "1.0.0"])
def test_release_accepts_semver_and_pep440_versions(version):
    """Stable and rc versions have canonical shared spellings."""
    assert release_contract.parse_version(version)


@pytest.mark.parametrize("version", ["0.8.0.post1", "0.8.0+local", "0.8.0.dev1", "0.08.0", "0.8.0-rc.1"])
def test_release_rejects_noncanonical_versions(version):
    """Reject ambiguous or unsupported prerelease versions."""
    with pytest.raises(ValueError):
        release_contract.parse_version(version)


def test_tag_cannot_name_another_package_version():
    """Do not release another package version under a tag."""
    with pytest.raises(ValueError, match="package"):
        release_contract.validate_identity("0.8.0", "0.7.1", "## [0.8.0] - 2026-11-22", tagged=True)


def test_tag_requires_dated_changelog_entry():
    """A pending changelog is not evidence of publication."""
    with pytest.raises(ValueError, match="dated changelog"):
        release_contract.validate_identity("0.8.0", "0.8.0", "## [0.8.0] - Pending release", tagged=True)


def test_release_candidate_allows_pending_changelog():
    """Release PRs can await their final publication date."""
    release_contract.validate_identity("0.8.0", "0.8.0", "## [0.8.0] - Pending release", tagged=False)


def test_missing_changelog_is_not_a_release():
    """Require a release-specific notes section."""
    with pytest.raises(ValueError, match="changelog"):
        release_contract.validate_identity("0.8.0", "0.8.0", "## [Unreleased]", tagged=False)


def test_wheel_metadata_cannot_disagree_with_tag(tmp_path):
    """Metadata, not filenames, defines a wheel's version."""
    wheel = tmp_path / "sirnaforge-0.7.1-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(
            "sirnaforge-0.7.1.dist-info/METADATA", "Metadata-Version: 2.3\nName: sirnaforge\nVersion: 0.7.1\n"
        )
    with pytest.raises(ValueError, match="wheel"):
        release_contract.validate_wheel(wheel, "0.8.0")


def test_matching_wheel_metadata_is_accepted(tmp_path):
    """Accept a wheel with the selected release version."""
    wheel = tmp_path / "sirnaforge-0.8.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr(
            "sirnaforge-0.8.0.dist-info/METADATA", "Metadata-Version: 2.3\nName: sirnaforge\nVersion: 0.8.0\n"
        )
    release_contract.validate_wheel(wheel, "0.8.0")


def test_published_wheel_must_match_build_bytes(tmp_path):
    """A same-name but different artifact must fail publication."""
    wheel = tmp_path / "sirnaforge-0.8.0-py3-none-any.whl"
    wheel.write_bytes(b"tested wheel")
    published = {"urls": [{"filename": wheel.name, "digests": {"sha256": "wrong"}}]}
    with pytest.raises(ValueError, match="digest"):
        release_contract.validate_published_wheel(wheel, published)


def test_published_wheel_accepts_identical_build(tmp_path):
    """Published bytes match the distribution that passed tests."""
    wheel = tmp_path / "sirnaforge-0.8.0-py3-none-any.whl"
    wheel.write_bytes(b"tested wheel")
    published = {
        "urls": [{"filename": wheel.name, "digests": {"sha256": hashlib.sha256(wheel.read_bytes()).hexdigest()}}]
    }
    release_contract.validate_published_wheel(wheel, published)
