"""Test Docker profiles functionality."""

import sys
from pathlib import Path

import pytest
from test_profiles import TEST_PROFILES, DockerTestProfile, get_available_profiles

# Add the docker test directory to the path for imports
sys.path.insert(0, str(Path(__file__).parent))


@pytest.mark.docker
def test_profile_definitions():
    """Test that all profiles are properly defined."""
    required_profiles = ["minimal", "bioinformatics", "workflow", "performance", "comprehensive"]

    for profile_name in required_profiles:
        assert profile_name in TEST_PROFILES
        profile = TEST_PROFILES[profile_name]
        assert isinstance(profile, DockerTestProfile)
        assert profile.name == profile_name
        assert profile.description
        assert isinstance(profile.required_tools, list)
        assert isinstance(profile.test_markers, list)
        assert isinstance(profile.environment_vars, dict)


@pytest.mark.docker
def test_profile_skip_logic():
    """Test that profile skip logic works correctly."""
    # Get a profile that should be available
    minimal_profile = TEST_PROFILES["minimal"]
    skip_reason = minimal_profile.should_skip()

    # Minimal profile should not be skipped (no required tools)
    assert skip_reason is None


@pytest.mark.docker
def test_get_available_profiles():
    """Test getting available profiles."""
    available = get_available_profiles()

    # Should always have at least minimal profile
    assert "minimal" in available
    assert isinstance(available["minimal"], DockerTestProfile)


@pytest.mark.docker
def test_profile_pytest_args():
    """Test that pytest args are generated correctly."""
    profile = TEST_PROFILES["minimal"]
    args = profile.get_pytest_args()

    assert "-v" in args
    assert len(args) >= 1
