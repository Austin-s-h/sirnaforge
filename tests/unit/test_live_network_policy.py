"""Exercise live-service opt-in through isolated pytest subprocesses, without HTTP."""

from pathlib import Path

import pytest

pytest_plugins = ["pytester"]
pytestmark = pytest.mark.unit
_CONFTEST = (Path(__file__).resolve().parents[1] / "conftest.py").read_text()


@pytest.fixture
def remote_suite(pytester, monkeypatch):
    """Build a standalone suite with the real hooks and no inherited CI state."""
    for variable in ("CI", "GITHUB_ACTIONS", "PYTEST_ADDOPTS"):
        monkeypatch.delenv(variable, raising=False)
    pytester.makeini("""
[pytest]
markers =
    requires_network: live service test
    dev: development tier
    ci: CI tier
    release: release tier
""")
    pytester.makepyfile("""
import pytest

def test_offline():
    assert True

@pytest.mark.requires_network
def test_remote():
    assert True
""")
    return pytester


def _hooks(suite, probe):
    suite.makeconftest(_CONFTEST + "\n\ndef _network_available():\n    " + probe + "\n")


@pytest.mark.parametrize("ci_variable", [None, "CI", "GITHUB_ACTIONS"])
def test_default_runs_never_probe_live_services(remote_suite, monkeypatch, ci_variable):
    """Normal local/CI runs skip remote bodies before any TLS probe."""
    if ci_variable:
        monkeypatch.setenv(ci_variable, "true")
    _hooks(remote_suite, 'raise AssertionError("unexpected network probe")')
    result = remote_suite.runpytest_subprocess("-q", "--strict-markers", timeout=30)
    result.assert_outcomes(passed=1, skipped=1)


def test_local_opt_in_runs_remote_tests(remote_suite):
    """A deliberate local invocation executes the retained remote functionality test."""
    _hooks(remote_suite, "return True")
    result = remote_suite.runpytest_subprocess("-q", "--run-network", timeout=30)
    result.assert_outcomes(passed=2)


def test_local_opt_in_does_not_hide_client_failures(remote_suite):
    """A genuine failure remains a failure when the developer enables live checks."""
    _hooks(remote_suite, "return True")
    remote_suite.makepyfile("""
import pytest
@pytest.mark.requires_network
def test_remote():
    assert False, 'client regression'
""")
    result = remote_suite.runpytest_subprocess("-q", "--run-network", timeout=30)
    result.assert_outcomes(failed=1)


def test_offline_local_opt_in_skips_without_reaching_test_body(remote_suite):
    """The existing verified-TLS capability check still applies to local opt-in."""
    _hooks(remote_suite, "return False")
    result = remote_suite.runpytest_subprocess("-q", "--run-network", timeout=30)
    result.assert_outcomes(passed=1, skipped=1)


@pytest.mark.parametrize("ci_variable", ["CI", "GITHUB_ACTIONS"])
def test_ci_cannot_enable_live_services(remote_suite, monkeypatch, ci_variable):
    """Even an explicit flag cannot introduce remote dependencies into CI."""
    monkeypatch.setenv(ci_variable, "true")
    _hooks(remote_suite, 'raise AssertionError("unexpected network probe")')
    result = remote_suite.runpytest_subprocess("-q", "--run-network", timeout=30)
    assert result.ret == pytest.ExitCode.USAGE_ERROR
    result.stderr.fnmatch_lines(["*--run-network is local-only*"])
