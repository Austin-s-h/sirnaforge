#!/usr/bin/env python3
"""Validation script for the fast CI/CD workflow setup."""

import subprocess
import sys
import time
from pathlib import Path

from tests.docker.test_config import smoke_profile


def check_toy_data() -> bool:
    """Validate toy data files exist and are appropriately sized."""
    print("🧪 Checking toy data files...")

    toy_files = ["tests/unit/data/toy_transcripts.fasta", "tests/unit/data/toy_candidates.fasta"]

    for file_path in toy_files:
        path = Path(file_path)
        if not path.exists():
            print(f"❌ Missing: {file_path}")
            return False

        size = path.stat().st_size
        if size > 500:
            print(f"⚠️  {file_path} is {size} bytes (should be < 500)")
        else:
            print(f"✅ {file_path}: {size} bytes")

    return True


def check_test_markers() -> bool:
    """Validate that smoke test markers work correctly."""
    print("🏷️ Checking test markers...")

    try:
        # Test smoke marker collection
        result = subprocess.run(
            ["python", "-m", "pytest", "--collect-only", "-m", "smoke", "-q"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )

        if result.returncode == 0:
            lines = result.stdout.count("\n")
            print(f"✅ Found {lines} smoke tests")
            return True
        print(f"❌ Pytest collection failed: {result.stderr}")
        return False

    except subprocess.TimeoutExpired:
        print("❌ Pytest collection timed out")
        return False
    except Exception as e:
        print(f"❌ Error checking markers: {e}")
        return False


def check_docker_config() -> bool:
    """Validate Docker test configuration."""
    print("🐳 Checking Docker test configuration...")

    try:
        project_root = Path(__file__).parent.parent
        sys.path.insert(0, str(project_root))

        profile = smoke_profile()
        print(f"✅ Smoke profile loaded: {profile.docker_memory}")

        # Check memory is minimal
        if profile.docker_memory in ["256m", "512m"]:
            print(f"✅ Memory allocation appropriate: {profile.docker_memory}")
        else:
            print(f"⚠️  Memory allocation high: {profile.docker_memory}")

        return True

    except Exception as e:
        print(f"❌ Error loading config: {e}")
        return False


def run_quick_smoke_test() -> bool:
    """Run a quick smoke test to validate functionality."""
    print("🔥 Running quick smoke test...")

    start_time = time.time()

    try:
        result = subprocess.run(
            ["python", "-m", "pytest", "tests/unit/test_toy_data.py", "-v", "-m", "smoke"],
            check=False,
            capture_output=True,
            text=True,
            timeout=60,
        )

        duration = time.time() - start_time

        if result.returncode == 0:
            print(f"✅ Smoke tests passed in {duration:.2f} seconds")
            return True
        print(f"❌ Smoke tests failed: {result.stderr}")
        return False

    except subprocess.TimeoutExpired:
        print("❌ Smoke tests timed out")
        return False
    except Exception as e:
        print(f"❌ Error running smoke tests: {e}")
        return False


def main() -> int:
    """Main validation function."""
    print("🚀 Validating Fast CI/CD Setup")
    print("=" * 50)

    checks = [check_toy_data, check_test_markers, check_docker_config, run_quick_smoke_test]

    passed = 0
    for check in checks:
        try:
            if check():
                passed += 1
            print()
        except Exception as e:
            print(f"❌ Check failed with exception: {e}")
            print()

    print(f"📊 Results: {passed}/{len(checks)} checks passed")

    if passed == len(checks):
        print("🎉 All checks passed! Fast CI/CD is ready.")
        return 0
    print("⚠️  Some checks failed. Review configuration.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
