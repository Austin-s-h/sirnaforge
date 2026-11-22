# 🧬 siRNAforge v0.8.0

**Comprehensive siRNA design toolkit with multi-species off-target analysis**

## [0.8.0] - 2026-11-22

### Added

- Report evidence schema 1.4.0: canonical availability states, stable generic gate-warning IDs,
  separate clean-evidence counts, unresolved guide views and explicit threshold units/provenance.
  Optional local Chromium acceptance checks cover all offline report modes.

- Pinned local mature miRNA reference manifests for workflow and off-target commands,
  with staged Nextflow inputs, checksum validation and task-produced reference evidence.
  See [local miRNA references](docs/mirna_references.md).


### Fixed

- `max_off_target_count` is now an opt-in query-species liability ceiling over `nm<=2`. The previous
  default of 15 was calibrated against an all-species, all-mismatch total and is not transferred to
  the newly scoped metric; all-species `off_target_count` remains descriptive.

- Recorded `unavailable`, `not_established` and legacy unknown gate states no longer fall through
  to numeric defaults. Threshold and guide panes scroll independently, and staged readers fit
  below their actual header height instead of overflowing the viewport.

- Native benchmark review preserves unpublished post-screen scores as missing and uses the
  supplied source weights for the current-profile comparison. Report register-neighbour lookup
  uses sorted position windows instead of repeatedly scanning entire transcripts. Nextflow
  version artifacts terminate their heredocs correctly and record the numeric BWA-MEM2 version.

- RefSeq gene lookup reads nested NCBI ELink `nuccore` groups and selects RefSeq RNA links,
  fixing empty transcript results for genes such as TP53. E-utilities calls are paced and
  retry HTTP 429 with bounded backoff. The README now documents existing
  local GENCODE and miRBase/combined-source support instead of marking them TODO.

- Host and container release gates exclude live-network tests. Remote functionality remains
  available locally through `make test-requires-network` / `pytest --run-network`; normal test
  runs skip it and CI refuses the opt-in, so public API outages cannot block releases.

- Automatic Nextflow CPU ceilings respect available processors, affinity and container quotas,
  preventing index/alignment tasks from requesting more CPUs than a small CI runner provides.
  Explicit production configuration overrides are accepted without duplicate-keyword errors.
- The container variant gate uses local transcripts and synthetic cached variants; live TP53/rsID
  checks remain in a separate integration suite instead of depending on public endpoint uptime
  during image publication. The default image gate excludes all `requires_network` checks;
  live container checks remain local-only and require `DOCKER_TEST_MARKERS=runs_in_container`
  plus `DOCKER_TEST_PYTEST_ARGS=--run-network`.

## Installation

### Docker (recommended)
```bash
# Pull this versioned release
docker pull ghcr.io/austin-s-h/sirnaforge:0.8.0

# Quick test - should complete in ~2 seconds
docker run --rm ghcr.io/austin-s-h/sirnaforge:0.8.0 sirnaforge --help
```

### Python package
```bash
# Via pip
pip install sirnaforge==0.8.0

# Via uv (recommended for speed)
uv add sirnaforge==0.8.0

# Verify installation
sirnaforge --help
```

### Development
```bash
# Clone repository
git clone https://github.com/austin-s-h/sirnaforge.git
cd sirnaforge

# Setup using the Makefile (uses uv under the hood)
make dev

# Verify with quick tests
make test-dev  # ~15 seconds
```

## Testing & quality

**Release validation is pending completion of `make test-release` and the release-candidate gate:**
- **Unit tests** - Core algorithm validation
- **Local Python tests** - Fastest development iteration
- **Docker smoke tests** - Critical functionality verification
- **Integration tests** - End-to-end workflow validation
- **Code quality** - Ruff formatting and MyPy type checking

## Resources

### Documentation
- [**Full documentation**](https://austin-s-h.github.io/sirnaforge) - Complete user guide
- [**Quick start**](https://github.com/austin-s-h/sirnaforge#-quick-start) - Get running in minutes
- [**API reference**](https://austin-s-h.github.io/sirnaforge/api) - Python API documentation
- [**Development guide**](https://github.com/austin-s-h/sirnaforge/blob/master/CONTRIBUTING.md) - Contributing instructions

### Container images
- **Versioned:** `ghcr.io/austin-s-h/sirnaforge:0.8.0`
- **Registry:** `ghcr.io/austin-s-h/sirnaforge`

### Support
- [**Source code**](https://github.com/austin-s-h/sirnaforge)
- [**Issues**](https://github.com/austin-s-h/sirnaforge/issues) - Bug reports & feature requests
- [**Discussions**](https://github.com/austin-s-h/sirnaforge/discussions) - Community support
- [**Changelog**](https://github.com/austin-s-h/sirnaforge/blob/master/CHANGELOG.md) - Version history
- [**All releases**](https://github.com/austin-s-h/sirnaforge/releases) - Previous versions

---

**Quick verification:**
```bash
# Test Docker image (should complete in ~2 seconds)
docker run --rm ghcr.io/austin-s-h/sirnaforge:0.8.0 sirnaforge version

# Expected output: siRNAforge v0.8.0
```


[Unreleased]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.2...HEAD
[0.2.2]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/austin-s-h/sirnaforge/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/austin-s-h/sirnaforge/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/austin-s-h/sirnaforge/releases/tag/v0.1.0
