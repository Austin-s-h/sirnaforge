# PRD: consolidate development and release environments with Pixi

Status: proposed; implementation follows v0.8.0.

Date: 2026-11-22.

Decision owner: repository maintainer.

## Decision and outcome

Evaluate promoting the already validated Linux Pixi prototype to the primary
environment manager for local development, CI and the Linux release image.
Cutover remains conditional on closing the outstanding adoption gaps below.
Contributors should install one manager, install the committed lock, and run the
same development and release commands used in CI. Docker remains a distribution
and isolation boundary; Python wheels remain installable through standard tools.
One manager may provide several named environments from one manifest and lock.

This is a new adoption proposal, not a claim that the previous evaluation approved
migration. The [performance review](performance_review.md) records adoption as
previously deferred. Its historical timings have different package closures and
some different versions; they do not establish a general performance win.
[Issue #124 and its discussion](https://github.com/Austin-s-h/sirnaforge/issues/124)
record substantially more than installation feasibility: full serial host suites,
a dedicated locked Pixi CI lane, an actual Pixi Docker image, native-tool imports,
offline Nextflow startup, and Conda-backed Nextflow processes under UID 12345
have passed. The earlier asyncio/provider hang was caused by the execution
sandbox denying a wakeup socket, not a Pixi defect. Earlier overwrite warnings
were resolved through explicit dependency ownership.

The matched-direct-version image measured 1,735,350,992 bytes versus
2,252,135,962 bytes for the micromamba baseline (23% smaller). Pixi container CLI
startup was slower (median 2.44 s versus 1.73 s); Nextflow startup was similar
(1.44 s versus 1.51 s). These one-host observations have different dependency
closures and are not universal performance claims. Minimal Python/pytest update
trials passed, but whole-application update safety and broader platform coverage
remain unproven. #121 closed with adoption deferred, not with adoption approved.

This proposal reopens the adoption decision around those remaining gaps. Reuse
recorded evidence; rerun it where the candidate changes its assumptions rather
than treating all prior work as incomplete.

## Problem and baseline

Dependency ownership currently spans `pyproject.toml`, `uv.lock`,
`environment-dev.yml`, `docker/environment-nextflow.yml`, and the experimental
`scripts/performance/pixi/pixi.toml` / `pixi.lock`, plus shell and CI entry points.
Maintaining them independently permits different tools and scientific libraries
for the same workflow.

Observed on the v0.8.0 development tree:

| Area | Current difference | Required resolution |
| --- | --- | --- |
| Nextflow | uv 26.4.6; Pixi 26.10.3; Docker uses a lower bound | One launcher/JAR owner and tested runtime version |
| fm-index | uv 3.0.2; Pixi 4.0.0 | Choose a baseline with backend correctness evidence |
| PyArrow | uv 25.0.1; Pixi 26.0.0 | Choose compatible cache/schema behavior |
| Alignment | Docker and Pixi include MEM3; local YAML omits it | Same executable inventory for supported full environments |
| Security tools | Local YAML includes safety; Python dev configuration documents pip-audit instead | One documented audit command |
| Platforms | Pixi trial declares linux-64 only | Explicit supported platform and fallback policy |

Lock versions are observations, not recommendations to select the highest version.
Python dependency metadata must continue to describe wheel consumers, even after
repository environment locking is consolidated.

## Scope and exclusions

In scope: local setup, editor interpreter guidance, Make targets, pre-commit,
Python and native dependencies, docs tooling, CI caches, release packaging, Docker,
and a documented dependency update procedure.

Excluded from this migration: scientific default changes, adding classic BWA,
changing the default aligner from MEM2, a simultaneous broad dependency upgrade,
or promising native support on platforms without tested package availability.
No environment migration is bundled into v0.8.0. Separate compatibility changes,
including Python support changes and major library upgrades, receive separate PRs.

## Requirements

1. Put the authoritative Pixi configuration at the repository root, preferably in
   `pyproject.toml` to reuse standard project metadata. Commit one root `pixi.lock`.
   Consume project dependencies and development groups where supported; enforce
   agreement automatically for any unavoidable duplicate constraints.
2. Define `runtime`, `dev`, and `docs` environments. Share runtime constraints and
   compatible resolved versions across development and release. Python compatibility
   jobs may need separate solves; do not force incompatible Python versions into
   one solve group. Runtime excludes test, notebook and docs tools.
3. Assign each installed distribution one owner. Conda owns Python and native tools
   including Java, Nextflow, BWA, samtools and RNAfold. Use Pixi PyPI dependencies
   for packages unavailable or unsuitable on Conda, documenting exceptions such as
   PyArrow/AWS CLI solver interactions. Never overlay a second Nextflow launcher
   or native Python package with an independent pip/uv installation.
4. Preserve the validated editable project installation using
   `pip install --no-deps --no-build-isolation` as the initial baseline. This is an
   explicit local-project installation step, not a second dependency resolver;
   Conda-owned dependencies must remain untouched. A manifest-managed editable
   dependency is an option only after proving it preserves that ownership,
   launcher/JAR pairing and offline behavior. Do not blindly import all project
   PyPI requirements into Pixi and recreate the resolved overwrite defect.
5. Normal setup and CI must use locked installs and fail on manifest drift.
   Dependency updates are explicit reviewed lock changes. A documented clean setup
   must work without an existing `.venv`, Conda prefix or ambient tool PATH.
6. Preserve public `make` entry points, including `test-release`, delegating to
   Pixi tasks without Make/Pixi recursion. Version lookup must not install or solve
   another environment. Pre-commit and editors must use the intended interpreter.
7. Build the Docker runtime from the same locked Linux runtime closure. Install
   the built wheel without re-resolving dependencies; inspect wheel requirements
   against the installed closure. Preserve offline Nextflow JAR lookup, entrypoint,
   CPU limits, image source fingerprints and artifact provenance.
8. Preserve release safeguards: dated changelog, package/tag identity, tested wheel
   bytes, release-candidate gate, master ancestry and versioned publication gates.
9. Start with Linux x86_64 and WSL2 full-tool support. Record the current Python
   3.10–3.12 compatibility promise and retain its tests until a separate support
   decision changes it. Evaluate macOS/ARM availability explicitly; document Docker
   fallback where native tools are missing. Do not silently drop a supported use case.
10. Retain standard PEP 517 wheel/sdist metadata. An isolated wheel installation
    outside Pixi must pass package smoke tests, with external tool requirements
    documented. A build/publish utility can remain pinned within Pixi; it must not
    become a second development resolver or require a second runtime lock.

## Delivery plan

### Phase 0: inventory and baseline

Record dependency owners, supported platforms, every setup/test/build entry point,
and current clean host/container release receipts at an exact commit. Reuse the
recorded Pixi image/provider/UID/host-suite evidence as the baseline. Resolve the
version differences above deliberately before attributing behavior changes to Pixi.
Capture workflow outputs and tool inventories, not just successful installation.

Exit: maintainer reviews the ownership/platform matrix and chosen version baseline.
The primary new experiment is a representative whole-application update: update a
throwaway copy of the full lock, review all transitive changes, install cleanly,
run the host and container release selections, validate native imports and offline
Nextflow/provider/UID behavior, then restore the prior lock and verify rollback.
The existing minimal pytest-update timing cannot substitute for this evidence.

### Phase 1: root environment and developer pilot

Promote the trial into a root manifest/lock on a dedicated branch. Add development,
docs, runtime and compatibility environments; convert Make and pre-commit paths.
Keep the current production path available during the pilot. Update contributor
and editor setup instructions and add manifest drift checks.

Exit: two independent clean Linux/WSL installs, with recorded commands and receipts,
pass lint/type checks, development tests, docs and both BWA smoke checks.

### Phase 2: CI and container parity

Move CI and Docker builds to the root lock. Run the candidate against the established
release path using matching direct versions and fixtures. Test a clean runtime
without development tools and an isolated wheel consumer outside Pixi.

Exit: all release and compatibility gates below pass on the candidate commit.

### Phase 3: cutover and removal

Make Pixi the documented default. Remove `uv.lock`, independently maintained Conda
YAML environments and the experimental manifest once all their consumers migrate.
If an external integration needs an export, generate it from the authoritative
configuration and check drift rather than maintaining another source by hand.
Archive historical performance evidence without rewriting historical measurements.

Exit: repository search finds no active independent setup path; onboarding, CI,
Docker and release documentation all use the same ownership model.

## Acceptance gates

| Gate | Required evidence |
| --- | --- |
| Reproducibility | Clean locked installs on two Linux/WSL machines; no tracked lock changes after setup/tests |
| Tool ownership | Inventory identifies one executable/provider for Nextflow and its JAR, Java, both BWA variants, samtools and RNAfold |
| Functional parity | Existing backend, repeated-hit alignment, index, workflow and report tests pass; expected scientific outputs remain equal under existing tolerances |
| Release | Full `make test-release` passes at the candidate commit, including docs, host/container coverage and offline container gates |
| Packaging | Wheel/sdist checks and isolated wheel smoke pass; container runs the same tested wheel bytes with satisfied metadata |
| Compatibility | Every retained Python support job passes; unsupported native platforms have an explicit tested container path |
| Offline behavior | After dependency/JAR/reference staging, representative offline workflow succeeds with network disabled |
| Maintenance | One documented install/update/audit procedure; CI detects unintended lock or duplicate-constraint drift |
| Performance | Matched-version cold/warm trials report downloads, disk and image size; median setup/startup regression over 20% requires maintainer acceptance with rationale |

An installation success or a fast startup alone does not satisfy the cutover gate.
Live public API tests remain separate local opt-ins; public endpoint uptime must
not decide release readiness.

## Update policy and rollback

After cutover, update security fixes promptly and review routine updates monthly.
Group Python libraries, workflow/JVM tools and alignment/native libraries into
separate reviewable changes. Attach the resolved version diff, audit output,
applicable correctness tests and any index/cache compatibility implications.
Change declared minimum versions only when justified by supported behavior.

Until Phase 3, retain the existing release path as rollback. If ownership,
scientific parity, platform support or release gates fail, stop cutover and record
the blocker. After cutover, revert the migration commit series and restore the
last tested environment definitions; do not silently reintroduce parallel locks.

## Open decisions before implementation

- Exact baseline versions for fm-index, PyArrow and Nextflow.
- Native macOS/ARM demand and package availability versus the Docker fallback.
- Whether Python 3.10 support ends in a separate release, and when 3.13/3.14 enter CI.
- Whether root Pixi integration can consume every existing dependency group without
  duplication, and how its drift check handles necessary native-package mappings.

## References

- [Issue #124: evaluation history and acceptance decisions](https://github.com/Austin-s-h/sirnaforge/issues/124).

- [Existing repository evaluation](performance_review.md).
- [Pixi migration from uv](https://pixi.prefix.dev/latest/switching_from/uv/):
  project metadata integration and locked installation instead of ad hoc pip changes.
- [Pixi manifest reference](https://pixi.prefix.dev/dev/reference/pixi_manifest/):
  Conda/PyPI resolution and ownership.
