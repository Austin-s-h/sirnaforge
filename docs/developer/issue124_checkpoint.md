# Issue #124 development checkpoint

ZFN output remains **EXPERIMENTAL** and requires independent validation; see
[#82](https://github.com/Austin-s-h/sirnaforge/issues/82).

Saved locally on 2026-10-24. This is an implementation checkpoint, not a claim
that issue #124 or its scientific validation milestones are complete.

## Accepted decisions

- The user reported 1,904 passing local tests, including network opt-in tests,
  before these changes. That run does not validate subsequent implementation.
- Prefer resilient pinned/local miRNA references; remote failure must remain
  failed/unknown evidence, never a successful zero-hit result.
- Extend the existing HTML report for evidence, identity and conservation.
- The user accepts **BWA-MEM3 plain mode** and explicitly waives further exact
  BWA-MEM2 parity proof. Continue plain-mode integration with parser, evidence
  completeness and gate checks. Compatibility measurements are historical only.
- Pixi remains an investigation: the application installation overwrites the
  Conda Nextflow launcher and breaks isolated offline startup.

## Saved implementation groups

- ZFN correctness: inward geometry, ambiguous-reference verification, shard
  ownership, cache identity and candidate summary schema v2 (`97fc086`).
- Native benchmark geometry, review-only unsupported records, four-species
  mapping capability, verified scientific references and honest GAA template
  (`b295b66`). Vendored target mappings are still human-only.
- Performance harness, Pixi Linux lock and real measurements (`3aadc80`).
- Pinned miRNA source manifests, checksums, atomic refresh and task-produced
  reference evidence; source/canonical/orthology provenance (`cb7af96`).
- Strict variant resolution and actual alternate-allele guide design (`f5198b1`).
- HTML evidence/identity/conservation, safe source links, projected artifacts,
  compact/paged/data-backed rendering and resource diagnostics (`792aa72`).
- Explicit target intent and separate complementary transcript-seed evidence
  (`f981db1`). Formatting checkpoint: `a89d119`.

Use `git log --oneline` for subsequent checkpoint commit identifiers. Changes
are local; no push, merge or default package-manager/backend migration is implied.

## Validation and commands

The final integrated offline host checkpoint passed **2,073 tests, 61 deselected
in 123.71 seconds**, after code commit `f08bf24`. Results are in
`/tmp/sirnaforge-124-checkpoint.xml` and `/tmp/sirnaforge-124-checkpoint.log`.
This excludes network/container tests; the user's earlier opt-in run predates
these changes. Ruff and formatting checks pass; mypy passes 106 source files.

The earlier full run exposed missing transcript-seed policy/schema registration.
Those integration gaps are fixed in `cba0599` and `f08bf24`; the policy/gate
focused run passed 152 tests and schema focused run passed 14 tests. The report
mock/documentation focused run passed 48 tests. No known failing test remains
in the final offline host checkpoint.

The sandbox hangs on asyncio thread shutdown; the same offline host command
works outside that restriction. A writable temporary cache is required here:

```bash
SIRNAFORGE_CACHE_DIR=/tmp/sirnaforge-goal-cache .venv/bin/python -m pytest \
  -m '(dev or ci or release) and not runs_in_container and not requires_network' \
  -n 0 --no-cov -o addopts='' --tb=short \
  --junitxml=/tmp/sirnaforge-124-checkpoint.xml
.venv/bin/ruff check src scripts/fetch_native_reference_release.py \
  scripts/native_mapping_generator.py scripts/performance tests/unit
.venv/bin/mypy src/sirnaforge
.venv/bin/python scripts/performance/test_review.py
```

The full-tree mypy check passed 106 source files after final agent handoff;
Ruff also passed the integrated source, scripts and unit tests.
The five standalone performance-harness tests passed. A partial commit's mypy
hook saw new untracked modules against temporarily stashed tracked dependencies;
that artificial snapshot failed. Run integrated checks against all saved groups.

## Resume priorities and closure limits

1. Finish and verify public CLI/workflow forwarding for target intent,
   complementary transcript-seed screening and all report modes. A seed screen
   must use the full resolved screening cDNA reference, independently of target
   annotation/enumeration. Missing 3'UTR annotation is unavailable evidence.
2. Check report-wide status/threshold behavior and genuinely large artifact
   performance. Data-backed reports require their local detail directory.
3. Integrate BWA-MEM3 **plain mode**, with backend-specific index identity and
   validated downstream parsing/completeness. Do not spend further time proving
   exact MEM2 parity. MEM3 local timing used unsupported GCC 13, recorded in the
   measurements; do not advertise a general performance gain from that fixture.
4. Generate actual native mappings from four pinned release references. The
   downloader and efficient mapping generator are saved; real four-species
   output is not yet proven. Any files in `/tmp/sirnaforge-ensembl-release116`
   are scratch artifacts; trust only verified manifests, never partial downloads.
5. Finish #82's public Nextflow ZFN route/flag forwarding and real-genome audit.
6. #46's actual GAA demonstration needs authoritative transcript and reviewed
   GRCh38 allele inputs; the bare c.428A>G identity is ambiguous. Do not invent
   a GAA variant or count the synthetic fixture as scientific validation.
7. #98/#102 require independent held-out biological validation before scientific
   defaults can be promoted. Existing synthetic tests cannot satisfy that.
8. Resolve Pixi dependency ownership, benchmark fair baselines and container/
   offline behavior before adoption. Keep its investigation separate from BWA.
9. Update #124 and parent issues with specific evidence. Keep partially addressed
   issues open; passing software tests alone does not satisfy scientific claims.

The large report-test stall was resolved in `3cb8c13`: a mistaken assertion
rejected intentionally inline gene identity labels, then pytest spent minutes
formatting the 820KB HTML string. The corrected assertion checks offloaded hit
transcripts. All three data-backed/paged tests passed in 3.39 seconds. This does
not establish compact-mode performance on large datasets.

The human GRCh38 Ensembl release-116 reference was downloaded and verified:
expanded SHA-256 `bd85494fb64766200ec18bae873f3196b451c6eecbcfbedadd72684f19692ff8`.
Its manifest and FASTA are in `/tmp/sirnaforge-ensembl-release116`; the other three
species and actual generated outputs remain pending.

Reference evidence smoke artifacts are under `/tmp/sirnaforge-pinned-smoke`.
Performance raw evidence and pinned workload preparation are tracked under
`scripts/performance/`; downloaded binaries and environments in `/tmp` are
recreatable and must not be treated as persistent deployment dependencies.

## Resume checkpoint — 2026-10-25

This section supersedes pending-work statements above where they conflict.

Committed and pushed on `dev`:

- `2fcfd5a`: live CI formatting/docs repairs and annotated installed-CLI variant fixtures.
- `dc09370`: optional explicit MEM3 backend integration and public workflow/report modes;
  compact Quilt registration, nullable mode metadata, backend-specific index identity.
- `8e766fb`: full release-116 four-species paired-target mapping tables and pinned provenance.
- `ebd68f6`: performance harness type annotations required by `make lint` (`src scripts`).

Issue #110 is closed with completion evidence. All 92,760 mapped transcript sites
were independently checked against the four checksummed full cDNA FASTAs and their
oriented query sequences. Public native import retained all 2,850 Ichihara rows.
Nonhuman mappings remain sequence projections of human assays; efficacy labels are
not transferred. Reference records are now tracked alongside the mapping tables.

Validation: 2,103 offline host tests passed in 127.20 seconds; two additional
compact/auto report CLI regressions passed. Ruff and format pass; mypy passes all
131 source/script files; five standalone performance harness tests pass. The final
warnings-as-errors offline docs build passes. Live CI on `ebd68f6` passed main
push/PR checks and documentation. Docker CI remains under observation.

A current Docker image built successfully, including offline Nextflow startup.
The full local network-disabled installed-wheel suite passed 45 tests with one
index-override failure that attempted a runtime download. The index-override case and both fresh-pin cases subsequently passed sequentially
in 48.13 seconds; do not treat the earlier full local run as passing. Earlier source-overlay/stale
image failures are not release evidence. A new fresh-pin container acceptance test
checks valid and bad-checksum miRBase inputs through the public CLI, Nextflow
producer and rendered report; both cases passed individually and together with the index-override rerun.

Priority remains: finish CI; #123/#112 reference availability and honest evidence;
#111–#114 report behavior across all modes; #101/#47 intent/seed/variant contracts;
#82 public ZFN route and #122 actual MEM3 runtime; #46 reviewed GAA and #121 Pixi;
#98/#102 independent biological validation. Keep partially completed issues open.
In particular, #123 has bounded socket attempts but no proven strict total download
deadline for a trickling server. No archive/mirror was automatically adopted.


## Reference resilience acceptance checkpoint — 2026-10-25

`ebd274f` replaces the per-socket-only acquisition behavior with a subprocess
worker that is killed and reaped on a total per-attempt deadline. This supersedes
the deadline gap identified above. The local-server regression sends bytes every
20 ms, exceeds a one-second deadline, and proves the old cached FASTA and metadata
survive unchanged. Normal downloads, gzip and decoding are checked through the
actual installed worker. No third-party mirror/archive has been adopted.

Paged and data-backed reports now consume the same source/guide evidence registry
as compact reports. They expose source provenance and remediation, an unresolved
view filter, separate gate-pass versus unresolved-evidence counts, and an
independently scrolling threshold inspector. Normalized guide tables have a fixed
Arrow schema; nullable fields in the first record no longer erase later values.
The canonical non-hit branch now checks the full transcript register rather than
the selected guide's own hit set.

Validation: 2,117 integrated offline host tests passed (67 deselected) in 141.22s;
all six valid/invalid pinned-reference cases passed through the installed public
CLI, real Nextflow, and each of the three report modes with networking disabled
in 80.80s. Ruff/format pass; mypy passes 132 source/script files; warnings-as-errors
docs pass. Focused reader tests execute the unresolved filter and stable warning
copy in Node without modifying recorded gate arrays. #123 is ready for closure
with this evidence after publication. #112 remains open for its broader vocabulary
and UI acceptance audit; #111/#113/#114 still need their remaining full criteria.
