# Performance investigations: Pixi and BWA-MEM3

Status on 2026-10-24: **Pixi deferred; BWA-MEM3 plain mode accepted by the user for integration. Investigations remain open.**
This implements the bounded F1/F2 investigation tooling from #124 for #121/#122.
Existing package, execution, index and reporting defaults are unchanged.

## Recorded local evidence

The machine-readable observations are in `scripts/performance/measurements/`.
These files contain the revision, machine, exact commands, timeouts, outcomes and
resource measurements where GNU time is available. Measurements refer to this
agent's environment, which differs from the user's comprehensive local test environment.

| Observation | Result | What it establishes |
| --- | --- | --- |
| Installed tools | Pixi 0.76.1, uv 0.12.2, micromamba 2.3.2, Nextflow 26.04.2, OpenJDK 21.0.12.1 | PATH inventory, not environment parity |
| `pixi info` on isolated prototype | Three successful runs; median approximately 0.041 s | Manifest parses; this is metadata/startup timing only |
| Pixi lock/install | Network-enabled solve produced the committed Linux lock; cached runtime/application installs took approximately 1.05/3.73 s (one sample each) | Reproducible prototype, neither cold installs nor baseline comparisons |
| BWA comparison | Official MEM2 2.2.1 and source-built MEM3 0.15.0 passed controlled synthetic and 320-guide chr21 SAM parity and independent exact-match checks | Limited compatibility evidence, with compiler caveat below |
| samtools and RNAfold | Unavailable on PATH | Container/runtime validation remains necessary |

The repository's existing uv cache was read-only in this session. Commands below
use the existing `.venv` Python directly; no global tools were installed.

## Reproduce inventory and bounded command measurements

Run from the repository root with Python 3.11+ (or 3.10 with `tomli`):

```bash
.venv/bin/python scripts/performance/review.py --timeout 5 \
  --output /tmp/performance-inventory.json inventory
.venv/bin/python scripts/performance/review.py --timeout 30 \
  --output /tmp/pixi-info.json measure --cache-state 'existing caches; metadata only' \
  --repetitions 3 -- pixi info --manifest-path scripts/performance/pixi/pixi.toml --json
.venv/bin/python scripts/performance/test_review.py
```

`measure` executes an argument vector without a shell. Failed or timed-out trials
are recorded and produce exit status 1; the median excludes those trials and
`complete` remains false. Timeouts kill the subprocess group. Linux GNU time
adds child CPU time and peak RSS; other platforms retain wall time only.
Cache-state labels describe operator preparation, not verified OS cache eviction.
The machine can be shared: these measurements do not imply exclusive CPU access.

## F1: isolated Pixi prototype

`scripts/performance/pixi/pixi.toml` declares a Linux runtime environment and an
optional editable application environment. It fixes the principal tool versions,
uses Python/Java ranges, and includes a generated Linux lockfile for
transitive/build reproducibility. Its versions are an investigation target, not
an assertion that they match this host. Only Linux is declared so an untested
platform is not presented as supported.

When channel access is available, solve without installation first:

```bash
.venv/bin/python scripts/performance/review.py --timeout 180 \
  --output /tmp/pixi-lock.json measure --cache-state 'isolated new cache; lock only' \
  --repetitions 1 -- env PIXI_CACHE_DIR=/tmp/sirnaforge-pixi-cold-cache \
  pixi lock --manifest-path scripts/performance/pixi/pixi.toml
```

Only after reviewing the lock, install with `pixi install --locked` using the same
manifest path. Record separate cold-cache, locked fresh-install, warm-install,
update, disk/image-size, Docker-build and first/warm application-startup results.
Use separate throwaway prototype copies and caches for each cold trial; do not
relabel a reused environment as cold. Hold resolved Python/Nextflow/Java/native
tool versions fixed in the uv/micromamba baseline. Use the generic `measure`
command for each exact installation or startup argv and record at least three
successful trials plus variability. Comparable baseline and genuinely cold install trials remain pending.

Dependency ownership requiring resolution before adoption:

- `pyproject.toml` and `uv.lock` own Python dependencies; `environment-dev.yml`
  duplicates much of the stack. The captured inventory embeds both Conda sources
  and lists Python/dev requirements rather than silently choosing an owner.
- Nextflow is supplied by both PyPI and Conda. The Dockerfile preserves the
  Conda launcher after pip installation and sets `NXF_DIST` to its framework JAR.
  The application install overwrote the Conda launcher with PyPI Nextflow. A fresh `NXF_HOME` proves runtime Nextflow 25.10.2 starts offline, but application Nextflow attempts to download 26.10.3 and fails. The isolated-offline measurement files record this adoption blocker.
- Native ViennaRNA/pysam packages overlap PyPI distributions. Resolve explicit
  ownership and ABI compatibility; installation success alone is insufficient.
- The prototype application feature does not reproduce the full development
  dependency group yet; add it after establishing runtime ownership.
- Pixi managing the parent environment is separate from Nextflow creating
  process environments. Preserve the current micromamba process-environment
  mechanism until a supported replacement is demonstrated.
- Host tests, actual container/Nextflow execution, arbitrary-UID writable caches,
  offline framework startup, CI and other platform checks remain outstanding.

Adopt only with a reviewed lock, comparable successful benchmarks, compatibility
checks, and a migration/rollback plan. A measured defer decision can close the
investigation, but these partial comparisons and launcher conflict do not satisfy #121.

Official references checked for this investigation:
[Pixi manifest](https://pixi.prefix.dev/latest/reference/pixi_manifest/) and
[Conda/PyPI integration](https://pixi.prefix.dev/latest/concepts/conda_pypi/).

## F2: BWA-MEM3 parity and performance harness

Supply installed binaries explicitly if they are outside PATH:

```bash
.venv/bin/python scripts/performance/review.py --timeout 120 \
  --output /tmp/bwa-plain.json bwa --work-dir /tmp/bwa-review \
  --mem2 /path/to/bwa-mem2 --mem3 /path/to/bwa-mem3 --threads 1 --repetitions 3
.venv/bin/python scripts/performance/review.py --timeout 120 \
  --output /tmp/bwa-compat.json bwa --work-dir /tmp/bwa-review \
  --mem2 /path/to/bwa-mem2 --mem3 /path/to/bwa-mem3 --mem3-mode compat
```

Each execution creates a new retained session directory, builds separate indices,
records index-member hashes/sizes and binary version/checksum, and runs the actual
current transcriptome flags `-a -k 12 -T 15 -w 100` with matching thread limits.
No index files are shared between backends. Results contain raw repetitions,
median/range, and strict SAM multiset differences, retaining flags, coordinates,
MAPQ, CIGAR, sequences and all alignment tags except program provenance `PG`.
Reference `@SQ` headers are checked separately; record/tag order is ignored.
Duplicate records are retained as multiset counts. Exit status 1 indicates
missing tools, execution failure, nondeterminism, or parity differences. Exit 0
establishes only the reported SAM comparison, never permission to adopt.

The deterministic seed-124 fixtures include 19–23 nt queries, reverse matches,
multiple exact loci, one/two mismatches, an ambiguous no-hit query, and potential
clipping/insertion/deletion cases. The latter names describe construction, not
asserted alignment truth. An independent exhaustive exact-match oracle checks
full-length ungapped matches on both strands, so identical backend omissions
can be detected. It does not establish mismatch/gap sensitivity or all biological
eligibility criteria. Oracle results are reported separately from SAM parity.

Official upstream [mode documentation](https://bwa-mem3.readthedocs.io/en/latest/whats-different/modes.html)
distinguishes plain, compatibility and fast modes; plain may add `MQ`/`HN` tags,
and compatibility targets BWA-MEM2 2.2.1. Consequently strict plain-mode tag
differences need review rather than being treated as missing hits. The harness
supports both plain and `--compat=bwa-mem2`; it does not enable fast mode.
[Installation documentation](https://bwa-mem3.readthedocs.io/en/latest/getting-started/installation.html)
now describes Bioconda distributions, while the documentation homepage still
says Bioconda is forthcoming. Verify package availability and pin actual release,
build/license/source before relying on either statement. The local results tested official MEM2 2.2.1 and source release MEM3 0.15.0; hashes/build commands are in `aligner-builds.json`. MEM3 required `ALLOW_UNSUPPORTED_COMPILER=1` with GCC 13, so timings are not representative upstream performance.

The controlled synthetic comparison produced 32 equal SAM records. The pinned
chr21 workload contains 320 deterministic exact antisense queries (19–23 nt).
Its alignment medians were approximately 1.127 s for MEM2 and 0.094 s for MEM3,
with one thread and `OMP_NUM_THREADS=1`. Both independent exact-match checks
passed. These workloads include startup overhead, reuse OS caches and use an
unsupported compiler for MEM3. They do not establish general speedup or gate
equivalence. Raw trials are in `bwa-mem3-v0.15.0-controlled.json` and
`bwa-mem3-v0.15.0-chr21.json`.

The user accepted **plain mode** on 2026-10-24 and explicitly removed further
exact BWA-MEM2 parity proof as a requirement. Compatibility-mode measurements
remain historical evidence; they do not define the acceptance criterion.
Proceed with plain-mode integration and verify parser compatibility, complete
screening evidence and gate behavior. Exact SAM/tag equality is not a blocker.

Outstanding before #122 can close:

1. Run both profiles against pinned binaries, synthetic truth and representative
   cDNA panels, including large/high-hit and multispecies workloads. Hash inputs.
2. Validate plain-mode outputs through the existing downstream parser: hit
   classes/counts, seed/mismatch coordinates, completion/censoring and gate
   verdicts must remain meaningful. Exact MEM2 parity is waived by the user.
3. Benchmark cold/warm indexing and alignment, small/large query batches, CPU/RSS,
   index size and end-to-end time. Fresh index construction does not make OS
   page caches cold. Never interpret fewer retained hits as a speedup.
4. Prove index-format compatibility; add backend/version identity, rebuild and
   cleanup/staging rules before changing cache reuse. Fixed assumptions reside
   in `core/off_target.py`, `data/transcriptome_manager.py`, `provenance.py`,
   `pipeline/nextflow_cli.py` and embedded Nextflow modules/subworkflows.
5. Validate host/container/Nextflow, binary/version checks, distribution support,
   fallback and rollback; then record adopt/defer with measured gains.

## Coordination with report performance (#111)

Share checksummed query/reference inputs and environment metadata. Report timing
should consume saved screening/evidence artifacts and exclude environment setup
and alignment. Evaluate manager changes with BWA-MEM2 held fixed; evaluate aligner
changes within a fixed environment. Measure a combined configuration only after
the independent effects are understood.
