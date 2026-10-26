# Performance investigations: Pixi and BWA-MEM3

Status on 2026-10-26: **Pixi evaluation remains open; BWA-MEM3 plain mode accepted by the user for integration.**
This implements the bounded F1/F2 investigation tooling from #124 for #121/#122.
Existing package, execution, index and reporting defaults are unchanged.

## Recorded local evidence

Machine-readable BWA observations are in `scripts/performance/measurements/`.
Pixi trial results are summarized here rather than checked in as raw,
machine-specific command logs. Measurements refer to this agent's environment,
which differs from the user's comprehensive local test environment.

| Observation | Result | What it establishes |
| --- | --- | --- |
| Installed tools | Pixi 0.76.1, uv 0.12.2, micromamba 2.3.2, Nextflow 26.04.2, OpenJDK 21.0.12.1 | PATH inventory, not environment parity |
| `pixi info` on isolated prototype | Three successful runs; median approximately 0.041 s | Manifest parses; this is metadata/startup timing only |
| Pixi runtime cold install | One isolated Linux x86_64 trial: 9.96 s, 246 MiB peak RSS, approximately 285 MB cache; OS cache was not evicted | A single observed cold environment install, not a stable benchmark |
| Micromamba runtime cold install | One isolated Linux x86_64 trial with matching direct tool versions: 35.87 s, 1,446 MiB peak RSS, approximately 374 MB downloads; installed environment approximately 1.1 GB | Direct versions match, but dependency closures differ and this is one trial |
| Warm environment install | Three idempotent trials each; Pixi median 0.060 s, micromamba median 5.862 s | Cached no-op install comparison only; not environment creation or update performance |
| Pixi application and offline Nextflow | Editable application install succeeded after pinning the PyPI launcher to 25.10.2 and pointing `NXF_DIST` at the Conda framework JAR; `sirnaforge --help` and three offline Nextflow launches succeeded | Resolves the observed launcher/framework mismatch for this locked Linux prototype |
| BWA comparison | Official MEM2 2.2.1 and source-built MEM3 0.15.0 passed controlled synthetic and 320-guide chr21 SAM parity and independent exact-match checks | Limited compatibility evidence, with compiler caveat below |
| samtools and RNAfold | Unavailable on PATH | Container/runtime validation remains necessary |

The repository's existing uv cache was read-only in this session. Commands below
use the existing `.venv` Python directly; no global tools were installed.

The cold install comparison used Pixi 0.76.1 against micromamba 2.3.2 on Linux
x86_64/WSL. Both runtime environments pinned Python 3.12.15, Nextflow 25.10.2,
OpenJDK 17.0.18, BWA-MEM2 2.2.1, samtools 1.24, pysam 0.24.1, and ViennaRNA
2.7.0. The Pixi lock resolved 94 packages; the micromamba baseline resolved 98,
including four additional Python packaging tools, so only the requested direct
versions match. Network downloads were enabled and the OS file cache was not
evicted. Cold results are single trials. Warm results reran each manager's
install/sync command against its already-created environment; they measure a
no-op check, not a warm environment creation.

Representative commands (use fresh roots/prefixes for cold trials and the same
roots for warm no-op trials):

```bash
PIXI_CACHE_DIR=/tmp/pixi-cold pixi install --locked \
  --manifest-path scripts/performance/pixi/pixi.toml --environment runtime
MAMBA_ROOT_PREFIX=/tmp/mamba-cold micromamba create -y \
  -n runtime -c conda-forge -c bioconda \
  python=3.12.15 nextflow=25.10.2 openjdk=17.0.18 bwa-mem2=2.2.1 \
  samtools=1.24 pysam=0.24.1 viennarna=2.7.0
PIXI_CACHE_DIR=/tmp/pixi-cold pixi install --locked \
  --manifest-path scripts/performance/pixi/pixi.toml --environment runtime
MAMBA_ROOT_PREFIX=/tmp/mamba-cold micromamba install -y -n runtime \
  -c conda-forge -c bioconda \
  python=3.12.15 nextflow=25.10.2 openjdk=17.0.18 bwa-mem2=2.2.1 \
  samtools=1.24 pysam=0.24.1 viennarna=2.7.0
```

Record timing with `scripts/performance/review.py measure` as described below;
for cold installation, set the manager-specific roots to new paths for every
repetition. The install commands above document direct requirements, not an
identical transitive closure.

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
successful trials plus variability. Three cold trials per manager, identical dependency closures, update behavior, and representative application startup benchmarks remain pending. Do not infer general Pixi speedup from warm no-op measurements.

Dependency ownership requiring resolution before adoption:

- `pyproject.toml` and `uv.lock` own Python dependencies; `environment-dev.yml`
  duplicates much of the stack. The captured inventory embeds both Conda sources
  and lists Python/dev requirements rather than silently choosing an owner.
- Nextflow is supplied by both PyPI and Conda. The Dockerfile preserves the
  Conda launcher after pip installation and sets `NXF_DIST` to its framework JAR.
  The application install initially overwrote the Conda launcher with PyPI
  Nextflow 26.10.3, which failed offline against the bundled 25.10.2 framework.
  The prototype now pins both to 25.10.2 and points `NXF_DIST` at the Conda
  framework; isolated offline launches succeed. Failed pre-fix observations
  remain only as historical evidence.
- Native ViennaRNA/pysam packages overlap PyPI distributions. Resolve explicit
  ownership and ABI compatibility; installation success alone is insufficient.
- The prototype test feature installs pytest dependencies. The focused host
  test `test_gencode_structure_versions_and_gzip` passes in Pixi Python 3.12.15
  (pytest 9.1.1, pytest-asyncio 1.4.1) when run outside the execution sandbox.
  Inside the sandbox, the executor thread's asyncio wakeup-socket write raises
  `EPERM`; asyncio suppresses that socket error, leaving the loop blocked in
  `selectors.select` after the worker future completes. This is a sandbox
  restriction rather than a Pixi/provider hang. The same sandbox condition can
  affect other `asyncio.to_thread` work that completes after the event loop has
  begun waiting. The full serial Pixi host selection
  `(dev or ci or release) and not runs_in_container and not requires_network`
  passed outside the sandbox: 2,195 passed and 75 were deselected. That run
  exposed a separate provenance bug: the current BWA executable version was
  reported as the build version for an unstamped index. The index provenance
  now leaves that version unknown unless its stamp records a usable build-time
  version; the current runtime executable is recorded separately. After the
  fix, the final full-host run passed 2,195 tests with 75 deselected. A fresh
  offline Nextflow 25.10.2 local-process smoke also passed from the Pixi test
  environment and wrote its expected output; `sirnaforge --help` passed from
  the editable install. A Nextflow process explicitly pointed at the Pixi
  prefix with `NXF_CONDA_CLI=micromamba` did not run: this host lacks the
  `conda` activation command/shim Nextflow invokes. The production Docker image
  supplies its Conda command path; this host experiment does not validate that
  image integration or justify replacing the existing process provider.
  Reproduce with `pixi run --manifest-path scripts/performance/pixi/pixi.toml
  --environment tests python -m pytest -m '(dev or ci or release) and not
  runs_in_container and not requires_network' -n 0`. Green GitHub CI and Docker
  runs also cover their configured environments, not the Pixi container; real
  Pixi container/provider/cache validation remains necessary before adoption.
- Pixi managing the parent environment is separate from Nextflow creating
  process environments. Preserve the current micromamba process-environment
  mechanism until a supported replacement is demonstrated.
- Docker is installed but the local Docker socket denies access. Actual
  container/Nextflow execution, arbitrary-UID writable caches, provider/index
  behavior, CI under Pixi, and other platform checks remain outstanding.

**Decision: defer Pixi adoption.** Keep uv/pip as the host Python owner and
micromamba/Conda as the existing Nextflow process-environment provider. Do not
replace either tool or the Docker setup with the prototype. Rollback is to
remove the isolated `scripts/performance/pixi/` prototype and its lockfile; no
production environment, cache, image, or user command needs migration because
none was switched. Reconsider only after the full Pixi host suite passes,
provider/index and arbitrary-UID cache behavior pass in a real container, and
comparable repeated cold/update/startup measurements justify the added
environment owner. Until then, retain this prototype as an isolated evaluation.

#121 remains open until host tests, container/Nextflow and cache/provider checks,
comparable benchmark repetitions, and a reviewed migration/rollback recommendation
are complete.

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
can be detected. The oracle uses a multi-pattern automaton over a streaming
reference scan; it does not need to load a full transcriptome into memory or run
one reference search per query. `--mem2-index` can resume a comparison with a
previously built MEM2 index after checking the exact FASTA checksum and required
index members; report that reuse because it does not measure a cold MEM2 index
build. The oracle does not establish mismatch/gap sensitivity or all biological
eligibility criteria. Oracle results are reported separately from SAM parity.

Official upstream [mode documentation](https://bwa-mem3.readthedocs.io/en/latest/whats-different/modes.html)
distinguishes plain, compatibility and fast modes; plain may add `MQ`/`HN` tags,
and compatibility targets BWA-MEM2 2.2.1. Consequently strict plain-mode tag
differences need review rather than being treated as missing hits. The harness
supports both plain and `--compat=bwa-mem2`; it does not enable fast mode.
[Installation documentation](https://bwa-mem3.readthedocs.io/en/latest/getting-started/installation.html)
describes Bioconda distributions. The [Bioconda bwa-mem3 recipe](https://bioconda.github.io/recipes/bwa-mem3/README.html)
currently lists release `0.15.0-0`, its mixed MIT/Apache-2.0/Zlib licensing,
and Linux aarch64 and macOS arm64 in addition to its default Linux x86-64
builds. This verifies package metadata, not a successful local install. An
attempt to add it to the isolated Pixi comparison environment could not resolve
the Bioconda index because network DNS is unavailable in this session; the
prototype lock remains unchanged. The local results tested official MEM2 2.2.1
and source release MEM3 0.15.0; hashes/build commands are in
`aligner-builds.json`. MEM3 required `ALLOW_UNSUPPORTED_COMPILER=1` with GCC 13,
so timings are not representative upstream performance.

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
The optional `--aligner bwa-mem3` route now builds and records a MEM3-specific
index, invokes plain `bwa-mem3 mem`, parses its SAM output, and records the
selected executable version and index identity. BWA-MEM2 remains the default.
A local MEM3 0.15.0 smoke built an index from a two-transcript FASTA and sent
plain-mode SAM records through the existing hit parser, producing six parsed
hits with transcript IDs, coordinates, strands, CIGAR, NM, and seed mismatch
fields. Focused backend, parser, and cache tests passed. This checks the software
path and output completeness for that fixture; it does not establish biological
eligibility, gate behavior across representative workloads, or a general speedup.
Exact SAM/tag equality is not a blocker.

### Release-116 cDNA panel results

To exercise representative references beyond the chr21 fixture, the same 3,395
unique full-length guides (20–21 nt; derived from 4,113 OligoGym observations)
were aligned once with one thread against Ensembl release-116 rat and macaque
cDNA. MEM2 2.2.1 and MEM3 0.15.0 plain mode each built their own index. Inputs
were hashed; the guide FASTA SHA-256 is
`5236565907d4940fa95b94fe0e1aae95122e7180775871c089a3dd36cea2884b`, rat
reference `14ad668039a84a3ffaf787944dd09194488392848574b20c484d50df56b07477`,
and macaque reference
`c702162e3053fb7bf4b0851f30b8b538d9885acfc737863bf2cc6e09d78db48d`.
There was one measurement per backend; wall-clock values include index load and
are not a multi-run performance claim.

| Reference | Exact placements (oracle) | SAM records | MEM2 / MEM3 alignment | Peak RSS, MEM2 / MEM3 |
| --- | ---: | ---: | ---: | ---: |
| Rat | 2,451 / 2,451 each | 17,025 each; SAM parity equal | 2.249 / 0.393 s | 909 / 666 MiB |
| Macaque | 4,701 expected; 4,697 observed each | 17,086 each; SAM parity equal | 2.062 / 0.324 s | 768 / 578 MiB |

The four macaque exact placements missed by both aligners and no unexpected
exact placements were reported. This shared omission needs interpretation in
the workload's biological context; equality between backends does not establish
completeness. For both species, the actual Python SAM parser produced identical
parsed hits across backends (16,323 rat; 16,573 macaque). A real local-profile
Nextflow run with MEM3 then completed the transcriptome tasks and evidence
reconciliation: 3,395/3,395 guides processed for each species, 16,323 rat sites
and 16,573 macaque sites, with both evidence statuses `complete`, no truncation,
and no unplanned evidence. The run disabled miRNA seed screening; that separate
database channel could not be validated offline.

The full Ensembl release-116 human cDNA contains 1.28 Gbp (FASTA SHA-256
`bd85494fb64766200ec18bae873f3196b451c6eecbcfbedadd72684f19692ff8`). MEM2
completed the full-reference alignment, but the exact-match oracle found 64,405
placements versus 62,247 reported (2,158 missing, all in two repetitive guides).
MEM3 0.15.0 refused to build its index: its suffix-array builder estimated a
28.5 GiB requirement, above this host's 24.8 GiB available-memory limit. This
is a concrete human-scale indexing constraint for this version and host; do not
infer that MEM3 is suitable for full human cDNA based on the smaller species.

Taken together, these results establish a working optional MEM3 host and
Nextflow path on these rat and macaque inputs, with matching SAM and parser
outputs, and expose an important large-reference limitation. They do not
establish exact-guide biological validation, gate-verdict equivalence, a
general speedup, or human-scale feasibility. MEM2 remains the default. The
temporary full-reference FASTAs, indices, SAM files, reports, and raw result
JSON were kept outside the repository and removed after this summary was
recorded.

The pinned Docker package now also passes a real Nextflow MEM3 smoke on a small
synthetic transcriptome: it builds a MEM3-specific index, emits and parses the
exact guide alignment, publishes complete one-guide evidence, and records the
selected MEM3 release in `transcriptome/versions.yml`. This ran in the Docker
workflow on commit `73463ef` (Actions run
https://github.com/Austin-s-h/sirnaforge/actions/runs/40271766698); Python CI
also passed (run
https://github.com/Austin-s-h/sirnaforge/actions/runs/40271766723). The test
caps the two synthetic-reference tasks at 2 GB for GitHub's 8 GB container;
production resource defaults are unchanged. This does not replace the larger
rat/macaque panel evidence above or establish full human cDNA feasibility.

Outstanding before #122 can close:

1. Run both profiles against pinned binaries, synthetic truth and representative
   cDNA panels, including large/high-hit and multispecies workloads. Hash inputs.
2. Exercise representative cDNA panels through the downstream parser and
   evidence reconciliation, including hit classes/counts, seed/mismatch
   coordinates, completion/censoring, and gate verdicts. Exact MEM2 parity is
   waived by the user.
3. Benchmark cold/warm indexing and alignment, small/large query batches, CPU/RSS,
   index size and end-to-end time. Fresh index construction does not make OS
   page caches cold. Never interpret fewer retained hits as a speedup.
4. The selected MEM3 executable, backend-specific index, plain SAM parser, and
   version output now pass the pinned Docker/Nextflow smoke above. Host
   packaging/distribution, index cache/rebuild behavior with the packaged
   binary, and other supported platforms remain outstanding before any default
   change.
5. Complete the remaining host checks, distribution support, fallback and
   rollback review, then record adopt/defer with measured gains.

## Coordination with report performance (#111)

Share checksummed query/reference inputs and environment metadata. Report timing
should consume saved screening/evidence artifacts and exclude environment setup
and alignment. Evaluate manager changes with BWA-MEM2 held fixed; evaluate aligner
changes within a fixed environment. Measure a combined configuration only after
the independent effects are understood.
