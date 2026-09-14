# Benchmark artifacts

The `sirnaforge benchmark` surface turns a published siRNA efficacy panel into a **versioned,
checksummed artifact** that the existing fixed-length design path consumes with no new execution
path, and records what every filter did to every measured observation under **two** policies at once.

Issue [#109](https://github.com/Austin-s-h/sirnaforge/issues/109). Variable-length and asymmetric
design are deliberately **not** here; they are [#110](https://github.com/Austin-s-h/sirnaforge/issues/110).

## Which panels are actually present

:::{warning}
**This repository vendors the bytes of exactly one benchmark panel.**

`huesken_subset` — 180 rows of measured knockdown across 8 transcripts, at
`tests/unit/data/sirna_efficacy_subset.csv`, with its transcripts beside it. Its primary citation,
the third-party redistribution it was obtained from, its SHA-256, its predeclared
development/held-out split and an explicit statement of what it does and does not prove are all in
[`tests/unit/data/README.md`](https://github.com/Austin-s-h/sirnaforge/blob/master/tests/unit/data/README.md).
Because the Huesken paper is not open access, those efficacy values were **not** verified against the
primary source.

**Ichihara, Martinelli, Shmushkovich and OligoGym are named in #109/#110 and are not vendored
anywhere in this repository.** The panel registry declares descriptors for them so a manifest can
name them and record `panel.data_present: false` honestly, but for those four the column mapping is a
*contract a future `--panel-csv` must satisfy*, not a description of a file anyone here has read, and
their citation strings say so where no primary source has been independently verified. `benchmark
prepare` refuses them outright without an explicit `--panel-csv`.

Every automated test over this surface therefore runs against `huesken_subset` or against synthetic
`tests/unit/data/benchmark/synthetic_*.csv` fixtures built to reproduce each declared architecture
(documented in
[`tests/unit/data/benchmark/README.md`](https://github.com/Austin-s-h/sirnaforge/blob/master/tests/unit/data/benchmark/README.md)).
**A passing test is evidence about the artifact contract, never about an unvendored panel's data.**

The PRD #109 cites as its source, `docs/prd_benchmark_artifacts_and_variable_length.md`, does not
exist in the working tree or in any branch's history. The issue body is the entire specification, and
nothing in this page summarises a document that was read.
:::

## The artifact

One directory per **(compatible panel × paired length)**, named `<panel_id>__len<paired_length>`,
with five fixed inner filenames:

| File                 | Written by | Holds                                                                       |
| -------------------- | ---------- | --------------------------------------------------------------------------- |
| `observations.csv`   | `prepare`  | One row per **source** row — compatible or not, never dropped               |
| `design_inputs.fasta`| `prepare`  | One record per **compatible** observation; a plain FASTA                     |
| `manifest.json`      | both       | Provenance, both policy blocks, counts, checksums                            |
| `candidates_all.csv` | `design`   | Every **enumerated** candidate, including enumeration-time rejects           |
| `accounting.csv`     | `design`   | The observation↔candidate join, with **both** filter verdict sets           |

Paired length is bounded to **19–23 nt** — the range `sirnaforge design` already accepts. An artifact
outside it could never be consumed by the path it exists to feed.

### Nothing measured is dropped, padded, trimmed or relabelled

An observation this build cannot pair to the requested fixed length still gets its
`observations.csv` row, with `compatibility_status: incompatible` and a stated
`compatibility_reason`. An asymmetric duplex (e.g. a 15/20 pair) is *never* sliced into a symmetric
core — silently truncating it would misrepresent a measured sequence, which #109 forbids outright.
Sequences are copied as measured and upper-cased only; a lower-case letter in a sequence column is
rejected by the schema, because it would mean some caller normalised rather than copied.

Two optional-column conventions are worth knowing before you parse a row. Every optional column uses
the ordinary "empty cell means `None`" convention, **except** `guide_3p_overhang` and
`passenger_3p_overhang`: there, an empty cell means *measured, and blunt* and the sentinel
`<not_stated>` means *the panel states no overhang*. Both are real, different facts, and CSV cannot
otherwise tell them apart.

## `benchmark prepare`

```bash
uv run sirnaforge benchmark prepare --panel huesken_subset --out-dir benchmark_artifacts
```

```text
┌──────────── Benchmark Summary ─────────────┐
│ 🧬 Benchmark artifact prepared             │
│ Panel: huesken_subset (data present: True) │
│ Paired length: 19 nt                       │
│ Observations kept: 180 (incompatible: 0)   │
└────────────────────────────────────────────┘
```

`--paired-length` defaults to the panel descriptor's declared length. A panel that ships vendored
bytes **refuses** an explicit `--panel-csv`, so a run cannot silently read different bytes than the
ones its own manifest names as vendored; a panel with no vendored bytes **requires** one.

### The prepared FASTA is the compatibility guarantee

`design_inputs.fasta` is a plain FASTA keyed by `design_context_id`, so the unmodified fixed-length
path consumes it with no benchmark code in the loop at all:

```bash
uv run sirnaforge design benchmark_artifacts/huesken_subset__len19/design_inputs.fasta --length 19
```

That is the executable form of #109's "consumable by the current fixed-length design/workflow path"
criterion, and it is pinned by a test rather than only documented here.

## `benchmark design`

```bash
uv run sirnaforge benchmark design --artifact benchmark_artifacts/huesken_subset__len19
```

```text
┌────────────────────── Benchmark Summary ──────────────────────┐
│ 🧬 Benchmark design complete                                  │
│ Entered design: 180  no candidate: 0                          │
│ Default pass: 111  Benchmark pass: 111                        │
└───────────────────────────────────────────────────────────────┘
```

### One run, two verdict sets

`accounting.csv` carries five accounting columns per row — `entered_design`,
`default_filter_status`, `default_filter_reasons`, `benchmark_filter_status`,
`benchmark_filter_reasons` — but **there is exactly one design pass**, under the benchmark policy.
The default set is *re-derived* by re-applying the default policy's filter descriptors to the values
that one run already recorded per candidate.

That re-derivation is exact rather than approximate only because widening is monotone: an observed
value like GC content or longest polynucleotide run is a property of the candidate, not of the
threshold, so a widened run enumerates a **superset** of what a default run would have. Which is why
a **narrowing** of either GC bound is refused before any file is opened:

```console
$ uv run sirnaforge benchmark design --artifact ... --gc-min 40
--gc-min 40.0 narrows the default floor of 30.0; the benchmark interface may only widen
gc_min, never raise it, because a narrowed run cannot stand in for the default run its
default_filter_status claims to re-derive
```

An observation with no candidate keeps its `accounting.csv` row with `entered_design: false`,
`candidate_id` empty and both statuses `not_evaluated`: **a gate cannot be blamed for a candidate
that was never built.** The schema enforces that pairing rather than trusting it.

### Per-run GC widening, with the default policy untouched

`--gc-min` / `--gc-max` widen the bounds **for that run only**. Nothing in this surface writes a
global default; the manifest carries *both* resolved policies so the separation is checkable rather
than asserted:

```bash
uv run sirnaforge benchmark design --artifact ... --gc-min 20 --gc-max 75
```

| Manifest field                                          | Value              |
| ------------------------------------------------------- | ------------------ |
| `gc_widening.gc_min` → `default` / `benchmark`           | `30.0` / `20.0`    |
| `gc_widening.gc_min.source`, `.widened`                 | `explicit`, `true` |
| `run_policy` resolved `gc_min`                          | `20.0` (explicit)  |
| `default_run_policy` resolved `gc_min`                  | `30.0` (builtin_profile) |
| `run_policy.profile.content_hash` vs `default_run_policy`'s | identical      |

`widened` is read off **provenance**, not off `benchmark != default`: a value test cannot tell an
omitted option from one typed with the default value.

### The polynucleotide-run requirement stays active

`max_poly_runs <= 3`, action `fail`, is resolved identically for both policies in every call — only
`gc_min`/`gc_max` are ever placed in the stated settings the benchmark policy is built from, so
nothing on this surface can relax it. It is recorded in **three** independent places:

1. per row, in `default_filter_reasons` **and** `benchmark_filter_reasons`;
2. per run, in `counts.excluded_by_filter["max_poly_runs"]`, split `default`/`benchmark`;
3. per artifact, in the `polynucleotide_run_requirement` block, which states the comparator,
   threshold, action, whether it was evaluated, and both exclusion counts.

The observed value itself is in `candidates_all.csv`'s `max_poly_runs_observed` column alongside its
`max_poly_runs_verdict`, so an exclusion can be checked, not just read.

### `candidates_all.csv` is a superset, on purpose

The designer drops a candidate outright when an enumeration-time gate rejects it with a `fail`
action; it survives only in the rejected list. A homopolymer observation that a default run would
also have rejected still has to appear in the accounting with a real verdict on both columns, so
`candidates_all.csv` is **survivors plus enumeration-time rejects** — a superset of what
`sirnaforge design` would write for the same input. On the vendored subset at length 19 that is
361 + 179 = 540 rows.

## Reproducibility

The manifest is designed so a second run over the same bytes can be *compared*, not merely trusted:

- every entry in `inputs` and `outputs` carries its own SHA-256, and `outputs` paths are the fixed
  inner filenames rather than wherever `--out-dir` happened to point;
- `invoked_command` is `sys.argv` verbatim, never reconstructed from options;
- `panel.descriptor_hash` covers the registry entry, so an edited descriptor is as visible as edited
  data, and both policy blocks carry the profile's content hash, covering thresholds and per-filter
  actions;
- neither CSV carries a timestamp or a path, and rows are written sorted by their join key.

`created_utc` is the single field the schema deliberately leaves nondeterministic, and it lives only
in `manifest.json`, never in a CSV. Two full `prepare` + `design` passes over the same panel bytes
with the same argv produce byte-identical `observations.csv`, `design_inputs.fasta`,
`candidates_all.csv` and `accounting.csv`, and manifests that differ in that one field alone.

## What this surface does not do

Out of scope for #109, and deliberately absent rather than stubbed:

- **variable-length or asymmetric design** — an asymmetric observation is recorded as incompatible,
  never coerced (#110);
- **native transcript mapping** — `target_transcript_id`, `target_start_1based`, `target_end_1based`
  and `target_strand` are placeholder fields the artifact contract needs; no registered panel maps a
  per-row position column, so `target_identity_status` is `unavailable` on every row and
  `counts.mapped_native` is pinned to `0` by the schema. A non-zero value there would be a claim this
  build cannot back up;
- **any change to global default GC thresholds** — widening is per run, through this interface only;
- `max_repeat_transcript_fraction` is excluded from both verdict sets, because the only code that
  ever measures it runs from the post-screening step against a real cDNA reference, never from the
  design-from-file path a benchmark run uses. Left in, it reported `unknown` on essentially every row
  and drowned the verdicts this surface exists to report.

See the `benchmark` section of the [CLI reference](cli_reference.md) for the generated option help.
