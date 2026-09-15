# Benchmark artifacts

The `sirnaforge benchmark` surface turns a published siRNA efficacy panel into a **versioned,
checksummed artifact** that the existing fixed-length design path consumes with no new execution
path, and records what every filter did to every measured observation under **two** policies at once.

Issue [#109](https://github.com/Austin-s-h/sirnaforge/issues/109). Variable-length and asymmetric
design are deliberately **not** here; they are [#110](https://github.com/Austin-s-h/sirnaforge/issues/110).

## `oligogym/` is a synthetic-flank benchmark

:::{warning}
**The 4,113 OligoGym-derived rows make design enumeration and parameter benchmarking runnable. They
do not validate native transcript accessibility.**

Every row of `tests/data/benchmarks/oligogym/records.csv` is `context_type =
"synthetic_neutral_flanks"` with `site_start = 71` — verified on all 4,113. The "transcript" each row
was written against is what `scripts/prepare_oligo_benchmarks.py` fabricated for it:

```text
70 × A  +  reverse-complement(guide)  +  70 × A
```

so the site span these rows carry (1-based 71..91 for the 21 nt guides, 71..90 for the 20 nt ones) is
an exact coordinate **in a fabricated context**. The guide/passenger sequences and the labels are as
measured; the coordinates are not native.

That is why the target-identity vocabulary gained a third value rather than stretching an existing
one: `synthetic_context_local`. `panel_local` would read as a coordinate in the panel's own measured
target, and there is no `confirmed` member at all — `panels.TargetIdentityStatus` has exactly
`unavailable | panel_local | synthetic_context_local`, so #109 structurally cannot emit a native
claim. Native mapping is [#110](https://github.com/Austin-s-h/sirnaforge/issues/110)'s, and #109
reads its Ensembl-mapped counterpart `tests/data/benchmarks/oligogym_native_design/` nowhere at all.

**Read every tally on this page accordingly.** A `default_pass` count over an OligoGym-derived panel
is evidence about design enumeration and filter behaviour. It is not an accessibility result, and
nothing on this surface scores an artifact against the measured efficacy the panel ships.

[`tests/data/benchmarks/README.md`](https://github.com/Austin-s-h/sirnaforge/blob/master/tests/data/benchmarks/README.md)
is the authority for these bytes and draws the same distinction. Read it before reading this page.
:::

## Which panels are present

Five of the six registered panels now ship vendored bytes. Commit `f4beab7` added
`tests/data/benchmarks/oligogym/records.csv` — 4,113 rows, 21 columns, SHA-256
`2fb449362985c52c89ef02d8c7606da641b8cb314d86dac1d7f0cba74611a9c1` — and four descriptors read it,
each selecting its own rows out of the one shared table by `dataset`.

| Panel id         | `data_present` | Bytes                                           | Rows it selects                                         | Declared geometry                            |
| ---------------- | -------------- | ----------------------------------------------- | ------------------------------------------------------- | -------------------------------------------- |
| `huesken_subset` | `true`         | `tests/unit/data/sirna_efficacy_subset.csv`     | 180                                                     | paired core 19 + overhang                    |
| `ichihara`       | `true`         | `tests/data/benchmarks/oligogym/records.csv`    | 2,850 (`ichihara_2007_1` 2,431 + `ichihara_2007_2` 419) | paired core 19 + overhang (21/19 nt)         |
| `martinelli`     | `true`         | same table                                      | 907 (`martinelli_2023_1`)                               | paired core 19 + overhang (21/21 nt)         |
| `shmushkovich`   | `true`         | same table                                      | 356 (`shmushkovich_2018_1`)                             | **asymmetric** (20/15 nt) — all incompatible |
| `oligogym`       | `true`         | same table                                      | all 4,113 — **refused**, `aggregate_of` the three above | mixed; see below                             |
| `huesken_full`   | `false`        | none — `work/sirna_bench.csv` is untracked here | —                                                       | paired core 19 + overhang                    |

`data_present` can no longer be a claim with nothing behind it: `PanelDescriptor` carries
`vendored_csv`, a validator enforces `data_present` ⇔ a named path, and a test asserts every claimed
path exists on disk. 2,850 + 907 + 356 = 4,113, and the three selectors are checked disjoint.

The geometries are re-derived from the bytes rather than from the issue prose, and in one case the
bytes contradicted it: #109 assumed Martinelli was a fully-complementary 21-mer, but of its 907 rows
858 are a 19 nt paired core with a 2 nt 3′ overhang on each strand, 18 are blunt over the full 21 nt,
and 31 match neither exactly. Declaring `FULLY_COMPLEMENTARY`/21 would have written
`guide_3p_overhang: ""` — "measured, and blunt" — for 889 rows whose bytes say otherwise.

`oligogym` is `data_present: true` and still not ingestible: its one table spans three geometries
(3,708 rows on the paired-core-19 geometry; the other 405 are Martinelli's 49 non-core rows plus all
356 asymmetric Shmushkovich rows), so `derive_observation` refuses the aggregate and names its
members rather than stamping one architecture on the file. Prepare a member panel instead.

### What the citations can and cannot show

Every placeholder citation is gone, and each one is now something this repository can be made to
show. What it can show is a DOI per row, a per-dataset row count, and the sequences and labels
themselves. What it cannot: `tests/data/external/oligogym/*.csv.gz`, the upstream extract every row
names in `source_file`, is **not** vendored here (verified absent), and no primary-source full text
has been read here either.

- **`ichihara`** — the Ichihara _et al._ 2007 _Nucleic Acids Research_ citation already vetted in
  `docs/models_and_scoring.md`, now corroborated by the bytes: every row of both datasets records
  `source_url` `https://doi.org/10.1093/nar/gkm699`.
- **`martinelli`** — the DOI its rows carry, `https://doi.org/10.1016/j.ygeno.2024.110815`, and
  nothing more. No author list, title or journal is recorded, because none has been read here. The
  dataset name says 2023 while the DOI string says 2024; **this repository holds nothing that
  resolves which is the publication year**, and the descriptor says so rather than picking one.
- **`shmushkovich`** — the same shape, over `https://doi.org/10.1093/nar/gky745`.
- **`oligogym`** — this repository holds **no citation for OligoGym itself**: no DOI, no URL, no
  upstream extract. What it holds is the adapter that wrote the vendored table,
  `scripts/prepare_oligo_benchmarks.py`, and the three member DOIs above.
- **`huesken_subset`** — unchanged, and still the only panel with a predeclared
  development/held-out split. Its citation, redistribution, SHA-256 and an explicit statement of what
  it does and does not prove are in
  [`tests/unit/data/README.md`](https://github.com/Austin-s-h/sirnaforge/blob/master/tests/unit/data/README.md).
  The Huesken paper is not open access, so those efficacy values were **not** verified against the
  primary source.

A test fails if the word "placeholder" reappears in any citation, endpoint or assay-label constant.

Measured endpoints stay per-panel and are deliberately unpoolable — spelled "percent" so nobody
averages them with `huesken_subset`'s 0–1 `inhibition_fraction`:
`percent_inhibition_relative_to_control_higher_is_better` (Ichihara, observed range −27.8..134.1),
`percent_knockdown_of_target_mrna_higher_is_better` (Martinelli, 0.0..98.0) and
`percent_target_mrna_remaining_lower_is_better` (Shmushkovich, 3.49..119.88). Shmushkovich reads
`label_processed`, its own direction, and not the `efficacy_higher_is_better` column beside it: for
that dataset alone that column is a derived `100 − label` flip, and reading it would record a
transformation as a measurement.

`split_rule_id` is `None` for all four OligoGym-derived panels. `records.csv` has no accession column
— `target` is a gene symbol, and empty on 1,275 of its 4,113 rows — so splitting on it would be a new
rule wearing the audited `sha256_accession_parity_v1` id.

### Shmushkovich's 356 rows are recorded and refused, not designed

All 356 are a 20 nt guide with a 15 nt passenger equal to `reverse-complement(guide[:15])`: a 15 nt
core with a 5 nt single-stranded guide tail. That is the asymmetric hsiRNA geometry #109 excludes by
construction, so the panel is declared `ASYMMETRIC` with `declared_paired_length: None` and every row
comes back `compatibility_status: incompatible` at any requested length, with a
`compatibility_reason` naming the asymmetry and #110.

Nothing is sliced, padded or relabelled to make them fit: `paired_guide_sequence ==
full_guide_sequence` on every row, with `guide_3p_overhang` and `passenger_3p_overhang` both `None`.
They are vendored, and ingested, **to be recorded and refused** — they reach `observations.csv` and
`accounting.csv`, and `design_inputs.fasta` for the panel is empty.

Because the descriptor declares no paired length, `prepare --panel shmushkovich` requires an explicit
`--paired-length`; inventing one that nothing measured is refused.

### What the tests are evidence about

`tests/unit/test_benchmark_panels.py` re-derives every geometry claim in the registry straight from
`records.csv`, and `tests/unit/test_benchmark_real_panels.py` carries those real rows through
`prepare` and `design` and reconciles the tallies against the vendored `manifest.json`. Tests over
the artifact _contract_ still run against `huesken_subset` and the synthetic
`tests/unit/data/benchmark/synthetic_*.csv` fixtures built to reproduce each declared architecture
(documented in
[`tests/unit/data/benchmark/README.md`](https://github.com/Austin-s-h/sirnaforge/blob/master/tests/unit/data/benchmark/README.md)).

A green run over the OligoGym panels is evidence that the surface carries measured rows through
without relabelling them. It is **not** evidence that a design is good, because nothing here compares
a design to the efficacy the panel measured.

The PRD #109 cites as its source, `docs/prd_benchmark_artifacts_and_variable_length.md`, does not
exist in the working tree or in any branch's history. The issue body is the entire specification, and
nothing in this page summarises a document that was read.

### `benchmark prepare` cannot read these four panels yet

:::{warning}
**At this commit the CLI refuses all four OligoGym-derived panels.** The registry knows the bytes are
vendored; `prepare.py` does not.

```console
$ uv run sirnaforge benchmark prepare --panel ichihara --paired-length 19
❌ Error: panel 'ichihara' is registered data_present=True but this module
declares no vendored path for it; that is a bug in prepare.py, not in your invocation
```

Four gaps remain, all in `prepare.py`/`artifact.py` and all outside the registry:

1. `_VENDORED_PANEL_CSV` still names only `huesken_subset`, so no OligoGym panel resolves a source
   table (the error above);
2. `PanelDescriptor.selects_row` is never called, so a reader handed the shared table would ingest
   all four datasets under one panel's architecture;
3. `target_identity_status` is hard-coded to `"unavailable"` with `target_start/end: None`, so the
   `71..91` span `derive_observation` computes is dropped, and `artifact.py`'s
   `TargetIdentityStatus = Literal["unavailable", "panel_local"]` cannot yet hold
   `"synthetic_context_local"`;
4. the aggregate `oligogym` raises a bare `ValueError` **after** `mkdir`, leaving an empty artifact
   directory behind.

**Gaps 1 and 2 must land together.** Wiring the vendored path while still ignoring the selector was
measured: `--panel ichihara` then ingests all 4,113 rows and reports `kept 4113 / incompatible 0`,
tallying `{(21 nt, compatible): 3757, (20 nt, compatible): 356}` — it stamps paired-core-19 on all
356 asymmetric Shmushkovich rows and writes them into `design_inputs.fasta`. That relabelling of a
measured sequence is precisely what #109 exists to prevent.

Gap 3 is lossy but not an overclaim: `unavailable` understates evidence the artifact could honestly
reproduce, and the manifest has no bucket for a `synthetic_context_local` row yet.
:::

## The artifact

One directory per **(compatible panel × paired length)**, named `<panel_id>__len<paired_length>`,
with five fixed inner filenames:

| File                  | Written by | Holds                                                              |
| --------------------- | ---------- | ------------------------------------------------------------------ |
| `observations.csv`    | `prepare`  | One row per **source** row — compatible or not, never dropped      |
| `design_inputs.fasta` | `prepare`  | One record per **compatible** observation; a plain FASTA           |
| `manifest.json`       | both       | Provenance, both policy blocks, counts, checksums                  |
| `candidates_all.csv`  | `design`   | Every **enumerated** candidate, including enumeration-time rejects |
| `accounting.csv`      | `design`   | The observation↔candidate join, with **both** filter verdict sets  |

Paired length is bounded to **19–23 nt** — the range `sirnaforge design` already accepts. An artifact
outside it could never be consumed by the path it exists to feed.

### Nothing measured is dropped, padded, trimmed or relabelled

An observation this build cannot pair to the requested fixed length still gets its
`observations.csv` row, with `compatibility_status: incompatible` and a stated
`compatibility_reason`. An asymmetric duplex is _never_ sliced into a symmetric core — silently
truncating it would misrepresent a measured sequence, which #109 forbids outright. That is no longer
a hypothetical: `shmushkovich`'s 356 vendored 20/15 rows all take this path.
Sequences are copied as measured and upper-cased only; a lower-case letter in a sequence column is
rejected by the schema, because it would mean some caller normalised rather than copied.

Two optional-column conventions are worth knowing before you parse a row. Every optional column uses
the ordinary "empty cell means `None`" convention, **except** `guide_3p_overhang` and
`passenger_3p_overhang`: there, an empty cell means _measured, and blunt_ and the sentinel
`<not_stated>` means _the panel states no overhang_. Both are real, different facts, and CSV cannot
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

`--paired-length` defaults to the panel descriptor's declared length, and is **required** for a panel
that declares none — only `ASYMMETRIC` does, i.e. `shmushkovich`. A panel that ships vendored bytes
**refuses** an explicit `--panel-csv`, so a run cannot silently read different bytes than the ones its
own manifest names as vendored; a panel with no vendored bytes **requires** one. With four more panels
now vendored, that first rule bites four more panel ids than it used to.

### The prepared FASTA is the compatibility guarantee

`design_inputs.fasta` is a plain FASTA keyed by `design_context_id`, so the unmodified fixed-length
path consumes it with no benchmark code in the loop at all:

```bash
uv run sirnaforge design benchmark_artifacts/huesken_subset__len19/design_inputs.fasta --length 19
```

That is the executable form of #109's "consumable by the current fixed-length design/workflow path"
criterion, and it is pinned by a test rather than only documented here.

Each record's `design_context_source` says what context it holds. `huesken_subset` has its
transcripts vendored beside it, so its records are `panel_transcript`. The OligoGym-derived panels do
not, so every record is `measured_target_site` — the 21 nt reverse complement of that row's measured
guide, and nothing wider. The synthetic transcripts in
`tests/data/benchmarks/oligogym/inputs/*.fasta` (161 nt, or 160 for the 20 nt guides) are **never
read** by `prepare` or `design`. A design
tally over those panels is therefore enumeration over a 21 nt window, which is another reason it
cannot be an accessibility result.

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
The default set is _re-derived_ by re-applying the default policy's filter descriptors to the values
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
global default; the manifest carries _both_ resolved policies so the separation is checkable rather
than asserted:

```bash
uv run sirnaforge benchmark design --artifact ... --gc-min 20 --gc-max 75
```

| Manifest field                                              | Value                    |
| ----------------------------------------------------------- | ------------------------ |
| `gc_widening.gc_min` → `default` / `benchmark`              | `30.0` / `20.0`          |
| `gc_widening.gc_min.source`, `.widened`                     | `explicit`, `true`       |
| `run_policy` resolved `gc_min`                              | `20.0` (explicit)        |
| `default_run_policy` resolved `gc_min`                      | `30.0` (builtin_profile) |
| `run_policy.profile.content_hash` vs `default_run_policy`'s | identical                |

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

## What the vendored OligoGym panels actually measure

Every number below was produced by running the real `prepare_artifact`/`design_artifact` code over
`tests/data/benchmarks/oligogym/records.csv` (SHA-256
`2fb449362985c52c89ef02d8c7606da641b8cb314d86dac1d7f0cba74611a9c1`) at paired length 19, under the
shipped default policy with no widening. Because `prepare.py` cannot resolve these panels yet (see
above), the two missing pieces — the vendored path and the row selector — were supplied in-process;
nothing else was substituted, and `tests/unit/test_benchmark_real_panels.py` reaches the same surfaces
the same way, with a probe descriptor.

| Panel / dataset   | in source | kept  | incompatible | entered design | no candidate | default pass | benchmark pass |
| ----------------- | --------- | ----- | ------------ | -------------- | ------------ | ------------ | -------------- |
| `ichihara_2007_1` | 2,431     | 2,431 | 0            | 2,431          | 0            | 1,450        | 1,450          |
| `ichihara_2007_2` | 419       | 419   | 0            | 419            | 0            | 232          | 232            |
| `ichihara` (both) | 2,850     | 2,850 | 0            | 2,850          | 0            | 1,682        | 1,682          |
| `martinelli`      | 907       | 907   | 0            | 907            | 0            | 683          | 683            |
| `shmushkovich`    | 356       | **0** | **356**      | 0              | 356          | 0            | 0              |
| `oligogym`        | 4,113     | —     | —            | —              | —            | —            | —              |

Every count reconciles exactly against the vendored `manifest.json`'s declared 2,431 / 419 / 907 /
356; there is no disagreement anywhere. `oligogym` is refused as an aggregate before a row is read.

`guide_match` is `paired_core_exact` on all 3,757 entered rows — the designer rebuilt each row's own
19 nt core from its own measured site — and `none` on all 356 Shmushkovich rows, whose FASTA is empty.
Martinelli's exclusions are `max_poly_runs` 179, `gc_content_max` 40, `max_paired_fraction` 17,
identical under both policies.

The default verdict set is invariant under a widened benchmark run, which is #109's core accounting
claim. Designing the same Martinelli artifact twice — once unwidened, once with `--gc-min 10
--gc-max 90` — leaves `default_pass` at 683 both times while `benchmark_pass` moves 683 → 711, and
`json.dumps(manifest["default_run_policy"], sort_keys=True)` is identical string-for-string across the
two runs, profile content hash included. The same holds on `ichihara_2007_2`: `default_pass` 232 → 232,
`benchmark_pass` 232 → 289. `gc_content_max` exclusions go `{default 40, benchmark 40}` →
`{default 40, benchmark 0}` while `max_poly_runs` stays `{default 179, benchmark 179}` in both runs.

Measured values survive the round trip **verbatim as strings**, not merely float-equal, on all 4,113
rows: the last `ichihara_2007_2` row's `74.00336740698765` comes out of `observations.csv` as the same
17-significant-digit string. `full_guide_sequence` equals `records.csv`'s `guide_sequence`
upper-cased on every row.

## Reproducibility

The manifest is designed so a second run over the same bytes can be _compared_, not merely trusted:

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
  never coerced (#110). `shmushkovich`'s 356 rows are the vendored proof of it;
- **native transcript mapping**, which is #110's. `counts.mapped_native` is pinned to `0` by the
  schema and `panels.TargetIdentityStatus` has no `confirmed` member, so a native claim is not
  something this build can express, let alone back up. The four OligoGym-derived panels declare
  `synthetic_context_local` and carry a real span in a fabricated context (71..91, or 71..90 for the
  20 nt guides); `prepare.py` still writes `unavailable` with `target_start/end: None` for every row,
  so that span is dropped rather than misreported. `tests/data/benchmarks/oligogym_native_design/` —
  the Ensembl-mapped table, its 11,391-row `mappings.csv`, its per-dataset orthology and its
  representative references — is vendored on this branch and **untouched by #109**;
- **any comparison of a design to the efficacy the panel measured.** Nothing here ranks, correlates or
  scores an artifact against `measured_value`. The measured labels are carried through unchanged so
  that such an evaluation becomes possible; performing one is not part of this surface, and no
  `default_pass` tally on this page is a statement about predictive accuracy;
- **any change to global default GC thresholds** — widening is per run, through this interface only;
- `max_repeat_transcript_fraction` is excluded from both verdict sets, because the only code that
  ever measures it runs from the post-screening step against a real cDNA reference, never from the
  design-from-file path a benchmark run uses. Left in, it reported `unknown` on essentially every row
  and drowned the verdicts this surface exists to report.

See the `benchmark` section of the [CLI reference](cli_reference.md) for the generated option help.
