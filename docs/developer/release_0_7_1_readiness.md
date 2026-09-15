# 0.7.1 release readiness

State of `integration/0.7.1` (PR #104) at the point work was deliberately stopped to stabilise the release.
Written so the scope decision and its honest limits survive outside a chat log.

## Verified at `57b38d2`

| check                                                | result                                                      |
| ---------------------------------------------------- | ----------------------------------------------------------- |
| `make test-dev`                                      | 1785 passed                                                 |
| `make lint`                                          | clean (ruff check, ruff format, mypy over 100 source files) |
| `SIRNAFORGE_DOCS_OFFLINE=1 make docs`                | exit 0 with `-W`, zero warnings                             |
| `tests/unit/test_no_private_programme_references.py` | 25 passed                                                   |

## `make test-release` end to end, at `41cfacd`

First full run of the target after #107 unblocked its `docs` step and after the tier fix pass
(`d126938`, `a57aaa3`). One invocation, exit 0, **9m38s** total on an idle-ish 10-core arm64 laptop
(load average 13 from a concurrent session in the same checkout). The image was stale, so this run
paid a rebuild; subtract it for the steady-state figure.

| stage                    | result | tests                                                                               | wall clock |
| ------------------------ | ------ | ----------------------------------------------------------------------------------- | ---------- |
| `docs`                   | pass   | `sphinx-build -W --keep-going`, 0 warnings                                          | 19s        |
| `test-release-host`      | pass   | 1822 selected: **1815 passed, 7 skipped**, 0 failed (46 container tests deselected) | 82.7s      |
| `docker-ensure` rebuild  | built  | fingerprint `d04d6683` → `496c3f66`                                                 | 3m11s      |
| `test-release-container` | pass   | 46 selected: **38 passed, 8 skipped**, 0 failed                                     | 4m39s      |
| `test-release-report`    | pass   | coverage rollup                                                                     | 6s         |

**All 15 skips are one reason**, reported by `tests/conftest.py:214`: `requires_network`, no verified
TLS route to `rest.ensembl.org`. Nothing skipped silently and nothing skipped for a reason the log
does not name.

| tier / marker                       | at task baseline | this run | why it moved                                      |
| ----------------------------------- | ---------------- | -------- | ------------------------------------------------- |
| `dev`                               | 1785             | 1788     |                                                   |
| `ci`                                | 85               | 86       |                                                   |
| `release`                           | 80               | 81       |                                                   |
| `release and not runs_in_container` | 37               | 37       |                                                   |
| `runs_in_container`                 | 45               | 46       | the Nextflow stub-run test                        |
| `requires_docker`                   | 0                | 0        | declared, deliberately unused                     |
| `requires_nextflow`                 | 1                | 2        | the stub run replaced a test that ran no Nextflow |
| `requires_tools`                    | 4                | 4        | all container-tier; all 4 **passed** in the image |
| `requires_network`                  | 13               | 15       | the two vacuous miRBase tests are now gated       |
| `slow`                              | 11               | 11       |                                                   |
| total                               | 1864             | 1868     |                                                   |

The two population changes that matter are the ones the fix pass aimed at, and both are visible here
rather than asserted: `requires_network` 13 → 15 is the two miRBase tests that used to pass after 108s
each having asserted nothing, and `requires_nextflow` 1 → 2 is a real `nextflow run -stub-run`
replacing a test that invoked a CLI option which does not exist. The host stage fell from 297s to 83s
as a direct consequence.

### Coverage, and what the number is worth

`TOTAL 16749 statements, 2866 missed, 83%`. Everything 0.7.1 added is genuinely exercised, not merely
imported:

| module                       | cover |
| ---------------------------- | ----- |
| `core/filtering.py`          | 100%  |
| `core/selection.py`          | 100%  |
| `models/evidence.py`         | 100%  |
| `models/policy.py`           | 100%  |
| `reporting/render.py`        | 100%  |
| `core/screening_evidence.py` | 99%   |
| `benchmark/design.py`        | 99%   |
| `reporting/quilt.py`         | 98%   |
| `benchmark/artifact.py`      | 97%   |
| `benchmark/panels.py`        | 97%   |
| `reporting/tracks.py`        | 97%   |
| `data/orthology.py`          | 96%   |
| `reporting/payload.py`       | 96%   |
| `benchmark/prepare.py`       | 93%   |
| `core/design.py`             | 91%   |
| `provenance.py`              | 90%   |
| `workflow.py`                | 87%   |
| `core/off_target.py`         | 75%   |

`off_target.py` is the one weak module and it is the pre-existing screening engine, not new code. Its
misses are one contiguous region — the `bwa-mem2`/`samtools` invocation and result-parsing paths — which
is exactly what the host cannot supply and what the container tier drives through `subprocess.run`,
where coverage cannot follow it. So 75% understates execution and 83% overstates validation.

⚠️ **The 83% is only readable at the moment the rollup runs.** Re-reading the _identical, untouched_
`.coverage` database twelve minutes later gave **78%, with `workflow.py` at 56%** — because a
concurrent session had edited `workflow.py` in between and `coverage report` re-parses source at report
time. That is the mechanism behind the drift recorded in `a57aaa3`: not a second process rewriting
`.coverage`, but source moving under a database that stores bare line numbers. Only `workflow.py`
drifted; every other module read identically on both passes.

## Size, and why that mattered to the scope decision

|                                      | files | insertions |
| ------------------------------------ | ----- | ---------- |
| 0.7.1 vs `dev`, all                  | 180   | 105,808    |
| — of which vendored benchmark panels | 26    | 52,779     |
| — `src/` code                        | 57    | 22,159     |
| 0.6.0 (`dev` → `master`), for scale  | 34    | 2,928      |

182 commits, and roughly **7.5× the previous release's `src` change**. That is the reason work stopped here
rather than continuing through the remaining issues: a release nobody can review is not safer for being
larger.

## Closed in this release

| issue | what landed                                                                                                                                         |
| ----- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| #100  | the screening evidence contract: plan/evidence producers, reconciliation, pure filtering and selection modules, plan threading through the pipeline |
| #103  | the single-candidate HTML report's remainder: client-side re-thresholding, preset views, `quilt_summarize.json`                                     |
| #105  | the isoform-coverage gate obeys its resolved action                                                                                                 |
| #106  | a completed zero-hit screen records measured passes rather than `not_evaluated`                                                                     |
| #107  | the docs release gate: 39 duplicate objects, then an unreachable inventory failing `-W`                                                             |
| #108  | the miRNA detail pointer names the aggregate it actually ingested                                                                                   |
| #109  | fixed-length benchmark artifacts and the `benchmark prepare`/`design` commands                                                                      |

Plus, not issue-numbered: the ten defects from a live in-browser audit of the report; removal of every
test dependency on `node`; and removal of an internal programme identifier from this public repository.

## Deferred out of 0.7.1

| issue     | state                                                                                                                                             |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| #101      | asymmetry floor **merged**; target intent and the transcript-seed channel built on `feat/issue-101-intent` at `8f1c923` and **RED (11 failures)** |
| #102      | not started — deterministic tied-score order and the unique-guide shortlist policy, both belonging in `core/selection.py`                         |
| #110      | not started; unblocked by #109, see its issue comment for two findings that change its scope                                                      |
| #111–#114 | not started — a chained report refactor, deliberately not stacked on a release that already absorbed a full report rewrite                        |

`feat/issue-101-intent` is pushed and its WIP commits name every failing test and why, so the checkpoint
cannot be mistaken for working code. Its 11 failures are not incidental: nine are the new gates declared in
the filter registry without discharging what that registry obliges — the counter on the row, an
independently clearable threshold, a verdict the report can re-derive — which is #101's own headline finding
turned back on itself.

## Honest limits of what is verified

Stated because a passing suite is not the same as a validated release.

- **Nothing was verified at production scale.** The tracked report fixture is 292 candidate rows / 28
  guides, about 0.8% of one internal run, and human-only. The guide cap was exercised by lowering the limit
  to 12 rather than by reaching 5,000, and the per-species liability split rendered with one species.
- **The container tier now runs, but only on toy references.** As of the run above it is green (38
  passed), and `bwa-mem2` really executes there — but against
  `tests/unit/data/toy_transcriptome_db.fasta`, not a full transcriptome. #100's headline magnitudes
  (`TRANSCRIPTOME_PERFECT_MATCH` 31,879 → 476, `REPEAT_ELEMENT` 0 → 10,460) are still **pre-fix** numbers
  and were not re-measured; direction is confirmed by execution, magnitude is not.
- **The end-to-end path from a gene symbol is the part that skipped.** All 8 container skips are
  `requires_network`: `test_docker_full_tp53_workflow`, the five `test_variant_integration` SNP modes,
  `test_full_workflow_with_gene_search` and `test_custom_transcriptome_offtarget`. Those are the most
  realistic tests in the suite, and behind a TLS-intercepting proxy none of them runs — so a green
  container tier does **not** say that gene name → Ensembl → transcripts → screen works.
- **`make docs` passes `-W` without resolving cross-project references.** On this host Sphinx dropped all
  four intersphinx inventories (python, numpy, pandas, biopython) to `CERTIFICATE_VERIFY_FAILED` and still
  exited 0, so a broken reference into another project's API is not caught locally.
- **No test executes the report's JavaScript.** The `node` parity harness was deleted by decision. A
  divergence between `payload.reevaluate_gates` and the template's evaluator would not fail the suite.
  `docs/html_report.md` carries a 22-property list of exactly what is no longer checked, mapped to the
  deleted test each property came from.
- **Browser-only report properties are asserted at the markup or source level, not driven**: layout fit,
  print output, refusal-banner visibility, and the per-keystroke timing.
- **`UNDETERMINED` and `undetermined_hits` remain largely unexercised on real data.** They only become
  non-trivial when a species has no index.
- **The internal identifier is removed from the working tree but remains in pushed history** — 12 commits,
  11 reachable from `origin/integration/0.7.1`, plus one on `origin/master`. Left by decision; a
  digest-based guard prevents recurrence.

## Known fragility, unrelated to any issue

**Fixed in `d126938`.** `tests/unit/test_repeat_detection.py::test_performance_guard_2mb_reference` used
to assert a 2.0 s _wall-clock_ ceiling, and failed at 2.845 s under concurrent load while passing at
0.96 s alone. It now measures `time.process_time()` against the same ceiling; `scan` is single-threaded
numpy, so the ceiling still discriminates (0.9 s now against 6.57 s pre-rewrite, both CPU-bound) but
scheduler contention cannot fail it. It passed in the run above, and separately at 0.92 s idle and
0.99 s under 12 competing CPU-bound processes.

It is worth knowing that this is the **only** timing assertion left in the suite — one `assert elapsed
< 2.0` across all 1868 tests. Everything else that could have been a stopwatch is a behavioural
assertion, which is why the tier is not load-sensitive.

## One thing worth promoting

The report audit built a 49-assertion verification harness that needs no JS engine — the shape asked for
after `node` was removed — but it is uncommitted scratch. Promoting it would restore real automated
coverage over the rendered document; the only cost is assembling the fixture run.
