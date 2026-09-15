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
- **The container tier was not re-run in this pass.** `bwa-mem2` is Docker-only on arm64, so no
  full-reference screen was executed. #100's headline magnitudes (`TRANSCRIPTOME_PERFECT_MATCH`
  31,879 → 476, `REPEAT_ELEMENT` 0 → 10,460) are **pre-fix** numbers and were not re-measured; direction is
  confirmed by execution, magnitude is not.
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

`tests/unit/test_repeat_detection.py::test_performance_guard_2mb_reference` asserts a 2.0 s wall-clock
ceiling. It failed at 2.845 s under concurrent load and passed at 0.96 s alone, so it will produce false
reds on a loaded machine or in CI.

## One thing worth promoting

The report audit built a 49-assertion verification harness that needs no JS engine — the shape asked for
after `node` was removed — but it is uncommitted scratch. Promoting it would restore real automated
coverage over the rendered document; the only cost is assembling the fixture run.
