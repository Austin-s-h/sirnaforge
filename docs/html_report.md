# The HTML report

`sirnaforge report <run_dir>` turns a finished run directory into **one self-contained HTML file**.

```bash
sirnaforge report results/TP53 -o results/TP53/report.html
```

It reads `candidates_all.csv`, the aggregated hit tables and `manifest.json`. It runs no alignment, no
folding and no network access, and it classifies nothing: `hit_class` is read from the table the run
published and every gate comes from the run's own resolved policy, so the report cannot disagree with
the run it describes.

The file works identically over `file://`, from a static web server, and inside the Quilt Catalog with
permissive HTML rendering disabled -- which is the default, and which withholds the
`allow-same-origin` sandbox token. There are no sidecar files, no external URLs, no `fetch`, no
cookies and no browser storage. Reader state lives in memory and in the URL fragment. A test in
`tests/unit/test_reporting_payload.py` enforces this by scanning the rendered document for any
`http(s)://` URL, any non-inline `src`/`href`, and any of `fetch(`, `XMLHttpRequest`, `localStorage`,
`sessionStorage` or `document.cookie` -- it is the report's core constraint and it does not get relaxed
to accommodate a feature.

## The unit is the guide sequence

A guide is enumerated once per transcript it was found on -- a median of 6 and up to 34 rows on the
public 0.7.1 baseline -- while off-target screening runs once per **distinct guide**. Joining evidence
on the candidate `id` would therefore leave most candidates falsely reading as having zero
off-targets, so the guide sequence is the join key and the per-transcript rows collapse into an
isoform sub-table. Two designs within `REGISTER_NEIGHBOUR_NT` (2) nt of each other on one transcript
share a window and can score identically; the isoform table flags them as register neighbours rather
than letting them look like unrelated duplicate rows.

## Five gate verdicts, and four guide statuses -- not the same thing

Each **gate** the report evaluates against one guide's best-scoring row reaches one of five codes:

| Verdict         | Meaning                                                                    |
| --------------- | --------------------------------------------------------------------------- |
| `pass`          | the value was available and met the threshold                              |
| `fail`          | the value was available, the gate's action is not `warn`, and it did not meet the threshold |
| `warn`          | the value was available, the gate's action is `warn`, and it did not meet the threshold -- a real finding, not a rejection |
| `unknown`       | the gate was in force but the run exports no value to compare, or the run itself recorded `UNKNOWN` for this row |
| `not evaluated` | the gate was not applied at all (its action is `off`, it has no threshold, or the run itself recorded `NOT_EVALUATED`), so it makes no claim either way |

Each **guide** then gets one of four statuses, derived from its gates: `fail` if it fails any gate,
else `unknown` if any gate is unknown (or the run rejected the guide by a name no declared gate
expresses -- see `GuideEntry.undeclared_run_rejection`), else `warn` if it exceeds any warn-action
gate, else `pass`. `warn` is kept apart from `fail` because folding a warn-action gate into `fail`
would make the report contradict a run `PASS` it actually agrees with; `unknown` is kept apart from
`pass` because a guide that fails nothing is only clean when every gate could actually be evaluated.
The report may report **less** than the run decided. It may never report more.

The header publishes the agreement between the report's own verdicts and the run's own
`passes_filters` column: how many run-`PASS` guides the report contradicts (should read 0), how many
run-`PASS` guides the report can only call `unknown` because it cannot re-derive them, and how many
run-`fail` guides the report calls clean (`overruled_run_fail` -- must stay 0, the fabricated-evidence
direction). The gate panel on every guide's page shows **every** gate it fails or is warned by, not
only the first the pipeline happened to report.

## Moving the thresholds

Each gate the run actually evaluated on this data gets a control seeded from the run's own resolved
threshold. Both the reader's current value and the run's own value stay on screen, so the two can
never be confused, and resetting returns every displayed verdict to exactly the payload's own.

**No threshold a reader can choose makes an unanswerable gate answerable.** Re-thresholding only ever
*selects among* candidates the run already judged; it does not re-judge them. Concretely: a gate's
triple `[value, verdict, reason]` is re-evaluated in the browser only when its `reason` is `OK` (the
run measured a value and compared it) or `EMPTY_VALUE` (the run measured nothing, which stays
`unknown` at every setting a reader tries). Every other reason -- the gate is off, undeclared, reads a
column the run never exported, or the run itself recorded `UNKNOWN`/`NOT_EVALUATED` for this row --
returns the frozen triple untouched, whatever threshold is asked for. That is the exact property a
reader cannot re-threshold their way out of an `unknown`, and it is derived entirely from data the
report already ships: no new per-filter field, no second copy of the threshold table.

The client-side evaluator restates exactly one thing -- the comparator table behind
`FilterComparator.passes` -- and contains no per-gate branch: it reads each gate from the `FILTERS`
payload the report already carries. A note on the report's own filter-picker UI, since an earlier
draft of this document claimed otherwise: the six reader-facing metric filters (`composite`,
`isoforms`, `gc`, `asym`, `offt`, `liab`) are a **separate, hand-written** row-selection layer over
plain candidate fields such as `gc_content`, `asymmetry_score` and `off_target_count` -- they narrow
which rows are shown, and they are not gate thresholds. Those columns do appear in the template. What
does not appear is a hand-written per-*gate* conditional: the gate panel and its re-thresholding read
every gate generically from the emitted descriptors.

## Preset views

Alongside "all", the index offers views computed from fields the payload already carries, not
re-derived logic: **passing** (`status == pass`), **near-miss** (fails exactly one gate and nothing is
unknown), **off-target-clean** (zero liability hits *and* the guide was actually screened -- a guide
never submitted to the aligner is not "clean", it is unscreened, which is why `off_target_screened`
rides along with `liability_count`), and **register-deduplicated** (one representative per cluster of
guides sharing a window, rather than every register neighbour shown as if independent).

## Registering with Quilt

`sirnaforge report` also writes `quilt_summarize.json` beside the HTML file by default (`--quilt-summarize`;
turn it off with `--no-quilt-summarize`). Quilt's package view renders exactly what a package's own
`quilt_summarize.json` names and nothing else, so a run published as a package without this file shows
no report at all in its package view.

The document is a JSON **array** of rows -- a `file` or an array of `file`, per
[quiltdata/quilt's schema](https://github.com/quiltdata/quilt/blob/main/shared/schemas/quilt_summarize.json)
-- not the dict shape some tooling hands back. `report.html` is always first, full-width and expanded,
titled and described from the run's own gene, guide count, verdict tally and embedded off-target
scope. After it: the passing and full candidate CSVs, the run manifest and workflow summary, the two
aggregated hit tables, and the ORF report and FASTAs -- each included only when it actually exists
under the run directory. A run built with `--skip-off-targets` leaves `off_target/` empty, and
`<GENE>_canonical.fasta` exists only on the gene-search path; either registered unconditionally would
give the catalog a path it can never render, which shows as a broken preview rather than an absent
one. `sirnaforge.reporting.quilt_summarize_entries(payload, report_path, run_dir)` builds the array
without writing it; `write_quilt_summarize(payload, report_path, run_dir)` writes it.

## The URL fragment

The fragment carries the selected guide and any moved thresholds:

```text
report.html#g=UUAUAGGAUUCAACCGGAGGA&t=min_asymmetry_score:0.2,gc_content_max:58
```

A bare fragment (`report.html#UUAUAGGAUUCAACCGGAGGA`) is still read as a guide, which is what earlier
reports wrote. Only re-thresholdable gates are honoured; anything else named in the fragment is listed
as ignored rather than silently dropped.

## Exporting a selection

`Download CSV` and `Download FASTA` write whatever is currently listed -- the search box, the status
toggles and the active preset decide that. Both are generated in the page from the payload already
loaded and handed to the browser as a `Blob`, so nothing is fetched and no sidecar is produced. A
sandboxed frame may withhold downloads entirely, so the same bytes are available as text to copy under
**Download blocked?**.

## Size, and the one figure

On an internal reference run -- 29,605 candidate rows, thousands of guides -- the report is about
**14.8 MB**: candidate rows, embedded off-target evidence at `human, nm<=2`, miRNA seed hits and the
count matrix, plus markup and script under 0.5 MB. Two decisions keep it that small: the filter
descriptors are emitted **once** and each guide carries only `[value, verdict, reason]` codes against
them, and row-level detail is embedded only for liability alignments within the embedded scope.
Everything outside that scope is still counted completely, and a guide whose only hits lie outside it
says so rather than rendering like a genuinely clean one. The index itself caps at `MAX_INDEX_GUIDES`
(5,000) guides, reported in the header so the cap is never silent.

The single figure is a hand-drawn inline SVG. `plotly` was adopted, implemented and then removed: its
4.29 MB bundle carries external URLs and browser-storage references in map traces the report never
invokes, and no static check can tell an inert string from a live call. One stacked bar did not justify
losing the zero-external-reach guarantee.

## Running the client-side parity check

The evaluator that ships inside the report is JavaScript, and the filter it must agree with is Python
(`sirnaforge.reporting.payload.reevaluate_gates`). Rather than a browser -- Playwright is not a
dependency of this repo, and its browser binary is a separate download that would make the test skip
instead of gate -- the parity fixture runs the shipped evaluator directly under Node, which this
project already requires (`.nvmrc` pins Node 20; Node 26 also works). The Jinja placeholder
substitution happens first, in Python, exactly as `render_html` does it, so the harness hands Node the
same JavaScript text a browser would receive rather than an unsubstituted template, and drives it
through `subprocess` the same way the test suite already shells out to other native tools (`cargo`,
`bwa`) rather than reimplementing them. Run it the same way as any other test in this repo:

```bash
uv run pytest tests/unit -k rethreshold -q -n 0
```

The fixture asserts the shipped evaluator and the Python filter agree at the run's own thresholds, and
that no threshold a reader supplies flips a gate whose `reason` is anything other than `OK` or
`EMPTY_VALUE` -- the property that a reader cannot re-threshold their way out of an `unknown`.
