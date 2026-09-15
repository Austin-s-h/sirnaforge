# The HTML report

`sirnaforge report <run_dir>` turns a finished run directory into **one self-contained HTML file**.

```bash
sirnaforge report results/TP53 -o results/TP53/report.html
```

It reads `candidates_all.csv`, the aggregated hit tables and `manifest.json`. It runs no alignment, no
folding and no network access, and it classifies nothing: `hit_class` is read from the table the run
published and every gate comes from the run's own resolved policy, so the report cannot disagree with
the run it describes. A hit table without the classification columns is **refused** rather than
rendered.

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

| Verdict         | Meaning                                                                                                                                                      |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `pass`          | the value was available and met the threshold                                                                                                                |
| `fail`          | the value was available, the gate's action is not `warn`, and it did not meet the threshold                                                                  |
| `warn`          | the value was available, the gate's action is `warn`, and it did not meet the threshold -- a real finding, not a rejection                                   |
| `unknown`       | the gate **was in force** and reached no decision: the run exports no value to compare, or the run itself recorded `UNKNOWN` or `NOT_EVALUATED` for this row |
| `not evaluated` | the gate was never in force -- its action is `off`, or it declares no threshold -- so it makes no claim either way                                           |

"Not applied" and "applied, evidence unavailable" are different claims, and only the second one bears
on whether a guide is clean, so they get different codes. A gate that is in force with a threshold and
that the **run itself** left undecided is therefore `unknown`, not `not evaluated`; the per-row reason
code is what keeps the two non-decisions apart, and the gate panel's Why column is what prints it
("the run recorded unknown for this gate; no value to re-compare" against "the run did not evaluate
this gate; value measured, no verdict applied"). `unknown` still nulls the value; the run's own
`not_evaluated` keeps the number the run recorded and applies no verdict to it.

Each **guide** then gets one of four statuses, derived from its gates: `fail` if it fails any gate,
else `unknown` if any gate is unknown (or the run rejected the guide by a name no declared gate
expresses -- see `GuideEntry.undeclared_run_rejection`), else `warn` if it exceeds any warn-action
gate, else `pass`. `warn` is kept apart from `fail` because folding a warn-action gate into `fail`
would make the report contradict a run `PASS` it actually agrees with; `unknown` is kept apart from
`pass` because a guide that fails nothing is only clean when every gate could actually be evaluated.
The report may report **less** than the run decided. It may never report more.

Because an in-force gate the run left undecided is `unknown`, **a guide carrying one is no longer
published as `pass`** -- it is `unknown`, and its panel carries the undecided-gate banner. Earlier
builds counted only the run's `UNKNOWN` towards the guide's unknown tally, so a gate the run recorded
`not_evaluated` was neither failed nor unknown and the status fell through to `pass`: a run whose
repeat scan never executed leaves `repeat_transcript_fraction` at its `0.0` default with a
`not_evaluated` verdict, from which the report published a confident `PASS` for a gate nobody had
evaluated. Runs missing evidence will therefore show **more `unknown` guides and fewer passes** than
they did before, on the same input.

The header publishes the agreement between the report's own verdicts and the run's own
`passes_filters` column: how many run-`PASS` guides the report contradicts (should read 0), how many
run-`PASS` guides the report can only call `unknown` because it cannot re-derive them, and how many
run-`fail` guides the report calls clean (`overruled_run_fail` -- must stay 0, the fabricated-evidence
direction). The gate panel on every guide's page shows **every** gate it fails or is warned by, not
only the first the pipeline happened to report.

## Moving the thresholds

Every gate this run actually evaluated on this data gets a **number box and a slider**, both seeded
from the run's own resolved threshold, with that threshold shown beside them as `run <value>` so
"yours" and "the run's" can never be confused. `Reset` restores every one of them, and the displayed
verdicts return to exactly the payload's own. Moving one re-derives that guide's gate triple, its
status pill, its gate panel, the index counts and the export -- all of them read one recomputation
(`recomputeLive`), so the pill and the panel underneath it cannot disagree.

The slider's bounds come from **the values this run produced** for that gate, not from a per-filter
display range invented in the template, and the run's own threshold is always reachable on it, so a
control can always be moved back to where the run left it. The bounds are also **snapped to the
slider's own step and clamped to the range the setting's own field declares** -- read off the
`FilterCriteria`/`OffTargetFilterCriteria` metadata rather than a hand-written table. Padding the
observed values by 5% instead published "at most -59 off-targets" and a GC floor of 27.804348 against a
step of 0.1, so 40 was unreachable and 39.704348 is what reached the URL and the cart TSV. The clamp is
also what makes `min_empirical_score` useful: its domain is 0.4-0.6, the two thresholds that answer
whether the gate does anything, rather than 0.39-0.61 either side of them.

The number box is not bounded by the slider, so the slider's step limits its resolution and never the
reachable thresholds -- which means the two can hold different numbers, and a range input cannot
represent that: it pins silently at its own edge. Typing 90 into a gate whose domain ends at 76.104348
used to leave the box at 90 and the slider at the pin, and one touch of the slider then dropped the real
threshold to 76.104348 without the reader asking. The domain is deliberately not widened -- the clamping
is what stops a count-valued gate offering a negative ceiling. Instead the threshold stays where it was
typed, the slider **says it cannot show it** and is disabled, so it cannot answer for a number it does
not hold.

### Whose threshold produced the verdict you are reading

A verdict re-derived at a reader's threshold is not the run's verdict, and the report says so in three
places at once rather than leaving the reader to infer it:

- the gate panel opens with `N thresholds moved from the run's` whenever any control has been touched,
  so a reader quoting a pill can see whose it is without reading every row;
- the row's Threshold cell reads `le 20 (yours; run 10)`; and
- the Why sentence states the comparison **that produced the verdict beside it**, then names the run's
  own: `29 le 30 at your threshold; the run used le 15`.

All three read one function (`liveThreshold`), which is the point: the Why column used to build its
sentence from the run's threshold while the pill beside it came from the live re-thresholded table, so
any reader move made the sentence arithmetically false in both directions -- `10 not le 10 so FAIL`
next to a `pass` pill, or a satisfied comparison next to a `fail`. A row whose stated arithmetic does
not produce its own verdict is publishing a verdict the pipeline would not reproduce.

The panel also carries `N of M gates undecided` when a guide has any -- the count is over gates the
report calls `unknown`, which after the vocabulary above includes the run's own `not_evaluated` on an
in-force gate, and deliberately excludes the rows whose pills read `not evaluated`.

### A number box refuses a non-number, visibly

Both kinds of box -- gate thresholds and reader-filter bounds -- are `type="number" step="any"` and go
through one boundary, `parseControlValue`: blank is _unset_, and anything non-finite (`abc`, `5px`,
`Infinity`, `1e999`, `NaN`) is **refused**. A refusal paints the control (`aria-invalid`, a red box) and
names it in a banner, so nothing is silently ignored. There is no `min`/`max` on the box on purpose:
the slider's domain bounds the slider, and the box is what keeps every threshold reachable, so what is
refused is not an unusual number but a non-number.

This matters because `Number('abc')` is `NaN` and every comparison against `NaN` is false: one junk
keystroke used to fail every re-decidable gate on every guide. `encodeHash` also filters on
`Number.isFinite`, so the fragment can never carry a value `applyHashFragment` would itself refuse --
including an empty one (`t=some_gate:`), which `Number('')` had been reading as a real threshold of `0`.

### A gate the run could not evaluate stays unevaluable

**No threshold a reader can choose makes an unanswerable gate answerable.** Re-thresholding only ever
_selects among_ candidates the run already judged; it does not re-judge them. A gate gets a control
only when all four of these hold:

1. the run applied it (its action is not `off`),
2. it has a declared threshold,
3. the run exports the column the report reads for it (`<filter_id>_observed`, or the descriptor's own
   column -- whichever `observed_column` resolved), and
4. the run **decided** it for at least one guide in this run.

The fourth condition used to read "at least one guide has a value in that column", which is a weaker
thing. A gate every row of which the run recorded `not_evaluated` still exports the number it measured,
so that test called it evaluable and the report shipped a live slider that could not change anything: a
reader moves it, every row stays frozen on its own reason code, nothing happens and nothing says why.
`evaluable` now means exactly one thing -- **moving this control can change a verdict** -- and such a
gate is published frozen with the honest reason instead. `n_values` still counts the measurements,
because they were measured.

Anything else is **frozen**: no control is rendered, the row states payload.py's own
`unevaluable_reason` rather than a re-derived guess at it, the evaluator returns the run's verdict
untouched whatever it is asked for, and a threshold arriving in the URL fragment for a frozen gate is
refused with the refusal shown on screen. On the tracked `baseline_0_7_1` slice, **7 of the 17
declared gates can be moved** (`gc_content_min`, `gc_content_max`, `max_repeat_transcript_fraction`,
`max_paired_fraction`, `min_asymmetry_score`, `min_empirical_score`, `max_off_target_count`) and 10
are frozen -- 6 for a column the run does not export and 4 because the run has the filter off.

Freezing is decided **per row**, not per filter, and it is keyed on that row's own `reason` code. A
gate's triple `[value, verdict, reason]` is re-evaluated only when its `reason` is `OK` (the run
measured a value and compared it) or `EMPTY_VALUE` (the run measured nothing, which stays `unknown` at
every setting a reader tries). Every other reason -- the gate is off, has no threshold, reads a column
the run never exported, or the run itself recorded `UNKNOWN`/`NOT_EVALUATED` _for this guide_ --
returns the frozen triple untouched. The per-row granularity matters: a filter can export its column
for most guides and still leave one row `NOT_EVALUATED` carrying a real measured number, and that
number is not a reader's to re-judge. This is the property that a reader cannot re-threshold their way
out of an `unknown`.

### The evaluator is generic

The client-side evaluator restates exactly one thing -- the comparator table behind
`FilterComparator.passes` -- and contains **no per-gate branch**: every gate, its comparator, its
action, its threshold and its control domain are read from the `FILTERS` payload the report already
carries. A unit test asserts that none of the 17 declared `filter_id`s appears as a literal string in
the template at all. (Column names such as `gc_content` do appear, in the cart's TSV header.)

The single Python implementation the browser has to agree with is
`sirnaforge.reporting.payload.reevaluate_gates`, and "Running the client-side parity check" below is
how that agreement is enforced.

### Seventeen sliders, ranked by the two that decided the run

Seventeen gates rendered as seventeen equal-looking controls is a panel with no shape: on one internal
run five of them rejected anything at all, two accounted for 1,765 rejections, and two could not reject
anybody because their threshold sat on the limit their own setting declares. So the payload publishes
what each gate did over the guides the report holds -- `rejects`, `sole_rejects`, `warns`, `unknowns` --
and the panel is a ranked list, ordered by the guides each gate is **solely** responsible for rejecting,
each row saying so (`9 rejected, 9 by this gate alone`). A gate that decided nothing is folded behind a
summary and labelled **inert in this run**, beside the payload's own reason why -- including the case
worth naming, a threshold sitting on the lowest value its own field permits, where "no permitted value
can fail it" is a property of the gate rather than of this run's data.

## Liability, and the species it is counted in

`max_off_target_count` compares the total over **every screened species**, and roughly half of that total
is non-human on a four-species screen (92,658 human against 91,536 mouse, 25,243 macaque and 23,427 rat
on one internal run). A guide that gate rejected may therefore carry no human liability at all, and the
number on screen could not be decomposed. The index cell now shows the total and its split into the query
species and the rest, the header states the run-wide split, and **one** selector re-scopes the column, its
sort and the liability bound together -- one control rather than two over the same number, which is how a
reader ends up believing the report says two things. The bound's own label names the species it counts,
because "at most 3 liabilities" means something different either side of that selector, and a species this
run did not screen is refused rather than adopted.

The selector re-scopes a **question**, never a verdict. Applying `max_off_target_count` to human-only
liabilities would be the report deciding something the run did not, so it is not offered: the payload
publishes `liability_by_species` per guide and `liability_rows_by_species` run-wide, and re-thresholding
still moves the all-species ceiling the run itself applied.

## Reader filters, which are not gates

Below the gate controls, in its own labelled group, sit three plain boxes over quantities **no gate
covers**:

| Filter            | Reads             | Direction |
| ----------------- | ----------------- | --------- |
| `composite score` | `composite_score` | at least  |
| `isoforms hit`    | `transcript_hits` | at least  |
| `liabilities`     | `liability_count` | at most   |

The distinction is the point, and the panel states it on both groups: **a gate control changes a
verdict the run computed; a reader filter only selects among rows.** Moving a reader filter recomputes
nothing -- no gate triple, no status pill, no counter -- it just narrows the list. None of the three is
relabelled as a gate, and none of them is reachable any other way: `composite_score` and
`transcript_hits` are gated by nothing at all (`min_isoform_coverage` is a _fraction of the run's
transcripts_, a different quantity), and `liability_count` is otherwise only expressible as _exactly
zero_, through the off-target-clean preset -- so "at most 3 liabilities" had no control until these
came back. The other three boxes the v1 report carried -- GC, asymmetry and off-target count -- are
deliberately **not** restored: all three are real movable gates now, and a second control over the same
number obeying a different rule is how a reader ends up believing the report says two things.

They obey the same rules as everything else in the panel. A blank box is _unset_, never `0`. **An
absent value cannot satisfy a threshold**: a guide with no composite score is excluded by a composite
floor rather than admitted because there was nothing to compare, exactly as an unevaluable gate is
never a pass. They compose (AND) with the presets, the status and conservation checkboxes, the search
box and the sort; `Reset` clears them along with every gate threshold and the preset; and they ride in
the URL fragment on the same footing as the thresholds, so a copied link restores the same rows.

## Preset views

Alongside "all", the index offers four views computed from fields the payload already carries, not
from logic re-derived in the browser:

- **Passing** -- live status `pass` or `warn`. It can never surface a guide the evaluator itself calls
  `unknown`, at any threshold.
- **Near miss** -- fails exactly one gate and is unknown on none. An unknown guide is never a near
  miss, because the distance to a pass is not known.
- **Off-target clean** -- zero liability hits **and** the guide was actually screened. A guide never
  submitted to the aligner is not clean, it is unscreened, and the two look like the same `0` unless
  `off_target_screened` rides along with `liability_count`; it does.
- **Register-deduplicated** -- one representative per cluster of guides sharing a window, rather than
  every register neighbour listed as if independent. The clusters are connected components over
  transcript positions within `REGISTER_NEIGHBOUR_NT`, computed once in Python over the whole
  candidate table (`register_cluster`, `register_representative`), keyed on each candidate's own
  score, so the keeper is the cluster's best-scoring member.

## Registering with Quilt

`sirnaforge report` also writes `quilt_summarize.json` beside the HTML file by default (`--quilt-summarize`;
turn it off with `--no-quilt-summarize`). Quilt's package view renders exactly what a package's own
`quilt_summarize.json` names and nothing else, so a run published as a package without this file shows
no report at all in its package view.

The document is a JSON **array** of rows -- a `file` or an array of `file`, per
[quiltdata/quilt's schema](https://github.com/quiltdata/quilt/blob/main/shared/schemas/quilt_summarize.json)
-- not the dict shape some tooling hands back. `report.html` is always first and **alone in its row**,
which is how a summarize file spans the package view's full width, and it is `expand`ed, titled and
described from the run's own gene, guide count, verdict tally, embedded off-target scope and
report-vs-run agreement. After it: the passing and full candidate CSVs, the run manifest and workflow
summary, the two aggregated hit tables, and the ORF report and FASTAs -- each included only when it
actually exists under the run directory. A run built with `--skip-off-targets` leaves `off_target/`
empty, and `<GENE>_canonical.fasta` exists only on the gene-search path; either registered
unconditionally would give the catalog a path it can never render, which shows as a broken preview
rather than an absent one. `sirnaforge.reporting.quilt_summarize_entries(payload, report_path, run_dir)`
builds the array without writing it; `write_quilt_summarize(payload, report_path, run_dir)` writes it.

A full `sirnaforge workflow` run writes both files itself, and puts the summarize file at the **run
root** rather than beside the report: the report lands in `sirnaforge/report.html`, but the run
directory is what gets published as a package, and Quilt reads a summarize file only at the package
root. Every path in the document is therefore relative to the summarize file's own directory, not the
report's -- `write_quilt_summarize(..., out_path=...)` is what tells the two apart. The `.nf` pipeline
itself renders nothing; a Nextflow-only run still needs `sirnaforge report` afterwards.

The workflow **re-issues** the registration once `logs/workflow_summary.json` has landed. Step 6 runs
before the summary is written, and `write_quilt_summarize` omits any artifact that does not exist yet, so
the run's own summary could never be registered by the step-6 write -- the one row a reader would go
looking for was the one row missing. Re-issued rather than reordered, because the summary's
`processing_time` has to keep measuring the whole run, and gated on the summary existing and the report
having rendered, so a `write_json_summary=False` run or a failed summary write registers nothing
dangling. Rendering the report and registering it are also separate `try`/`except` regions with separate
log messages now: one failure used to be reported as the other, and the step's console artifact list
announced files it had not written -- `candidates_pass.fasta` is deliberately deleted when nothing
passes, and was listed anyway. Every line in that list is existence-checked.

## The URL fragment

The fragment is the whole view, not a part of it. Nine keys, written in one canonical order and each
omitted at its default, so `Reset` leaves an empty fragment: the selected guide (`g=`), moved thresholds
(`t=`), reader-filter bounds (`r=`), the liability scope (`sp=`), the active preset (`preset=`), the
search box (`q=`), the status checkboxes (`s=`), the conservation species and seed-tolerance toggle
(`c=`) and the cart (`k=`).

```text
report.html#g=UUAUAGGAUUCAACCGGAGGA&t=min_asymmetry_score:0.2,gc_content_max:58&r=composite:60,liab:3&sp=human&preset=near_miss&q=CCGGA&s=pass,warn&k=UUAUAGGAUUCAACCGGAGGA
```

The four state keys were added because half the visible state used to stay on screen. The default status
set hides every `fail` and `unknown` row -- 3,617 of them on one internal run -- so a URL carrying moved
thresholds without `s=` carried the verdicts and not the rows they were chosen to show, and the search box
never wrote to the fragment at all. `s=` with an **empty** value is meaningful and is not a missing key: it
is a reader who unchecked all four statuses, which is a real view of every guide. Sets **replace** rather
than union, because a re-applied fragment must not merge with whatever happened to be on screen, and
`applyFragmentState` resets the controls first for the same reason -- what a URL omits is a default, not
the last thing the page was showing.

The cart travels as guide **sequences**, never as indices: an index would resolve to a different guide in a
differently built report, whereas a sequence this file does not hold can be refused. Browser storage is not
an option -- the report reaches nothing outside itself, and the sandbox it is usually read in withholds
same-origin anyway -- so the address is the cart, which the cart card says on screen. The first
`CART_IN_URL_MAX` (500) guides in the index's own order travel; above that the card says the address
carries fewer guides than the cart holds, and the TSV export is the way out.

A bare fragment (`report.html#UUAUAGGAUUCAACCGGAGGA`) is still read as a guide, which is what earlier
reports wrote. Only re-thresholdable gates are honoured; a frozen or unknown `filter_id`, an `r=` name that
is not one of the three reader filters, a species this run did not screen, or a `k=` sequence this file does
not hold is **listed as refused** in the panel rather than silently dropped, and an unknown preset name
leaves the view at "all". Because that refusal banner lives inside the thresholds `<details>`, which ships
closed, anything the reader has to see opens the panel -- a refused name, a refused keystroke, or a fragment
that set any control at all.

Writing the fragment uses `history.replaceState`, not `location.hash =`: the latter is a navigation, so it
pushed one history entry per keystroke on every threshold, Back then changed the URL without changing the
view, and leaving the report took dozens of presses. Reading it is not a load-time-only affair either -- a
`hashchange` listener means Back, Forward, an edited address bar and a pasted link all reach the code a
fresh load runs. A reader who reloads the URL sees the same rows: `encodeHash` and `applyHashFragment`
round-trip byte for byte, and the set of keys the encoder writes is asserted equal to the set the parser
reads back, derived from the source rather than listed, so the next control encoded and never parsed fails
a test instead of shipping a URL that restores less than it carries.

## Exporting a selection

The cart is an explicit pick list: click the `+` column in the index, or `Add top <n> to cart` to take
the highest-scoring guides currently listed. `Export TSV` downloads exactly what the cart holds --
guide, passenger, live status, both scores, isoform coverage, the metric columns the thresholds were
applied to, the structure string and per-species conservation -- generated in the page from the
payload already loaded and handed to the browser as a `Blob`, so nothing is fetched and no sidecar is
produced. A sandboxed frame may withhold downloads entirely (a Quilt iframe without `allow-downloads`
does), and **nothing in the page can see that**: the embedder refuses the synthetic `a.click()` itself,
with no exception, no callback and nothing observable, so a fallback conditional on catching the refusal
never runs and the reader was told the file had been written. The claim is therefore not made. The same
bytes sit in a text box above the button on **every** export, one keystroke from the data, and the note
says both outcomes rather than picking one. "Download refused by this viewer" survives only in the branch
where a refusal really was observed -- a throw out of `Blob` or `createObjectURL`.

The `status` column is the reader's **live** status, so the file opens with `#`-prefixed, tab-delimited
`key=value` comment lines saying which thresholds produced it. Without them the export was a table of
verdicts with no way to tell whose they were once it left the page:

```text
#sirnaforge_cart	schema=1	gene=TP53	guides=2	preset=near_miss
#status_basis=reader_rethresholded	moved_gates=1
#moved_gate=max_off_target_count	comparator=le	run_threshold=15	reader_threshold=30
#reader_filter=composite	direction=min	bound=30
```

`status_basis` is `run_thresholds` when nothing has been moved and `reader_rethresholded` otherwise,
with one `#moved_gate=` line per moved gate. The `#reader_filter=` lines are recorded for a different
reason: a reader filter can never touch the `status` column, but `Add top <n> to cart` picks from the
_filtered_ view, so the bounds decide which guides are in the file. The comment says exactly that, so
the two kinds of provenance are not confused. The columns and their order are unchanged, and are frozen
in the parity harness.

## The candidate map re-colours with the reader

The design map's backdrop -- axes, gridlines, the CDS/UTR region bar -- is drawn in Python, because one
function draws both a notebook figure and this card. Its **points and its legend are not**: they used to
be, and a threshold a reader moves cannot reach into a server-rendered SVG. Setting `gc_content_min` to 60
on one internal run made the live counts 23 pass / 16 warn / 4,961 fail while the legend still read
`PASS, all gates evaluated (394)`, byte-identical, over 394 green dots -- this report's own failure mode,
in the largest visual on the page.

So the windows ride as numbers with the guide behind each one, and both the dots and the legend are
painted from each guide's live gate table, in draw order, with a class carrying no window on the current
transcript omitted rather than offered as a class the reader failed to find. A window whose guide the cap
dropped has no live verdict to read, so it keeps its own key and says whose verdict that colour is. When
any threshold has moved the legend stamps itself -- "Re-coloured at your thresholds" -- and says that no
colour or count beside it is the run's own. The isoform picker is ordered canonical-first (by length after
that), since window count is the wrong first key: it opened one report on a 7,988 nt transcript with the
canonical one tenth in the list. Canonical status comes from the transcript FASTA the run writes and is
**unknown**, never `false`, when the run recorded none. A reader who picks an isoform keeps it: the
selection following each new guide to an isoform that carried it silently re-pointed the map away from the
transcript they had chosen, so a pinned isoform stays, the note says the selection is not enumerated
there, and leaving it is a button.

## The page as a document, and the index as a control

The two panes size themselves. `main{height:calc(100vh - 86px)}` was arithmetic on a header the stylesheet
cannot see -- it renders 129.5px on one internal run and grows with the gene query, the provenance line,
the cap statement and the off-target scope -- so the panes ran 43.5px past the fold, the document itself
became scrollable and each pane then had a scrollbar inside a page that already had one. Nothing holds a
copy of that number now: the header takes its own height, the panes take what is left, and
`grid-template-rows:minmax(0,1fr)` is the other half of it, because an auto row sizes to its content and a
5,000-row pane would overflow the box it was told to fit. A narrow window stacks the panes rather than
overflowing a 430px minimum inside a 40% track, and **printing** is a supported view rather than an
accident: the nesting becomes ordinary flow, the controls are dropped -- a slider printed at some position
invites a reader to believe it -- and a print-only line says that every verdict on the paper is the one in
force when it was printed, and that the index printed the rows that were drawn.

The index is a control, not a table. Six columns sort, from the keyboard as well as the pointer, and say
which way in the accessibility tree (`aria-sort`) and in a glyph; rows and pick cells answer Enter and
Space; the arrows walk the index; and focus survives the rebuild a keypress itself causes, matched by guide
rather than by row position. The seventh column is the blank pick header, and it used to carry the same
sort handler as the six real ones -- clicking it set the sort key to `undefined`, every row's key became 0
and the sort was silently destroyed.

It also draws a **window** of `INDEX_DRAW_MAX` (400) rows in the current sort rather than every matching
row. A threshold keystroke cost 148 ms at the audited scale, almost all of it building 5,000 rows of seven
cells and wiring 5,000 click handlers; the guide's index is stashed once, the handlers are delegated to the
`tbody`, and the keystroke path measures 19 ms against 61 ms for the old index build alone. This is not
virtualisation and does not pretend to be: the last row of the table says how many rows matched and says
that the counts, the cart, `Add top <n>` and the export are computed over every one of them. A reader who
wants row 3,000 sorts or filters to it. For the same reason the detail pane says where the open guide sits
relative to the index beside it -- a full guide, every card real, next to an index reading "0 of 5,000" is
a screen whose honest reading is that the index is broken -- and names which of the thresholds or the
search box excluded it.

`Add top <n>` is a labelled number box of its own, read through the same boundary as every threshold. It
shipped **inside** its own button: invalid markup that read that way, with the accessible name "Add top to
cart" and the value nowhere in it, a click in the field landing on the button until a handler guessed
otherwise, and `parseInt('ten',10)||0` adding the top zero guides so junk looked like a broken button.
Asking for more than the filters admit is not an error, and it now says what it did.

## Size, and the one figure

The tracked `baseline_0_7_1` fixture slice -- 292 candidate rows, 28 guides, 2,173 classified
alignments and 268 miRNA seed hits -- renders through `sirnaforge report` to **585,987 bytes** (572 KiB),
of which markup, CSS and script are under 0.5 MB in total and effectively constant. That figure moves
with the template rather than with the run: it was 519,084 bytes before the reader filters, the
threshold-provenance banner, the refusal painting and the cart's comment header were added, and 554,923
bytes before the reader-facing repairs described above -- the cap statement, the client-painted
map legend and pills, the ranked gate panel, the per-species liability cell, the two `@media` blocks and
the index's own controls -- none of which scale with the number of candidates. A full internal reference run of tens of
thousands of candidate rows lands near 15 MB; the growth is candidate rows, embedded off-target
evidence at `human, nm<=2`, miRNA seed hits and the count matrix. Two decisions keep it that small:
the filter descriptors are emitted **once** and each guide carries only `[value, verdict, reason]`
codes against them, and row-level detail is embedded only for liability alignments within the embedded
scope. Everything outside that scope is still counted completely, and a guide whose only hits lie
outside it says so rather than rendering like a genuinely clean one.

The **size** of that scope is a run-level fact and is stated as one, in the header, with the arithmetic
that closes it and the breakdown by mismatch count and by species. It has to be, because a per-guide
banner is worthless once it is universal: 197,727 of one run's 232,864 liability alignments sit at
`nm >= 3`, so 4,616 of 5,000 guides raised "not clean, but nothing in the embedded scope" and a reader
learned to skip it. The per-guide banner is kept only for the case that could change a decision -- a guide
nothing here rejects, carrying liability this file cannot itemise -- and a guide already being rejected has
the same fact stated plainly instead.

The one cap on how many guides a report holds is `DEFAULT_MAX_EMBEDDED_GUIDES` (5,000), and it lives in
`build_payload` **above** every count the report prints: `guides` is what the file holds, and
`guides_total`, `guides_dropped` and `guides_dropped_by_status` say what the cap cost. When it bites, the
header states it -- how many guides are missing, the limit that dropped them, the breakdown by verdict
with its zeros, and how many of the dropped were not rejections -- so no number on the page describes a
set the file does not contain. The dropped set is chosen by verdict, worst first, never by truncating a
score sort. Separately and for a different reason, the index **draws** a window of `INDEX_DRAW_MAX` (400)
rows in the current sort: the last row of the table says how many rows matched and says that the pills,
the agreement line, the cart, `Add top <n>` and the export are all computed over every one of them, not
over the window. Neither number is silent, and they are not the same number.

The single figure is a hand-drawn inline SVG. `plotly` was adopted, implemented and then removed: its
4.29 MB bundle carries external URLs and browser-storage references in map traces the report never
invokes, and no static check can tell an inert string from a live call. One stacked bar did not justify
losing the zero-external-reach guarantee.

## Running the client-side parity check

The evaluator that ships inside the report is JavaScript, and the filter it must agree with is Python
(`sirnaforge.reporting.payload.reevaluate_gates`). Rather than a browser -- Playwright is not a
dependency of this repo, and its browser binary is a separate download that would make the test skip
instead of gate -- the parity fixture runs the shipped evaluator directly under **Node**, which this
project already expects (`.nvmrc` pins Node 20; Node 26 also works). The placeholder substitution
happens first, in Python, by calling `render_html` itself, so Node is handed the same JavaScript text a
browser would receive rather than an unsubstituted template, and it is driven through `subprocess` the
same way the suite already shells out to other native tools. Every DOM-touching statement in the
template lives inside a function or behind one `typeof document !== 'undefined'` guard, which is what
makes that possible. Missing Node **fails** the test with the install hint; it does not skip.

```bash
uv run pytest tests/unit/test_report_client_evaluator_parity.py tests/unit/test_reporting_rethreshold.py -q -n 0
```

Both files run in `make test-dev`. The Node harness compares the shipped evaluator against the Python
filter over **every guide x every filter at five threshold sets** -- the run's own, every control at
zero, one below each threshold, far above each threshold, and a set naming both a frozen filter and a
filter id that does not exist -- and asserts the gate triples, all three counters and the derived
status match on all of them. It also pins the three presets that are not a plain status filter to
their trap cases, and asserts that a fragment naming a frozen or unknown filter is refused visibly. A
companion assertion guards the fixture itself against going vacuous: every reason code in `payload.py`
must still fire on at least one guide.

The same harness drives the three reader filters through the functions the report ships --
`passesReaderFilters`, `passesFilters`, `encodeHash`, `applyHashFragment` and the `resetControls` the
`Reset` button itself calls: each filter's selection is compared against the payload's own numbers, a
guide with **no** composite score is asserted to be excluded by a composite floor, the three are shown
to narrow together with a preset and the status set (each cutting rows the other two keep), and the
fragment is round-tripped and then reset. The fixture carries a guide with 4 liabilities and one with 2
so that "at most 3" has a real answer, and one guide with no composite score at all.

Four more drivers cover the parts of the report that are text rather than arithmetic, because a report
whose words disagree with its own numbers is wrong in the way that matters:

- **the gate panel.** `gatesCard` touches no DOM -- it returns HTML as a string -- so the harness builds
  the real panel at the run's thresholds and again with one gate moved, parses every row out of it, and
  asserts as a property that the threshold each Why sentence cites is the one in force, that the
  sentence's polarity is arithmetically true, and that it agrees with the pill in the same row. Two rows
  cover the run's own non-decisions: reason 5 and reason 6 must each state _which_ non-decision it was
  and never a comparison, and the panel's banner must name the set it counts.
- **the cart export.** The TSV is generated twice, and the `#` provenance lines, the `status_basis` and
  the frozen column order are all asserted against the payload's own thresholds.
- **the numeric boundary.** `parseControlValue` is driven over thirteen strings, and `encodeHash` is
  shown to omit a poisoned `T`/`R` rather than write a fragment `applyHashFragment` refuses. The DOM
  half -- `type="number"`, `validity.badInput`, the refusal banner -- is asserted structurally, because
  only a browser sets that flag.
- **the slice itself.** `_slice_script` ends where a browser ends a script (the **first** `</script>`)
  and refuses a document that closes its own script twice, and a companion test renders a payload
  spelling `</script>` to prove the escape in `render._embed` holds. Without both, a truncated document
  could pass parity.
