"""Render a :class:`ReportPayload` to one self-contained HTML file.

Issue #103. The file must work identically over ``file://``, served as a static file, and inside the
Quilt Catalog with permissive HTML rendering **disabled** -- which is the default, and which withholds
the ``allow-same-origin`` sandbox token. That rules out ``fetch``, ``LocalStorage`` and cookies, so
everything lives in the one file and reader state lives in memory and the URL fragment.

The one figure is a hand-drawn inline SVG. plotly was tried and removed (D16): its bundle is 4.29 MB
and its map traces carry external URLs and browser-storage/network API references this report never
invokes, which no static check can tell apart from live ones. One stacked bar does not justify that.
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import asdict
from pathlib import Path
from typing import Any

from jinja2 import Environment, select_autoescape

from sirnaforge.reporting.payload import (
    EMBED_MAX_NM,
    EMBED_SPECIES,
    LIABILITY,
    MIN_UNCOVERED_NT,
    REGISTER_NEIGHBOUR_NT,
    ReportPayload,
)
from sirnaforge.reporting.tracks import SERIES_FILL, PointSeries, TranscriptRegions, transcript_map_svg

#: Hit classes that are a liability, as the strings the payload's own matrix carries.
_LIABILITY_CLASSES = frozenset(c.value for c in LIABILITY)

#: Mismatch bands ``payload._offtarget_views`` emits, in reading order. ``unknown`` is a row whose
#: ``nm`` the run left blank, which is neither inside the embedded scope nor countable against a band.
_NM_BANDS = ("0", "1", "2", ">=3", "unknown")


def _liability_scope(payload: ReportPayload) -> dict[str, Any]:
    """What the embedded off-target scope shows, and how much of the run's liability sits outside it.

    The report embeds per-row detail for human liabilities at nm <= 2 only; everything else rides as
    counts. That is a deliberate scope decision, but it was stated only once per *guide*, as a warning
    -- and on one internal run 197,727 of 232,864 liability alignments sit at nm >= 3, so 4,616 of
    5,000 guides raised it. A warning on 92% of rows is not a warning, so the size of the scope is
    stated once here, at run level, and the per-guide note is reserved for the guides a reader might
    act on (see ``offtargetCard``).

    Every number is summed over the guides this report embedded, so the arithmetic closes against the
    pills and the index rather than against a run total the file does not contain.
    """
    bands: Counter[str] = Counter()
    species: Counter[str] = Counter()
    for guide in payload.guides:
        for cell in guide.offtarget_matrix:
            if cell["hit_class"] in _LIABILITY_CLASSES:
                bands[str(cell["nm"])] += int(cell["n"])
                species[str(cell["species"])] += int(cell["n"])
    liabilities = sum(g.liability_count for g in payload.guides)
    in_scope = sum(len(g.offtarget_rows) for g in payload.guides)
    return {
        "liabilities": liabilities,
        "in_scope": in_scope,
        "outside": liabilities - in_scope,
        "bands": [(band, bands.get(band, 0)) for band in _NM_BANDS if bands.get(band, 0)],
        "species": sorted(species.items()),
        "guides_with_liabilities": sum(1 for g in payload.guides if g.liability_count),
        # The guides the per-guide note is about: they carry liability the reader cannot open here.
        "guides_with_liability_outside": sum(1 for g in payload.guides if g.liability_count > len(g.offtarget_rows)),
        "max_nm": EMBED_MAX_NM,
    }


def _embed(value: Any) -> str:
    r"""JSON for a ``<script>`` body: no payload string can close the block or open a tag.

    An HTML parser ends a script at the first literal ``</script>`` inside it, whatever the JavaScript
    means -- so a gene query, transcript id or gene symbol spelling that would truncate the document at
    that point and take every panel below it with it. ``json.dumps`` does not escape ``<``, and there is
    no sanitising layer between the payload and this file, so the escape belongs here. ``\\u003c`` is
    ordinary JSON and parses back to the same string, so nothing downstream sees a difference.
    """
    text = json.dumps(value, separators=(",", ":"))
    return text.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")


_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>siRNAforge report — {{ p.run.gene_query }}</title>
<style>
:root{--bg:#fbfbfc;--fg:#1a1d21;--mut:#6b7280;--line:#e5e7eb;--card:#fff;
--pass:#15803d;--fail:#b91c1c;--unk:#b45309;--off:#9ca3af;--acc:#1d4ed8}
*{box-sizing:border-box}
/* The two panes fill whatever the header leaves, measured by the layout rather than guessed at. The
   height was `calc(100vh - 86px)` against a header that renders 129.5px on a real run -- gene query,
   provenance, four pills, the agreement line, the dropped-guide banner and the off-target scope, every
   one of which grows with the run -- so `main` ran 43.5px past the fold, the document itself became
   scrollable and both panes got a second scrollbar inside a page that already had one. A flex column
   with `min-height:0` cannot be wrong about a number it never holds. */
html,body{height:100%}
body{margin:0;font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;background:var(--bg);
color:var(--fg);display:flex;flex-direction:column}
header{flex:0 0 auto;padding:18px 24px;border-bottom:1px solid var(--line);background:var(--card)}
h1{margin:0 0 4px;font-size:17px;font-weight:650}
.sub{color:var(--mut);font-size:12.5px}
/* `grid-template-rows:minmax(0,1fr)` as well as `min-height:0`: an auto-sized row takes its height from
   its content, which is how a pane 5,000 rows tall overflows a box that was told to be short. */
main{flex:1 1 auto;min-height:0;display:grid;grid-template-columns:minmax(430px,40%) 1fr;
grid-template-rows:minmax(0,1fr);gap:0}
/* One header cell carries no label a reader needs, and every label a screen reader needs. */
.vh{position:absolute;width:1px;height:1px;overflow:hidden;clip-path:inset(50%);white-space:nowrap}
.printonly{display:none}
/* The index carries six columns including a 23-nt monospace sequence. Without tightened padding and
   short headers it overflows its pane and the last column -- isoform coverage -- is the one lost. */
#idx th,#idx td{padding:6px 7px}
#idx td:first-child{white-space:nowrap}
#idx td.pick,#idx th.pick{width:24px;padding:6px 2px 6px 8px;text-align:center}
.filters{padding:0 12px 10px;border-bottom:1px solid var(--line);background:var(--card)}
.filters summary{cursor:pointer;font-size:11.5px;text-transform:uppercase;letter-spacing:.04em;
color:var(--mut);font-weight:600;padding:8px 0}
.frow{display:grid;grid-template-columns:1fr auto auto;gap:6px 8px;align-items:center;font-size:12px;
margin-bottom:5px}
.frow input{width:62px;padding:2px 5px;border:1px solid var(--line);border-radius:4px;font:inherit;font-size:12px}
/* Every threshold box is type=number, so the browser itself refuses text that is not a number and
   leaves the only trace of it in validity.badInput -- which fnum() reads. The spinner is hidden
   because a threshold is typed here and stepped on the slider, and aria-invalid is painted red so a
   refused keystroke is visible at the control rather than silently ignored. */
.frow input[type=number]{appearance:textfield;-moz-appearance:textfield}
.frow input[type=number]::-webkit-outer-spin-button,
.frow input[type=number]::-webkit-inner-spin-button{-webkit-appearance:none;margin:0}
.frow input[aria-invalid="true"]{border-color:var(--fail);background:#fee2e2}
/* A gate control changes a verdict; a reader filter only selects rows. The two groups are separated
   and each says which it is, because mistaking one for the other misreads what the report claims. */
.fsec{font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:var(--mut);font-weight:600;
margin:9px 0 5px;padding-top:7px;border-top:1px solid var(--line)}
.fsec i{text-transform:none;letter-spacing:0;font-weight:400;font-style:normal}
/* What each gate did to the guides in this file, under its own control: the panel is a ranked list,
   and a slider that decided nothing has to say so where the reader is about to reach for it. */
.geffect{font-size:11.5px;color:var(--mut);margin:-1px 0 8px}
/* A threshold the slider's domain cannot represent. Painted in the refusal colour and sitting between
   the two controls, because the whole point is that the box and the slider no longer agree and the box
   is the one in force. Empty when they do agree, so the row keeps its height. */
.offscale{font-size:11.5px;color:var(--fail);margin:-1px 0 7px}
.offscale:empty{display:none}
input[type=range]:disabled{opacity:.45}
.fstat label{margin-right:9px;font-size:12px;white-space:nowrap}
.cartbtn{cursor:pointer;border:1px solid var(--line);background:var(--card);border-radius:5px;
padding:3px 9px;font:inherit;font-size:12px;color:var(--fg)}
.cartbtn:hover{background:#f3f4f6}
.picked{color:var(--acc);font-weight:700}
#cartlist{max-height:210px;overflow:auto;font-size:12.5px}
/* Tab-separated columns only line up if the text is not wrapped. */
#carttsv{width:100%;height:96px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;
border:1px solid var(--line);border-radius:5px;padding:7px;resize:vertical;white-space:pre;overflow:auto}
#left{border-right:1px solid var(--line);overflow:auto;min-height:0;background:var(--card)}
#right{overflow:auto;min-height:0;padding:20px 24px}
#q{width:100%;padding:9px 12px;border:1px solid var(--line);border-radius:6px;font:inherit;font-family:ui-monospace,monospace}
.searchbar{padding:12px;position:sticky;top:0;background:var(--card);border-bottom:1px solid var(--line);z-index:2}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th,td{text-align:left;padding:6px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{font-weight:600;color:var(--mut);font-size:11px;text-transform:uppercase;letter-spacing:.04em;
position:sticky;top:57px;background:var(--card);cursor:pointer;user-select:none}
#right th{position:static;top:auto}
tbody tr{cursor:pointer}
tbody tr:hover{background:#f3f4f6}
tbody tr.sel{background:#dbeafe}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:16px;margin-bottom:16px}
.card h2{margin:0 0 12px;font-size:13px;text-transform:uppercase;letter-spacing:.05em;color:var(--mut)}
.kv{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px 18px}
.kv div{font-size:12.5px}.kv b{display:block;color:var(--mut);font-weight:600;font-size:11px}
.pill{display:inline-block;padding:1px 7px;border-radius:10px;font-size:11px;font-weight:650}
.v-pass{color:var(--pass);background:#dcfce7}.v-fail{color:var(--fail);background:#fee2e2}
.v-unknown{color:var(--unk);background:#fef3c7}.v-not_evaluated{color:var(--off);background:#f3f4f6}
.v-warn{color:#92400e;background:#fef9c3}
#mapwrap{position:relative;overflow-x:auto}
/* Same max-width as the base map: without it a wide viewport stretches the overlay past the
   base, preserveAspectRatio centres it, and the selection marker points at the wrong position. */
#mapover{position:absolute;left:0;top:0;width:100%;height:100%;max-width:1120px;pointer-events:none}
select{font:inherit;font-size:12px;padding:2px 6px;border:1px solid var(--line);border-radius:5px;
background:var(--card);color:var(--fg);text-transform:none;letter-spacing:0}
.struct{display:flex;gap:20px;flex-wrap:wrap;align-items:flex-start}
.struct>div:first-child{flex:0 0 380px;max-width:380px}
.struct>div:last-child{flex:1 1 260px;min-width:240px}
.struct .db{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px;
letter-spacing:.5px;line-height:1.5;word-break:break-all}
.warn{background:#fffbeb;border-left:3px solid var(--unk);padding:9px 12px;font-size:12.5px;margin-bottom:12px;border-radius:0 5px 5px 0}
.empty{color:var(--mut);font-style:italic;font-size:12.5px}
.big{font-size:22px;font-weight:680}
.liab{color:var(--fail);font-weight:650}.nonliab{color:var(--mut)}
code{background:#f3f4f6;padding:1px 4px;border-radius:3px;font-size:12px}
/* A sort a reader can see, on a column they can reach: the arrow is the visible half of `aria-sort`,
   and the focus ring is what makes a `<th>` a control rather than a word that happens to react. */
th .si{margin-left:4px;font-size:9px;color:var(--acc)}
th[aria-sort="none"] .si{color:var(--line)}
#idx th:focus-visible,#idx tbody tr:focus-visible,#idx td.pick:focus-visible{outline:2px solid var(--acc);outline-offset:-2px}
#idx th.pick{cursor:default}
/* The index draws a window of the matching rows and says so in the last row rather than in a tooltip. */
tr.drawcap td{color:var(--mut);font-style:italic;background:#f9fafb;cursor:default}
tr.drawcap:hover td{background:#f9fafb}
/* Same treatment as every threshold box, because it is refused the same way rather than degraded to 0. */
#topn{padding:2px 5px;border:1px solid var(--line);border-radius:4px;font:inherit;font-size:12px;
appearance:textfield;-moz-appearance:textfield}
#topn::-webkit-outer-spin-button,#topn::-webkit-inner-spin-button{-webkit-appearance:none;margin:0}
#topn[aria-invalid="true"]{border-color:var(--fail);background:#fee2e2}
/* Narrow viewports: the left pane's own 430px minimum plus a 40% track is wider than a phone or a
   half-screen window, and the grid does not wrap. Stacked, both panes keep their own scroll. */
@media (max-width:820px){
  main{grid-template-columns:minmax(0,1fr);grid-template-rows:minmax(0,45%) minmax(0,55%)}
  #left{border-right:0;border-bottom:1px solid var(--line)}
  #right{padding:16px}
}
/* Paper. Nothing here scrolls, so every nested scroller has to become ordinary flow or its content is
   cut at the first page: a printed copy of this report used to be the header and the top 40% of two
   panes, once. The controls are dropped because they cannot be operated on paper, and the print-only
   line says what the reader is therefore not seeing. */
@media print{
  html,body{height:auto}
  body{display:block;background:#fff}
  main{display:block;min-height:0}
  #left,#right{overflow:visible;min-height:0;border:0;padding:0}
  #left{margin-bottom:14px}
  .searchbar,.filters,.cartbtn,#carttsv{display:none}
  #idx th,.searchbar{position:static}
  .printonly{display:block}
  .card,tr,svg{break-inside:avoid}
  #idx{font-size:10.5px}
  @page{margin:14mm}
}
</style></head><body>
<header>
  {# Every number in this header counts the guides the FILE holds. `p.run.guides` is that count
     (payload._select_embedded); the run's own total is named beside it whenever the two differ, and
     never substituted for it. A header that headlined the run's 5,706 while the file held 5,000 was
     the sharpest form of that defect: the agreement line then read "0/1,424 run PASS contradicted"
     as a coverage claim over 41 guides the reader could not open. #}
  <h1>siRNAforge — {{ p.run.gene_query }} ·
    {{ '{:,}'.format(p.run.guides) }} guides{% if p.run.guides_dropped %} of the run's
    {{ '{:,}'.format(p.run.guides_total) }}{% endif %}</h1>
  <div class="sub">
    {{ '{:,}'.format(p.run.candidate_rows) }}{% if p.run.candidate_rows_total != p.run.candidate_rows %}
    of the run's {{ '{:,}'.format(p.run.candidate_rows_total) }}{% endif %}
    candidate rows collapsed onto guide sequence ·
    {{ '{:,}'.format(p.run.hit_rows) }} alignments
    (<span class="liab">{{ '{:,}'.format(p.run.liability_rows) }} liabilities</span>,
     <span class="nonliab">{{ '{:,}'.format(p.run.non_liability_rows) }} not liabilities</span>) ·
    payload schema {{ p.schema_version }} · profile {{ p.provenance.policy_profile }} ·
    run mode {{ p.provenance.run_mode }} · gates from {{ p.provenance.policy_source }} ·
    sirnaforge {{ p.provenance.tool_version }}
  </div>
  {# The four counts are Jinja-seeded with the run's own and re-rendered by renderStatusPills() the
     moment a threshold moves: frozen pills beside a live map is the same defect as a frozen legend. The
     agreement line is NOT re-derived, because comparing the run to a threshold the run never used is not
     an agreement -- so it says which thresholds it was computed at, and says when they are no longer
     the ones in force. #}
  <div class="sub" style="margin-top:6px">
    <span class="pill v-pass"><span id="n_pass">{{ p.run.status_counts["pass"] }}</span> pass</span>
    <span class="pill v-warn"><span id="n_warn">{{ p.run.status_counts["warn"] }}</span> pass, warned</span>
    <span class="pill v-unknown"><span id="n_unknown">{{ p.run.status_counts["unknown"] }}</span> not established</span>
    <span class="pill v-fail"><span id="n_fail">{{ p.run.status_counts["fail"] }}</span> fail</span>
    &nbsp;· over the {{ '{:,}'.format(p.run.guides) }} guides in this file<span id="cbasis"></span>
    {% if p.run.agreement.comparable %}
    <div style="margin-top:3px">Against the run's own verdicts:
    <b>{{ p.run.agreement.contradicted_run_pass }}</b>/{{ p.run.agreement.run_pass_guides }} run PASS
    contradicted, <b>{{ p.run.agreement.overruled_run_fail }}</b> run rejections overruled,
    <b>{{ p.run.agreement.run_failed_not_rederivable }}</b> run rejections not re-derivable here
    — all four counted over those same {{ '{:,}'.format(p.run.guides) }}, at the run's own
    thresholds.<span id="agreestale"></span></div>
    {% endif %}
  </div>
  {% if p.run.guides_dropped %}
  <div class="warn" style="margin-top:8px">
    <b>{{ '{:,}'.format(p.run.guides_dropped) }} of the run's {{ '{:,}'.format(p.run.guides_total) }}
    guides are not in this file.</b> The embed limit is
    {{ '{:,}'.format(p.run.guide_embed_limit) }} guides, and they cannot be searched, carted or
    exported here. Not embedded, by verdict:
    {% for status, n in p.run.guides_dropped_by_status.items() %}{{ '{:,}'.format(n) }} {{ status }}{{ ", " if not loop.last }}{% endfor %}.
    {% if dropped_shippable %}<b>{{ '{:,}'.format(dropped_shippable) }} of them are not rejections</b> —
    re-run the report with a higher guide limit to see them.{% else %}All of them are guides the run
    and this report agree to reject.{% endif %}
  </div>
  {% endif %}
  {# Finding 8: the embedded off-target scope, stated once, at run level, with its own arithmetic. #}
  <div class="sub" style="margin-top:6px">
    <b>Off-target scope.</b> Per-alignment detail is embedded for {{ p.run.embed_scope }} only:
    {{ '{:,}'.format(scope.in_scope) }} of the {{ '{:,}'.format(scope.liabilities) }} liability
    alignments these {{ '{:,}'.format(p.run.guides) }} guides carry
    ({{ '%.0f'|format(100 * scope.in_scope / scope.liabilities if scope.liabilities else 0) }}%).
    The other {{ '{:,}'.format(scope.outside) }} are counted, charted and gated but have no row of
    their own here, which is why
    {{ '{:,}'.format(scope.guides_with_liability_outside) }} of
    {{ '{:,}'.format(p.run.guides) }} guides carry liability this file cannot itemise.
    {% if scope.bands or scope.species %}
    <details style="margin-top:4px"><summary style="cursor:pointer">what sits outside it</summary>
      <div style="margin:4px 0 0">Liability alignments by mismatch count:
        {% for band, n in scope.bands %}<code>nm {{ band }}</code> {{ '{:,}'.format(n) }}{{ ", " if not loop.last }}{% endfor %}.
        Only <code>nm 0</code>–<code>nm {{ scope.max_nm }}</code> can be in scope.</div>
      <div>By species:
        {% for name, n in scope.species %}{{ name }} {{ '{:,}'.format(n) }}{{ ", " if not loop.last }}{% endfor %}.
        Only {{ p.run.embed_scope.split(',')[0] }} can be in scope, so every alignment in another
        species is outside it by definition — a rejection driven by those is a non-human-index
        finding, not a human liability.</div>
      {% if p.run.screened_species %}
      <div>Species screened: {{ p.run.screened_species|join(', ') }}. Run-wide liability alignments by
        species (all guides, not only the embedded ones):
        {% for name, n in p.run.liability_rows_by_species.items() %}{{ name }} {{ '{:,}'.format(n) }}{{ ", " if not loop.last }}{% endfor %}.</div>
      {% endif %}
    </details>
    {% endif %}
  </div>
  {# Only ever seen on paper (`.printonly`), where the controls are dropped because they cannot be
     operated there. What a printed copy cannot show, it says. #}
  <div class="printonly sub" style="margin-top:6px"><b>Printed copy.</b> The thresholds and filters are
  controls, so they are not on paper: every verdict printed below is the one in force when this was
  printed, and the legend under the candidate map says whether those thresholds are the run's own. The
  index prints the rows that were drawn on screen, in the sort order they were drawn in.</div>
</header>
<main>
 <div id="left">
  <div class="searchbar"><input id="q" placeholder="Search guide sequence, candidate id or transcript…" autocomplete="off"></div>
  <details class="filters" id="filters">
    <summary>Thresholds and filters — <span id="fcount"></span></summary>
    <div id="presets" style="margin-bottom:8px"></div>
    <div id="refused" class="warn" style="display:none"></div>
    <div id="badnum" class="warn" style="display:none"></div>
    <div class="fstat" id="fstat"></div>
    <div class="fsec">Gate thresholds <i>— moving one re-derives the verdict this run computed.
    Ranked by the guides each gate is solely responsible for rejecting in this file.</i></div>
    <div id="frows"></div>
    <div class="fsec">Reader filters <i>— these only select among rows; no verdict changes</i></div>
    <div id="rrows"></div>
    <div class="fstat" id="fcons" style="margin-top:6px;padding-top:6px;border-top:1px solid var(--line)">
    </div>
    {# The number box was INSIDE the button: invalid markup, an accessible name that read "Add top  to
       cart" with the value nowhere in it, and a click on the field counted as a click on the button
       until a handler was written to guess otherwise. It is its own labelled control now, beside a
       button whose name is what the button does, and junk in it is refused where it is typed rather
       than degraded to 0 by `parseInt(...)||0` and silently adding nothing. #}
    <div style="margin-top:8px">
      <label for="topn">Add top</label>
      <input id="topn" type="number" step="1" min="1" value="10" style="width:52px"
        aria-label="how many of the highest-scoring matching guides to add to the cart">
      <button class="cartbtn" id="addtop">Add to cart</button>
      <button class="cartbtn" id="freset">Reset</button>
      <div id="topnnote" class="empty" style="margin-top:4px"></div>
    </div>
  </details>
  {# Seven columns, six of them sortable, none of them reachable before: `<th onclick>` with no tabindex
     and no key handler, `aria-sort` unset on all seven so nothing said which column was sorted or which
     way, and the blank pick header carried the same handler as the rest -- clicking it set `sortK` to
     undefined and every row's key to 0, silently destroying the sort. The pick column has no key, so
     `sortBy` refuses it; the six that do carry one are buttons in every sense a keyboard can tell. #}
  <table id="idx"><thead><tr>
    <th class="pick" scope="col"><span class="vh">In cart</span></th>
    <th scope="col" tabindex="0" aria-sort="none" data-k="status">Status<span class="si" aria-hidden="true"></span></th>
    <th scope="col" tabindex="0" aria-sort="none" data-k="guide">Guide<span class="si" aria-hidden="true"></span></th>
    <th scope="col" tabindex="0" aria-sort="none" data-k="composite">Score<span class="si" aria-hidden="true"></span></th>
    <th scope="col" tabindex="0" aria-sort="none" data-k="failed" title="gates failed / gates undecided">Gates<span class="si" aria-hidden="true"></span></th>
    {# One liability column, decomposed in the cell rather than summarised in it. The bare total is the
       number max_off_target_count gates, and it is also the number that cannot tell a human liability
       from a non-human-index artefact: on one internal run roughly half of every liability total is
       non-human. `#liabscope` re-scopes the sort, the emphasis and the liability bound together. #}
    <th scope="col" tabindex="0" aria-sort="none" data-k="liab" id="thliab"><span
      id="liablbl">Liab.</span><span class="si" aria-hidden="true"></span></th>
    <th scope="col" tabindex="0" aria-sort="none" data-k="rows"
      title="distinct isoforms carrying this guide, of all in the run">Isoforms<span class="si" aria-hidden="true"></span></th>
  </tr></thead><tbody></tbody></table>
 </div>
 <div id="right">
   <div class="card" id="cartcard" style="display:none">
     <h2>Cart &nbsp;<span id="cartn"></span></h2>
     <div id="cartlist"></div>
     <div style="margin:10px 0 8px">
       <button class="cartbtn" id="cartdl">Export TSV</button>
       <button class="cartbtn" id="cartcopy">Select all text</button>
       <button class="cartbtn" id="cartclear">Clear</button>
       <span id="dlnote" class="empty"></span>
     </div>
     <textarea id="carttsv" readonly aria-label="cart as TSV"></textarea>
     <p class="empty" style="margin:6px 0 0">Tab-separated, one row per guide, with the values the
     thresholds above were applied to. The <code>#</code> lines above the header name the thresholds the
     <code>status</code> column was decided at, and say whether they are the run's own or yours: a
     status pasted into a ticket has to carry what produced it.
     <b>Export TSV</b> downloads it. Some viewers -- including a
     Quilt iframe without <code>allow-downloads</code> -- block that, and they do it without telling this
     page: no error reaches the report, so it cannot say whether your file arrived. The box above is
     therefore selected every time you export, so a blocked download is one copy away from the data.</p>
     {# Where the cart lives between reloads, said on screen. It cannot be browser storage: this file
        reaches nothing outside itself, and the sandbox it is read in withholds storage anyway. The URL
        is the only place left, so the reader is told that the address IS the cart. #}
     <div id="carturl" class="empty" style="margin:6px 0 0"></div>
   </div>
   <div class="card" id="mapcard">
     <h2>Candidate positions</h2>
     <div style="font-size:12px;color:#6b7280;margin:-4px 0 10px">Every candidate the run enumerated
     on one isoform, coloured by the verdict in force <b>now</b> — move a gate threshold and the dots
     and the legend below move with it. <label for="tx">Isoform</label>
     <select id="tx" aria-label="isoform"></select>
     <span id="txnote"></span>
     {# Only shown once the reader has chosen an isoform and the selected guide is not on it: the map used
        to leave a chosen transcript silently, and now that it stays, leaving it has to be an action. #}
     <button class="cartbtn" id="txfollow" style="display:none">Show an isoform carrying it</button>
     {% if not p.run.canonical_source %}
     <div>This run recorded no canonical transcript, so none of these is marked as one.</div>
     {% endif %}</div>
     <div id="maplegend"></div>
     <div id="mapwrap"><div id="mapbase"></div><svg id="mapover"></svg></div>
     <div id="mapnote" class="empty"></div>
   </div>
   <div id="detail"><p class="empty">Select a guide.</p></div>
 </div>
</main>
<script>
const G = GUIDES_JSON_PLACEHOLDER;
const FILTERS = FILTERS_JSON_PLACEHOLDER;
const MAPS = MAPS_JSON_PLACEHOLDER;
const MAP_LEGEND = MAP_LEGEND_JSON_PLACEHOLDER;
const NOT_EMBEDDED = NOT_EMBEDDED_JSON_PLACEHOLDER;
const GAPS_KEY = GAPS_KEY_JSON_PLACEHOLDER;
const LAYOUTS = LAYOUTS_JSON_PLACEHOLDER;
const TX_IDS = TX_IDS_JSON_PLACEHOLDER;
const SCREENED_SPECIES = SCREENED_SPECIES_JSON_PLACEHOLDER;
const EMBED_SPECIES = EMBED_SPECIES_JSON_PLACEHOLDER;
const EMBED_SCOPE = EMBED_SCOPE_JSON_PLACEHOLDER;
const GENE = GENE_JSON_PLACEHOLDER;
const REGISTER_NT = REGISTER_NT_PLACEHOLDER;
const fmt = n => n===null||n===undefined ? '—' : (typeof n==='number' ? (Number.isInteger(n)?n.toLocaleString():n.toFixed(3)) : n);
const esc = s => String(s??'').replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
let view = G.slice(), sortK='composite', sortAsc=false, selected=null;

// Each guide's own place in this file, stashed once. `renderIndex` read `G.indexOf(g)` INSIDE the row
// map, so one keystroke on a threshold cost 5,000 scans of a 5,000-element array -- 25 million
// comparisons, and most of the 148 ms every keystroke took on one internal run. A guide's index never
// changes, and computing it once is the whole fix.
G.forEach((g,i)=>{ g._i=i; });

// How many of the matching rows the index draws. Not virtualisation: 5,000 rows of 7 cells is 35,000
// nodes, each with a decomposed liability cell and a title built per row, and building all of them is
// what remains of the keystroke cost once the O(n^2) above is gone. A window is proportionate and can be
// honest about itself -- the last row says how many matched, and every count, the cart and the export
// are computed over all of them and never over the window. A reader who needs a row outside it sorts or
// filters, which is what the panel above the table is for.
const INDEX_DRAW_MAX = 400;

// ---- liabilities, decomposed by species ----------------------------------------------------------
// A single all-species liability count is what `max_off_target_count` compares, and it is also the one
// number that cannot tell a human liability from a non-human-index artefact. On one internal run about
// half of every liability total is non-human, so a guide rejected on that ceiling may be carrying no
// human liability at all -- and the reader had no way to see which. `liability_by_species` sums to
// `liability_count` (payload.GuideEntry), so the decomposition is exact rather than an estimate.
// Species names are read from the run (SCREENED_SPECIES); a spelling hardcoded here would outlive the
// screen that produced it. EMBED_SPECIES is the query species, the one the per-row detail is embedded
// for, so it is named first and is what "the rest" is the rest of.
const OTHER_SPECIES = SCREENED_SPECIES.filter(s => s !== EMBED_SPECIES);
function liabOf(g, sp){ return (g.liability_by_species||{})[sp] || 0; }
//: null means every screened species -- the gate's own scope, and the default, so nothing about the
//: index or the liability bound changes meaning until a reader deliberately narrows it.
let liabScope = null;
function liabIn(g){ return liabScope===null ? g.liability_count : liabOf(g, liabScope); }
function scopeLabel(){ return liabScope===null ? 'all screened species' : liabScope; }

// Total first, because that is the gated number, then the split -- with the scoped part emphasised so
// the column the reader chose to rank by is the one their eye lands on.
function liabCell(g){
  const total=g.liability_count;
  if(!SCREENED_SPECIES.length || !total) return String(total);
  const parts=[[EMBED_SPECIES, liabOf(g,EMBED_SPECIES)]]
    .concat(OTHER_SPECIES.length ? [['other', OTHER_SPECIES.reduce((n,sp)=>n+liabOf(g,sp),0)]] : []);
  const split=parts.map(([name,n])=>{
    const scoped = liabScope===null ? false
      : (name==='other' ? OTHER_SPECIES.includes(liabScope) : liabScope===EMBED_SPECIES);
    return `<span${scoped?' class="picked"':''}>${n} ${esc(name)}</span>`;
  }).join(', ');
  return `${total} <span class="empty">(${split})</span>`;
}

function liabTitle(g){
  if(!SCREENED_SPECIES.length) return 'no hit table in this run, so no species is recorded';
  return SCREENED_SPECIES.map(sp=>`${sp} ${liabOf(g,sp)}`).join(', ')
    + ` — ${g.liability_count} in total, the count max_off_target_count compares`;
}

const STATUS_RANK={pass:0,warn:1,unknown:2,fail:3};
// Every key reads g._live -- the client-recomputed gate table -- not the frozen payload snapshot.
// Sorting or painting the index off the run's own g.status/g.n_gates_failed would let the pill and the
// gate panel disagree the moment a threshold moves (render.py:751-755 injects those as a display
// convenience for the render's own defaults, not as the live truth once a control exists to move it).
function rowKey(g,k){
  if(k==='guide')return g.guide; if(k==='composite')return g.composite_score??-1;
  if(k==='design')return g.design_score??-1;
  if(k==='failed')return g._live.n_gates_failed*100+g._live.n_gates_unknown;
  if(k==='status')return -STATUS_RANK[g._live.status];
  // Sorts on the SCOPED count, so re-scoping the column re-ranks it: "worst human liabilities first"
  // and "worst total first" are different questions and the column answers whichever is selected.
  if(k==='liab')return liabIn(g);
  if(k==='rows')return g.n_rows; return 0;
}
// Which column is sorted and which way, said in the markup a screen reader reads and in a glyph a reader
// sees. `aria-sort` was unset on all seven columns and there was no visual indicator at all, so the sort
// was a fact only the person who clicked knew.
function paintSortHeaders(){
  document.querySelectorAll('#idx th[data-k]').forEach(th=>{
    const on = th.dataset.k===sortK;
    th.setAttribute('aria-sort', on ? (sortAsc?'ascending':'descending') : 'none');
    const si=th.querySelector('.si'); if(si) si.textContent = on ? (sortAsc?'\u25b2':'\u25bc') : '\u25c6';
  });
}

// One place for a sort, reached by click and by Enter/Space. The guard is the defect: the blank pick
// header carried the same handler as the six real ones, so clicking it set `sortK` to undefined, rowKey
// returned 0 for every guide and the index kept whatever order the last comparison happened to leave.
function sortBy(k){
  if(!k) return;
  if(k===sortK) sortAsc=!sortAsc; else { sortK=k; sortAsc=(k==='guide'); }
  renderIndex();
}

function renderIndex(){
  const tb=document.querySelector('#idx tbody');
  // The tbody is replaced wholesale on every keystroke, which throws focus back to the document. A reader
  // who opens a row with Enter would get one row and then have to tab in from the top again, so the
  // keyboard path has to survive the rebuild that the keypress itself caused.
  const act=document.activeElement, inRow=act&&act.closest?act.closest('#idx tbody tr'):null;
  const heldPick=inRow?!!act.closest('td.pick'):false;
  const held=inRow&&inRow.querySelector('td.pick')?inRow.querySelector('td.pick').dataset.pick:null;
  view.sort((a,b)=>{const x=rowKey(a,sortK),y=rowKey(b,sortK);
    const c = typeof x==='string' ? x.localeCompare(y) : x-y; return sortAsc?c:-c;});
  const V={pass:'v-pass',warn:'v-warn',unknown:'v-unknown',fail:'v-fail'};
  const drawn = view.slice(0, INDEX_DRAW_MAX);
  // `data-i` is the guide's stashed index, and `tabindex`/`aria-current`/`aria-pressed` are what make a
  // row and its pick cell operable without a mouse: the row was a `<tr onclick>` and the pick cell a bare
  // `<td>` with no role, so neither existed for a keyboard or a screen reader.
  tb.innerHTML = drawn.map(g=>`<tr data-i="${g._i}" tabindex="0"${
    g.guide===selected?' class="sel" aria-current="true"':''}>
    <td class="pick" data-pick="${esc(g.guide)}" role="button" tabindex="0"
      aria-pressed="${CART.has(g.guide)?'true':'false'}" aria-label="${CART.has(g.guide)?'remove':'add'} ${esc(g.guide)} ${CART.has(g.guide)?'from':'to'} cart"
      >${CART.has(g.guide)?'<span class="picked">\u2713</span>':'<span style="color:#d1d5db">+</span>'}</td>
    <td><span class="pill ${V[g._live.status]}">${g._live.status}</span></td>
    <td class="mono">${esc(g.guide)}</td><td>${fmt(g.composite_score)}</td>
    <td>${g._live.n_gates_failed} / ${g._live.n_gates_unknown}</td>
    <td class="${liabIn(g)?'liab':'nonliab'}" title="${esc(liabTitle(g))}">${liabCell(g)}</td>
    <td title="${g.n_rows} enumeration${g.n_rows===1?'':'s'}">${isoformFrac(g)}</td></tr>`).join('')
    + (view.length>drawn.length ? `<tr class="drawcap"><td colspan="7">Showing the first ${
      drawn.length.toLocaleString()} of ${view.length.toLocaleString()} matching rows, in this sort
      order. Every count above, the cart, <b>Add top</b> and the export are computed over all ${
      view.length.toLocaleString()}; sort or filter to bring another row into view.</td></tr>` : '');
  if(held){
    const cells=[...tb.querySelectorAll('td.pick')].filter(c=>c.dataset.pick===held);
    if(cells.length) (heldPick?cells[0]:cells[0].parentElement).focus();
  }
  paintSortHeaders();
  renderPlacement();     // the open guide may be outside the filter, or outside the drawn window
}
// Both wired in wireUI(), at the bottom: every DOM-touching statement in this file is either inside a
// function or behind that one guard, so the substituted script also loads under node with nothing
// stubbed but `document` itself (tests/unit/test_report_client_evaluator_parity.py).

// Hand-drawn inline SVG, deliberately: plotly.min.js is 4.29 MB and its map traces carry external
// URLs and browser-storage/network API references this report never invokes, which no static check
// can tell apart from live ones (D16). One stacked bar does not justify that.
const CLASS_FILL={off_target:'#b91c1c',undetermined:'#b45309',on_target:'#9ca3af',ortholog:'#6b7280',repeat:'#d1d5db'};
const LIABILITY_CLASSES=new Set(['off_target','undetermined']);
function drawMatrix(id, matrix){
  const host=document.getElementById(id); if(!host) return;
  const bands=['0','1','2','>=3','unknown'].filter(b=>matrix.some(m=>m.nm===b));
  const cls=[...new Set(matrix.map(m=>m.hit_class))].sort(
    (a,b)=>(LIABILITY_CLASSES.has(b)?1:0)-(LIABILITY_CLASSES.has(a)?1:0)||a.localeCompare(b));
  const at=(b,c)=>matrix.filter(m=>m.nm===b&&m.hit_class===c).reduce((a,m)=>a+m.n,0);
  const totals=bands.map(b=>cls.reduce((a,c)=>a+at(b,c),0));
  const max=Math.max(1,...totals);
  const W=Math.max(320,host.clientWidth||520), H=200, L=46, R=8, T=10, B=30;
  const pw=W-L-R, ph=H-T-B, bw=Math.min(64, pw/bands.length*0.62), step=pw/bands.length;
  const y=v=>T+ph-(v/max)*ph;
  // Two gridlines and the max, so a bar height is readable without a tooltip.
  const ticks=[0,Math.round(max/2),max].filter((v,i,a)=>a.indexOf(v)===i);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" role="img"
    aria-label="alignments by mismatch count, stacked by hit class">`;
  s+=ticks.map(v=>`<g><line x1="${L}" x2="${W-R}" y1="${y(v)}" y2="${y(v)}" stroke="#e5e7eb"/>
    <text x="${L-6}" y="${y(v)+3.5}" text-anchor="end" font-size="10" fill="#6b7280">${v.toLocaleString()}</text></g>`).join('');
  bands.forEach((b,bi)=>{
    const x=L+bi*step+(step-bw)/2; let acc=0;
    cls.forEach(c=>{const n=at(b,c); if(!n) return;
      const h=(n/max)*ph; acc+=n;
      s+=`<rect x="${x}" y="${y(acc)}" width="${bw}" height="${Math.max(1,h)}" fill="${CLASS_FILL[c]||'#9ca3af'}">
        <title>${c}, ${b===' >=3'?b:b} mismatch: ${n.toLocaleString()} alignment(s)</title></rect>`;});
    s+=`<text x="${x+bw/2}" y="${H-B+14}" text-anchor="middle" font-size="10" fill="#6b7280">${b}</text>`;
    if(totals[bi]) s+=`<text x="${x+bw/2}" y="${y(totals[bi])-4}" text-anchor="middle" font-size="10"
      fill="#1a1d21" font-weight="650">${totals[bi].toLocaleString()}</text>`;
  });
  s+=`<text x="${L+pw/2}" y="${H-2}" text-anchor="middle" font-size="10" fill="#6b7280">mismatches</text></svg>`;
  s+=`<div style="font-size:11px;color:#6b7280;margin-top:2px">`+cls.map(c=>
    `<span style="margin-right:12px"><span style="display:inline-block;width:9px;height:9px;border-radius:2px;
     background:${CLASS_FILL[c]||'#9ca3af'};margin-right:4px"></span>${c}${LIABILITY_CLASSES.has(c)?' (liability)':''}</span>`).join('')+`</div>`;
  host.innerHTML=s;
}
const VERDICT=['pass','fail','unknown','not_evaluated','warn'];
const V_PASS=0, V_FAIL=1, V_UNKNOWN=2, V_NOT_EVALUATED=3, V_WARN=4;
//: payload.py's seven reason codes, restated so the client evaluator can read them -- only REASON_OK is
//: ever re-thresholded; every other reason means the run itself never decided. The two the evaluator
//: names are the two it branches on; gateReason below spells out all seven for the reader.
const REASON_OK=0, REASON_EMPTY_VALUE=2;
// The threshold a verdict shown beside this gate was ACTUALLY decided against: the reader's, once they
// have moved it, and the run's otherwise. Every part of the gate panel reads this one function, because
// the panel's Why sentence read `f.threshold` -- the run's -- while the verdict beside it came from the
// live re-thresholded table, so any reader move made the sentence arithmetically false in both
// directions: "10 > 10 so FAIL" next to a PASS pill, or a passing comparison next to a FAIL (#103).
// A report that states a comparison its own verdict did not come from is publishing a verdict the
// pipeline would not reproduce, which is the one failure this report exists to prevent.
function liveThreshold(f){
  return (evaluable(f) && Object.prototype.hasOwnProperty.call(T, f.filter_id)) ? T[f.filter_id] : f.threshold;
}
function gateReason(f,value,verdict,reason,threshold){
  const measured = (value!==null&&value!==undefined);
  if(reason===1) return 'the run exports neither '+f.filter_id+'_observed nor '+f.column;
  if(reason===2) return (f.read_column||f.column)+' is empty';
  // An off gate still measured something on most runs. Saying only "filter is off" hid a real
  // per-guide number -- and max_mirna_1mm_seed's own definition promises the number is reported.
  if(reason===3) return measured ? 'filter is off; value measured, nothing acts on it'
                                 : 'filter is off, and the run exports no value for '+f.column;
  if(reason===4) return measured ? 'no threshold declared; value measured, nothing acts on it'
                                 : 'no threshold declared, and no value exported';
  // The run's own two non-decisions on a gate it DID hold in force. Both ship under the unknown verdict
  // code, because "applied, evidence unavailable" is not "not applied" (payload.py:_evaluate, #103), and
  // only the reason code tells them apart -- so only this column can. Reasons 5 and 6 used to fall
  // through to the comparison line below, which printed a comparison that was never made and, for
  // reason 6, one whose arithmetic contradicts the pill beside it: the run records a real measured value
  // with no verdict, so `5 not ge 1` appeared next to `unknown` on a gate 5 does satisfy. A row that
  // states a comparison the run did not make is the fabricated-evidence direction, in words.
  if(reason===5) return 'the run recorded unknown for this gate; no value to re-compare';
  if(reason===6) return measured ? 'the run did not evaluate this gate; value measured, no verdict applied'
                                 : 'the run did not evaluate this gate';
  // States the comparison that produced `verdict`, against the threshold that produced it -- and names
  // the run's own value whenever that is not the same number, so the reader's threshold is never
  // silently substituted for the run's in a sentence the reader will quote.
  const t = threshold===undefined ? liveThreshold(f) : threshold;
  const cmp=fmt(value)+(verdict===0?' ':' not ')+f.comparator+' '+fmt(t)
    +(t!==f.threshold ? ' at your threshold; the run used '+f.comparator+' '+fmt(f.threshold) : '');
  return verdict===4 ? cmp+'; action=warn, not a rejection' : cmp;
}
// Reads the guide's LIVE gate table (g._live, kept current by recomputeLive()), not the frozen
// payload snapshot -- otherwise this panel would go stale the moment a reader moves a threshold and
// contradict the status pill sitting right above it.
function gatesCard(g){
  const live=g._live, order=[3,0,1,4,2];  // fail, unknown, warn, pass, not_evaluated -- worst news first
  const idx=live.gates.map((t,i)=>i).sort((a,b)=>order[live.gates[a][1]]-order[live.gates[b][1]]);
  const rows=idx.map(i=>{const [value,verdict,reason]=live.gates[i], f=FILTERS[i], v=VERDICT[verdict];
    const t=liveThreshold(f), moved=t!==f.threshold;
    return `<tr><td class="mono">${esc(f.filter_id)}</td>
    <td><span class="pill v-${v}">${v.replace('_',' ')}</span></td>
    <td>${fmt(value)}</td><td class="mono">${esc(f.comparator)} ${fmt(t)}${
      moved?` <span class="empty">(yours; run ${fmt(f.threshold)})</span>`:''}</td>
    <td>${esc(f.scope_label)}</td><td>${esc(f.stage)}</td><td>${esc(gateReason(f,value,verdict,reason,t))}</td></tr>`;}).join('');
  // The banner below counts UNKNOWN verdicts, and after #103 that set includes an in-force gate the run
  // itself recorded `not_evaluated`. It said "N gates not evaluated" over rows whose pills read `unknown`,
  // while the rows whose pills literally read `not evaluated` -- off, or no declared threshold -- are the
  // ones it does NOT count. "Undecided" is the one word for the set the number is actually over.
  const nu=live.n_gates_unknown;
  // Whose thresholds these verdicts came from, said once at the top of the panel: a reader quoting a
  // pill has to be able to see that it is theirs and not the run's without reading every row.
  const nm=FILTERS.filter(f=>liveThreshold(f)!==f.threshold).length;
  return `<div class="card"><h2>Gates — all ${live.gates.length}, independently evaluated</h2>
    ${nm?`<div class="warn"><b>${nm} threshold${nm===1?'':'s'} moved from the run's.</b> Every verdict
      below is re-derived at your value; the run's own threshold is shown beside it and named again in
      the Why column, so nothing here is the run's verdict unless it says so.</div>`:''}
    ${nu?`<div class="warn"><b>${nu} of ${live.gates.length} gates undecided.</b> An undecided gate is
      not a pass. See the Why column for which non-decision it was.</div>`:''}
    ${(!nu&&live.status==='unknown')?`<div class="warn"><b>Run verdict: ${esc(g.run_verdict)}.</b>
      No declared filter covers this rejection, so it cannot be re-derived from the ${live.gates.length}
      gates below, all of which are satisfied. Reported as not established.</div>`:''}
    <table><thead><tr><th>Filter</th><th>Verdict</th><th>Value</th><th>Threshold</th>
    <th>Scope</th><th>Stage</th><th>Why</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}
// Distinct isoforms carrying the guide, over every isoform the run enumerated. Deliberately not the
// enumeration count: a site can occur twice in one transcript, so the two differ.
function isoformFrac(g){
  const n=g.transcript_hits, d=TX_IDS.length;
  return (n===null||n===undefined||!d) ? '—' : `${n}/${d}`;
}

function isoformStrip(g){
  const hit=new Set(g.isoforms.map(i=>i.transcript));
  return TX_IDS.map(t=>`<span title="${esc(t)}${hit.has(t)?' — carries this guide':''}"
    style="display:inline-block;width:13px;height:13px;margin-right:3px;border-radius:2px;
    background:${hit.has(t)?'#15803d':'#e5e7eb'}"></span>`).join('');
}

function isoformCard(g){
  const any=g.isoforms.some(i=>i.register_neighbours.length);
  const n=g.transcript_hits, d=TX_IDS.length;
  const pct=(n!==null&&n!==undefined&&d)?` (${(100*n/d).toFixed(0)}%)`:'';
  const extra=g.isoforms.length-(n||0);
  return `<div class="card"><h2>Isoform coverage</h2>
   <div class="kv" style="margin-bottom:12px">
     <div><b>Isoforms carrying this guide</b><span class="big" style="font-size:17px">${isoformFrac(g)}</span>${pct}</div>
     <div><b>Enumerations</b>${g.isoforms.length}${extra>0?` <span class="pill v-not_evaluated">${extra} repeat site${extra===1?'':'s'}</span>`:''}</div>
     <div><b>Isoforms in this run</b>${d}</div>
   </div>
   <div style="margin-bottom:12px">${isoformStrip(g)}
     <div style="font-size:11px;color:#6b7280;margin-top:5px">One square per isoform in the run;
     filled where this guide's site occurs. Hover for the transcript id.</div></div>
   ${any?`<div class="warn"><b>Register neighbour.</b> Another design starts within
     ${REGISTER_NT} nt on the same transcript. The two overlap and can score identically; they are
     still distinct designs and are not interchangeable.</div>`:''}
   <table><thead><tr><th>Candidate id</th><th>Transcript</th><th>Position</th><th>Register neighbours</th></tr></thead>
   <tbody>${g.isoforms.map(i=>`<tr><td class="mono">${esc(i.candidate_id)}</td><td class="mono">${esc(i.transcript)}</td>
     <td>${fmt(i.position)}</td><td>${i.register_neighbours.length?`<span class="pill v-unknown">${i.register_neighbours.join(', ')}</span>`:'—'}</td></tr>`).join('')}</tbody></table></div>`;
}
function offtargetCard(g,idx){
  if(!g.offtarget_by_symbol.length) return `<div class="card"><h2>Off-targets</h2>
    <p class="empty">No alignments for this guide in the published hit table.</p></div>`;
  const rows=g.offtarget_by_symbol.map(s=>`<tr><td class="mono">${esc(s.symbol)}</td>
    <td><span class="pill ${s.is_liability?'v-fail':'v-not_evaluated'}">${esc(s.hit_class)}</span></td>
    <td>${esc(s.species)}</td><td>${s.n}</td><td>${fmt(s.min_nm)}</td></tr>`).join('');
  // Finding 8. The payload's `offtarget_embedded_scope_empty_but_counts_exist` flag is true for every
  // guide whose liabilities all sit outside the embedded scope -- 4,616 of 5,000 guides on one internal
  // run, because 197,727 of that run's 232,864 liability alignments are at nm >= 3. A warning on 92% of
  // rows is not a warning: it trains a reader to skip it and buries the guides it matters for. The size
  // of the scope is now stated once, at run level, in the header. Here the banner is reserved for the
  // case that could change a decision -- a guide this report is NOT rejecting whose liability the file
  // cannot itemise -- and for a guide already being rejected the same fact is stated plainly instead.
  const outside = g.liability_count - g.offtarget_rows.length;
  const shortlistable = g._live.status==='pass' || g._live.status==='warn';
  const note = outside <= 0 ? ''
    : shortlistable
      ? `<div class="warn"><b>Nothing here rejects this guide, and ${outside} of its
         ${g.liability_count} liability alignments lie outside the embedded scope of ${esc(EMBED_SCOPE)}.</b> They are counted,
         charted and gated, but this file holds no row for them — so the table below is not the whole
         picture for a guide you could shortlist.</div>`
      : `<p class="empty">${outside} of this guide's ${g.liability_count} liability alignments lie outside
         the embedded scope of ${esc(EMBED_SCOPE)}: counted, charted and gated, but not itemised here. The header states how
         much of this run sits outside that scope.</p>`;
  const split = SCREENED_SPECIES.filter(sp=>liabOf(g,sp)).map(sp=>`${esc(sp)} ${liabOf(g,sp)}`).join(', ');
  return `<div class="card"><h2>Off-targets by gene — ${g.liability_count} liabilit${g.liability_count===1?'y':'ies'}
      of ${g.offtarget_by_symbol.reduce((a,s)=>a+s.n,0)} alignments</h2>
    ${split?`<p class="empty" style="margin:-6px 0 10px">By species: ${split}. Only ${esc(EMBED_SPECIES)}
      is a human liability; the rest is a non-human-index finding, and
      <code>max_off_target_count</code> compares the total of both.</p>`:''}
    ${note}<div id="mx${idx}" style="height:230px"></div>
    <table><thead><tr><th>Gene</th><th>Class</th><th>Species</th><th>Alignments</th><th>Best nm</th></tr></thead>
    <tbody>${rows}</tbody></table>
    ${g.offtarget_rows.length?`<h2 style="margin-top:18px">Liability alignments — ${g.offtarget_rows.length}
      of ${g.liability_count}, the ones inside the embedded scope of ${esc(EMBED_SCOPE)}</h2>
      <table><thead><tr><th>Transcript</th><th>Gene</th><th>Class</th><th>nm</th><th>Seed mm</th><th>CIGAR</th><th>Strand</th></tr></thead>
      <tbody>${g.offtarget_rows.map(r=>`<tr><td class="mono">${esc(r.t)}</td><td class="mono">${esc(r.s)}</td>
        <td><span class="pill v-fail">${esc(r.c)}</span></td><td>${fmt(r.nm)}</td><td>${fmt(r.sm)}</td>
        <td class="mono">${esc(r.cig)}</td><td>${esc(r.st)}</td></tr>`).join('')}</tbody></table>`
      :`<p class="empty" style="margin-top:14px">No liability alignment inside the embedded scope of ${esc(EMBED_SCOPE)}.
        Non-liability classes are grouped above and counted completely in the chart.</p>`}</div>`;
}
function mirnaCard(g){
  if(!g.mirna.length) return `<div class="card"><h2>miRNA seed resemblance</h2>
    <p class="empty">No named miRNA seed hit at 0 or 1 seed mismatch.</p></div>`;
  const perfect=g.mirna.filter(m=>m.seed_mismatches===0);
  const offset=g.mirna.filter(m=>m.coord!==null&&m.coord!==undefined&&m.coord!==1);
  return `<div class="card"><h2>miRNA seed resemblance — ${perfect.length} perfect,
    ${g.mirna.length-perfect.length} at 1 mismatch</h2>
   ${offset.length?`<div class="warn"><b>${offset.length} hit(s) do not start at seed position 1.</b>
     A guide-seed motif matching elsewhere on a miRNA is not a seed match; that was a real defect in
     this scanner once, so the coordinate is shown rather than assumed.</div>`:''}
   <table><thead><tr><th>miRNA</th><th>Seed mm</th><th>nm</th><th>Coord</th><th>Species</th>
   <th>Database</th></tr></thead>
   <tbody>${g.mirna.map(m=>`<tr><td class="mono">${esc(m.mirna)||'<span class="empty">unnamed</span>'}</td>
     <td>${m.seed_mismatches===0?'<span class="pill v-fail">0</span>':fmt(m.seed_mismatches)}</td>
     <td>${fmt(m.nm)}</td><td>${fmt(m.coord)}</td><td>${esc(m.source)}</td>
     <td>${esc(m.database)}</td></tr>`).join('')}</tbody></table></div>`;
}

// ---- cross-species conservation: the ortholog class, read as evidence rather than discarded ----
const CONS_SPECIES = ['mouse','macaque','rat'];

function consOf(g, sp){ return (g.ortholog||{})[sp] || null; }

// "Usable in this species" is a reader's call, not the tool's: nm=0, or nm<=1 with the seed intact.
function consTier(d){
  if(!d || d.nm===null || d.nm===undefined) return null;
  if(d.nm===0) return 'perfect';
  if(d.nm<=1 && d.seed_mismatches===0) return 'seed-intact 1 mm';
  return `${d.nm} mm`;
}

function conservationCard(g){
  const seen=CONS_SPECIES.filter(sp=>consOf(g,sp));
  if(!seen.length) return `<div class="card"><h2>Cross-species conservation</h2>
    <p class="empty">No ortholog alignment for this guide in any screened species. Either the site is
    human-specific, or no ortholog index was available — the screen cannot tell those apart.</p></div>`;
  const cell=sp=>{const d=consOf(g,sp), t=consTier(d);
    if(!d) return `<div><b>${sp}</b><span class="empty">no ortholog hit</span></div>`;
    const cls = t==='perfect' ? 'v-pass' : (t&&t.startsWith('seed-intact') ? 'v-warn' : 'v-fail');
    return `<div><b>${sp}</b><span class="pill ${cls}">${esc(t)}</span>
      <span class="empty"> nm ${fmt(d.nm)}, seed ${fmt(d.seed_mismatches)}</span></div>`;};
  return `<div class="card"><h2>Cross-species conservation</h2>
    <div class="kv">${CONS_SPECIES.map(cell).join('')}</div>
    <p class="empty" style="margin:10px 0 0">Best ortholog alignment per species, from the off-target
    screen's own <code>ortholog</code> class. This is transcript-level conservation of the target site,
    not a genomic claim, and it depends on the ortholog gene mapping. A mismatch outside guide
    positions 2&ndash;8 leaves the seed intact; whether that is acceptable is a programme decision, so
    it is offered as a threshold above rather than enforced as a gate.</p></div>`;
}
// ---- design map: a Python-drawn backdrop, and dots the reader's own thresholds colour ------------
// The backdrop -- regions, gridlines, axis, uncovered shading -- is drawn in Python by reporting.tracks,
// so a notebook figure and this card cannot diverge. The POINTS are not: a dot's colour is a verdict,
// and the largest visual on the page was the last thing on it still frozen at the run's thresholds. Set
// one gate so 39 of 5,000 guides pass and the legend still read `PASS, all gates evaluated (394)`,
// byte-identical, over 394 green dots -- #103's own failure mode, in the one panel a reader looks at
// first. Both the dots and the legend are now painted from each guide's live gate table.
const TXS = Object.keys(MAPS);
//: `txPinned` is whether the READER chose this isoform. The picker defaults to whatever `_picker_order`
//: puts first -- the canonical one, where the run recorded it -- and `showOn` then followed each new
//: selection to an isoform carrying it, which silently re-pointed the map away from the transcript the
//: reader had deliberately picked. A pinned isoform stays: the note says the selection is not enumerated
//: there, and the button beside the picker is how you leave it.
let curTx = TXS[0] || null, txPinned = false;
const MAP_FILL = {};
MAP_LEGEND.forEach(e => { MAP_FILL[e.key] = e.fill; });

// A window whose guide is not embedded has no live verdict to read, so it keeps its own key rather than
// borrowing a colour: it is the one class on this map a moved threshold cannot answer for.
function mapClass(pt){ const i=pt[2]; return i<0 ? NOT_EMBEDDED : G[i]._live.status; }

// The primitive's own two scales, against the plot area and the value bounds it rounded its gridlines
// to and handed over in `geometry`. Re-deriving either here would put a second copy of a scale in this
// file, which is exactly how a dot and its own gridline come to disagree.
function mapX(geo,p){
  return geo.x0+(geo.x1-geo.x0)*(Math.min(Math.max(p,1)-1,geo.length-1))/Math.max(geo.length-1,1);
}
function mapY(geo,v){
  const lo=geo.value_min, hi=geo.value_max;
  return geo.y1-(geo.y1-geo.y0)*(Math.min(Math.max(v,lo),hi)-lo)/((hi-lo)||1);
}

// `canonical` is true / false / null, and null means this run recorded canonical status NOWHERE. An
// explicit === null check, because `!canonical` would label an unknown transcript as not canonical --
// inventing the fact the run declined to state.
// `mp` and never `m` for a map entry, deliberately: every property the client reads off one is spelled
// `mp.<key>`, so tests/unit/test_reporting_rethreshold.py can check each against the keys _design_maps
// actually emits. A renamed or dropped key used to reach the browser as `undefined` painted into the
// page -- `MAPS[t].legend` outlived the Python that built it and printed the word "undefined" where the
// legend belongs, with nothing failing.
function txLabel(t){
  const mp=MAPS[t];
  return esc(mp.label) + (mp.canonical===true ? ' \u2022 canonical' : (mp.canonical===null ? '' : ' \u2022 not canonical'));
}

function renderMap(){
  const sel=document.getElementById('tx');
  if(!TXS.length){ document.getElementById('mapcard').style.display='none'; return; }
  const carries=new Set((G.find(x=>x.guide===selected)||{isoforms:[]}).isoforms.map(i=>i.transcript));
  sel.innerHTML = TXS.map(t=>`<option value="${esc(t)}"${t===curTx?' selected':''}>${
    txLabel(t)}${carries.has(t)?' \u25cf carries this guide':''}</option>`).join('');
  const stayed = txPinned && selected && !carries.has(curTx);
  document.getElementById('txnote').textContent = selected
    ? ` \u2014 ${carries.size} of ${TXS.length} shown isoforms carry the selected guide.`
      + (stayed ? ' The map stayed on the isoform you chose, which is not one of them.' : '') : '';
  // The escape hatch from a pinned isoform, offered only when the pin is costing the reader the marker.
  document.getElementById('txfollow').style.display = stayed ? '' : 'none';
  const mp=MAPS[curTx];
  document.getElementById('mapbase').innerHTML=mp.svg;
  document.getElementById('mapover').setAttribute('viewBox', mp.viewBox);
  drawOverlay();
}

// Painted from the classes just drawn, in draw order, so the legend cannot describe a colouring the map
// no longer has. A class with no window on this transcript is omitted -- a key for an absent class
// reads as a class the reader failed to find -- and the gap key is a property of the transcript, not a
// count of points, so it is offered on its own terms.
function renderMapLegend(byClass, mp){
  const swatch = fill => `<span style="display:inline-block;width:9px;height:9px;border-radius:2px;`
    + (fill ? `background:${fill}` : `background:#f8fafc;border:1px solid var(--line)`)
    + `;margin-right:5px"></span>`;
  const items = MAP_LEGEND
    .filter(e => e.key===GAPS_KEY ? mp.gaps : (byClass[e.key]||[]).length)
    .map(e => `<span style="margin-right:14px;white-space:nowrap">${swatch(e.fill)}${esc(e.label)}${
      e.key===GAPS_KEY ? '' : ` (${byClass[e.key].length.toLocaleString()})`}</span>`).join('');
  const moved=FILTERS.filter(f=>liveThreshold(f)!==f.threshold).length;
  document.getElementById('maplegend').innerHTML =
    `<div style="font-size:11.5px;color:var(--mut);margin:6px 0 2px">${items}</div>`
    + (moved ? `<div style="font-size:11.5px;color:var(--acc);margin:0 0 2px"><b>Re-coloured at your
        thresholds</b> \u2014 ${moved} moved from the run's, so no colour or count here is the run's own.
        Every number counts candidate windows on this transcript, not guides.</div>`
             : `<div style="font-size:11.5px;color:var(--mut);margin:0 0 2px">Counts are candidate
        windows on this transcript, not guides, at the run's own thresholds.</div>`);
}

// Where this guide sits on the transcript now shown -- a guide is enumerated once per transcript, so
// a guide absent from the current one is a fact worth saying rather than an empty overlay.
function selectedPositions(){
  const g=G.find(x=>x.guide===selected);
  if(!g||!curTx) return [];
  return g.isoforms.filter(i=>i.transcript===curTx && i.position!==null).map(i=>i.position);
}

function cartPositions(){
  const out=[];
  for(const k of CART){ const g=G.find(x=>x.guide===k); if(!g) continue;
    for(const i of g.isoforms) if(i.transcript===curTx && i.position!==null) out.push(i.position); }
  return out;
}

function drawOverlay(){
  const over=document.getElementById('mapover'); if(!over||!curTx) return;
  const mp=MAPS[curTx], geo=mp.geometry, ps=selectedPositions(), cart=cartPositions();
  // One path per class rather than one element per point: 40,000 <circle> nodes cost megabytes and
  // seconds of layout, while the same marks as subpaths of one path cost neither. Legend order is draw
  // order, so a passing window is never hidden underneath a rejected one.
  const byClass={};
  for(const pt of mp.points){ const k=mapClass(pt); (byClass[k]=byClass[k]||[]).push(pt); }
  const dots = MAP_LEGEND.filter(e=>e.key!==GAPS_KEY && (byClass[e.key]||[]).length).map(e=>
    `<path d="${byClass[e.key].map(pt=>
      `M${mapX(geo,pt[0]).toFixed(1)} ${mapY(geo,pt[1]).toFixed(1)}h2v2h-2z`).join('')}"
      fill="${e.fill}" fill-opacity="0.85"/>`).join('');
  // Cart members next, so the selection marker draws over them rather than under.
  over.innerHTML = dots
    + cart.map(p=>`<path d="M${mapX(geo,p).toFixed(1)} ${geo.y1}l-3.5 7h7z" fill="#15803d"/>`).join('')
    + ps.map(p=>{const px=mapX(geo,p).toFixed(1);
    return `<line x1="${px}" y1="${geo.y0}" x2="${px}" y2="${geo.y1}" stroke="#1d4ed8" stroke-width="1.4" stroke-dasharray="3 2"/>`
         + `<path d="M${px} ${geo.y0-2}l-4.5-7h9z" fill="#1d4ed8"/>`;}).join('');
  renderMapLegend(byClass, mp);
  const bits=[];
  if(selected) bits.push(ps.length
    ? `${selected} is enumerated at ${ps.map(p=>p.toLocaleString()).join(', ')} on this transcript.`
    : `${selected} is not enumerated on this transcript.`);
  if(cart.length) bits.push(`${cart.length} cart position${cart.length===1?'':'s'} marked below the axis.`);
  if(mp.gap_note) bits.push(mp.gap_note);
  document.getElementById('mapnote').textContent = bits.join(' ');
}

function showOn(g){
  // Follow the selection to an isoform that actually carries it, so the marker is never off-screen,
  // then re-label the picker: which isoforms carry the guide is part of the selection, not the map.
  // Unless the reader chose the isoform, in which case following would take the map off a transcript they
  // picked on purpose -- and since the run's canonical is only the picker's DEFAULT, that made the
  // canonical one impossible to hold on to across two selections.
  if(!curTx) return;
  const here=g.isoforms.some(i=>i.transcript===curTx);
  if(!here && !txPinned){ const first=g.isoforms.find(i=>MAPS[i.transcript]); if(first) curTx=first.transcript; }
  renderMap();
}

// ---- secondary structure: the run's own dot-bracket, laid out, never refolded ----
function balanced(db){
  let depth=0;
  for(const c of db){ if(c==='(')depth++; else if(c===')'&&--depth<0) return false; }
  return depth===0;
}

function pairsOf(db){
  const p=new Array(db.length).fill(-1), st=[];
  for(let i=0;i<db.length;i++){
    if(db[i]==='(') st.push(i);
    else if(db[i]===')'){ const j=st.pop(); if(j!==undefined){ p[i]=j; p[j]=i; } }
  }
  return p;
}

function structureSvg(seq, db){
  const W=380, H=230, pad=14, foot=14, SEED=[2,8];
  const p=pairsOf(db), xy=LAYOUTS[db];
  let body='';
  if(xy && xy.length===db.length){
    const xs=xy.map(a=>a[0]), ys=xy.map(a=>a[1]);
    const xl=Math.min(...xs), xh=Math.max(...xs), yl=Math.min(...ys), yh=Math.max(...ys);
    const sc=Math.min((W-2*pad)/Math.max(xh-xl,1e-6),(H-foot-pad)/Math.max(yh-yl,1e-6));
    const dx=pad+((W-2*pad)-(xh-xl)*sc)/2, dy=pad+((H-foot-pad)-(yh-yl)*sc)/2;
    const X=i=>(dx+(xy[i][0]-xl)*sc), Y=i=>(dy+(yh-xy[i][1])*sc);
    body+=`<path d="${xy.map((_,i)=>(i?'L':'M')+X(i).toFixed(1)+' '+Y(i).toFixed(1)).join('')}" fill="none" stroke="#9ca3af" stroke-width="1.6"/>`;
    for(let i=0;i<p.length;i++) if(p[i]>i)
      body+=`<line x1="${X(i).toFixed(1)}" y1="${Y(i).toFixed(1)}" x2="${X(p[i]).toFixed(1)}" y2="${Y(p[i]).toFixed(1)}" stroke="#1d4ed8" stroke-width="1.2" stroke-opacity="0.55"/>`;
    for(let i=0;i<seq.length;i++){
      if(i+1>=SEED[0]&&i+1<=SEED[1]) body+=`<circle cx="${X(i).toFixed(1)}" cy="${Y(i).toFixed(1)}" r="7" fill="#fef3c7"/>`;
      body+=`<text x="${X(i).toFixed(1)}" y="${(Y(i)+3.5).toFixed(1)}" text-anchor="middle" font-size="9.5" fill="#1a1d21">${esc(seq[i]||'')}</text>`;
    }
  } else {
    const n=seq.length, base=H-foot-24, step=(W-2*pad)/Math.max(n-1,1), px=i=>pad+i*step;
    body+=`<rect x="${(px(SEED[0]-1)-step/2).toFixed(1)}" y="${base-12}" width="${((SEED[1]-SEED[0]+1)*step).toFixed(1)}" height="20" fill="#fef3c7"/>`;
    body+=`<line x1="${pad}" y1="${base}" x2="${W-pad}" y2="${base}" stroke="#9ca3af" stroke-width="1.2"/>`;
    for(let i=0;i<p.length;i++) if(p[i]>i){
      const r=(px(p[i])-px(i))/2, top=Math.max(base-r,14);
      body+=`<path d="M${px(i).toFixed(1)} ${base}Q${((px(i)+px(p[i]))/2).toFixed(1)} ${(2*top-base).toFixed(1)} ${px(p[i]).toFixed(1)} ${base}" fill="none" stroke="#1d4ed8" stroke-width="1.2" stroke-opacity="0.55"/>`;
    }
    for(let i=0;i<n;i++) body+=`<text x="${px(i).toFixed(1)}" y="${base+14}" text-anchor="middle" font-size="9.5" fill="#1a1d21">${esc(seq[i])}</text>`;
  }
  return `<svg viewBox="0 0 ${W} ${H}" width="100%" style="max-width:${W}px" role="img">${body}</svg>`;
}

function structureCard(g){
  const db=g.structure;
  if(!db) return '';
  // structure.py refuses a dot-bracket that cannot describe this guide's structure; refuse the same
  // here rather than indexing past the end of a coordinate array, which threw inside the template
  // literal and took the entire detail pane -- gates, isoforms, off-targets -- down with it.
  if(db.length!==g.guide.length || /[^.()]/.test(db) || !balanced(db))
    return `<div class="card"><h2>Guide secondary structure</h2>
      <div class="warn"><b>Not drawn.</b> The published <code>structure</code> is ${db.length} nt against
      a ${g.guide.length} nt guide, or is not a balanced dot-bracket, so it cannot be laid out against
      this sequence. Value as published: <code class="mono">${esc(db)}</code></div></div>`;
  const pairCount=(db.match(/\(/g)||[]).length, degenerate=/^\.+$/.test(db);
  const laid=!!(LAYOUTS[db] && LAYOUTS[db].length===db.length);
  // No run-specific figures in this text: the card ships with the tool and is read against every target.
  const caption = degenerate
    ? `Open chain over ${db.length} nt: mfe 0 is the minimum-free-energy floor for an unstructured oligo of this length.`
    : `${pairCount} base pair${pairCount===1?'':'s'}, ${2*pairCount} of ${db.length} nt paired.`;
  return `<div class="card"><h2>Guide secondary structure</h2>
    ${degenerate?`<div class="warn"><b>No base pairs predicted.</b> ${caption}
      <code>paired_fraction</code> and <code>EXCESS_PAIRING</code> are both derived from this
      structure.</div>`:''}
    <div class="struct">
      <div>${structureSvg(g.guide, db)}</div>
      <div>
        <div class="db">${[...db].map((c,i)=>`<span style="color:${c==='.'?'#9ca3af':'#1d4ed8'}">${c}</span>`).join('')}</div>
        <div class="kv" style="margin-top:10px">
          <div><b>Pairs</b>${pairCount}</div>
          <div><b>paired_fraction</b>${fmt(g.metrics.paired_fraction)}</div>
          <div><b>mfe</b>${fmt(g.metrics.mfe)}</div>
          <div><b>Layout</b>${laid?'ViennaRNA naview':'arc diagram (no layout in payload)'}</div>
        </div>
        <p class="empty" style="margin:10px 0 0">Laid out from the <code>structure</code> column as
        published by the run; not refolded, so the drawing and <code>paired_fraction</code> describe the
        same structure. Guide positions 2&ndash;8 (seed) shaded.</p>
      </div>
    </div></div>`;
}

// Wired in wireUI(), at the bottom, with every other DOM-touching statement in this file.

// One card raising must cost that card, not the pane. Every panel below is independent evidence, so
// losing the off-target table because a structure string was malformed is never the right trade.
function card(render, g){
  try { return render(g); }
  catch(e){ return `<div class="card"><h2>Panel failed to render</h2>
    <div class="warn"><b>${esc(e && e.message || String(e))}</b> The other panels are unaffected.</div></div>`; }
}

// ---- thresholds and cart -------------------------------------------------------------------------
// Every evaluable filter (#103 follow-up) gets a live control: moving it re-derives that guide's
// gate, its status, and the counts the index and gate panel show. There is no frozen "selects among
// the run's own verdicts" reading left -- a filter with no control (action=off, or no declared
// threshold) cannot be asked to move at all, and the evaluator below refuses it whatever a URL
// fragment names, with the refusal shown rather than swallowed.
const CART = new Set();
const STATUS_KEYS = ['pass','warn','unknown','fail'];
//: The status set a report opens on. Named rather than spelled twice, because encodeHash has to know
//: what the default IS: it omits the status set when it is this one, so a URL carries a status filter
//: only when the reader chose it -- and `Reset` therefore leaves an empty fragment.
const DEFAULT_STATUS = ['pass','warn'];
// Every reader-side selection that is not a threshold: the search text, the status set, the
// conservation species and the seed tolerance. `q` lives here rather than being read off the input on
// demand, so there is one source of truth for it and the fragment can restore it -- reading the DOM
// meant a URL could carry a search the box did not show, or a box the URL did not carry.
const F = {status:new Set(DEFAULT_STATUS), cons:new Set(), consSeedIntact:true, q:''};

// A filter a reader may move: payload.py's own `evaluable` (payload.py's _filter_view), decided from
// this run's output -- applied, thresholded, and this run measured at least one value for it. Read
// as data rather than re-derived, so a filter whose column is exported but empty for every guide
// freezes here exactly as it does server-side, with no second copy of that rule to drift from it.
function evaluable(f){ return !!f.evaluable; }

// The ONE comparator table, restated once: FilterComparator.passes (models/policy.py:216-227). Every
// gate here is one lookup into this, never a per-filter branch.
function CMP(comparator, v, t){
  if(comparator==='le') return v<=t;
  if(comparator==='lt') return v<t;
  if(comparator==='ge') return v>=t;
  if(comparator==='gt') return v>t;
  return false;
}

//: filter_id -> the threshold a reader has moved it to. Seeded from the run's own threshold for
//: every evaluable filter, so "yours" and "run" start identical and a reset is just restoring this.
const T = {};
FILTERS.forEach(f => { if(evaluable(f)) T[f.filter_id] = f.threshold; });
let refusedFilters = [];
let activePreset = 'all';

// Re-thresholds ONE gate triple -- the same two rules as payload.py's own `reevaluate_gates`, restated
// here because a client-side evaluator has to restate it, not re-derived some other way. `reason`
// (payload.py:187-196) already recorded *why* this guide's row landed where it did:
// REASON_EMPTY_VALUE has nothing to compare and stays unknown; every other non-OK reason -- an off
// filter, an undeclared threshold, a missing column, or the run's own recorded UNKNOWN/NOT_EVALUATED
// -- is returned byte-for-byte untouched, because a row the run itself did not decide is not this
// evaluator's to decide either. REASON_OK is the one case a comparison actually runs.
function reevaluateGate(f, triple, t){
  const [value, verdict, reason] = triple;
  if(reason === REASON_EMPTY_VALUE) return [value, V_UNKNOWN, REASON_EMPTY_VALUE];
  if(reason !== REASON_OK) return triple;
  const passed = CMP(f.comparator, value, t===undefined?f.threshold:t);
  const code = passed ? V_PASS : (f.action==='warn' ? V_WARN : V_FAIL);
  return [value, code, REASON_OK];
}

// Re-derives one guide's whole gate table plus the counters and status render.py:751-755 otherwise
// freezes at the run's own thresholds. Mirrors GuideEntry.status and .undeclared_run_rejection
// (payload.py:87-134) exactly: a rejection the run made that no gate here expresses stays `unknown`
// no matter how every movable gate is set, because the evaluable gates all being satisfied says
// nothing about the one the run used that this report cannot re-derive.
function reevaluateGates(g, overrides){
  const ov = overrides || {};
  const gates = FILTERS.map((f,i) => reevaluateGate(f, g.gates[i],
    Object.prototype.hasOwnProperty.call(ov, f.filter_id) ? ov[f.filter_id] : f.threshold));
  const n_gates_failed = gates.reduce((n,t)=>n+(t[1]===V_FAIL?1:0), 0);
  const n_gates_unknown = gates.reduce((n,t)=>n+(t[1]===V_UNKNOWN?1:0), 0);
  const n_gates_warned = gates.reduce((n,t)=>n+(t[1]===V_WARN?1:0), 0);
  const undeclared = g.run_verdict!=null && g.run_verdict!=='PASS' && !(n_gates_failed || n_gates_unknown);
  let status;
  if(n_gates_failed) status='fail';
  else if(n_gates_unknown || undeclared) status='unknown';
  else status = n_gates_warned ? 'warn' : 'pass';
  return {gates, n_gates_failed, n_gates_unknown, n_gates_warned, status};
}

// Kept current on every guide after a threshold moves, so the index, the sort, the gate panel and the
// status pill all read the same recomputation and none of them can disagree.
function recomputeLive(){ for(const g of G) g._live = reevaluateGates(g, T); }
recomputeLive();

// ---- off-target-clean: liability_count===0 alone cannot mean "clean" ------------------------------
// `off_target_screened` rides along on every guide (payload.py's GuideEntry, from candidates_all.csv's
// own column) precisely so a never-screened guide and a genuinely clean one are never the same 0.
// liability_count already sums off_target + undetermined hits (payload.py:51), so once screening is
// confirmed, liability_count===0 is the whole answer -- no second copy of that sum needed here.
function offTargetClean(g){ return g.off_target_screened===true && g.liability_count===0; }

// ---- register-deduplicated: read the cluster payload.py already computed cross-guide ---------------
// register_cluster/register_representative (payload.py's _register_clusters) are Python-computed
// once over the whole candidate table, keyed on each candidate row's own score -- not approximated
// again here from register_neighbours, which names positions with no guide behind them. A candidate
// with no numeric position is trivially its own cluster and always carries register_representative
// true, so a guide is a cluster's keeper if any of its isoform rows says so.
function isRegisterRepresentative(g){
  return g.isoforms.length===0 || g.isoforms.some(i=>i.register_representative);
}

// ---- five preset views, each a pure predicate over (g, live) --------------------------------------
// 'passing' filters on live.status alone -- derived fail > unknown > warn > pass, mirroring
// GuideEntry.status -- so it can never surface a `live.status==='unknown'` guide, whatever threshold
// produced it: the property a reader cannot re-threshold their way out of.
const PRESETS = {
  all:              {label:'All',                  test:(g,live)=>true},
  passing:          {label:'Passing',               test:(g,live)=>live.status==='pass'||live.status==='warn'},
  near_miss:        {label:'Near miss',             test:(g,live)=>live.n_gates_failed===1 && live.n_gates_unknown===0},
  off_target_clean: {label:'Off-target clean',      test:(g,live)=>offTargetClean(g)},
  register_dedup:   {label:'Register-deduplicated',  test:(g,live)=>isRegisterRepresentative(g)},
};

// ---- reader filters: select among rows, and never touch a verdict --------------------------------
// Beside the gate controls, deliberately not inside them. A gate control moves a threshold this run
// applied and re-derives that guide's verdict; these three are no gate's quantity, so a bound over one
// can only narrow the list. All three were reader boxes in the first report and are restored as such:
// composite_score and transcript_hits are gated by nothing (min_isoform_coverage is a *fraction of the
// run's transcripts*, a different quantity), and liability_count is otherwise reachable only as
// "exactly zero", through the off-target-clean preset. The other three v1 boxes are gone for good --
// GC, asymmetry and off-target count are real movable gates now, and a second control over the same
// number, obeying a different rule, is how a reader ends up believing the report says two things.
const READER_FILTERS = [
  {k:'composite', label:'composite score', dir:'min', get:g=>g.composite_score},
  {k:'isoforms',  label:'isoforms hit',    dir:'min', get:g=>g.transcript_hits},
  //: Read off the guide, like the index's own Liab. column: moving a gate re-decides verdicts and the
  //: counters over them, never this count, so _live carries no copy of it to prefer instead. `liabIn`
  //: and not `liability_count`, so the species selector beside the box re-scopes the bound too: "at
  //: most 0 human liabilities" is the question a single all-species ceiling cannot ask, and it is the
  //: one that separates a real liability from a non-human-index artefact. Unscoped -- the default --
  //: liabIn IS liability_count, so the bound's meaning does not change until a reader narrows it.
  {k:'liab',      label:'liabilities',     dir:'max', get:g=>liabIn(g)},
];

//: reader filter key -> the bound a reader typed, or null for "not set". A blank box is null and never
//: 0 (fnum's rule), because a 0 floor on liabilities is a real and very different request.
const R = {};
READER_FILTERS.forEach(f => { R[f.k] = null; });

// The same rule the gate evaluator holds to: an absent value cannot satisfy a threshold. A guide with
// no composite score is EXCLUDED by a composite floor rather than admitted because there was nothing
// to compare -- silently keeping it would let a floor report guides that never cleared it.
function passesReaderFilters(g, bounds){
  const b = bounds || R;
  for(const f of READER_FILTERS){
    const t = b[f.k];
    if(t===null||t===undefined) continue;
    const v = f.get(g);
    if(v===null||v===undefined) return false;
    if(f.dir==='min' ? v<t : v>t) return false;
  }
  return true;
}

// Conservation is a selection criterion, not a verdict: requiring mouse and macaque is one
// programme's requirement and would be wrong baked into the tool, so it lives here as a threshold.
function passesConservation(g){
  if(!F.cons.size) return true;
  for(const sp of F.cons){
    const d=consOf(g,sp);
    if(!d || d.nm===null || d.nm===undefined) return false;
    if(d.nm>1) return false;
    if(d.nm===1 && (!F.consSeedIntact || d.seed_mismatches!==0)) return false;
  }
  return true;
}

// The ONE numeric boundary for every control in the panel, and for the URL fragment (#103). Blank is
// `null` -- "no bound", the reader's own way of clearing one. Anything that is not a finite number --
// `abc`, `Infinity`, `1e999` -- returns `undefined`, and every caller REFUSES it rather than passing it
// on: `Number('abc')` is NaN and every comparison against NaN is false, so one junk keystroke would
// fail every re-decidable gate on every guide, and encodeHash would then write a fragment
// applyHashFragment itself refuses -- the reader's own reload could not reproduce what they were
// looking at. The old boxes carried only a soft-keyboard hint, which constrains nothing at all.
function parseControlValue(text){
  const v=String(text===null||text===undefined?'':text).trim();
  if(v==='') return null;
  const n=Number(v);
  return Number.isFinite(n) ? n : undefined;
}
function fnum(id){ const el=document.getElementById(id); if(!el) return null;
  // type=number blanks `value` for text it cannot parse, so validity.badInput is the only trace left
  // of a reader having typed something; without this branch that junk reads as an empty box.
  if(el.validity && el.validity.badInput) return undefined;
  return parseControlValue(el.value); }

// Both refusal banners live inside the thresholds `<details>`, and that element ships closed. On a real
// load of `#t=<gate>:<value>` the banner was in the DOM with `display` cleared and `checkVisibility()`
// still returned false: `elementFromPoint` at its own box returned the index's table header behind it.
// Both banners exist so that a refusal is visible rather than swallowed, and the container was doing the
// swallowing. Anything the reader has to see opens it -- a refused filter, a refused keystroke, or a
// fragment that set any control at all.
function openFilters(){ const d=document.getElementById('filters'); if(d) d.open=true; }

// A refused keystroke has to be visible at the control. Silently keeping the last good value would
// leave the reader looking at a threshold they did not type and believing they had moved it.
function markControl(el, ok){
  el.setAttribute('aria-invalid', ok?'false':'true');
  const note=document.getElementById('badnum'); if(!note) return;
  const bad=[...document.querySelectorAll('#frows input[aria-invalid="true"], #rrows input[aria-invalid="true"]')];
  note.style.display = bad.length ? '' : 'none';
  if(bad.length) openFilters();
  note.innerHTML = bad.length ? `<b>Ignored.</b> ${bad.map(b=>esc(b.getAttribute('aria-label')||b.id)).join(', ')}
    ${bad.length===1?'is':'are'} not a finite number, so nothing moved: no verdict here was decided
    against it, and the URL still describes the thresholds actually in force.` : '';
}

function passesFilters(g){
  const live = g._live;
  if(F.status.size && !F.status.has(live.status)) return false;
  if(!passesConservation(g)) return false;
  if(!passesReaderFilters(g)) return false;
  if(!PRESETS[activePreset].test(g, live)) return false;
  return true;
}

function buildPresetButtons(){
  document.getElementById('presets').innerHTML = Object.entries(PRESETS).map(([k,p])=>
    `<button class="cartbtn" data-preset="${k}"${k===activePreset?' style="border-color:var(--acc);color:var(--acc)"':''}>${esc(p.label)}</button>`).join('');
  document.querySelectorAll('#presets button').forEach(b=>b.onclick=()=>{
    activePreset=b.dataset.preset; buildPresetButtons(); applyFilters(); syncHash();
  });
}

// What one gate actually did to the guides in this file, in words, from the counts payload.py published
// (rejects / sole_rejects / warns / unknowns / inert / inert_reason) rather than a re-derived guess.
// Seventeen identical sliders was the panel the audit found: on one internal run five of the seventeen
// reject anything at all, two of them account for 1,765 of the rejections, and two cannot reject anybody
// because their threshold sits on the limit their own setting declares. Presenting all seventeen as
// equals invited a reader to spend their attention on the ones that decide nothing.
function gateEffect(f){
  if(!evaluable(f))
    return `<b>frozen</b> — ${esc(f.unevaluable_reason||'not re-thresholdable')}`;
  if(f.inert)
    return `<b>inert in this run</b> — ${esc(f.inert_reason||'no guide in this file is failed or flagged by it')}`;
  const bits=[`<b>${(f.rejects||0).toLocaleString()} rejected</b>`];
  if(f.sole_rejects) bits.push(`${f.sole_rejects.toLocaleString()} by this gate alone`);
  if(f.warns) bits.push(`${f.warns.toLocaleString()} warned`);
  if(f.unknowns) bits.push(`${f.unknowns.toLocaleString()} left undecided`);
  return bits.join(', ');
}

// 0 = decided something, 1 = evaluable but decided nothing, 2 = cannot be re-thresholded at all.
function gateTier(f){ return !evaluable(f) ? 2 : (f.inert ? 1 : 0); }

// `sole_rejects` is the ranking key: a guide rejected by one gate alone is a guide that gate is solely
// responsible for losing, so moving that gate is the move that changes the shortlist. Filter order is
// the final tiebreak, so the ranking is stable and reads the same on every reload.
function gateOrder(){
  return FILTERS.map((f,i)=>[f,i]).sort((a,b)=>
       (gateTier(a[0])-gateTier(b[0]))
    || ((b[0].sole_rejects||0)-(a[0].sole_rejects||0))
    || ((b[0].rejects||0)-(a[0].rejects||0))
    || ((b[0].warns||0)-(a[0].warns||0))
    || (a[1]-b[1]));
}

// ---- the box and the slider are ONE threshold ----------------------------------------------------
// The slider's domain is this run's own observed values, snapped to the step and clamped to the
// setting's declared range (payload._control_domain) -- deliberately narrower than what the number box
// can reach, because a domain wider than the run's values spends its travel where no verdict changes.
// A range input cannot represent a value outside its own min/max: it pins silently at the edge. So
// typing 90 into a gate whose domain ends at 76.104348 left the box reading 90, T at 90, and the slider
// sitting at 76.104348 -- and one touch of that slider then dropped the real threshold to the pin
// without the reader asking for it. The code asserted the two could not disagree by keeping them in
// sync on every keystroke; the pin is the case where the sync silently fails.
//
// The domain is NOT re-derived here to make the disagreement go away: widening it would undo the
// clamping that stopped this report publishing "at most -59 off-targets" on a count-valued gate. The
// threshold stays where the reader typed it, the slider says it cannot show it, and it is disabled so it
// cannot answer for a number it does not hold.
function offScale(f){
  const c=f.control||{}, v=T[f.filter_id];
  return Number.isFinite(v) && Number.isFinite(c.min) && Number.isFinite(c.max) && (v<c.min || v>c.max);
}
function offScaleNote(f){
  if(!offScale(f)) return '';
  const c=f.control||{};
  return `<b>Off this slider's scale.</b> ${fmt(T[f.filter_id])} is outside ${fmt(c.min)} to ${fmt(c.max)},
    the span of the values this run measured for this gate, so the slider cannot show it: the box above is
    the threshold in force and every verdict here is decided at <span class="mono">${esc(f.comparator)}
    ${fmt(T[f.filter_id])}</span>. The slider is disabled until the threshold is back inside its scale.`;
}

// Keeps one gate row's two controls telling the same story after a move: whichever the reader touched
// stays as they left it, the other follows, and the slider goes disabled-and-labelled the moment it
// cannot represent T rather than pinning at its own edge.
function paintControl(f, moved){
  const id=f.filter_id, off=offScale(f);
  const box=document.getElementById('ctl_'+id), rng=document.getElementById('rng_'+id);
  if(box && box!==moved) box.value = T[id];
  if(rng){ rng.disabled = off; if(!off && rng!==moved) rng.value = T[id]; }
  const note=document.getElementById('os_'+id);
  if(note) note.innerHTML = offScaleNote(f);
}

// One generic control per filter -- no filter_id and no filter column is ever a branch in this function,
// only data read off FILTERS. A frozen filter gets a read-only row, no input, and payload.py's own
// `unevaluable_reason` (not a re-derived guess at why): the evaluator above already refuses to move it,
// and this is where that refusal is visible rather than a control nobody could use. An INERT gate keeps
// its control, deliberately: it decided nothing at the run's threshold, which is not the same as being
// unable to decide anything at any threshold.
function gateRow(f){
  if(!evaluable(f)) return `<div class="frow" data-filter="${esc(f.filter_id)}">
      <span class="mono">${esc(f.filter_id)}</span>
      <span class="empty" style="grid-column:span 2">${gateEffect(f)}</span></div>`;
  const c = f.control || {};
  return `<div class="frow" data-filter="${esc(f.filter_id)}">
      <span title="${esc(f.definition||'')}"><span class="mono">${esc(f.filter_id)}</span> ${esc(f.comparator)}</span>
      <input id="ctl_${esc(f.filter_id)}" type="number" step="any" value="${T[f.filter_id]}"
        aria-label="${esc(f.filter_id)} threshold">
      <span class="empty">run ${fmt(f.threshold)}</span></div>
      <input type="range" id="rng_${esc(f.filter_id)}" min="${c.min}" max="${c.max}" step="${c.step}"
        value="${T[f.filter_id]}"${offScale(f)?' disabled':''} style="grid-column:1/-1;width:100%"
        aria-label="${esc(f.filter_id)} slider">
      <div class="offscale" id="os_${esc(f.filter_id)}">${offScaleNote(f)}</div>
      <div class="geffect">${gateEffect(f)}</div>`;
}

function buildGateControls(){
  const ordered=gateOrder();
  const loud=ordered.filter(e=>gateTier(e[0])===0), quiet=ordered.filter(e=>gateTier(e[0])>0);
  const inert=quiet.filter(e=>evaluable(e[0])).length;
  document.getElementById('frows').innerHTML = loud.map(e=>gateRow(e[0])).join('')
    + (quiet.length ? `<details><summary style="cursor:pointer;font-size:11.5px;color:var(--mut);
        padding:5px 0">${quiet.length} of ${FILTERS.length} gates decided nothing for the guides in this
        file — ${inert} inert, ${quiet.length-inert} frozen. Each row says why.</summary>
        ${quiet.map(e=>gateRow(e[0])).join('')}</details>` : '');
  // The number box and the slider are the same threshold; typing in one moves the other, so "yours"
  // never has two disagreeing readouts. The box carries `type=number` and `step=any` but deliberately
  // no min/max: the slider's domain bounds the slider, and payload.py's own contract is that the number
  // box can still reach any threshold (payload._control_domain). What is rejected is not an unusual
  // number, it is a value that is not a number at all -- and a number the slider's domain cannot hold is
  // kept, with the slider disabled and saying so, rather than quietly pinned (paintControl).
  document.querySelectorAll('#frows input[id^="ctl_"], #frows input[id^="rng_"]').forEach(el=>el.oninput=()=>{
    const id=el.id.slice(el.id.indexOf('_')+1), f=FILTERS.find(x=>x.filter_id===id), v=fnum(el.id);
    if(v===undefined){ markControl(el,false); return; }   // refused here, so T can never hold a NaN
    markControl(el,true);
    // While the threshold is off the slider's scale the slider's value is a pin, not a reading, and it
    // is disabled so a reader cannot touch it. `disabled` alone is the browser's promise; this is the
    // report's: an input arriving from that control anyway cannot drop the threshold to the pin.
    if(el.id.startsWith('rng_') && offScale(f)){ paintControl(f, null); return; }
    T[id] = v===null ? f.threshold : v;
    paintControl(f, el);
    onThresholdsChanged();
  });
}

// One box per reader filter, reusing the gate rows' own layout so the two groups read as siblings --
// but with no slider and no "run <value>", because there is no run threshold over these to be seeded
// from or reset to. Moving one deliberately does NOT call onThresholdsChanged: nothing is re-derived,
// so recomputing every guide's gate table would be a claim that something changed when nothing did.
function buildReaderFilters(){
  document.getElementById('rrows').innerHTML = READER_FILTERS.map(f=>
    `<div class="frow" data-reader="${esc(f.k)}"><span>${esc(readerLabel(f))}</span>
      <input id="rf_${esc(f.k)}" type="number" step="any" placeholder="${f.dir}"
        value="${R[f.k]===null||R[f.k]===undefined?'':R[f.k]}"
        aria-label="${esc(readerLabel(f))}, ${f.dir==='min'?'at least':'at most'}">
      <span class="empty">${f.dir==='min'?'at least':'at most'}</span></div>`).join('')
    + liabScopeControl();
  document.querySelectorAll('#rrows input').forEach(el=>el.oninput=()=>{
    const v=fnum(el.id);
    if(v===undefined){ markControl(el,false); return; }   // a bound is a number or it is not a bound
    markControl(el,true);
    R[el.id.slice(3)] = v;                 // blank clears the bound; it does not set it to 0
    applyFilters(); syncHash();
  });
  const sel=document.getElementById('liabscope');
  if(sel) sel.onchange = e => {
    liabScope = e.target.value || null;
    // Re-scoping is not a re-derivation: no verdict moves, so onThresholdsChanged is deliberately not
    // called. What does change is the column, its sort, and which guides the liability bound admits.
    labelLiabColumn(); buildReaderFilters(); applyFilters(); syncHash();
  };
}

// The bound's own label says which species it counts, because "at most 3 liabilities" means two
// different things before and after the selector moves and a box that did not say so would be lying.
function readerLabel(f){ return f.k==='liab' ? `liabilities (${scopeLabel()})` : f.label; }

// One selector, three effects: the index column's emphasis and sort, and the liability bound above.
// Offered only when the run has a hit table to decompose -- otherwise there is no species to choose.
function liabScopeControl(){
  if(!SCREENED_SPECIES.length) return '';
  const opt=(v,label)=>`<option value="${esc(v)}"${(liabScope||'')===v?' selected':''}>${esc(label)}</option>`;
  return `<div class="frow" data-reader="liabscope"><span>count liabilities in</span>
    <select id="liabscope" style="grid-column:span 2" aria-label="species the liability count is over">
      ${opt('','all screened species')}${SCREENED_SPECIES.map(sp=>opt(sp,sp)).join('')}</select></div>
    <div class="empty" style="font-size:11.5px;margin:-2px 0 5px">Every liability total on this page is
    over all ${SCREENED_SPECIES.length} screened species, which is what <code>max_off_target_count</code>
    compares. Narrowing to ${esc(EMBED_SPECIES)} asks the different question of whether a rejection is a
    human liability or a non-human-index artefact; it re-scopes the column and the bound, never a verdict.</div>`;
}

// Set from JS, not the template, because the scope is a reader's choice and the species is run data.
function labelLiabColumn(){
  // The label is a child span, not the cell's own text: the cell also holds the sort indicator, and
  // `th.textContent =` would delete it -- a header that stops saying which way it sorts as soon as the
  // reader narrows the species.
  const th=document.getElementById('thliab'), lbl=document.getElementById('liablbl');
  if(!th || !lbl) return;
  lbl.textContent = SCREENED_SPECIES.length ? `Liab. (${scopeLabel()})` : 'Liab.';
  th.title = SCREENED_SPECIES.length
    ? `off-target liabilities: the total over all ${SCREENED_SPECIES.length} screened species -- the count `
      + `max_off_target_count compares -- then the split into ${EMBED_SPECIES} and the rest. Sorted on `
      + scopeLabel() + '.'
    : 'off-target liabilities in every screened species';
}

function buildFilterUI(){
  // The search box is a control like any other, so it is restored from F here rather than left holding
  // whatever the reader last typed: a fragment that carries `q=` has to put the text back in the box, and
  // Reset has to clear it.
  const qbox=document.getElementById('q');
  if(qbox) qbox.value = F.q;
  document.getElementById('fstat').innerHTML = STATUS_KEYS.map(k=>
    `<label><input type="checkbox" data-status="${k}"${F.status.has(k)?' checked':''}> ${k}</label>`).join('');
  // syncHash, because this is the single biggest thing the URL used to leave out: the default set hides
  // every fail and unknown row -- 3,617 of them on one internal run -- so a URL that did not carry it did
  // not carry the rows the reader was looking at.
  document.querySelectorAll('#fstat input').forEach(cb=>cb.onchange=()=>{
    cb.checked ? F.status.add(cb.dataset.status) : F.status.delete(cb.dataset.status);
    applyFilters(); syncHash(); });
  buildGateControls();
  buildReaderFilters();
  labelLiabColumn();
  document.getElementById('fcons').innerHTML =
    '<span style="color:var(--mut)">conserved in</span> ' +
    CONS_SPECIES.map(sp=>`<label><input type="checkbox" data-cons="${sp}"${F.cons.has(sp)?' checked':''}> ${sp}</label>`).join('') +
    `<label title="accept one mismatch outside guide positions 2-8">
      <input type="checkbox" id="consseed"${F.consSeedIntact?' checked':''}> allow 1 mm outside the seed</label>`;
  document.querySelectorAll('#fcons input[data-cons]').forEach(cb=>cb.onchange=()=>{
    cb.checked ? F.cons.add(cb.dataset.cons) : F.cons.delete(cb.dataset.cons); applyFilters(); syncHash(); });
  document.getElementById('consseed').onchange = e => {
    F.consSeedIntact=e.target.checked; applyFilters(); syncHash(); };
  buildPresetButtons();
  renderRefused();
  // Rebuilding the boxes discards every aria-invalid with them, so the refusal banner would otherwise
  // outlive the input it was about -- most visibly after Reset, which replaces all of them.
  const badnum=document.getElementById('badnum');
  if(badnum){ badnum.style.display='none'; badnum.innerHTML=''; }
}

// Everything #freset does, as a function rather than inline in its handler: "Reset" has to mean every
// control -- the status and conservation checkboxes, every gate threshold, every reader filter and the
// preset -- and the parity harness drives this same function rather than a restatement of it.
// The cart is deliberately NOT cleared here: it is a pick list a reader built, not a filter setting, and
// wiping it from a Reset button inside the thresholds panel would destroy work nothing else on this page
// can rebuild. A fragment is a different matter -- see applyFragmentState.
function resetControls(){
  F.status = new Set(DEFAULT_STATUS); F.cons = new Set(); F.consSeedIntact = true; F.q = '';
  FILTERS.forEach(f => { if(evaluable(f)) T[f.filter_id] = f.threshold; });
  READER_FILTERS.forEach(f => { R[f.k] = null; });
  liabScope = null;                        // back to the gate's own all-species scope
  activePreset = 'all'; refusedFilters = [];
}

// The four header pills, over the guides in this file, at the thresholds in force NOW. Seeded from the
// payload by the template so the numbers are right before any script runs, then re-rendered from the
// live gate tables: leaving them at the run's own values beside a re-coloured map would be the frozen
// legend again with the numbers swapped. The agreement clause beside them is deliberately not
// re-derived -- it is a comparison WITH the run, so it belongs at the run's own thresholds -- and says
// so, plus says when those are no longer the thresholds the pills came from.
function renderStatusPills(){
  const n={pass:0,warn:0,unknown:0,fail:0};
  for(const g of G) n[g._live.status]++;
  for(const k of STATUS_KEYS){
    const el=document.getElementById('n_'+k); if(el) el.textContent=n[k].toLocaleString();
  }
  const moved=FILTERS.filter(f=>liveThreshold(f)!==f.threshold).length;
  const basis=document.getElementById('cbasis');
  if(basis) basis.textContent = moved
    ? `, at your thresholds — ${moved} moved from the run's`
    : ", at the run's own thresholds";
  const stale=document.getElementById('agreestale');
  if(stale) stale.innerHTML = moved
    ? ` <b>Not re-derived at your thresholds:</b> an agreement with the run has to be measured against
      the run's own, so these four no longer describe the pills above them.` : '';
}

function matching(){ return G.filter(passesFilters); }

// Reads F.q, not the input: the search text is reader state that the fragment carries and Reset clears,
// so the box is a view of it rather than its home.
function applyFilters(){
  const raw=(F.q||'').trim();
  const seq=raw.toUpperCase().replace(/T/g,'U'), id=raw.toLowerCase();
  view = matching().filter(g => !raw || g.guide.includes(seq)
    || g.isoforms.some(i=>String(i.candidate_id).toLowerCase().includes(id)
                        || String(i.transcript).toLowerCase().includes(id)));
  document.getElementById('fcount').textContent = `${view.length.toLocaleString()} of ${G.length.toLocaleString()} guides`;
  renderIndex();
}

function togglePick(guide){
  CART.has(guide) ? CART.delete(guide) : CART.add(guide);
  renderIndex(); renderCart(); drawOverlay(); syncHash();
}

// ---- where the cart lives between reloads ---------------------------------------------------------
// Not browser storage -- the obvious answer, and the one this report may not use: it reaches nothing
// outside itself, and the sandbox it is read in withholds same-origin, so no storage API is even
// reachable there. The URL fragment is the only in-page place a pick list can persist, and it has the
// property storage does not: a reader can send it. So the cart rides in `k=`, as guide sequences and
// never as indices into G -- a sequence not in this file is refused, where an index would silently
// resolve to a different guide in a differently built
// report. The cap is the cost of that choice: 5,706 guides at 22 characters is a 125 KB address that no
// mail client keeps intact, so the fragment carries the first CART_IN_URL_MAX in the index's own order
// and the cart card says on screen when it is carrying fewer than the cart holds.
const CART_IN_URL_MAX = 500;
function cartInUrl(){ return G.filter(g=>CART.has(g.guide)).map(g=>g.guide).slice(0, CART_IN_URL_MAX); }

function cartRows(){
  return [...CART].map(k=>G.find(g=>g.guide===k)).filter(Boolean)
    .sort((a,b)=>(b.composite_score??-1)-(a.composite_score??-1));
}

// ---- the cart as TSV, with the thresholds its `status` column came from -------------------------
// The columns and their order are the export's contract and are unchanged. What was missing is above
// them: the `status` column carries the reader's LIVE, possibly re-thresholded status, and the file
// recorded nothing about the thresholds that produced it -- so a TSV pasted into a ticket claimed a
// verdict nobody could reproduce from the run (#103). The provenance rides as `#`-prefixed,
// tab-delimited `key=value` comment lines, the one shape that is both readable in a ticket and
// skippable by every TSV reader, so a consumer that ignores comments still gets exactly the old file.
const CART_COLUMNS=['guide','passenger','status','composite_score','design_score','isoforms_hit',
  'isoforms_in_run','gc_content','asymmetry_score','off_target_count','liabilities','structure',
  ...CONS_SPECIES.flatMap(sp=>[sp+'_nm', sp+'_seed_mm'])];

// A comment field is one tab-delimited cell, so a gene query carrying a tab or a newline would split
// it into a line the reader would mis-parse. Filter ids and numbers cannot; free text can.
function tsvField(s){ return String(s===null||s===undefined?'':s).replace(/[\t\r\n]+/g,' '); }

// Also records the reader filters and the preset, not because they touch `status` -- they cannot -- but
// because `Add top n to cart` picks from the filtered view, so they decide which guides are in the file.
function cartProvenance(rows){
  const moved=FILTERS.filter(f=>evaluable(f) && liveThreshold(f)!==f.threshold);
  const lines=['#sirnaforge_cart\tschema=1\tgene='+tsvField(GENE)+'\tguides='+rows.length
    +'\tpreset='+activePreset,
    '#status_basis='+(moved.length?'reader_rethresholded':'run_thresholds')+'\tmoved_gates='+moved.length];
  for(const f of moved) lines.push('#moved_gate='+f.filter_id+'\tcomparator='+f.comparator
    +'\trun_threshold='+f.threshold+'\treader_threshold='+liveThreshold(f));
  for(const f of READER_FILTERS) if(Number.isFinite(R[f.k]))
    lines.push('#reader_filter='+f.k+'\tdirection='+f.dir+'\tbound='+R[f.k]);
  return lines;
}

function cartTsv(rows){
  return [...cartProvenance(rows), CART_COLUMNS.join('\t'), ...rows.map(g=>[
    g.guide, g.passenger??'', g._live.status, g.composite_score??'', g.design_score??'',
    g.transcript_hits??'', TX_IDS.length, g.metrics.gc_content??'', g.metrics.asymmetry_score??'',
    g.metrics.off_target_count??'', g.liability_count, g.structure??'',
    ...CONS_SPECIES.flatMap(sp=>{const d=consOf(g,sp); return [d?.nm??'', d?.seed_mismatches??''];})
    ].join('\t'))].join('\n');
}

function renderCart(){
  const rows=cartRows(), card=document.getElementById('cartcard');
  card.style.display = rows.length ? '' : 'none';
  if(!rows.length){ return; }
  document.getElementById('cartn').textContent = `${rows.length} guide${rows.length===1?'':'s'}`;
  document.getElementById('cartlist').innerHTML =
    `<table><thead><tr><th></th><th>Guide</th><th>Status</th><th>Score</th><th>Isoforms</th>
      <th>Liab.</th></tr></thead><tbody>${rows.map(g=>
      `<tr data-cart="${esc(g.guide)}"><td class="pick"><span class="picked">×</span></td>
       <td class="mono">${esc(g.guide)}</td>
       <td><span class="pill ${({pass:'v-pass',warn:'v-warn',unknown:'v-unknown',fail:'v-fail'})[g._live.status]}">${g._live.status}</span></td>
       <td>${fmt(g.composite_score)}</td><td>${isoformFrac(g)}</td>
       <td class="${g.liability_count?'liab':'nonliab'}">${g.liability_count}</td></tr>`).join('')}</tbody></table>`;
  document.querySelectorAll('#cartlist tr[data-cart]').forEach(tr=>
    tr.querySelector('td.pick').onclick=()=>togglePick(tr.dataset.cart));
  document.getElementById('carttsv').value = cartTsv(rows);
  const carried=cartInUrl().length, note=document.getElementById('carturl');
  if(note) note.textContent = carried < rows.length
    ? `This page's URL carries ${carried} of these ${rows.length} guides. The cart is kept in the address —
       nothing is stored in this browser — and ${rows.length-carried} of them do not fit, so a copied URL
       restores ${carried}. Export the TSV to keep all ${rows.length}.`
    : `The cart is in this page's URL: copy the address and these ${rows.length}
       guide${rows.length===1?'':'s'} travel with your thresholds, filters and search. Nothing is stored in
       this browser, so closing the page without the URL loses the cart.`;
}

// ---- URL fragment: the whole reader state, and the only place it lives -----------------------------
// A bare fragment with none of the '=' syntax below is still read as a guide -- what every earlier
// report in this run wrote. Reader filters ride in their own `r=` list, on the same footing as `t=`:
// a reader who sends the URL sends the rows they were looking at, not just the verdicts.
//
// "The rows they were looking at" was false as shipped: `t=`, `r=`, `sp=`, `preset=` and `g=` were
// encoded, and the search box, the four status checkboxes, the conservation species, the seed tolerance
// and the cart were not -- and the status set alone hides every fail and unknown row by default, 3,617
// of them on one internal run. So `q=`, `s=`, `c=` and `k=` are here too, and every control that changes
// any of them calls syncHash. Each is omitted when it holds its default, which is what keeps `Reset` on
// an empty fragment and keeps a URL down to what the reader actually chose.
//
// The token for "1 mismatch outside the seed is NOT acceptable". It rides inside `c=` because it is part
// of one conservation question, and it is spelled out rather than a bare flag so an unknown token in that
// list can be refused like any other name the URL invents.
const SEED_STRICT = 'no1mm';
function statusIsDefault(){
  return F.status.size===DEFAULT_STATUS.length && DEFAULT_STATUS.every(k=>F.status.has(k));
}
function encodeHash(){
  const parts=[];
  if(selected) parts.push('g='+encodeURIComponent(selected));
  // Number.isFinite, not merely "different from the run's": the fragment must never carry a value
  // applyHashFragment would itself refuse, or the reader's own reload cannot reproduce the page the URL
  // was copied from. The controls reject non-finite input before it reaches T or R (fnum), so this is
  // the second half of one rule rather than a new one.
  const moved=Object.keys(T).filter(id=>{
    const f=FILTERS.find(x=>x.filter_id===id); return f && Number.isFinite(T[id]) && T[id]!==f.threshold;
  });
  if(moved.length) parts.push('t='+moved.map(id=>id+':'+T[id]).join(','));
  const bounded=READER_FILTERS.filter(f=>Number.isFinite(R[f.k]));
  if(bounded.length) parts.push('r='+bounded.map(f=>f.k+':'+R[f.k]).join(','));
  // The liability scope rides too: it decides which guides the `liab` bound admits, so a URL without it
  // would restore the bound and not the question it was asked about.
  if(liabScope!==null) parts.push('sp='+encodeURIComponent(liabScope));
  if(activePreset!=='all') parts.push('preset='+activePreset);
  // Canonical order for each set -- the status keys' own order, the species' own order, the index's own
  // order -- so two readers looking at the same view copy the same URL, and a fragment round-trips.
  if(F.q) parts.push('q='+encodeURIComponent(F.q));
  if(!statusIsDefault()) parts.push('s='+STATUS_KEYS.filter(k=>F.status.has(k)).join(','));
  const cons=CONS_SPECIES.filter(sp=>F.cons.has(sp));
  if(cons.length || !F.consSeedIntact)
    parts.push('c='+cons.concat(F.consSeedIntact?[]:[SEED_STRICT]).join(','));
  const carted=cartInUrl();
  if(carted.length) parts.push('k='+carted.join(','));
  return parts.join('&');
}

// `location.hash = ...` is a navigation: it pushed one history entry per keystroke on every threshold, so
// `history.length` climbed 5 -> 7 during light use, Back then changed the URL without changing the view --
// the page showing one state while its URL claimed another -- and escaping the report took dozens of
// presses. replaceState keeps the URL current without owning the reader's Back button, and does not fire
// `hashchange`, so the listener below only ever sees changes the reader made.
// The catch is not defensive noise: pushState and replaceState throw a SecurityError on an opaque origin,
// and the sandbox this report is read in withholds allow-same-origin. Self-navigation to a fragment is
// still allowed there, so the fallback is the old behaviour rather than losing the URL altogether.
function syncHash(){
  if(typeof location==='undefined') return;
  const h=encodeHash();
  try { history.replaceState(null, '', '#'+h); }
  catch(e){ if(location.hash.slice(1)!==h) location.hash = h; }
}

// Reads a fragment written by encodeHash (or a bare guide sequence from an earlier report) and
// returns the names it had to refuse -- an unknown filter id, one with no control, a species this run
// never screened, a guide sequence this file does not hold -- so the caller can show that refusal rather
// than silently drop it. `controls` says whether the fragment set any control at all, which is what tells
// the caller to open the panel those controls and their refusal banners live in.
function applyHashFragment(raw){
  const refused=[];
  if(raw.indexOf('=')<0) return {guide: raw ? decodeURIComponent(raw) : null, refused, controls:false};
  let guide=null, controls=false;
  for(const part of raw.split('&')){
    const eq=part.indexOf('='); if(eq<0) continue;
    const k=part.slice(0,eq), v=part.slice(eq+1);
    if(k==='g'){ guide=decodeURIComponent(v); }
    else if(k==='t'){
      controls=true;
      for(const pair of v.split(',')){
        if(!pair) continue;
        const ci=pair.indexOf(':'); if(ci<0) continue;
        // parseControlValue, the same boundary the controls use: a non-finite value, and an empty one
        // (`t=some_gate:`, which Number() reads as 0 and would silently set a real threshold), are
        // refused rather than applied.
        const id=pair.slice(0,ci), num=parseControlValue(pair.slice(ci+1));
        const f=FILTERS.find(x=>x.filter_id===id);
        if(f && evaluable(f) && typeof num==='number') T[id]=num;
        else refused.push(id);
      }
    } else if(k==='r'){
      controls=true;
      for(const pair of v.split(',')){
        if(!pair) continue;
        const ci=pair.indexOf(':'); if(ci<0) continue;
        const key=pair.slice(0,ci), num=parseControlValue(pair.slice(ci+1));
        if(READER_FILTERS.some(f=>f.k===key) && typeof num==='number') R[key]=num;
        else refused.push(key);
      }
    } else if(k==='sp'){
      // A species this run did not screen is refused, not adopted: scoping the column to a name with no
      // alignments behind it would read as "no liabilities here" for every guide in the file.
      controls=true;
      const sp=decodeURIComponent(v);
      if(SCREENED_SPECIES.includes(sp)) liabScope=sp; else refused.push(sp);
    } else if(k==='preset'){ controls=true; if(PRESETS[v]) activePreset=v; }
    else if(k==='q'){ controls=true; F.q=decodeURIComponent(v); }
    else if(k==='s'){
      // A set, so the key REPLACES it rather than adding to it -- and `s=` with nothing after it is a
      // reader who unchecked all four, which is a real view (every guide) and not a missing key.
      controls=true;
      F.status=new Set();
      for(const key of v.split(',')){
        if(!key) continue;
        if(STATUS_KEYS.includes(key)) F.status.add(key); else refused.push(key);
      }
    } else if(k==='c'){
      // Conservation is one question -- which species, and whether a mismatch outside the seed is
      // acceptable -- so both halves ride in one key and both are reset by it.
      controls=true;
      F.cons=new Set(); F.consSeedIntact=true;
      for(const token of v.split(',')){
        if(!token) continue;
        if(CONS_SPECIES.includes(token)) F.cons.add(token);
        else if(token===SEED_STRICT) F.consSeedIntact=false;
        else refused.push(token);
      }
    } else if(k==='k'){
      // Cart members as sequences, resolved against this file. A sequence this report does not hold is
      // refused rather than dropped: a URL from a differently built report names guides that are not
      // here, and a cart quietly missing them would read as a shortlist the reader never made.
      controls=true;
      for(const seq of v.split(',')){
        if(!seq) continue;
        const key=decodeURIComponent(seq);
        if(G.some(g=>g.guide===key)) CART.add(key); else refused.push(key);
      }
    }
  }
  return {guide, refused, controls};
}
function renderRefused(){
  const el=document.getElementById('refused');
  if(!refusedFilters.length){ el.style.display='none'; return; }
  el.style.display='';
  openFilters();          // a banner inside a closed <details> is a refusal that was swallowed
  el.innerHTML=`<b>Refused.</b> The URL asked for ${refusedFilters.map(esc).join(', ')} -- ${
    refusedFilters.length===1?'that name has':'those names have'} no control here and no match in this
    file (a frozen gate, an unscreened species, a guide this report does not hold, or no such name at
    all), so nothing was set from ${refusedFilters.length===1?'it':'them'}.`;
}

// Where the open guide sits relative to the index beside it. A guide arrives here by URL (`g=`) or stays
// here while a threshold moves under it, and the pane then showed a full guide -- every card, every
// metric -- while the index read "0 of 5,000". Nothing said the two were describing different sets, so
// the index looked broken and the guide looked like a result. The pane says which it is.
function placementNote(){
  const g=G.find(x=>x.guide===selected);
  if(!g) return '';
  const at=view.indexOf(g);
  if(at<0) return `<div class="warn"><b>Outside your current filter.</b> This guide is open because you
    asked for it — by link, or before the ${passesFilters(g)?'search box':'thresholds and filters'} above
    excluded it — and the index counts only what they admit. Nothing below is filtered: it is this
    guide's own evidence, at the thresholds in force.</div>`;
  if(at>=INDEX_DRAW_MAX) return `<div class="warn"><b>Not among the drawn rows.</b> This guide is row ${
    (at+1).toLocaleString()} of the ${view.length.toLocaleString()} your filter admits, and the index
    draws the first ${INDEX_DRAW_MAX.toLocaleString()} in this sort order, so it is not highlighted
    there. Everything below is this guide's own evidence.</div>`;
  return '';
}
function renderPlacement(){
  const el=document.getElementById('placement');
  if(el) el.innerHTML = placementNote();
}

function renderDetail(g){
  const i=g._i, live=g._live;
  document.getElementById('detail').innerHTML = `
   <div id="placement"></div>
   <div class="card"><h2>Guide</h2>
     <div class="mono big">${esc(g.guide)} <span class="pill ${({pass:'v-pass',warn:'v-warn',unknown:'v-unknown',fail:'v-fail'})[live.status]}" style="font-size:12px;vertical-align:middle">${live.status==='unknown'?'not established':live.status}</span></div>
     <div class="kv" style="margin-top:12px">
       <div><b>Passenger</b><span class="mono">${esc(g.passenger||'—')}</span></div>
       <div><b>Overhang</b>${esc(g.overhang||'—')}</div>
       <div><b>Modifications</b>${esc(g.modifications||'—')}</div>
       <div><b>Composite</b>${fmt(g.composite_score)}${g.weight_vector?` <span class="pill v-not_evaluated">${esc(g.weight_vector)}</span>`:''}</div>
       <div><b>Design score</b>${fmt(g.design_score)}</div>
       ${Object.entries(g.metrics).map(([k,v])=>`<div><b>${esc(k)}</b>${fmt(v)}</div>`).join('')}
     </div></div>
   ${card(gatesCard,g)}${card(structureCard,g)}${card(isoformCard,g)}${card(conservationCard,g)}${card(()=>offtargetCard(g,i),g)}${card(mirnaCard,g)}`;
  renderPlacement();
  if(g.offtarget_matrix.length) drawMatrix('mx'+i, g.offtarget_matrix);
}

function show(g){
  selected=g.guide; renderIndex();
  renderDetail(g);
  showOn(g);
  syncHash();
}

// Fired by every gate-threshold input. Deliberately does not rebuild the controls (that would steal
// focus mid-keystroke) -- only the derived state: live gates, the view, the cart, the map markers,
// the open detail pane and the fragment.
function onThresholdsChanged(){
  recomputeLive();
  renderStatusPills();
  applyFilters();
  renderCart();
  drawOverlay();
  if(selected){ const g=G.find(x=>x.guide===selected); if(g) renderDetail(g); }
  syncHash();
}

// One path for every way a fragment can arrive: the first load, Back, Forward, an edited address bar, a
// pasted link. The fragment was read exactly once, at load, with no `hashchange` listener anywhere -- so
// after Back the page showed one state while its URL claimed another. A fragment describes the WHOLE
// reader state, which is why resetControls runs first: what a URL omits is a default, not whatever the
// page happened to be showing. The cart is cleared here and not in resetControls for the same reason --
// `k=` is part of that state, so a URL without one is an empty cart, while the Reset button beside the
// thresholds must never destroy a pick list.
function applyFragmentState(raw){
  resetControls();
  CART.clear();
  const applied = applyHashFragment(raw);
  refusedFilters = applied.refused;
  recomputeLive();        // thresholds named in the fragment must be live before anything paints
  buildFilterUI();        // controls, checkboxes and the search box reflect T, R and F
  renderStatusPills();    // and so do the pills: a fragment can arrive with a threshold already moved
  applyFilters();
  renderCart();
  renderMap();
  // Whatever the fragment set is inside a container that ships closed, and so is any refusal of it.
  if(applied.controls || applied.refused.length) openFilters();
  const start = G.find(g=>g.guide===applied.guide) || G[0];
  if(start) show(start);
}

// One bootstrap block, guarded: the substituted <script> body is also handed directly to node for the
// parity harness (tests/unit/test_report_client_evaluator_parity.py), which has no `document` and no
// `location`. Every DOM-touching statement in the file lives inside a function or behind this guard,
// so the pure evaluator above (CMP, reevaluateGate, reevaluateGates, the presets, the register
// clustering) loads and runs there with nothing stubbed beyond the check itself.
function wireUI(){
  // Sorting, by mouse and by keyboard. `<th onclick>` with the default `tabIndex` of -1 is not a control:
  // it could not be focused, Enter and Space did nothing, and the sort was unreachable without a pointer.
  document.querySelectorAll('#idx th[data-k]').forEach(th=>{
    th.onclick=()=>sortBy(th.dataset.k);
    th.onkeydown=e=>{ if(e.key==='Enter'||e.key===' '){ e.preventDefault(); sortBy(th.dataset.k); } };
  });
  // Delegated once onto the tbody rather than assigned per row: the rows are rebuilt on every keystroke,
  // so 5,000 handler assignments were part of what a keystroke cost, and a row drawn after wiring used to
  // have to be re-wired to work at all. Enter/Space opens a row, the arrows walk the index, and the pick
  // cell is a control of its own -- it was a bare `<td>` with no role, invisible to a keyboard.
  const tb=document.querySelector('#idx tbody');
  const rowOf=e=>{ const tr=e.target.closest('tr'); return tr && tr.dataset.i!==undefined ? tr : null; };
  tb.onclick=e=>{
    const tr=rowOf(e); if(!tr) return;
    const cell=e.target.closest('td.pick');
    if(cell){ togglePick(cell.dataset.pick); return; }   // the pick column selects, it does not navigate
    show(G[+tr.dataset.i]);
  };
  tb.onkeydown=e=>{
    const tr=rowOf(e); if(!tr) return;
    if(e.key==='Enter'||e.key===' '){
      e.preventDefault();
      const cell=e.target.closest('td.pick');
      if(cell) togglePick(cell.dataset.pick); else show(G[+tr.dataset.i]);
      return;
    }
    if(e.key==='ArrowDown'||e.key==='ArrowUp'){
      const next = e.key==='ArrowDown' ? tr.nextElementSibling : tr.previousElementSibling;
      if(next && next.dataset.i!==undefined){ e.preventDefault(); next.focus(); }
    }
  };
  document.getElementById('q').oninput = e => { F.q = e.target.value; applyFilters(); syncHash(); };
  // A deliberate choice of isoform is pinned, so the map stops following the selection off it.
  document.getElementById('tx').onchange = e => { curTx = e.target.value; txPinned = true; renderMap(); };
  document.getElementById('txfollow').onclick = () => {
    txPinned = false;
    const g=G.find(x=>x.guide===selected); if(g) showOn(g); else renderMap();
  };
  // Read through `fnum`, the same boundary every threshold box uses, and refused at the control rather
  // than degraded: `parseInt('ten',10)||0` is 0, and adding the top 0 guides is indistinguishable from a
  // button that does not work. A count is also a whole number of guides, so a fraction is refused too.
  document.getElementById('addtop').onclick = () => {
    const box=document.getElementById('topn'), note=document.getElementById('topnnote');
    const n=fnum('topn'), ok = typeof n==='number' && Number.isInteger(n) && n>0;
    box.setAttribute('aria-invalid', ok?'false':'true');
    if(!ok){
      note.className='offscale';
      note.textContent = `That is not a whole number of guides above zero, so nothing was added to the cart.`;
      return;
    }
    const pool=matching().slice().sort((a,b)=>(b.composite_score??-1)-(a.composite_score??-1));
    const take=pool.slice(0,n);
    take.forEach(g=>CART.add(g.guide));
    note.className='empty';
    note.textContent = take.length < n
      ? `Added ${take.length.toLocaleString()} — only ${take.length.toLocaleString()} guides match the filters in force, of the ${n.toLocaleString()} asked for.`
      : `Added the top ${take.length.toLocaleString()} of ${pool.length.toLocaleString()} matching guides, by composite score.`;
    renderIndex(); renderCart(); drawOverlay(); syncHash();
  };
  document.getElementById('cartclear').onclick = () => {
    CART.clear(); renderIndex(); renderCart(); drawOverlay(); syncHash(); };
  document.getElementById('cartcopy').onclick = () => { const t=document.getElementById('carttsv'); t.focus(); t.select(); };
  // A Blob download, with the textarea as the declared fallback -- and the fallback CANNOT be reached by
  // a failure, which was the defect. The case the old catch named is an iframe without `allow-downloads`,
  // exactly the Quilt one: there the embedder refuses the synthetic `a.click()` itself, silently, with no
  // exception, no callback and no observable state, so the catch never ran and the reader was told the
  // file had been written. Nothing in this page can see that refusal, so it is not claimed: the box is
  // selected on every export, one keystroke from the data, and the note says both outcomes. The catch
  // stays for what does throw -- URL.createObjectURL and Blob are refused outright by some CSPs.
  document.getElementById('cartdl').onclick = () => {
    const note=document.getElementById('dlnote'), n=CART.size;
    const stamp=new Date().toISOString().slice(0,10).replace(/-/g,'');
    const name=`${GENE.replace(/[^A-Za-z0-9_.-]/g,'_')}_cart_${n}guides_${stamp}.tsv`;
    const box=document.getElementById('carttsv');
    box.focus(); box.select();
    let asked=false;
    try{
      const url=URL.createObjectURL(new Blob([box.value], {type:'text/tab-separated-values'}));
      const a=document.createElement('a');
      a.href=url; a.download=name; a.style.display='none';
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(()=>URL.revokeObjectURL(url), 0);
      asked=true;
    }catch(e){ asked=false; }
    note.textContent = asked
      ? ` ${name} requested. A viewer that blocks downloads refuses it without telling this page, so if no`
        + ` file arrived, the box below is already selected — copy it.`
      : ` download refused by this viewer — it could not build the file at all, which is a refusal this`
        + ` page CAN see. The box below is selected instead: copy it.`;
  };
  document.getElementById('freset').onclick = () => { resetControls(); buildFilterUI(); onThresholdsChanged(); };

  // The reader's own navigation, back through the report's own URL. Nothing here fires on the report's
  // own writes: replaceState does not raise hashchange, and the guard covers the location.hash fallback
  // syncHash uses in a sandbox that refuses the History API -- re-applying a fragment we just wrote would
  // fight the reader's typing.
  if(typeof window!=='undefined') window.addEventListener('hashchange', () => {
    if(location.hash.slice(1) === encodeHash()) return;
    applyFragmentState(location.hash.slice(1));
  });
  applyFragmentState(location.hash.slice(1));
}
if (typeof document !== 'undefined') { wireUI(); }
</script></body></html>
"""


def render_html(payload: ReportPayload) -> str:
    """Render the payload to one self-contained HTML document.

    Args:
        payload: Built by :func:`sirnaforge.reporting.payload.build_payload`.

    Returns:
        The whole document. No sidecar files, no external URLs.
    """
    # Every guide the payload holds, deliberately with no slice of its own. The cap lives in
    # build_payload, above every count the report prints (payload.DEFAULT_MAX_EMBEDDED_GUIDES); a second
    # one here would truncate `max_guides=None` again downstream and make the header over-count.
    guides = []
    for g in payload.guides:
        d = asdict(g)
        d["n_gates_failed"] = g.n_gates_failed
        d["n_gates_unknown"] = g.n_gates_unknown
        d["status"] = g.status
        guides.append(d)

    dropped = payload.run.get("guides_dropped_by_status") or {}
    env = Environment(autoescape=select_autoescape(default=True), trim_blocks=True, lstrip_blocks=True)
    # Rendered before the payload is substituted, so no JSON string can be parsed as template syntax.
    html = env.from_string(_TEMPLATE).render(
        p=payload,
        scope=_liability_scope(payload),
        # Dropped guides a reader might have shipped: the cap's real cost, as opposed to guides the run
        # and this report agree to reject.
        dropped_shippable=sum(n for status, n in dropped.items() if status != "fail"),
    )
    maps = _design_maps(payload)
    for placeholder, value in (
        ("GUIDES_JSON_PLACEHOLDER", guides),
        ("FILTERS_JSON_PLACEHOLDER", payload.filters),
        ("MAPS_JSON_PLACEHOLDER", maps),
        ("MAP_LEGEND_JSON_PLACEHOLDER", _MAP_LEGEND),
        # The two legend keys the client branches on: one is not a verdict (no guide to re-decide), the
        # other is not a point at all. Substituted rather than spelled in the JS so there is one name.
        ("NOT_EMBEDDED_JSON_PLACEHOLDER", _NOT_EMBEDDED),
        ("GAPS_KEY_JSON_PLACEHOLDER", _GAPS_KEY),
        ("LAYOUTS_JSON_PLACEHOLDER", payload.run.get("structure_layouts") or {}),
        ("TX_IDS_JSON_PLACEHOLDER", payload.run.get("transcript_ids") or []),
        # The species split the index column and its bound read (#101's decomposition), and which of
        # them the embedded per-row detail is for. Both as data: a species name hardcoded in the JS
        # would outlive the screen that produced it.
        ("SCREENED_SPECIES_JSON_PLACEHOLDER", payload.run.get("screened_species") or []),
        ("EMBED_SPECIES_JSON_PLACEHOLDER", EMBED_SPECIES),
        ("EMBED_SCOPE_JSON_PLACEHOLDER", str(payload.run.get("embed_scope") or "")),
        ("GENE_JSON_PLACEHOLDER", str(payload.run.get("gene_query") or "run")),
        ("REGISTER_NT_PLACEHOLDER", REGISTER_NEIGHBOUR_NT),
    ):
        html = html.replace(placeholder, _embed(value))
    return html


#: Point classes on the design map, drawn in this order so a passing window is never hidden under a
#: rejected one. The keys are the report's own verdicts, so the labels are the words the status pills
#: use rather than a second vocabulary: the dots are re-coloured from each guide's LIVE verdict, and a
#: legend reading `run PASS, warn threshold exceeded` beside a dot the reader's own threshold decided
#: would be describing a comparison that no longer happened.
_MAP_SERIES = (
    ("fail", "fail"),
    ("unknown", "not established"),
    ("warn", "pass, warned"),
    ("pass", "pass, every gate evaluated"),
)

#: Windows whose guide is not in this file, so no live verdict exists for them. They keep the colour the
#: run's own verdict gave them and are drawn underneath everything else, because they are the one class
#: on the map a moved threshold cannot answer for. Only a capped report has any (payload's
#: ``_transcript_maps`` plots every candidate row either way, deliberately).
_NOT_EMBEDDED = "not_embedded"

#: The gap key. Outlined rather than filled: ``transcript_map_svg`` shades an uncovered stretch almost
#: white, which is invisible as a 9 px swatch, and the fill itself is that primitive's private business.
_GAPS_KEY = "gaps"

#: The legend, as data, because the client paints it from the counts it just drew. Pre-rendering it in
#: Python is what froze it at the run's thresholds while every other number on the page went live: on
#: one internal run, moving one gate to 60 left 39 guides passing and the legend still read
#: `PASS, all gates evaluated (394)`, byte-identical, over 394 green dots (#103's own failure mode, in
#: the largest visual on the page). Order is draw order, so a passing dot is never hidden underneath a
#: rejected one and the windows no live verdict can answer for sit at the bottom.
_MAP_LEGEND = [
    {"key": _NOT_EMBEDDED, "label": "not in this file — the run's own verdict", "fill": SERIES_FILL["reference"]},
    *({"key": key, "label": label, "fill": SERIES_FILL[key]} for key, label in _MAP_SERIES),
    {"key": _GAPS_KEY, "label": "no candidate in this table", "fill": None},
]


def _guide_index_by_position(payload: ReportPayload) -> dict[str, dict[int, int]]:
    """``transcript -> position -> index into the embedded guide list``.

    The link the map needs to go live: the payload's point series carry ``(position, value)`` and the
    class the run's thresholds gave them, but not which guide each dot is. A position on a transcript is
    one enumerated window, so the guide's own isoform rows -- which carry both -- resolve it exactly.
    """
    out: dict[str, dict[int, int]] = {}
    for i, guide in enumerate(payload.guides):
        for row in guide.isoforms:
            position = row.get("position")
            if position is None:
                continue
            out.setdefault(str(row["transcript"]), {})[int(position)] = i
    return out


def _picker_order(transcripts: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Transcripts in the order the isoform picker should offer them: canonical first.

    A stable sort on one key only. The payload orders these by window count, which is very nearly
    length, and that opened one internal report on a 7,988 nt transcript with the canonical one tenth --
    so the map a reader sees first was of an isoform nobody designs against. ``canonical`` is ``None``
    when the run recorded canonical status nowhere: unknown sorts between known-canonical and
    known-not, because promoting an unknown over a stated non-canonical would be inventing the fact.
    """
    rank = {True: 0, None: 1, False: 2}
    return sorted(transcripts, key=lambda t: rank[t.get("canonical")])


def _design_maps(payload: ReportPayload) -> dict[str, dict[str, Any]]:
    """Pre-render one transcript backdrop per transcript, and hand the client the dots.

    The backdrop -- regions, gridlines, axis, not-enumerated shading -- is drawn server-side by the
    shared track primitive, so a notebook figure and this card still come out of one function. The
    **points** are handed over as numbers instead of pixels, because their colour is a verdict: it has
    to follow the reader's thresholds like every other verdict on the page, and a pre-rendered dot
    cannot. The client re-draws them into the marker overlay in the geometry returned here, so the two
    share one coordinate system rather than two copies of a scale.
    """
    out: dict[str, dict[str, Any]] = {}
    guide_at = _guide_index_by_position(payload)
    for entry in _picker_order(payload.run.get("transcripts") or []):
        regions = TranscriptRegions(
            transcript_id=entry["transcript_id"],
            length=entry["length"],
            cds_start=entry.get("cds_start"),
            cds_end=entry.get("cds_end"),
        )
        series = [
            PointSeries(label, key, [(int(p), float(v)) for p, v in entry["series"].get(key, [])])
            for key, label in _MAP_SERIES
        ]
        gaps = [(int(a), int(b)) for a, b in entry.get("gaps") or []]
        title = f"{regions.transcript_id} - {regions.length:,} nt, {entry['windows']:,} candidate windows"
        value_label = entry.get("value_column") or "composite score"
        # Drawn twice on purpose: the first call is the only public way to learn the value bounds the
        # primitive rounds its own gridlines to, and the second draws the same backdrop with those
        # bounds and no points on it. Re-deriving the rounding here instead would put a second copy of
        # the y-scale in this file, which is exactly how a dot and its gridline come to disagree.
        # Both calls pass the same title and value_label: they shift the plot area down, so a geometry
        # measured without them would put every client-drawn dot 24 px above the axes it belongs to.
        _, bounds = transcript_map_svg(
            regions, series, title=title, value_label=value_label, gaps=gaps, standalone=False
        )
        svg, geometry = transcript_map_svg(
            regions,
            [],
            title=title,
            value_label=value_label,
            value_range=(bounds.value_min, bounds.value_max),
            gaps=gaps,
            standalone=False,
        )
        here = guide_at.get(regions.transcript_id, {})
        points = [
            [int(p), round(float(v), 3), here.get(int(p), -1)]
            for key, _label in _MAP_SERIES
            for p, v in entry["series"].get(key, [])
        ]
        out[regions.transcript_id] = {
            "label": f"{regions.transcript_id} - {regions.length:,} nt, {entry['windows']:,} candidates",
            "canonical": entry.get("canonical"),
            "svg": svg,
            # ``[position, value, guide index or -1]``. The run's own class is deliberately absent: it
            # would be a second, frozen verdict for the same dot, and only the guide's live one paints.
            "points": points,
            "geometry": geometry.as_dict(),
            "viewBox": _view_box(svg),
            "gaps": len(gaps),
            "gap_note": (
                f"{len(gaps)} stretch(es) of {MIN_UNCOVERED_NT} nt or more carry no candidate in this table."
                if gaps
                else ""
            ),
        }
    return out


def _view_box(svg: str) -> str:
    """The base map's viewBox, so the marker overlay shares its coordinate system exactly."""
    start = svg.index('viewBox="') + len('viewBox="')
    return svg[start : svg.index('"', start)]


def write_report(payload: ReportPayload, out_path: Path | str) -> Path:
    """Render and write the report, returning the path written."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_html(payload), encoding="utf-8")
    return out
