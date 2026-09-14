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
from dataclasses import asdict
from pathlib import Path
from typing import Any

from jinja2 import Environment, select_autoescape

from sirnaforge.reporting.payload import MIN_UNCOVERED_NT, REGISTER_NEIGHBOUR_NT, ReportPayload
from sirnaforge.reporting.tracks import PointSeries, TranscriptRegions, legend_html, transcript_map_svg

#: Guides embedded in the index. The per-guide detail is rendered for all of them; this bounds only
#: how many rows the client holds, and the number is reported in the header so it is never a silent cap.
MAX_INDEX_GUIDES = 5000


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
body{margin:0;font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;background:var(--bg);color:var(--fg)}
header{padding:18px 24px;border-bottom:1px solid var(--line);background:var(--card)}
h1{margin:0 0 4px;font-size:17px;font-weight:650}
.sub{color:var(--mut);font-size:12.5px}
main{display:grid;grid-template-columns:minmax(430px,40%) 1fr;gap:0;height:calc(100vh - 86px)}
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
.fstat label{margin-right:9px;font-size:12px;white-space:nowrap}
.cartbtn{cursor:pointer;border:1px solid var(--line);background:var(--card);border-radius:5px;
padding:3px 9px;font:inherit;font-size:12px;color:var(--fg)}
.cartbtn:hover{background:#f3f4f6}
.picked{color:var(--acc);font-weight:700}
#cartlist{max-height:210px;overflow:auto;font-size:12.5px}
/* Tab-separated columns only line up if the text is not wrapped. */
#carttsv{width:100%;height:96px;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:11px;
border:1px solid var(--line);border-radius:5px;padding:7px;resize:vertical;white-space:pre;overflow:auto}
#left{border-right:1px solid var(--line);overflow:auto;background:var(--card)}
#right{overflow:auto;padding:20px 24px}
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
</style></head><body>
<header>
  <h1>siRNAforge — {{ p.run.gene_query }} · {{ '{:,}'.format(p.run.guides) }} guides</h1>
  <div class="sub">
    {{ '{:,}'.format(p.run.candidate_rows) }} candidate rows collapsed onto guide sequence ·
    {{ '{:,}'.format(p.run.hit_rows) }} alignments
    (<span class="liab">{{ '{:,}'.format(p.run.liability_rows) }} liabilities</span>,
     <span class="nonliab">{{ '{:,}'.format(p.run.non_liability_rows) }} not liabilities</span>) ·
    liability rows embedded for {{ p.run.embed_scope }} ·
    payload schema {{ p.schema_version }} · profile {{ p.provenance.policy_profile }} ·
    run mode {{ p.provenance.run_mode }} · gates from {{ p.provenance.policy_source }} ·
    sirnaforge {{ p.provenance.tool_version }}
  </div>
  <div class="sub" style="margin-top:6px">
    <span class="pill v-pass">{{ p.run.status_counts["pass"] }} pass</span>
    <span class="pill v-warn">{{ p.run.status_counts["warn"] }} pass, warned</span>
    <span class="pill v-unknown">{{ p.run.status_counts["unknown"] }} not established</span>
    <span class="pill v-fail">{{ p.run.status_counts["fail"] }} fail</span>
    {% if p.run.agreement.comparable %}
    &nbsp;· against the run's own verdicts:
    <b>{{ p.run.agreement.contradicted_run_pass }}</b>/{{ p.run.agreement.run_pass_guides }} run PASS
    contradicted, <b>{{ p.run.agreement.overruled_run_fail }}</b> run rejections overruled,
    <b>{{ p.run.agreement.run_failed_not_rederivable }}</b> run rejections not re-derivable here.
    {% endif %}
  </div>
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
    <div class="fsec">Gate thresholds <i>— moving one re-derives the verdict this run computed</i></div>
    <div id="frows"></div>
    <div class="fsec">Reader filters <i>— these only select among rows; no verdict changes</i></div>
    <div id="rrows"></div>
    <div class="fstat" id="fcons" style="margin-top:6px;padding-top:6px;border-top:1px solid var(--line)">
    </div>
    <div style="margin-top:8px">
      <button class="cartbtn" id="addtop">Add top <input id="topn" value="10" style="width:42px"> to cart</button>
      <button class="cartbtn" id="freset">Reset</button>
    </div>
  </details>
  <table id="idx"><thead><tr>
    <th class="pick" title="in cart"></th>
    <th data-k="status">Status</th><th data-k="guide">Guide</th><th data-k="composite">Score</th>
    <th data-k="failed" title="gates failed / gates not evaluated">Gates</th>
    <th data-k="liab" title="off-target liabilities">Liab.</th>
    <th data-k="rows" title="distinct isoforms carrying this guide, of all in the run">Isoforms</th>
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
     Quilt iframe without <code>allow-downloads</code> -- block that; the box above is then the way
     out, and says so if the download is refused.</p>
   </div>
   <div class="card" id="mapcard">
     <h2>Candidate positions</h2>
     <div style="font-size:12px;color:#6b7280;margin:-4px 0 10px">Every candidate the run enumerated
     on one isoform. <label for="tx">Isoform</label>
     <select id="tx" aria-label="isoform"></select>
     <span id="txnote"></span></div>
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
const LAYOUTS = LAYOUTS_JSON_PLACEHOLDER;
const TX_IDS = TX_IDS_JSON_PLACEHOLDER;
const GENE = GENE_JSON_PLACEHOLDER;
const REGISTER_NT = REGISTER_NT_PLACEHOLDER;
const fmt = n => n===null||n===undefined ? '—' : (typeof n==='number' ? (Number.isInteger(n)?n.toLocaleString():n.toFixed(3)) : n);
const esc = s => String(s??'').replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
let view = G.slice(), sortK='composite', sortAsc=false, selected=null;

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
  if(k==='liab')return g.liability_count; if(k==='rows')return g.n_rows; return 0;
}
function renderIndex(){
  const tb=document.querySelector('#idx tbody');
  view.sort((a,b)=>{const x=rowKey(a,sortK),y=rowKey(b,sortK);
    const c = typeof x==='string' ? x.localeCompare(y) : x-y; return sortAsc?c:-c;});
  const V={pass:'v-pass',warn:'v-warn',unknown:'v-unknown',fail:'v-fail'};
  tb.innerHTML = view.map(g=>`<tr data-i="${G.indexOf(g)}" class="${g.guide===selected?'sel':''}">
    <td class="pick" data-pick="${esc(g.guide)}" title="add to / remove from cart"
      >${CART.has(g.guide)?'<span class="picked">\u2713</span>':'<span style="color:#d1d5db">+</span>'}</td>
    <td><span class="pill ${V[g._live.status]}">${g._live.status}</span></td>
    <td class="mono">${esc(g.guide)}</td><td>${fmt(g.composite_score)}</td>
    <td>${g._live.n_gates_failed} / ${g._live.n_gates_unknown}</td>
    <td class="${g.liability_count?'liab':'nonliab'}">${g.liability_count}</td>
    <td title="${g.n_rows} enumeration${g.n_rows===1?'':'s'}">${isoformFrac(g)}</td></tr>`).join('');
  tb.querySelectorAll('tr').forEach(tr=>tr.onclick=e=>{
    const cell=e.target.closest('td.pick');
    if(cell){ togglePick(cell.dataset.pick); return; }   // the pick column selects, it does not navigate
    show(G[+tr.dataset.i]);
  });
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
//: payload.py's reason codes (payload.py:187-196), restated so the client evaluator can read them --
//: only REASON_OK is ever re-thresholded; every other reason means the run itself never decided.
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
  const nu=live.n_gates_unknown;
  // Whose thresholds these verdicts came from, said once at the top of the panel: a reader quoting a
  // pill has to be able to see that it is theirs and not the run's without reading every row.
  const nm=FILTERS.filter(f=>liveThreshold(f)!==f.threshold).length;
  return `<div class="card"><h2>Gates — all ${live.gates.length}, independently evaluated</h2>
    ${nm?`<div class="warn"><b>${nm} threshold${nm===1?'':'s'} moved from the run's.</b> Every verdict
      below is re-derived at your value; the run's own threshold is shown beside it and named again in
      the Why column, so nothing here is the run's verdict unless it says so.</div>`:''}
    ${nu?`<div class="warn"><b>${nu} of ${live.gates.length} gates not evaluated.</b> An unevaluated gate is
      not a pass. See the Why column for the missing input.</div>`:''}
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
  const note = g.offtarget_embedded_scope_empty_but_counts_exist
    ? `<div class="warn"><b>Nothing in the embedded scope, but this guide is not clean.</b>
       Its alignments all lie outside ${esc('the embedded scope')} and are shown as counts only.
       A guide with hits at nm&nbsp;&ge;&nbsp;3 must never render like a genuinely clean one.</div>` : '';
  return `<div class="card"><h2>Off-targets by gene — ${g.liability_count} liabilit${g.liability_count===1?'y':'ies'}
      of ${g.offtarget_by_symbol.reduce((a,s)=>a+s.n,0)} alignments</h2>
    ${note}<div id="mx${idx}" style="height:230px"></div>
    <table><thead><tr><th>Gene</th><th>Class</th><th>Species</th><th>Alignments</th><th>Best nm</th></tr></thead>
    <tbody>${rows}</tbody></table>
    ${g.offtarget_rows.length?`<h2 style="margin-top:18px">Liability alignments — ${g.offtarget_rows.length} in the embedded scope</h2>
      <table><thead><tr><th>Transcript</th><th>Gene</th><th>Class</th><th>nm</th><th>Seed mm</th><th>CIGAR</th><th>Strand</th></tr></thead>
      <tbody>${g.offtarget_rows.map(r=>`<tr><td class="mono">${esc(r.t)}</td><td class="mono">${esc(r.s)}</td>
        <td><span class="pill v-fail">${esc(r.c)}</span></td><td>${fmt(r.nm)}</td><td>${fmt(r.sm)}</td>
        <td class="mono">${esc(r.cig)}</td><td>${esc(r.st)}</td></tr>`).join('')}</tbody></table>`
      :`<p class="empty" style="margin-top:14px">No liability alignment in the embedded scope. Non-liability
        classes are grouped above and counted completely in the chart.</p>`}</div>`;
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
// ---- design map: one pre-rendered transcript at a time, plus a live marker for the selection ----
// The map itself is drawn in Python by reporting.tracks, so a notebook figure and this card cannot
// diverge. Only the selection marker is client-side, because only it depends on what is selected.
const TXS = Object.keys(MAPS);
let curTx = TXS[0] || null;

function renderMap(){
  const sel=document.getElementById('tx');
  if(!TXS.length){ document.getElementById('mapcard').style.display='none'; return; }
  const carries=new Set((G.find(x=>x.guide===selected)||{isoforms:[]}).isoforms.map(i=>i.transcript));
  sel.innerHTML = TXS.map(t=>`<option value="${esc(t)}"${t===curTx?' selected':''}>${
    esc(MAPS[t].label)}${carries.has(t)?' \u25cf carries this guide':''}</option>`).join('');
  document.getElementById('txnote').textContent = selected
    ? ` \u2014 ${carries.size} of ${TXS.length} shown isoforms carry the selected guide.` : '';
  const m=MAPS[curTx];
  document.getElementById('mapbase').innerHTML=m.svg;
  document.getElementById('maplegend').innerHTML=m.legend;
  document.getElementById('mapover').setAttribute('viewBox', m.viewBox);
  document.getElementById('mapnote').textContent=m.note;
  drawOverlay();
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
  const geo=MAPS[curTx].geometry, ps=selectedPositions();
  const x=p=>geo.x0+(geo.x1-geo.x0)*(Math.min(Math.max(p,1)-1,geo.length-1))/Math.max(geo.length-1,1);
  // Cart members first, so the selection marker draws over them rather than under.
  over.innerHTML = cartPositions().map(p=>
      `<path d="M${x(p).toFixed(1)} ${geo.y1}l-3.5 7h7z" fill="#15803d"/>`).join('')
    + ps.map(p=>{const px=x(p).toFixed(1);
    return `<line x1="${px}" y1="${geo.y0}" x2="${px}" y2="${geo.y1}" stroke="#1d4ed8" stroke-width="1.4" stroke-dasharray="3 2"/>`
         + `<path d="M${px} ${geo.y0-2}l-4.5-7h9z" fill="#1d4ed8"/>`;}).join('');
  const note=document.getElementById('mapnote');
  if(selected) note.textContent = ps.length
    ? `${esc(selected)} is enumerated at ${ps.map(p=>p.toLocaleString()).join(', ')} on this transcript.`
    : `${esc(selected)} is not enumerated on this transcript.`;
  else note.textContent = MAPS[curTx].note;
  const nc=cartPositions().length;
  if(nc) note.textContent += ` ${nc} cart position${nc===1?'':'s'} marked below the axis.`;
}

function showOn(g){
  // Follow the selection to an isoform that actually carries it, so the marker is never off-screen,
  // then re-label the picker: which isoforms carry the guide is part of the selection, not the map.
  if(!curTx) return;
  const here=g.isoforms.some(i=>i.transcript===curTx);
  if(!here){ const first=g.isoforms.find(i=>MAPS[i.transcript]); if(first) curTx=first.transcript; }
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
const F = {status:new Set(['pass','warn']), cons:new Set(), consSeedIntact:true};

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
  //: counters over them, never this count, so _live carries no copy of it to prefer instead.
  {k:'liab',      label:'liabilities',     dir:'max', get:g=>g.liability_count},
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

// A refused keystroke has to be visible at the control. Silently keeping the last good value would
// leave the reader looking at a threshold they did not type and believing they had moved it.
function markControl(el, ok){
  el.setAttribute('aria-invalid', ok?'false':'true');
  const note=document.getElementById('badnum'); if(!note) return;
  const bad=[...document.querySelectorAll('#frows input[aria-invalid="true"], #rrows input[aria-invalid="true"]')];
  note.style.display = bad.length ? '' : 'none';
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

// One generic control per evaluable filter -- no filter_id and no filter column is ever a branch in
// this function, only data read off FILTERS. A frozen filter gets a read-only row, no input, and
// payload.py's own `unevaluable_reason` (not a re-derived guess at why): the evaluator above already
// refuses to move it, and this is where that refusal is visible rather than a control nobody could use.
function buildGateControls(){
  document.getElementById('frows').innerHTML = FILTERS.map(f=>{
    if(!evaluable(f)) return `<div class="frow" data-filter="${esc(f.filter_id)}">
      <span class="mono">${esc(f.filter_id)}</span>
      <span class="empty" style="grid-column:span 2">frozen — ${esc(f.unevaluable_reason||'not re-thresholdable')}</span></div>`;
    const c = f.control || {};
    return `<div class="frow" data-filter="${esc(f.filter_id)}">
      <span title="${esc(f.definition||'')}"><span class="mono">${esc(f.filter_id)}</span> ${esc(f.comparator)}</span>
      <input id="ctl_${esc(f.filter_id)}" type="number" step="any" value="${T[f.filter_id]}"
        aria-label="${esc(f.filter_id)} threshold">
      <span class="empty">run ${fmt(f.threshold)}</span></div>
      <input type="range" id="rng_${esc(f.filter_id)}" min="${c.min}" max="${c.max}" step="${c.step}"
        value="${T[f.filter_id]}" style="grid-column:1/-1;width:100%" aria-label="${esc(f.filter_id)} slider">`;
  }).join('');
  // The number box and the slider are the same threshold; typing in one moves the other, so "yours"
  // never has two disagreeing readouts. The box carries `type=number` and `step=any` but deliberately
  // no min/max: the slider's domain bounds the slider, and payload.py's own contract is that the number
  // box can still reach any threshold (payload._control_domain). What is rejected is not an unusual
  // number, it is a value that is not a number at all.
  document.querySelectorAll('#frows input[id^="ctl_"], #frows input[id^="rng_"]').forEach(el=>el.oninput=()=>{
    const id=el.id.slice(el.id.indexOf('_')+1), f=FILTERS.find(x=>x.filter_id===id), v=fnum(el.id);
    if(v===undefined){ markControl(el,false); return; }   // refused here, so T can never hold a NaN
    markControl(el,true);
    T[id] = v===null ? f.threshold : v;
    const other=document.getElementById((el.id.startsWith('ctl_')?'rng_':'ctl_')+id);
    if(other) other.value = T[id];
    onThresholdsChanged();
  });
}

// One box per reader filter, reusing the gate rows' own layout so the two groups read as siblings --
// but with no slider and no "run <value>", because there is no run threshold over these to be seeded
// from or reset to. Moving one deliberately does NOT call onThresholdsChanged: nothing is re-derived,
// so recomputing every guide's gate table would be a claim that something changed when nothing did.
function buildReaderFilters(){
  document.getElementById('rrows').innerHTML = READER_FILTERS.map(f=>
    `<div class="frow" data-reader="${esc(f.k)}"><span>${esc(f.label)}</span>
      <input id="rf_${esc(f.k)}" type="number" step="any" placeholder="${f.dir}"
        value="${R[f.k]===null||R[f.k]===undefined?'':R[f.k]}"
        aria-label="${esc(f.label)}, ${f.dir==='min'?'at least':'at most'}">
      <span class="empty">${f.dir==='min'?'at least':'at most'}</span></div>`).join('');
  document.querySelectorAll('#rrows input').forEach(el=>el.oninput=()=>{
    const v=fnum(el.id);
    if(v===undefined){ markControl(el,false); return; }   // a bound is a number or it is not a bound
    markControl(el,true);
    R[el.id.slice(3)] = v;                 // blank clears the bound; it does not set it to 0
    applyFilters(); syncHash();
  });
}

function buildFilterUI(){
  document.getElementById('fstat').innerHTML = STATUS_KEYS.map(k=>
    `<label><input type="checkbox" data-status="${k}"${F.status.has(k)?' checked':''}> ${k}</label>`).join('');
  document.querySelectorAll('#fstat input').forEach(cb=>cb.onchange=()=>{
    cb.checked ? F.status.add(cb.dataset.status) : F.status.delete(cb.dataset.status); applyFilters(); });
  buildGateControls();
  buildReaderFilters();
  document.getElementById('fcons').innerHTML =
    '<span style="color:var(--mut)">conserved in</span> ' +
    CONS_SPECIES.map(sp=>`<label><input type="checkbox" data-cons="${sp}"${F.cons.has(sp)?' checked':''}> ${sp}</label>`).join('') +
    `<label title="accept one mismatch outside guide positions 2-8">
      <input type="checkbox" id="consseed"${F.consSeedIntact?' checked':''}> allow 1 mm outside the seed</label>`;
  document.querySelectorAll('#fcons input[data-cons]').forEach(cb=>cb.onchange=()=>{
    cb.checked ? F.cons.add(cb.dataset.cons) : F.cons.delete(cb.dataset.cons); applyFilters(); });
  document.getElementById('consseed').onchange = e => { F.consSeedIntact=e.target.checked; applyFilters(); };
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
function resetControls(){
  F.status = new Set(['pass','warn']); F.cons = new Set(); F.consSeedIntact = true;
  FILTERS.forEach(f => { if(evaluable(f)) T[f.filter_id] = f.threshold; });
  READER_FILTERS.forEach(f => { R[f.k] = null; });
  activePreset = 'all'; refusedFilters = [];
}

function matching(){ return G.filter(passesFilters); }

function applyFilters(){
  const raw=(document.getElementById('q').value||'').trim();
  const seq=raw.toUpperCase().replace(/T/g,'U'), id=raw.toLowerCase();
  view = matching().filter(g => !raw || g.guide.includes(seq)
    || g.isoforms.some(i=>String(i.candidate_id).toLowerCase().includes(id)
                        || String(i.transcript).toLowerCase().includes(id)));
  document.getElementById('fcount').textContent = `${view.length.toLocaleString()} of ${G.length.toLocaleString()} guides`;
  renderIndex();
}

function togglePick(guide){
  CART.has(guide) ? CART.delete(guide) : CART.add(guide);
  renderIndex(); renderCart(); drawOverlay();
}

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
}

// ---- URL fragment: selected guide, moved thresholds, reader filters, active preset ----------------
// A bare fragment with none of the '=' syntax below is still read as a guide -- what every earlier
// report in this run wrote. Reader filters ride in their own `r=` list, on the same footing as `t=`:
// a reader who sends the URL sends the rows they were looking at, not just the verdicts.
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
  if(activePreset!=='all') parts.push('preset='+activePreset);
  return parts.join('&');
}
function syncHash(){ location.hash = encodeHash(); }

// Reads a fragment written by encodeHash (or a bare guide sequence from an earlier report) and
// returns the filter ids it had to refuse -- an unknown id, or one with no control -- so the caller
// can show that refusal rather than silently drop it.
function applyHashFragment(raw){
  const refused=[];
  if(raw.indexOf('=')<0) return {guide: raw ? decodeURIComponent(raw) : null, refused};
  let guide=null;
  for(const part of raw.split('&')){
    const eq=part.indexOf('='); if(eq<0) continue;
    const k=part.slice(0,eq), v=part.slice(eq+1);
    if(k==='g') guide=decodeURIComponent(v);
    else if(k==='t'){
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
      for(const pair of v.split(',')){
        if(!pair) continue;
        const ci=pair.indexOf(':'); if(ci<0) continue;
        const key=pair.slice(0,ci), num=parseControlValue(pair.slice(ci+1));
        if(READER_FILTERS.some(f=>f.k===key) && typeof num==='number') R[key]=num;
        else refused.push(key);
      }
    } else if(k==='preset'){ if(PRESETS[v]) activePreset=v; }
  }
  return {guide, refused};
}
function renderRefused(){
  const el=document.getElementById('refused');
  if(!refusedFilters.length){ el.style.display='none'; return; }
  el.style.display='';
  el.innerHTML=`<b>Refused.</b> The URL asked to set ${refusedFilters.map(esc).join(', ')} -- ${
    refusedFilters.length===1?'that name has':'those names have'} no control here (a frozen gate, or no
    such control at all), so nothing moved.`;
}

function renderDetail(g){
  const i=G.indexOf(g), live=g._live;
  document.getElementById('detail').innerHTML = `
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
  applyFilters();
  renderCart();
  drawOverlay();
  if(selected){ const g=G.find(x=>x.guide===selected); if(g) renderDetail(g); }
  syncHash();
}

// One bootstrap block, guarded: the substituted <script> body is also handed directly to node for the
// parity harness (tests/unit/test_report_client_evaluator_parity.py), which has no `document` and no
// `location`. Every DOM-touching statement in the file lives inside a function or behind this guard,
// so the pure evaluator above (CMP, reevaluateGate, reevaluateGates, the presets, the register
// clustering) loads and runs there with nothing stubbed beyond the check itself.
function wireUI(){
  document.querySelectorAll('#idx th').forEach(th=>th.onclick=()=>{
    const k=th.dataset.k; if(k===sortK) sortAsc=!sortAsc; else {sortK=k;sortAsc=(k==='guide');} renderIndex();});
  document.getElementById('q').oninput = () => applyFilters();
  document.getElementById('tx').onchange = e => { curTx = e.target.value; renderMap(); };
  document.getElementById('addtop').onclick = e => {
    if(e.target.id==='topn') return;                     // typing in the field is not a click on the button
    const n=Math.max(0, parseInt(document.getElementById('topn').value,10)||0);
    matching().slice().sort((a,b)=>(b.composite_score??-1)-(a.composite_score??-1))
      .slice(0,n).forEach(g=>CART.add(g.guide));
    renderIndex(); renderCart(); drawOverlay();
  };
  document.getElementById('cartclear').onclick = () => { CART.clear(); renderIndex(); renderCart(); drawOverlay(); };
  document.getElementById('cartcopy').onclick = () => { const t=document.getElementById('carttsv'); t.focus(); t.select(); };
  // A Blob download, because that is what "export" means, with the textarea as the declared fallback:
  // the report's own sandbox may withhold allow-downloads, and a button that silently does nothing is
  // worse than one that says why.
  document.getElementById('cartdl').onclick = () => {
    const note=document.getElementById('dlnote'), n=CART.size;
    const stamp=new Date().toISOString().slice(0,10).replace(/-/g,'');
    const name=`${GENE.replace(/[^A-Za-z0-9_.-]/g,'_')}_cart_${n}guides_${stamp}.tsv`;
    try{
      const url=URL.createObjectURL(new Blob([document.getElementById('carttsv').value],
        {type:'text/tab-separated-values'}));
      const a=document.createElement('a');
      a.href=url; a.download=name; a.style.display='none';
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(()=>URL.revokeObjectURL(url), 0);
      note.textContent=` ${name}`;
    }catch(e){
      note.textContent=' download refused by this viewer — copy the box below instead';
      const t=document.getElementById('carttsv'); t.focus(); t.select();
    }
  };
  document.getElementById('freset').onclick = () => { resetControls(); buildFilterUI(); onThresholdsChanged(); };

  const {guide, refused} = applyHashFragment(location.hash.slice(1));
  refusedFilters = refused;
  recomputeLive();       // thresholds named in the fragment must be live before anything paints
  buildFilterUI();        // gate controls reflect T, including any fragment overrides
  applyFilters();
  renderCart();
  renderMap();
  const start = G.find(g=>g.guide===guide) || G[0];
  if(start) show(start);
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
    guides = []
    for g in payload.guides[:MAX_INDEX_GUIDES]:
        d = asdict(g)
        d["n_gates_failed"] = g.n_gates_failed
        d["n_gates_unknown"] = g.n_gates_unknown
        d["status"] = g.status
        guides.append(d)

    env = Environment(autoescape=select_autoescape(default=True), trim_blocks=True, lstrip_blocks=True)
    # Rendered before the payload is substituted, so no JSON string can be parsed as template syntax.
    html = env.from_string(_TEMPLATE).render(p=payload)
    for placeholder, value in (
        ("GUIDES_JSON_PLACEHOLDER", guides),
        ("FILTERS_JSON_PLACEHOLDER", payload.filters),
        ("MAPS_JSON_PLACEHOLDER", _design_maps(payload)),
        ("LAYOUTS_JSON_PLACEHOLDER", payload.run.get("structure_layouts") or {}),
        ("TX_IDS_JSON_PLACEHOLDER", payload.run.get("transcript_ids") or []),
        ("GENE_JSON_PLACEHOLDER", str(payload.run.get("gene_query") or "run")),
        ("REGISTER_NT_PLACEHOLDER", REGISTER_NEIGHBOUR_NT),
    ):
        html = html.replace(placeholder, _embed(value))
    return html


#: Point classes on the design map, drawn in this order so a passing window is never hidden under a
#: rejected one. Labels are the legend text.
_MAP_SERIES = (
    ("fail", "rejected"),
    ("unknown", "run PASS, gate evidence incomplete"),
    ("warn", "run PASS, warn threshold exceeded"),
    ("pass", "PASS, all gates evaluated"),
)


def _design_maps(payload: ReportPayload) -> dict[str, dict[str, Any]]:
    """Pre-render one transcript map per transcript, using the shared track primitive.

    Server-side because the drawing belongs in one place: a notebook figure and this card come out of
    the same function, so they cannot drift. Only the selection marker is left to the client, because
    only it depends on which guide is selected.
    """
    out: dict[str, dict[str, Any]] = {}
    for entry in payload.run.get("transcripts") or []:
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
        svg, geometry = transcript_map_svg(
            regions,
            series,
            title=f"{regions.transcript_id} - {regions.length:,} nt, {entry['windows']:,} candidate windows",
            value_label=entry.get("value_column") or "composite score",
            gaps=gaps,
            standalone=False,
        )
        counts = ", ".join(f"{len(s.points):,} {s.label}" for s in reversed(series))
        note = f"{counts}." + (
            f" {len(gaps)} stretch(es) of {MIN_UNCOVERED_NT} nt or more carry no candidate in this table."
            if gaps
            else ""
        )
        out[regions.transcript_id] = {
            "label": f"{regions.transcript_id} - {regions.length:,} nt, {entry['windows']:,} candidates",
            "svg": svg,
            "legend": legend_html(series, gaps=bool(gaps)),
            "geometry": geometry.as_dict(),
            "viewBox": _view_box(svg),
            "note": note,
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
