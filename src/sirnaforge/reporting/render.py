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
    <summary>Thresholds — <span id="fcount"></span></summary>
    <div class="fstat" id="fstat"></div>
    <div id="frows"></div>
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
     thresholds above were applied to. <b>Export TSV</b> downloads it. Some viewers -- including a
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
function rowKey(g,k){
  if(k==='guide')return g.guide; if(k==='composite')return g.composite_score??-1;
  if(k==='design')return g.design_score??-1; if(k==='failed')return g.n_gates_failed*100+g.n_gates_unknown;
  if(k==='status')return -STATUS_RANK[g.status];
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
    <td><span class="pill ${V[g.status]}">${g.status}</span></td>
    <td class="mono">${esc(g.guide)}</td><td>${fmt(g.composite_score)}</td>
    <td>${g.n_gates_failed} / ${g.n_gates_unknown}</td>
    <td class="${g.liability_count?'liab':'nonliab'}">${g.liability_count}</td>
    <td title="${g.n_rows} enumeration${g.n_rows===1?'':'s'}">${isoformFrac(g)}</td></tr>`).join('');
  tb.querySelectorAll('tr').forEach(tr=>tr.onclick=e=>{
    const cell=e.target.closest('td.pick');
    if(cell){ togglePick(cell.dataset.pick); return; }   // the pick column selects, it does not navigate
    show(G[+tr.dataset.i]);
  });
}
document.querySelectorAll('#idx th').forEach(th=>th.onclick=()=>{
  const k=th.dataset.k; if(k===sortK) sortAsc=!sortAsc; else {sortK=k;sortAsc=(k==='guide');} renderIndex();});
document.getElementById('q').oninput = () => applyFilters();

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
function gateReason(f,value,verdict,reason){
  const measured = (value!==null&&value!==undefined);
  if(reason===1) return 'the run exports neither '+f.filter_id+'_observed nor '+f.column;
  if(reason===2) return (f.read_column||f.column)+' is empty';
  // An off gate still measured something on most runs. Saying only "filter is off" hid a real
  // per-guide number -- and max_mirna_1mm_seed's own definition promises the number is reported.
  if(reason===3) return measured ? 'filter is off; value measured, nothing acts on it'
                                 : 'filter is off, and the run exports no value for '+f.column;
  if(reason===4) return measured ? 'no threshold declared; value measured, nothing acts on it'
                                 : 'no threshold declared, and no value exported';
  const cmp=fmt(value)+(verdict===0?' ':' not ')+f.comparator+' '+fmt(f.threshold);
  return verdict===4 ? cmp+'; action=warn, not a rejection' : cmp;
}
function gatesCard(g){
  const order=[3,0,1,4,2];  // fail, unknown, warn, pass, not_evaluated -- worst news first
  const idx=g.gates.map((t,i)=>i).sort((a,b)=>order[g.gates[a][1]]-order[g.gates[b][1]]);
  const rows=idx.map(i=>{const [value,verdict,reason]=g.gates[i], f=FILTERS[i], v=VERDICT[verdict];
    return `<tr><td class="mono">${esc(f.filter_id)}</td>
    <td><span class="pill v-${v}">${v.replace('_',' ')}</span></td>
    <td>${fmt(value)}</td><td class="mono">${esc(f.comparator)} ${fmt(f.threshold)}</td>
    <td>${esc(f.scope_label)}</td><td>${esc(f.stage)}</td><td>${esc(gateReason(f,value,verdict,reason))}</td></tr>`;}).join('');
  const nu=g.n_gates_unknown;
  return `<div class="card"><h2>Gates — all ${g.gates.length}, independently evaluated</h2>
    ${nu?`<div class="warn"><b>${nu} of ${g.gates.length} gates not evaluated.</b> An unevaluated gate is
      not a pass. See the Why column for the missing input.</div>`:''}
    ${(!nu&&g.status==='unknown')?`<div class="warn"><b>Run verdict: ${esc(g.run_verdict)}.</b>
      No declared filter covers this rejection, so it cannot be re-derived from the ${g.gates.length}
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

document.getElementById('tx').onchange = e => { curTx = e.target.value; renderMap(); };

// One card raising must cost that card, not the pane. Every panel below is independent evidence, so
// losing the off-target table because a structure string was malformed is never the right trade.
function card(render, g){
  try { return render(g); }
  catch(e){ return `<div class="card"><h2>Panel failed to render</h2>
    <div class="warn"><b>${esc(e && e.message || String(e))}</b> The other panels are unaffected.</div></div>`; }
}

// ---- thresholds and cart -------------------------------------------------------------------------
// Client-side only: the run's verdicts are fixed and are never recomputed here. A threshold below
// *selects among* candidates the run already judged; it does not re-judge them, and the pill on each
// row still shows what the run and the gate panel concluded.
const CART = new Set();
const STATUS_KEYS = ['pass','warn','unknown','fail'];

//: label, how to read the value off a guide, direction, and the bound the field starts at.
const FILTERS_UI = [
  {k:'composite', label:'composite score', get:g=>g.composite_score, dir:'min'},
  {k:'isoforms',  label:'isoforms hit',    get:g=>g.transcript_hits,  dir:'min'},
  {k:'gc',        label:'GC %',            get:g=>g.metrics.gc_content, dir:'range'},
  {k:'asym',      label:'asymmetry',       get:g=>g.metrics.asymmetry_score, dir:'min'},
  {k:'offt',      label:'off-target count',get:g=>g.metrics.off_target_count, dir:'max'},
  {k:'liab',      label:'liabilities',     get:g=>g.liability_count, dir:'max'},
];
const F = {status:new Set(['pass','warn']), cons:new Set(), consSeedIntact:true};

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

function fnum(id){ const el=document.getElementById(id); if(!el) return null;
  const v=el.value.trim(); return v===''?null:Number(v); }

function passesFilters(g){
  if(F.status.size && !F.status.has(g.status)) return false;
  if(!passesConservation(g)) return false;
  for(const f of FILTERS_UI){
    const v=f.get(g);
    const lo=fnum('f_'+f.k+'_lo'), hi=fnum('f_'+f.k+'_hi');
    if(lo!==null || hi!==null){
      if(v===null||v===undefined) return false;   // an absent value cannot satisfy a threshold
      if(lo!==null && v<lo) return false;
      if(hi!==null && v>hi) return false;
    }
  }
  return true;
}

function buildFilterUI(){
  document.getElementById('fstat').innerHTML = STATUS_KEYS.map(k=>
    `<label><input type="checkbox" data-status="${k}"${F.status.has(k)?' checked':''}> ${k}</label>`).join('');
  document.getElementById('frows').innerHTML = FILTERS_UI.map(f=>{
    const lo = f.dir==='max' ? '' : `<input id="f_${f.k}_lo" placeholder="min" inputmode="decimal">`;
    const hi = f.dir==='min' ? '' : `<input id="f_${f.k}_hi" placeholder="max" inputmode="decimal">`;
    return `<div class="frow"><span>${esc(f.label)}</span>${lo||'<span></span>'}${hi||'<span></span>'}</div>`;
  }).join('');
  document.querySelectorAll('#fstat input').forEach(cb=>cb.onchange=()=>{
    cb.checked ? F.status.add(cb.dataset.status) : F.status.delete(cb.dataset.status); applyFilters(); });
  document.querySelectorAll('#frows input').forEach(el=>el.oninput=applyFilters);
  document.getElementById('fcons').innerHTML =
    '<span style="color:var(--mut)">conserved in</span> ' +
    CONS_SPECIES.map(sp=>`<label><input type="checkbox" data-cons="${sp}"${F.cons.has(sp)?' checked':''}> ${sp}</label>`).join('') +
    `<label title="accept one mismatch outside guide positions 2-8">
      <input type="checkbox" id="consseed"${F.consSeedIntact?' checked':''}> allow 1 mm outside the seed</label>`;
  document.querySelectorAll('#fcons input[data-cons]').forEach(cb=>cb.onchange=()=>{
    cb.checked ? F.cons.add(cb.dataset.cons) : F.cons.delete(cb.dataset.cons); applyFilters(); });
  document.getElementById('consseed').onchange = e => { F.consSeedIntact=e.target.checked; applyFilters(); };
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
       <td><span class="pill ${({pass:'v-pass',warn:'v-warn',unknown:'v-unknown',fail:'v-fail'})[g.status]}">${g.status}</span></td>
       <td>${fmt(g.composite_score)}</td><td>${isoformFrac(g)}</td>
       <td class="${g.liability_count?'liab':'nonliab'}">${g.liability_count}</td></tr>`).join('')}</tbody></table>`;
  document.querySelectorAll('#cartlist tr[data-cart]').forEach(tr=>
    tr.querySelector('td.pick').onclick=()=>togglePick(tr.dataset.cart));
  const cols=['guide','passenger','status','composite_score','design_score','isoforms_hit','isoforms_in_run',
              'gc_content','asymmetry_score','off_target_count','liabilities','structure',
              ...CONS_SPECIES.flatMap(sp=>[sp+'_nm', sp+'_seed_mm'])];
  document.getElementById('carttsv').value = [cols.join('\t'), ...rows.map(g=>[
    g.guide, g.passenger??'', g.status, g.composite_score??'', g.design_score??'',
    g.transcript_hits??'', TX_IDS.length, g.metrics.gc_content??'', g.metrics.asymmetry_score??'',
    g.metrics.off_target_count??'', g.liability_count, g.structure??'',
    ...CONS_SPECIES.flatMap(sp=>{const d=consOf(g,sp); return [d?.nm??'', d?.seed_mismatches??''];})
    ].join('\t'))].join('\n');
}

document.getElementById('addtop').onclick = e => {
  if(e.target.id==='topn') return;                       // typing in the field is not a click on the button
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
document.getElementById('freset').onclick = () => {
  F.status = new Set(['pass','warn']); F.cons = new Set(); F.consSeedIntact = true;
  document.querySelectorAll('#frows input').forEach(el=>{ el.value=''; });
  buildFilterUI(); applyFilters();
};

function show(g){
  selected=g.guide; renderIndex();
  const i=G.indexOf(g);
  document.getElementById('detail').innerHTML = `
   <div class="card"><h2>Guide</h2>
     <div class="mono big">${esc(g.guide)} <span class="pill ${({pass:'v-pass',warn:'v-warn',unknown:'v-unknown',fail:'v-fail'})[g.status]}" style="font-size:12px;vertical-align:middle">${g.status==='unknown'?'not established':g.status}</span></div>
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
  showOn(g);
  location.hash = encodeURIComponent(g.guide);
}
// One bootstrap block, at the end: renderMap reaches CART through drawOverlay, and a `const` is
// hoisted but not initialised, so calling it from the map section threw a TDZ ReferenceError.
buildFilterUI();
applyFilters();
renderCart();
renderMap();
const initial = decodeURIComponent(location.hash.slice(1));
const start = G.find(g=>g.guide===initial) || G[0];
if(start) show(start);
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
        html = html.replace(placeholder, json.dumps(value, separators=(",", ":")))
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
