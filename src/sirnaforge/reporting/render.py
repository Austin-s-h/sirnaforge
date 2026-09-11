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

from jinja2 import Environment, select_autoescape

from sirnaforge.reporting.payload import ReportPayload

#: Guides embedded in the index. The per-guide detail is rendered for all of them; this bounds only
#: how many rows the client holds, and the number is reported in the header so it is never a silent cap.
MAX_INDEX_GUIDES = 5000


_TEMPLATE = """<!DOCTYPE html>
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
main{display:grid;grid-template-columns:minmax(340px,38%) 1fr;gap:0;height:calc(100vh - 86px)}
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
    <span class="pill v-warn">{{ p.run.status_counts["warn"] }} pass with a warning</span>
    <span class="pill v-unknown">{{ p.run.status_counts["unknown"] }} not established</span>
    <span class="pill v-fail">{{ p.run.status_counts["fail"] }} fail</span>
    {% if p.run.agreement.comparable %}
    &nbsp;· vs the run's own verdict: <b>{{ p.run.agreement.contradicted_run_pass }}</b> of
    {{ p.run.agreement.run_pass_guides }} run-PASS guides contradicted;
    <b>{{ p.run.agreement.run_failed_not_rederivable }}</b> the run failed cannot be re-derived here
    — reported <em>not established</em>, never flipped to pass;
    <b>{{ p.run.agreement.overruled_run_fail }}</b> run rejections overruled.
    {% endif %}
  </div>
</header>
<main>
 <div id="left">
  <div class="searchbar"><input id="q" placeholder="Search guide sequence, candidate id or transcript…" autocomplete="off"></div>
  <table id="idx"><thead><tr>
    <th data-k="status">Status</th><th data-k="guide">Guide</th><th data-k="composite">Composite</th>
    <th data-k="failed">Fail / unknown</th><th data-k="liab">Liabilities</th><th data-k="rows">Isoforms</th>
  </tr></thead><tbody></tbody></table>
 </div>
 <div id="right"><p class="empty">Select a guide.</p></div>
</main>
<script>
const G = GUIDES_JSON_PLACEHOLDER;
const FILTERS = FILTERS_JSON_PLACEHOLDER;
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
    <td><span class="pill ${V[g.status]}">${g.status}</span></td>
    <td class="mono">${esc(g.guide)}</td><td>${fmt(g.composite_score)}</td>
    <td>${g.n_gates_failed} / ${g.n_gates_unknown}</td>
    <td class="${g.liability_count?'liab':'nonliab'}">${g.liability_count}</td><td>${g.n_rows}</td></tr>`).join('');
  tb.querySelectorAll('tr').forEach(tr=>tr.onclick=()=>show(G[+tr.dataset.i]));
}
document.querySelectorAll('#idx th').forEach(th=>th.onclick=()=>{
  const k=th.dataset.k; if(k===sortK) sortAsc=!sortAsc; else {sortK=k;sortAsc=(k==='guide');} renderIndex();});
document.getElementById('q').oninput = e => {
  const t=e.target.value.trim().toUpperCase().replace(/T/g,'U');
  const raw=e.target.value.trim().toLowerCase();
  view = !t ? G.slice() : G.filter(g => g.guide.includes(t)
      || g.isoforms.some(i => String(i.candidate_id).toLowerCase().includes(raw)
                           || String(i.transcript).toLowerCase().includes(raw)));
  renderIndex();
};

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
  if(reason===1) return 'the run exports neither '+f.filter_id+'_observed nor '+f.column;
  if(reason===2) return (f.read_column||f.column)+' is empty';
  if(reason===3) return 'filter is off';
  if(reason===4) return 'no threshold declared';
  const cmp=fmt(value)+(verdict===0?' ':' not ')+f.comparator+' '+fmt(f.threshold);
  return verdict===4 ? cmp+' — action is warn, so the run did not reject it' : cmp;
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
  return `<div class="card"><h2>Gates — all ${g.gates.length}, not only the first to fire</h2>
    ${nu?`<div class="warn"><b>${nu} gate${nu>1?'s':''} could not be evaluated.</b> An unknown is not a pass.
      The Why column names what each one was missing.</div>`:''}
    ${(!nu&&g.status==='unknown')?`<div class="warn"><b>The run rejected this guide as
      ${esc(g.run_verdict)}, and no declared gate expresses that.</b> Every gate below is satisfied, so
      this report cannot re-derive the rejection — and will not overrule it.</div>`:''}
    <table><thead><tr><th>Filter</th><th>Verdict</th><th>Value</th><th>Threshold</th>
    <th>Scope</th><th>Stage</th><th>Why</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}
function isoformCard(g){
  const any=g.isoforms.some(i=>i.register_neighbours.length);
  return `<div class="card"><h2>Isoforms — ${g.isoforms.length} enumeration${g.isoforms.length===1?'':'s'} of this one guide</h2>
   ${any?`<div class="warn"><b>Register neighbour.</b> Another design sits within 2 nt on the same transcript,
     so the two share a window and can score identically. On the reference panel two such designs
     differed 1.4&times; in measured knockdown — this is not a duplicate row.</div>`:''}
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
  return `<div class="card"><h2>miRNA seed resemblance — ${g.mirna.length} named hit(s)</h2>
   <table><thead><tr><th>miRNA</th><th>Source</th><th>Seed mismatches</th><th>nm</th></tr></thead>
   <tbody>${g.mirna.map(m=>`<tr><td class="mono">${esc(m.mirna)}</td><td>${esc(m.source)}</td>
     <td>${fmt(m.seed_mismatches)}</td><td>${fmt(m.nm)}</td></tr>`).join('')}</tbody></table></div>`;
}
function show(g){
  selected=g.guide; renderIndex();
  const i=G.indexOf(g);
  document.getElementById('right').innerHTML = `
   <div class="card"><h2>Guide</h2>
     <div class="mono big">${esc(g.guide)} <span class="pill ${({pass:'v-pass',unknown:'v-unknown',fail:'v-fail'})[g.status]}" style="font-size:12px;vertical-align:middle">${g.status==='unknown'?'not established':g.status}</span></div>
     <div class="kv" style="margin-top:12px">
       <div><b>Passenger</b><span class="mono">${esc(g.passenger||'—')}</span></div>
       <div><b>Overhang</b>${esc(g.overhang||'—')}</div>
       <div><b>Modifications</b>${esc(g.modifications||'—')}</div>
       <div><b>Composite</b>${fmt(g.composite_score)}${g.weight_vector?` <span class="pill v-not_evaluated">${esc(g.weight_vector)}</span>`:''}</div>
       <div><b>Design score</b>${fmt(g.design_score)}</div>
       ${Object.entries(g.metrics).map(([k,v])=>`<div><b>${esc(k)}</b>${fmt(v)}</div>`).join('')}
     </div></div>
   ${gatesCard(g)}${isoformCard(g)}${offtargetCard(g,i)}${mirnaCard(g)}`;
  if(g.offtarget_matrix.length) drawMatrix('mx'+i, g.offtarget_matrix);
  location.hash = encodeURIComponent(g.guide);
}
renderIndex();
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
    return html.replace("GUIDES_JSON_PLACEHOLDER", json.dumps(guides, separators=(",", ":"))).replace(
        "FILTERS_JSON_PLACEHOLDER", json.dumps(payload.filters, separators=(",", ":"))
    )


def write_report(payload: ReportPayload, out_path: Path | str) -> Path:
    """Render and write the report, returning the path written."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_html(payload), encoding="utf-8")
    return out
