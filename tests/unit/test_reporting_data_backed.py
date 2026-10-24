"""Vertical contract tests for the offline, shard-loaded report mode."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from sirnaforge.reporting import data_backed as data_backed_module
from sirnaforge.reporting import write_report_from_run
from tests.unit.test_reporting_payload import GUIDE, _candidate_row, _hit_row, _write_run


def _base4(value: int, length: int = 21) -> str:
    letters = "ACGU"
    chars = []
    for _ in range(length):
        chars.append(letters[value & 3])
        value >>= 2
    return "".join(reversed(chars))


@pytest.mark.unit
def test_data_backed_report_keeps_every_guide_and_loads_local_detail_shard(tmp_path: Path) -> None:
    """A large fixture stays queryable, reports exact status totals, and loads detail from local JS."""
    n_guides = 1_200
    guides = [_base4(index) for index in range(n_guides)]
    candidates = [_candidate_row(f"c{index}", guide, f"ENST{index}", index + 1) for index, guide in enumerate(guides)]
    hits = [
        _hit_row(guide, "human", f"ENSTH{index}", 1, "off_target", f"GENE{index}") for index, guide in enumerate(guides)
    ]
    run = _write_run(tmp_path, candidates, hits)
    manifest_path = run / "sirnaforge/manifest.json"
    run_manifest = json.loads(manifest_path.read_text())
    run_manifest["target_identity"] = {
        "schema_version": 1,
        "transcripts": [
            {
                "transcript_id": f"ENST{index}",
                "gene_id": f"ENSG{index}",
                "gene_symbol": f"GENE{index}",
                "transcript_biotype": "protein_coding",
                "is_canonical": index == 0,
                "canonical_status": "established",
                "canonical_source": "test authority",
                "canonical_release": "110",
            }
            for index in range(n_guides)
        ],
        "resolvers": {"ensembl": {"template": "https://www.ensembl.org/id/{identifier}"}},
    }
    manifest_path.write_text(json.dumps(run_manifest))
    output = tmp_path / "large-report.html"

    result = write_report_from_run(run, output, mode="data_backed")

    assert result.output == output
    assert result.mode == "data_backed"
    assert result.manifest_path.exists() and result.diagnostics_path.exists()
    assert result.diagnostics["mode_selection"]["mode"] == "data_backed"
    artifact_root = output.with_suffix(".data")
    manifest = json.loads((artifact_root / "report_manifest.json").read_text())
    assert manifest["report_mode"] == "data_backed"
    assert {"guides", "isoforms", "hit_summary"} <= set(manifest["tables"])
    assert manifest["summary"]["guides"] == n_guides
    # Every fixture candidate passes its recorded PASS run verdict, while the report's active
    # evidence gates observe no corresponding metrics in the fixture. The report must therefore
    # retain every guide as unknown rather than infer a pass from the run label.
    assert manifest["summary"]["status_counts"] == {"fail": 0, "unknown": n_guides, "warn": 0, "pass": 0}
    assert manifest["summary"]["candidate_rows"] == n_guides
    assert pq.read_metadata(artifact_root / "hit_summary.parquet").num_rows == n_guides
    isoforms = pq.read_table(artifact_root / "isoforms.parquet").to_pylist()
    canonical = next(row for row in isoforms if row["transcript_id"] == "ENST0")
    assert canonical["canonical"] == "canonical"
    assert canonical["canonical_source"] == "test authority"
    guide_rows = pq.read_table(artifact_root / "guides.parquet").to_pylist()
    assert next(row for row in guide_rows if row["guide"] == guides[0])["canonical_hit"] == "hit"
    assert manifest["reader"]["offline_file_origin_supported"] is True
    diagnostics = json.loads(output.with_suffix(".diagnostics.json").read_text())
    assert diagnostics["peak_memory_scope"] == "process_lifetime"
    assert diagnostics["phases_seconds"]["normalized_table_build"] >= 0
    assert diagnostics["output_html_bytes"] == output.stat().st_size
    html = output.read_text()
    # Gene labels belong to the summary/identity metadata; detailed hit rows belong in shards.
    detailed_hit_is_inline = "ENSTH1199" in html
    assert not detailed_hit_is_inline
    assert len(manifest["detail_shards"]) == 10

    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for the offline detail-shard execution check")
    shard_path = output.with_name("large-report.details") / manifest["detail_shards"][-1]["file"]
    script = (
        "const fs=require('fs'),vm=require('vm');global.window={};"
        f"vm.runInThisContext(fs.readFileSync({json.dumps(str(shard_path))},'utf8'));"
        f"const d=window.SF_DETAIL_SHARDS[9][{json.dumps(guides[-1])}];"
        "if(!d||d.hits.length!==1||d.hits[0].hit_symbol!=='GENE1199')process.exit(2);"
        "process.stdout.write(String(Object.keys(window.SF_DETAIL_SHARDS[9]).length));"
    )
    result = subprocess.run([node, "-e", script], text=True, capture_output=True, check=True, timeout=10)
    assert result.stdout == "48"


@pytest.mark.unit
def test_data_backed_writer_does_not_call_monolithic_payload_builder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The large mode is a staged path, not a wrapper around the full in-memory HTML payload."""
    run = _write_run(
        tmp_path,
        [_candidate_row("c1", GUIDE, "ENST1", 10)],
        [_hit_row(GUIDE, "human", "ENST9", 1, "off_target", "GENE9")],
    )

    def fail_if_called(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("data_backed called build_payload")

    monkeypatch.setattr(data_backed_module, "build_payload", fail_if_called)
    output = tmp_path / "report.html"
    result = write_report_from_run(run, output, mode="data_backed")
    assert result.mode == "data_backed"
    manifest = json.loads((output.with_suffix(".data") / "report_manifest.json").read_text())
    assert manifest["report_mode"] == "data_backed"


@pytest.mark.unit
def test_paged_single_file_navigates_all_filtered_guides_and_opens_inline_detail(tmp_path: Path) -> None:
    """Page navigation and filters span the full guide set while detail remains selectable offline."""
    n_guides = 120
    guides = [_base4(index) for index in range(n_guides)]
    run = _write_run(
        tmp_path,
        [_candidate_row(f"c{index}", guide, f"ENST{index}", index + 1) for index, guide in enumerate(guides)],
        [
            _hit_row(guide, "human", f"ENSTH{index}", 1, "off_target", f"GENE{index}")
            for index, guide in enumerate(guides)
        ],
    )
    output = tmp_path / "paged.html"
    result = write_report_from_run(run, output, mode="paged_single_file")

    assert result.mode == "paged_single_file"
    assert result.manifest["report_mode"] == "paged_single_file"
    assert result.diagnostics["mode_selection"]["thresholds"]
    assert "Page ${currentPage+1} of ${pageCount}" in output.read_text()
    assert not output.with_name("paged.details").exists()
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for the page-control execution check")
    scripts = re.findall(r"<script>(.*?)</script>", output.read_text(), flags=re.DOTALL)
    bootstrap = r"""
const vm=require('vm');global.window={};
const elements={};function el(){return {value:'',selectedIndex:0,textContent:'',innerHTML:'',disabled:false,handlers:{},addEventListener(k,fn){this.handlers[k]=fn;}}}
for(const id of ['#counts','#gatecounts','#search','#status','#gate','#gate-result','#visible-count','#page-info','#previous','#next','#rows','#details'])elements[id]=el();
elements['#gate'].appendChild=()=>{};global.document={querySelector:s=>elements[s],querySelectorAll:()=>[],createElement:()=>el(),head:{appendChild:()=>{}}};
"""
    code = (
        bootstrap
        + "\n".join(scripts)
        + "\nif(!elements['#page-info'].textContent.endsWith('of 3'))process.exit(2);"
        + "elements['#next'].handlers.click();if(!elements['#page-info'].textContent.startsWith('Page 2'))process.exit(3);"
        + f"elements['#search'].value={json.dumps(guides[-1])};elements['#search'].handlers.input();"
        + "if(!elements['#visible-count'].textContent.startsWith('1 guides'))process.exit(4);"
        + "elements['#search'].value='';elements['#gate'].selectedIndex=1;const target=G.find(g=>g.guide==="
        + f"{json.dumps(guides[-1])});const outcome=JSON.parse(target.gates)[0][1];elements['#gate-result'].value=String(outcome);elements['#gate-result'].handlers.change();"
        + "const expected=G.filter(g=>JSON.parse(g.gates)[0][1]===outcome).length;if(!elements['#visible-count'].textContent.startsWith(expected+' guides'))process.exit(6);"
        + f"showGuide({json.dumps(guides[-1])}).then(()=>{{if(!elements['#details'].innerHTML.includes('GENE119'))process.exit(5);process.stdout.write('ok');}});"
    )
    script_path = tmp_path / "check-paged.js"
    script_path.write_text(code)
    completed = subprocess.run([node, str(script_path)], text=True, capture_output=True, check=True, timeout=10)
    assert completed.stdout == "ok"
