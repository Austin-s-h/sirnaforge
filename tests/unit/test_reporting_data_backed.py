"""Vertical contract tests for the offline, shard-loaded report mode."""

from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pyarrow.parquet as pq
import pytest
from typer.testing import CliRunner

from sirnaforge.cli import app
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
    candidates.append(_candidate_row("c1199-alt", guides[-1], "ENST_ALT", n_guides + 1))
    hits = [
        _hit_row(guide, "human", f"ENSTH{index}", 1, "off_target", f"GENE{index}") for index, guide in enumerate(guides)
    ]
    run = _write_run(tmp_path, candidates, hits)
    mirna_table = run / "off_target/results/aggregated/combined_mirna_hits.tsv"
    mirna_table.write_text(
        f"qseq\tmirna_id\tdatabase\tcoord\tseed_mismatches\n{guides[-1]}\tHsa-Mir-24-P2_3p\tmirgenedb\t1\t0\n"
    )
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
        ]
        + [
            {
                "transcript_id": "ENST_ALT",
                "gene_id": "ENSG1199",
                "gene_symbol": "GENE1199",
                "transcript_biotype": "protein_coding",
                "is_canonical": False,
                "canonical_status": "established",
                "canonical_source": "test authority",
                "canonical_release": "110",
            }
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
    assert manifest["summary"]["candidate_rows"] == n_guides + 1
    assert manifest["summary"]["gene_filter_available"] is True
    assert len(manifest["summary"]["target_gene_options"]) == n_guides
    assert manifest["guide_index"]["rows"] == n_guides
    assert manifest["reader"]["guide_index_file"] == "guide-index.js"
    guide_index_path = output.with_name("large-report.details") / manifest["guide_index"]["file"]
    guide_index = json.loads(guide_index_path.read_text().removeprefix("window.SF_GUIDE_INDEX=").removesuffix(";\n"))
    assert len(guide_index) == n_guides
    assert guide_index[0]["guide"] == guides[0]
    assert pq.read_metadata(artifact_root / "hit_summary.parquet").num_rows == n_guides
    isoforms = pq.read_table(artifact_root / "isoforms.parquet").to_pylist()
    canonical = next(row for row in isoforms if row["transcript_id"] == "ENST0")
    assert canonical["canonical"] == "canonical"
    assert canonical["canonical_source"] == "test authority"
    guide_rows = pq.read_table(artifact_root / "guides.parquet").to_pylist()
    assert next(row for row in guide_rows if row["guide"] == guides[0])["canonical_hit"] == "hit"
    assert next(row for row in guide_rows if row["guide"] == guides[1])["canonical_hit"] == "not_hit"
    assert manifest["reader"]["offline_file_origin_supported"] is True
    diagnostics = json.loads(output.with_suffix(".diagnostics.json").read_text())
    assert diagnostics["peak_memory_scope"] == "process_lifetime"
    assert diagnostics["phases_seconds"]["normalized_table_build"] >= 0
    assert diagnostics["output_html_bytes"] == output.stat().st_size
    html = output.read_text()
    assert '<script src="large-report.details/guide-index.js"></script>' in html
    assert "window.SF_GUIDE_INDEX=[" not in html
    # Gene labels belong to the summary/identity metadata; detailed hit rows belong in shards.
    detailed_hit_is_inline = "ENSTH1199" in html
    assert not detailed_hit_is_inline
    assert len(manifest["detail_shards"]) == 10
    assert "Canonical transcript" in html
    assert "Target transcript identity" in html
    assert "function sourceAnchor" in html
    assert "Named miRNA matches" in html

    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for the offline detail-shard execution check")
    shard_path = output.with_name("large-report.details") / manifest["detail_shards"][-1]["file"]
    script = (
        "const fs=require('fs'),vm=require('vm');global.window={};"
        f"vm.runInThisContext(fs.readFileSync({json.dumps(str(shard_path))},'utf8'));"
        f"const d=window.SF_DETAIL_SHARDS[9][{json.dumps(guides[-1])}];"
        "if(!d||d.hits.length!==1||d.hits[0].hit_symbol!=='GENE1199')process.exit(2);"
        "if(d.isoforms.length!==2)process.exit(3);"
        "const iso=d.isoforms.find(i=>i.transcript_id==='ENST1199');if(!iso||iso.canonical!=='other_isoform')process.exit(13);"
        "const tx=JSON.parse(iso.transcript_link);"
        "if(tx.url!=='https://e110.ensembl.org/id/ENST1199')process.exit(4);"
        "const mirna=d.mirna[0].source_link;"
        "if(mirna.url!=='https://mirgenedb.org/show/hsa/Mir-24-P2')process.exit(5);"
        "process.stdout.write(String(Object.keys(window.SF_DETAIL_SHARDS[9]).length));"
    )
    result = subprocess.run([node, "-e", script], text=True, capture_output=True, check=True, timeout=10)
    assert result.stdout == "48"

    html_scripts = re.findall(r"<script>(.*?)</script>", html, flags=re.DOTALL)
    assert len(html_scripts) == 1
    ui_script = tmp_path / "check-data-backed-ui.js"
    ui_script.write_text(
        r"""
const fs=require('fs'),vm=require('vm');global.window={};
const elements={};function el(){return {value:'',selectedIndex:0,checked:false,textContent:'',innerHTML:'',disabled:false,attrs:{},handlers:{},addEventListener(k,fn){this.handlers[k]=fn;},appendChild(){},setAttribute(k,v){this.attrs[k]=v;}}}
for(const id of ['#counts','#gatecounts','#unresolved-only','#geneenable','#targetgene','#gene-state','#search','#status','#gate','#gate-result','#visible-count','#page-info','#previous','#next','#rows','#details'])elements[id]=el();
elements['#gate'].appendChild=()=>{};
global.window.location={hash:''};window.addEventListener=()=>{};
global.document={querySelector:s=>elements[s],querySelectorAll:()=>[],createElement:()=>el(),head:{appendChild:()=>{}}};
"""
        + f"vm.runInThisContext(fs.readFileSync({json.dumps(str(guide_index_path))},'utf8'));\n"
        + f"vm.runInThisContext({json.dumps(html_scripts[0])});\n"
        + f"vm.runInThisContext(fs.readFileSync({json.dumps(str(shard_path))},'utf8'));\n"
        + f"const detailData=window.SF_DETAIL_SHARDS[9][{json.dumps(guides[-1])}];"
        + "const txInfo=detailData.isoforms.find(i=>i.transcript_id==='ENST1199');"
        + "const map=candidateMap(detailData.candidates,'ENST1199',txInfo);"
        + "if(!map.includes('role=\"img\"')||!map.includes('transcript position')||!map.includes('GENE1199 · ENST1199 · Other isoform'))process.exit(11);"
        + "const context=candidateContext(detailData.candidates,detailData.isoforms);"
        + "if(!context.includes('GENE1199 · ENST1199 · Other isoform')||!context.includes('GENE1199 · ENST_ALT · Other isoform')||!context.includes('candidate-transcript'))process.exit(12);"
        + "elements['#geneenable'].checked=true;elements['#targetgene'].value='id:ENSG1199';elements['#geneenable'].handlers.change();"
        + "if(elements['#visible-count'].textContent!=='1 guides match current filters')process.exit(6);"
        + "if(!window.location.hash.includes('gene=id%3AENSG1199'))process.exit(7);"
        + "if(elements['#gene-state'].textContent!=='Target gene filter active: id:ENSG1199')process.exit(8);"
        + "elements['#geneenable'].checked=false;elements['#targetgene'].value='';restoreFilterState();renderRows();"
        + "if(!elements['#geneenable'].checked||elements['#targetgene'].value!=='id:ENSG1199'||elements['#targetgene'].disabled)process.exit(9);"
        + "if(elements['#visible-count'].textContent!=='1 guides match current filters')process.exit(10);"
        + f"showGuide({json.dumps(guides[-1])}).then(()=>{{"
        + "const detail=elements['#details'].innerHTML;"
        + "for(const value of ['Canonical transcript','Other isoform','https://e110.ensembl.org/id/ENST1199','https://mirgenedb.org/show/hsa/Mir-24-P2'])"
        + "if(!detail.includes(value))process.exit(2);"
        + "process.stdout.write('ok');});\n"
    )
    ui_result = subprocess.run([node, str(ui_script)], text=True, capture_output=True, check=True, timeout=10)
    assert ui_result.stdout == "ok"


@pytest.mark.unit
def test_oligogym_scale_report_keeps_details_out_of_the_html_tail(tmp_path: Path) -> None:
    """Native OligoGym guide observations complete in data-backed mode with synthetic hit stress rows."""
    observations = Path("tests/data/benchmarks/oligogym_native_design/observations.csv")
    with observations.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    assert len(source_rows) == 4_113
    guides = [str(row["full_guide_sequence"]) for row in source_rows]
    candidates = [
        _candidate_row(str(row["benchmark_id"]), guide, "", index + 1)
        for index, (row, guide) in enumerate(zip(source_rows, guides, strict=True))
    ]
    # These alignments are intentionally synthetic and exist only to stress staged hit/detail handling.
    hits = [
        _hit_row(guide, "human", f"SYNTHETIC_TARGET_{index:06d}", 1, "off_target", "SYNTHETIC")
        for index, guide in enumerate(guides)
        for _ in range(5)
    ]
    run = _write_run(tmp_path, candidates, hits)
    candidate_path = run / "sirnaforge/candidates_all.csv"
    candidate_path.write_text(
        "id,guide_sequence,transcript_id,position\n"
        + "\n".join(f"oligogym-{index},{guide},,{index + 1}" for index, guide in enumerate(guides))
        + "\n"
    )
    output = tmp_path / "oligogym-report.html"

    result = write_report_from_run(run, output, mode="data_backed")

    assert result.artifact_dir is not None
    manifest = json.loads((result.artifact_dir / "report_manifest.json").read_text())
    assert result.mode == "data_backed"
    assert manifest["summary"]["candidate_rows"] == 4_113
    assert manifest["summary"]["guides"] == len(set(guides))
    assert manifest["tables"]["hits"]["rows"] == 4_113 * 5
    assert len(manifest["detail_shards"]) == (len(set(guides)) + 127) // 128
    assert manifest["diagnostics"]["peak_rss_bytes"] > 0
    assert manifest["diagnostics"]["output_html_bytes"] == output.stat().st_size
    detail_bytes = sum(item["bytes"] for item in manifest["detail_shards"])
    assert output.stat().st_size < 10_000_000
    assert detail_bytes > output.stat().st_size
    assert "SYNTHETIC_TARGET_000000" not in output.read_text()
    first_guide = pq.read_table(result.artifact_dir / "guides.parquet", columns=["guide"]).column(0)[0].as_py()
    first_detail = result.output.with_name(f"{result.output.stem}.details") / manifest["detail_shards"][0]["file"]
    assert first_guide in first_detail.read_text()


@pytest.mark.unit
def test_data_backed_gene_filter_stays_hidden_for_single_gene_with_many_isoforms(tmp_path: Path) -> None:
    """Transcript count alone does not turn a confirmed single-gene report into a multi-gene report."""
    run = _write_run(
        tmp_path,
        [
            _candidate_row("c1", GUIDE, "ENST1", 1),
            _candidate_row("c2", GUIDE, "ENST2", 2),
        ],
        [],
    )
    manifest_path = run / "sirnaforge/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["target_identity"] = {
        "schema_version": 1,
        "transcripts": [
            {
                "transcript_id": transcript,
                "gene_id": "ENSG1",
                "gene_symbol": "TP53",
                "is_canonical": transcript == "ENST1",
                "canonical_status": "established",
                "canonical_source": "test authority",
                "canonical_release": "110",
            }
            for transcript in ("ENST1", "ENST2")
        ],
        "resolvers": {"ensembl": {"template": "https://www.ensembl.org/id/{identifier}"}},
    }
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / "single-gene.html"

    write_report_from_run(run, output, mode="data_backed")

    summary = json.loads((output.with_suffix(".data") / "report_manifest.json").read_text())["summary"]
    assert summary["gene_filter_available"] is False
    assert summary["target_gene_options"] == [{"key": "id:ENSG1", "label": "TP53 · ENSG1"}]
    assert 'id="geneenable"' not in output.read_text()


@pytest.mark.unit
@pytest.mark.parametrize("mode", ["paged_single_file", "data_backed"])
def test_staged_gencode_provenance_does_not_pin_ensembl_links_to_gencode_release(tmp_path: Path, mode: str) -> None:
    """Staged isoform links retain GENCODE authority versions separately from Ensembl URLs."""
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST1.2", 1)], [])
    manifest_path = run / "sirnaforge/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["target_identity"] = {
        "schema_version": 1,
        "transcripts": [
            {
                "transcript_id": "ENST1.2",
                "gene_id": "ENSG1.3",
                "gene_symbol": "TP53",
                "is_canonical": True,
                "canonical_status": "established",
                "canonical_source": "GENCODE explicit canonical/MANE tags (annotation 49)",
                "canonical_release": "49",
            }
        ],
    }
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / f"gencode-{mode}.html"

    result = write_report_from_run(run, output, mode=mode)

    isoforms = pq.read_table(result.artifact_dir / "isoforms.parquet").to_pylist()
    assert len(isoforms) == 1
    assert isoforms[0]["canonical_release"] == "49"
    assert json.loads(isoforms[0]["transcript_link"])["url"] == "https://www.ensembl.org/id/ENST1.2"
    assert json.loads(isoforms[0]["gene_link"])["url"] == "https://www.ensembl.org/id/ENSG1.3"


@pytest.mark.unit
@pytest.mark.parametrize("mode", ["paged_single_file", "data_backed"])
def test_staged_conservation_filters_distinguish_matches_and_unavailable_evidence(tmp_path: Path, mode: str) -> None:
    """Staged UI filters recorded ortholog evidence and explains unavailable zero-match cases."""
    guides = [_base4(index) for index in range(4)]
    run = _write_run(
        tmp_path,
        [_candidate_row(f"c{index}", guide, f"ENST{index}", index + 1) for index, guide in enumerate(guides)],
        [
            _hit_row(guides[0], "mouse", "ENSMUST1", 1, "ortholog", "GeneA"),
            _hit_row(guides[1], "mouse", "ENSMUST2", 1, "ortholog", "GeneB"),
        ],
    )
    candidate_path = run / "sirnaforge/candidates_all.csv"
    candidate_lines = candidate_path.read_text().splitlines()
    candidate_path.write_text(
        candidate_lines[0]
        + ",off_target_screened\n"
        + "\n".join(line + "," + ("False" if index == 2 else "True") for index, line in enumerate(candidate_lines[1:]))
        + "\n"
    )
    hit_path = run / "off_target/results/aggregated/combined_offtargets.tsv"
    hit_lines = hit_path.read_text().splitlines()
    hit_path.write_text(hit_lines[0] + "\n" + hit_lines[1] + "\n" + hit_lines[2].replace("\t0\t", "\t1\t", 1) + "\n")
    manifest_path = run / "sirnaforge/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["provenance"] = {
        "databases": {
            "orthology": {
                "source": "Ensembl Compara",
                "attempted": True,
                "compara_release": {"value": "110"},
                "conservation_scope": {
                    "denominator_species": ["mouse"],
                    "unscreened_in_denominator": {"value": []},
                },
                "per_species": [{"species": "mouse", "status": "resolved", "gene_ids": ["ENSG1"]}],
            }
        }
    }
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / f"conservation-{mode}.html"
    write_report_from_run(run, output, mode=mode)

    artifact_root = output.with_suffix(".data")
    staged = json.loads((artifact_root / "report_manifest.json").read_text())
    assert staged["conservation_evidence"]["per_species"]["mouse"]["screening_complete"] is True
    guide_rows = pq.read_table(artifact_root / "guides.parquet").to_pylist()
    by_guide = {row["guide"]: row for row in guide_rows}
    assert json.loads(by_guide[guides[0]]["ortholog"])["mouse"] == {
        "nm": 1,
        "seed_mismatches": 0,
        "transcript_id": "ENSMUST1",
        "gene_symbol": "GeneA",
    }
    assert json.loads(by_guide[guides[1]]["ortholog"])["mouse"] == {
        "nm": 1,
        "seed_mismatches": 1,
        "transcript_id": "ENSMUST2",
        "gene_symbol": "GeneB",
    }
    assert by_guide[guides[2]]["off_target_screened"] is False
    assert by_guide[guides[3]]["off_target_screened"] is True
    html = output.read_text()
    assert '<fieldset id="conservation-filters"' in html
    assert 'role="status" aria-live="polite"' in html
    guide_index_path = output.with_name(f"conservation-{mode}.details") / "guide-index.js"
    scripts = re.findall(r"<script>(.*?)</script>", html, flags=re.DOTALL)
    bootstrap = r"""
const vm=require('vm');global.window={location:{hash:'#c=mouse'},addEventListener(){}};
const elements={};function el(){return {value:'',selectedIndex:0,checked:false,textContent:'',innerHTML:'',disabled:false,handlers:{},addEventListener(k,fn){this.handlers[k]=fn;},appendChild(){},setAttribute(k,v){this[k]=v;}}}
for(const id of ['#counts','#gatecounts','#unresolved-only','#search','#status','#gate','#gate-result','#visible-count','#page-info','#previous','#next','#rows','#consdiagnostic','#cons-state','#cons-seed','#details'])elements[id]=el();
const mouse=el();mouse.value='mouse';global.document={querySelector:s=>elements[s],querySelectorAll:s=>s==='[data-conservation]'?[mouse]:[],createElement:()=>el(),head:{appendChild(){}}};
"""
    guide_index_loader = (
        f"\nvm.runInThisContext(fs.readFileSync({json.dumps(str(guide_index_path))},'utf8'));\n"
        if mode == "data_backed"
        else ""
    )
    code = bootstrap + guide_index_loader + "\n".join(scripts)
    code += "if(!mouse.checked||elements['#cons-seed'].checked!==true)process.exit(2);"
    code += "if(elements['#visible-count'].textContent!=='1 guides match current filters')process.exit(3);"
    code += "if(!elements['#cons-state'].textContent.startsWith('Conservation filter active: mouse'))process.exit(10);"
    code += "mouse.checked=false;mouse.handlers.change();if(window.location.hash.includes('c='))process.exit(4);"
    code += "mouse.checked=true;mouse.handlers.change();elements['#cons-seed'].checked=false;elements['#cons-seed'].handlers.change();"
    code += "if(elements['#visible-count'].textContent!=='2 guides match current filters')process.exit(5);"
    code += "G[0].ortholog='{}';elements['#cons-seed'].checked=true;elements['#cons-seed'].handlers.change();"
    code += "if(elements['#visible-count'].textContent!=='0 guides match current filters')process.exit(6);"
    code += (
        "if(!elements['#consdiagnostic'].textContent.includes('does not establish biological absence'))process.exit(7);"
    )
    code += (
        "if(conservationDiagnostic(G[3],'mouse').label!=='no conserved site found in completed screen')process.exit(9);"
    )
    code += (
        "M.conservation_evidence.per_species.mouse.state='legacy';"
        "const legacy=conservationDiagnostic(G[3],'mouse');"
        "if(!legacy.unavailable||legacy.label!=='orthology evidence not established')process.exit(11);"
    )
    code += (
        f"G[0].ortholog={json.dumps(json.dumps({'mouse': {'nm': 1, 'seed_mismatches': 0, 'transcript_id': 'ENSMUST1', 'gene_symbol': 'GeneA'}}))};"
        f"window.SF_DETAIL_SHARDS={{0:{{{json.dumps(guides[0])}:{{candidates:[],hits:[],mirna:[],isoforms:[]}}}}}};"
        f"showGuide({json.dumps(guides[0])}).then(()=>{{const detail=elements['#details'].innerHTML;"
        "for(const value of ['Cross-species conservation','Ensembl Compara','ENSMUST1','GeneA'])"
        "if(!detail.includes(value))process.exit(8);});"
    )
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required for the staged conservation UI test")
    result = subprocess.run([node, "-e", code], check=False, text=True, capture_output=True, timeout=10)
    assert result.returncode == 0, result.stderr


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
    manifest_path = run / "sirnaforge/manifest.json"
    run_manifest = json.loads(manifest_path.read_text())
    run_manifest["target_identity"] = {
        "schema_version": 1,
        "transcripts": [
            {"transcript_id": f"ENST{index}", "gene_id": f"ENSG{index}", "gene_symbol": f"GENE{index}"}
            for index in range(n_guides)
        ],
    }
    manifest_path.write_text(json.dumps(run_manifest))
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
const elements={};function el(){return {value:'',selectedIndex:0,textContent:'',innerHTML:'',disabled:false,checked:false,attrs:{},handlers:{},addEventListener(k,fn){this.handlers[k]=fn;},setAttribute(k,v){this.attrs[k]=v;}}}
for(const id of ['#counts','#gatecounts','#unresolved-only','#geneenable','#targetgene','#gene-state','#search','#status','#gate','#gate-result','#visible-count','#page-info','#previous','#next','#rows','#details'])elements[id]=el();
elements['#gate'].appendChild=()=>{};global.window.location={hash:''};window.addEventListener=()=>{};global.document={querySelector:s=>elements[s],querySelectorAll:()=>[],createElement:()=>el(),head:{appendChild:()=>{}}};
"""
    code = (
        bootstrap
        + "\n".join(scripts)
        + "\nif(!elements['#page-info'].textContent.endsWith('of 3'))process.exit(2);"
        + "elements['#next'].handlers.click();if(!elements['#page-info'].textContent.startsWith('Page 2'))process.exit(3);"
        + "elements['#geneenable'].checked=true;elements['#targetgene'].value='id:ENSG119';elements['#geneenable'].handlers.change();"
        + "if(elements['#visible-count'].textContent!=='1 guides match current filters'||!window.location.hash.includes('gene=id%3AENSG119'))process.exit(7);"
        + "elements['#geneenable'].checked=false;elements['#geneenable'].handlers.change();"
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


@pytest.mark.unit
def test_cli_data_backed_quilt_does_not_build_monolithic_payload(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CLI report generation and Quilt registration stay on the summary-only path."""
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST1", 1)], [])
    output = tmp_path / "cli-report.html"

    def forbidden_payload(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("data-backed CLI must not build the monolithic payload")

    monkeypatch.setattr("sirnaforge.cli.build_payload", forbidden_payload)
    result = CliRunner().invoke(
        app, ["report", str(run), "-o", str(output), "--mode", "data_backed", "--quilt-summarize"]
    )

    assert result.exit_code == 0, result.output
    assert output.exists() and output.with_name("cli-report.details").is_dir()
    assert (tmp_path / "quilt_summarize.json").exists()
    entries = json.loads((tmp_path / "quilt_summarize.json").read_text())
    assert entries[0]["description"].startswith("1 guides:")


@pytest.mark.unit
def test_report_cli_forwards_mode_artifact_and_threshold_options(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The public CLI passes user-selected report budgets into the writer API."""
    calls: list[dict[str, object]] = []
    run = tmp_path / "run"
    run.mkdir()
    output = tmp_path / "custom.html"
    artifact_dir = tmp_path / "report-data"

    def fake_writer(_run: Path, out: Path, **kwargs: object) -> object:
        calls.append(kwargs)
        out.write_text("<html></html>")
        diagnostics_path = out.with_suffix(".diagnostics.json")
        diagnostics_path.write_text("{}")
        return SimpleNamespace(
            output=out,
            mode="paged_single_file",
            manifest={"summary": {"guides": 0, "candidate_rows": 0, "status_counts": {}}},
            diagnostics={},
            diagnostics_path=diagnostics_path,
        )

    monkeypatch.setattr("sirnaforge.cli.write_report_from_run", fake_writer)
    result = CliRunner().invoke(
        app,
        [
            "report",
            str(run),
            "-o",
            str(output),
            "--mode",
            "paged_single_file",
            "--artifact-dir",
            str(artifact_dir),
            "--compact-max-input-bytes",
            "111",
            "--compact-max-candidate-rows",
            "222",
            "--compact-max-guides",
            "3",
            "--paged-max-input-bytes",
            "444",
            "--paged-max-candidate-rows",
            "555",
            "--paged-max-guides",
            "6",
            "--diagnostics-json",
            str(tmp_path / "custom-diagnostics.json"),
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls[0]["mode"] == "paged_single_file"
    assert calls[0]["artifact_dir"] == artifact_dir
    assert calls[0]["thresholds"].compact_max_input_bytes == 111  # type: ignore[union-attr]
    assert calls[0]["thresholds"].paged_max_guides == 6  # type: ignore[union-attr]
    assert (tmp_path / "custom-diagnostics.json").exists()


@pytest.mark.unit
@pytest.mark.parametrize("mode", ["compact_single_file", "auto"])
def test_compact_cli_handles_absent_mode_facts(tmp_path: Path, mode: str) -> None:
    """Compact reports can print their summary and register Quilt with nullable mode facts."""
    run = _write_run(tmp_path, [_candidate_row("c1", GUIDE, "ENST1", 1)], [])
    output = tmp_path / "compact.html"
    result = CliRunner().invoke(app, ["report", str(run), "-o", str(output), "--mode", mode, "--quilt-summarize"])
    assert result.exit_code == 0, result.output
    assert "1 guides from" in result.output
    assert (tmp_path / "quilt_summarize.json").is_file()
