"""Evidence state and threshold provenance survive paged and sharded reporting."""

import json
import re
import shutil
import subprocess
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from sirnaforge.reporting import write_report_from_run
from tests.unit.test_reporting_payload import GUIDE, OTHER, _write_fully_evidenced_run

pytestmark = pytest.mark.unit


@pytest.mark.parametrize("mode", ["paged_single_file", "data_backed"])
@pytest.mark.parametrize(
    ("acquisition", "screening", "warning_id"),
    [
        ("failed", "failed", "MIRNA_ACQUISITION_FAILED"),
        ("complete", "failed", "MIRNA_SCREENING_FAILED"),
        ("skipped", "skipped", "MIRNA_SCREENING_SKIPPED"),
        ("complete", "complete", None),
    ],
)
def test_staged_report_retains_evidence_warnings_and_policy(
    tmp_path: Path, mode: str, acquisition: str, screening: str, warning_id: str | None
) -> None:
    """Missing matches stay distinct from completed zero hits, with stable copy and inspectable provenance."""
    run = _write_fully_evidenced_run(tmp_path)
    manifest_path = run / "sirnaforge/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    reference = {
        "source": "mirbase",
        "requested_species": ["human"],
        "acquisition": {"status": acquisition},
        "screening": {"status": screening, "guide_sequences": [GUIDE]},
        "sources": [
            {
                "source": "mirbase",
                "species": "human",
                "acquisition": {"status": acquisition},
                "screening": {"status": screening},
                "references": {"release": "22.1", "acquisition_mode": "pinned_import", "raw_sha256": "a" * 64},
            }
        ],
    }
    manifest["reference_evidence"] = {"schema_version": 1, "mirna": reference}
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / "report.html"
    result = write_report_from_run(run, output, mode=mode)
    assert result.artifact_dir is not None
    guides = pq.read_table(result.artifact_dir / "guides.parquet").to_pylist()
    assert guides[0]["mirna_warning_id"] == warning_id
    assert guides[0]["status"] == "pass", "evidence inspection must not reinterpret the run's gate results"
    assert result.manifest["summary"]["status_counts"]["pass"] == 1
    assert result.manifest["summary"]["passes_without_unresolved_evidence"] == (1 if warning_id is None else 0)
    html = output.read_text()
    assert 'aria-label="Evidence and thresholds"' in html
    assert "<dt>Stage</dt>" in html and "<dt>Scope</dt>" in html
    assert "22.1" in html and "pinned_import" in html
    if warning_id:
        assert f'data-warning-id="{warning_id}"' in html
        assert "Show unresolved evidence (1)" in html
        assert result.manifest["summary"]["evidence_warning_counts"] == {warning_id: 1}
    else:
        assert '<details id="evidence-notices">' not in html
        assert '<label hidden><input id="unresolved-only"' in html


def test_unresolved_view_and_guide_detail_use_registry_without_changing_gates(tmp_path: Path) -> None:
    """Execute filtering and a registry-driven guide notice with one submitted and one unsubmitted guide."""
    node = shutil.which("node")
    if node is None:
        pytest.skip("Node.js is required to execute the offline reader controls")
    run = _write_fully_evidenced_run(tmp_path, guides=(GUIDE, OTHER))
    manifest_path = run / "sirnaforge/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["reference_evidence"] = {
        "schema_version": 1,
        "mirna": {
            "source": "mirbase",
            "acquisition": {"status": "complete"},
            "screening": {"status": "complete", "guide_sequences": [GUIDE]},
        },
    }
    manifest_path.write_text(json.dumps(manifest))
    output = tmp_path / "paged.html"
    result = write_report_from_run(run, output, mode="paged_single_file")
    assert result.manifest["summary"]["unresolved_guides"] == 1
    scripts = re.findall(r"<script>(.*?)</script>", output.read_text(), flags=re.DOTALL)
    bootstrap = r"""
global.window={};const elements={};
function el(){return {checked:false,value:'',selectedIndex:0,textContent:'',innerHTML:'',handlers:{},addEventListener(k,fn){this.handlers[k]=fn;},appendChild(){},setAttribute(k,v){this[k]=v;}};}
global.document={querySelector:s=>elements[s]||(elements[s]=el()),querySelectorAll:()=>[],createElement:()=>el()};
"""
    code = (
        bootstrap
        + "\n".join(scripts)
        + "\nconst before=JSON.stringify(G.map(g=>g.gates));"
        + "elements['#unresolved-only'].checked=true;elements['#unresolved-only'].handlers.change();"
        + "if(!elements['#visible-count'].textContent.startsWith('1 guides'))throw Error('wrong unresolved count');"
        + f"showGuide({json.dumps(OTHER)}).then(()=>{{"
        + "if(!elements['#details'].innerHTML.includes('MIRNA_EVIDENCE_UNKNOWN'))throw Error('missing warning identity');"
        + "if(!elements['#details'].innerHTML.includes('Repeat screening with provenance enabled'))throw Error('missing remedy');"
        + "if(before!==JSON.stringify(G.map(g=>g.gates)))throw Error('gates changed');process.stdout.write('ok');});"
    )
    script = tmp_path / "reader.js"
    script.write_text(code)
    completed = subprocess.run([node, str(script)], capture_output=True, text=True, timeout=10, check=True)
    assert completed.stdout == "ok"
