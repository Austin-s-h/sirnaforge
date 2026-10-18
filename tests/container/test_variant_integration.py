"""Installed-CLI variant integration, independent of live gene/variant services.

Synthetic cached variants exercise report generation and repeated --snp/coordinate
parsing. They do not establish allele-specific design accuracy. Live TP53 retrieval
and rsID resolution remain in tests/integration/test_variant_live.py.
"""

import json
import os
import subprocess

import pytest

from sirnaforge.data.variant_resolver import VariantResolver
from sirnaforge.models.variant import ClinVarSignificance, VariantRecord, VariantSource

pytestmark = [pytest.mark.integration, pytest.mark.runs_in_container]


@pytest.mark.parametrize("mode", ["avoid", "target", "both"])
@pytest.mark.parametrize("multiple", [False, True])
def test_cached_variant_workflow(tmp_path, realistic_transcripts_fasta, mode, multiple):
    """All variant modes retain resolved records and render reports through the installed CLI."""
    cache = tmp_path / "cache"
    output = tmp_path / "output"
    min_af = 0.01 if mode == "avoid" else 0.001
    resolver = VariantResolver(
        cache_dir=cache / "variants",
        variant_mode=mode,
        min_af=min_af,
        clinvar_filters=[ClinVarSignificance.PATHOGENIC, ClinVarSignificance.LIKELY_PATHOGENIC],
    )
    coordinates = ["chr1:100:C:T", "chr1:200:G:A"] if multiple else ["chr1:100:C:T"]
    for coordinate in coordinates:
        query = resolver.parse_identifier(coordinate)
        variant = VariantRecord(
            chr=query.chr,
            pos=query.pos,
            ref=query.ref,
            alt=query.alt,
            sources=[VariantSource.LOCAL_VCF],
            af=0.05,
            clinvar_significance=ClinVarSignificance.PATHOGENIC,
        )
        resolver.cache.put(resolver._get_cache_key(query), variant)

    command = [
        "sirnaforge",
        "workflow",
        "SYNTHETIC_VARIANT_TEST",
        "--input-fasta",
        str(realistic_transcripts_fasta),
        "--output-dir",
        str(output),
        "--variant-mode",
        mode,
        "--min-af",
        str(min_af),
        "--species",
        "human",
        "--skip-off-targets",
        "--clinvar-filter-levels",
        "Pathogenic,Likely pathogenic",
    ]
    for coordinate in coordinates:
        command.extend(["--snp", coordinate])
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
        env={**os.environ, "SIRNAFORGE_CACHE_DIR": str(cache)},
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"
    report = json.loads((output / "logs/resolved_variants.json").read_text())
    assert report["gene"] == "SYNTHETIC_VARIANT_TEST"
    assert report["variant_mode"] == mode
    assert report["filters"]["min_af"] == min_af
    assert report["filters"]["clinvar_significance"] == ["Pathogenic", "Likely pathogenic"]
    assert report["summary"]["total_variants"] == len(coordinates)
    assert {variant["vcf_style"] for variant in report["variants"]} == set(coordinates)
    assert (output / "sirnaforge/candidates_all.csv").is_file()
    assert (output / "sirnaforge/report.html").is_file()
