"""Installed-CLI variant integration, independent of live gene/variant services.

Synthetic cached variants exercise report generation and repeated --snp/coordinate
parsing. They do not establish allele-specific design accuracy. Live TP53 retrieval
and rsID resolution remain in tests/integration/test_variant_live.py.
"""

import json
import os
import subprocess

import pandas as pd
import pytest
from Bio import SeqIO

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
    # A synthetic genomic exon maps this actual transcript sequence explicitly.
    # Keep one transcript so unrelated, unannotated input cannot mask mapping failures.
    record = next(SeqIO.parse(realistic_transcripts_fasta, "fasta"))
    fasta = tmp_path / "transcript.fa"
    SeqIO.write([record], fasta, "fasta")
    annotation = tmp_path / "annotation.gtf"
    attributes = f'gene_id "SYNTHETIC_GENE"; transcript_id "{record.id}";'
    annotation.write_text(
        f"1\tfixture\ttranscript\t1001\t{1000 + len(record)}\t.\t+\t.\t{attributes}\n"
        f"1\tfixture\texon\t1001\t{1000 + len(record)}\t.\t+\t.\t{attributes}\n"
    )
    offsets = [100, 200] if multiple else [100]
    coordinates = [
        f"chr1:{1001 + offset}:{record.seq[offset]}:{'A' if record.seq[offset] != 'A' else 'G'}" for offset in offsets
    ]
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
        str(fasta),
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
        env={
            **os.environ,
            "SIRNAFORGE_CACHE_DIR": str(cache),
            "SIRNAFORGE_ANNOTATION_PROVIDER": "gencode",
            "SIRNAFORGE_GENCODE_GTF": str(annotation),
            "SIRNAFORGE_GENCODE_ASSEMBLY": "GRCh38",
        },
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

    candidates = pd.read_csv(output / "sirnaforge/candidates_all.csv")
    assert not candidates.empty
    if mode == "target":
        assert candidates["targeted_alleles"].map(json.loads).map(lambda alleles: alleles == ["alt"]).all()
    elif mode == "both":
        alleles = candidates["targeted_alleles"].map(json.loads)
        assert any("alt" in value for value in alleles)
        assert any("ref" in value for value in alleles)
    else:
        assert candidates["overlapped_variants"].map(json.loads).map(len).eq(0).all()
