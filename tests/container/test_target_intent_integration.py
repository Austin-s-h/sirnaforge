"""Installed CLI coverage for target intent across isoforms."""

import hashlib
import json
import os
import subprocess

import pandas as pd
import pytest

from sirnaforge.data.variant_resolver import VariantResolver
from sirnaforge.models.variant import ClinVarSignificance, VariantRecord, VariantSource

pytestmark = [pytest.mark.integration, pytest.mark.runs_in_container]


@pytest.mark.parametrize(
    ("selectivity", "required_ids", "excluded_ids", "seed_scope", "expected_selective_pass", "expected_qualified"),
    [
        ("pan_isoform", [], [], "full_cdna", False, False),
        ("isoform_selective", ["ENST00000000001"], ["ENST00000000002"], "full_cdna", True, True),
        ("isoform_selective", ["ENST00000000001"], ["ENST00000000002"], "three_prime_utr", True, False),
    ],
)
def test_cli_target_intent_changes_candidate_eligibility(
    tmp_path, selectivity, required_ids, excluded_ids, seed_scope, expected_selective_pass, expected_qualified
):
    """Intent, seed liability, excluded isoforms and missing UTRs stay distinct in the public workflow."""
    guide = "AACCTGTAAGCTGATCGTACG"
    target = guide.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    first_sequence = f"{'G' * 35}{target}{'C' * 130}"
    # The sibling has the reverse complement of guide positions 2–8, but not the full target site.
    # Pan-isoform intent must therefore fail full-site coverage while seed screening still reports
    # a real sibling-transcript liability.
    seed_target = guide[1:8].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    second_sequence = f"{'A' * 35}{seed_target}{'C' * 130}"

    fasta = tmp_path / "transcripts.fa"
    fasta.write_text(f">ENST00000000001.1\n{first_sequence}\n>ENST00000000002.1\n{second_sequence}\n")
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        "".join(
            f"1\tfixture\t{feature}\t1\t{len(sequence)}\t.\t+\t.\t"
            f'gene_id "ENSG00000000001"; transcript_id "{transcript_id}"; gene_name "SYNTHETIC_TARGET_INTENT";\n'
            for transcript_id, sequence in (
                ("ENST00000000001.1", first_sequence),
                ("ENST00000000002.1", second_sequence),
            )
            for feature in ("transcript", "exon")
        )
    )
    intent = tmp_path / "target-intent.json"
    intent.write_text(
        json.dumps(
            {
                "selectivity": selectivity,
                "target_species": ["human"],
                "offtarget_screen_species": [],
                "required_transcript_ids": required_ids,
                "excluded_transcript_ids": excluded_ids,
            }
        )
    )
    mirna_reference = tmp_path / "mature.fa"
    mirna_reference.write_text(">hsa-synthetic-fixture\nACGUACGUACGUACGUACGUA\n")
    mirna_manifest = tmp_path / "mirna-references.json"
    mirna_manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "references": [
                    {
                        "source": "mirbase",
                        "species": "human",
                        "release": "target-intent-offline-fixture",
                        "sha256": hashlib.sha256(mirna_reference.read_bytes()).hexdigest(),
                        "path": mirna_reference.name,
                    }
                ],
            }
        )
    )
    output = tmp_path / "output"
    result = subprocess.run(
        [
            "sirnaforge",
            "workflow",
            "SYNTHETIC_TARGET_INTENT",
            "--input-fasta",
            str(fasta),
            "--output-dir",
            str(output),
            "--target-intent-config",
            str(intent),
            "--run-mode",
            "qualified",
            "--species",
            "human",
            "--query-species",
            "human",
            "--transcriptome-fasta",
            f"human:{fasta}",
            "--mirna-db",
            "mirbase",
            "--mirna-reference-manifest",
            str(mirna_manifest),
            "--transcript-seed-species",
            "human",
            "--transcript-seed-required-species",
            "human",
            "--max-transcript-seed-sites",
            "1000",
            "--transcript-seed-scope",
            seed_scope,
            # This intent-focused fixture has only two synthetic isoforms, so the on-target sequence
            # itself exceeds the default repeat fraction. Keep that gate visible as a warning while
            # allowing the test to isolate the isoform-intent eligibility decision.
            "--filter-action",
            "max_repeat_transcript_fraction=warn",
        ],
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
        env={
            **os.environ,
            "SIRNAFORGE_ANNOTATION_PROVIDER": "gencode",
            "SIRNAFORGE_GENCODE_GTF": str(gtf),
            "SIRNAFORGE_GENCODE_ASSEMBLY": "GRCh38",
        },
    )
    expected_returncode = 2 if seed_scope == "three_prime_utr" else 0
    assert result.returncode == expected_returncode, f"{result.stdout}\n{result.stderr}"

    resolved = json.loads((output / "target_intent/resolved_intent.json").read_text())
    assert sorted(resolved["intent"]["annotation_universe"]) == ["ENST00000000001", "ENST00000000002"]
    assert sorted(resolved["intent"]["enumeration_inputs"]) == ["ENST00000000001", "ENST00000000002"]
    expected_denominator = (
        ["ENST00000000001"]
        if expected_selective_pass
        else [
            "ENST00000000001",
            "ENST00000000002",
        ]
    )
    assert sorted(resolved["intent"]["coverage_denominator"]) == expected_denominator

    observations = json.loads((output / "target_intent/guide_observations.json").read_text())["guides"]
    candidates = pd.read_csv(output / "sirnaforge/candidates_all.csv")
    assert not candidates.empty
    target_candidates = candidates[candidates["guide_sequence"] == guide]
    assert not target_candidates.empty
    assert set(candidates["id"]) == set(observations)
    coverage = target_candidates["isoform_coverage"].dropna().unique().tolist()
    assert coverage == ([1.0] if expected_selective_pass else [0.5])
    target_verdicts = [observations[candidate_id]["verdict"] for candidate_id in target_candidates["id"]]
    assert ("pass" in target_verdicts) is expected_selective_pass
    assert ("fail" in target_verdicts) is not expected_selective_pass
    assert all(
        observations[candidate_id]["coverage_numerator_ids"] == ["ENST00000000001"]
        and observations[candidate_id]["coverage_denominator_ids"] == expected_denominator
        for candidate_id in target_candidates["id"]
    )

    assert sorted(resolved["intent"]["excluded_transcript_ids"]) == sorted(excluded_ids)

    workflow_summary = json.loads((output / "logs/workflow_summary.json").read_text())
    seed_summary = workflow_summary["transcript_seed"]
    assert seed_summary["scope"] == seed_scope
    assert seed_summary["required_species"] == ["human"]
    seed_observations = seed_summary["per_candidate"]
    if seed_scope == "full_cdna":
        assert seed_summary["status"] == "complete"
        assert all(
            seed_observations[candidate_id]["by_species"]["human"]["status"] == "complete"
            for candidate_id in target_candidates["id"]
        )
        assert all(
            seed_observations[candidate_id]["by_species"]["human"]["site_count"] > 0
            for candidate_id in target_candidates["id"]
        )
    else:
        assert all(
            seed_observations[candidate_id]["by_species"]["human"]["status"] == "unknown"
            and seed_observations[candidate_id]["required_species_site_count"] is None
            for candidate_id in target_candidates["id"]
        )
    assert "max_transcript_seed_sites_verdict" in candidates
    expected_seed_verdicts = {"pass"} if seed_scope == "full_cdna" else {"unknown"}
    assert set(candidates.loc[target_candidates.index, "max_transcript_seed_sites_verdict"]) <= expected_seed_verdicts

    qualified = pd.read_csv(output / "sirnaforge/candidates_qualified.csv")
    qualified_fasta = output / "sirnaforge/candidates_qualified.fasta"
    if expected_qualified:
        assert guide in set(qualified["guide_sequence"]), target_candidates[
            ["id", "selection_state", "passes_filters", "quality_issues"]
        ].to_dict(orient="records")
        assert qualified_fasta.is_file()
        assert guide in qualified_fasta.read_text()
    else:
        assert guide not in set(qualified["guide_sequence"])
        expected_state = "not_eligible" if not expected_selective_pass else "withheld_incomplete_evidence"
        assert set(target_candidates["selection_state"]) == {expected_state}


def test_alt_variant_target_composes_with_isoform_intent_and_seed_liability(tmp_path):
    """A GRCh38 alt guide covers its selected isoform while the excluded sibling keeps a seed site."""
    guide_ref = "AACCTGTAAGCTGATCGTACG"
    target_ref = guide_ref.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    variant_offset = 10
    ref = target_ref[variant_offset]
    alt = next(base for base in "ACGT" if base != ref)
    target_alt = target_ref[:variant_offset] + alt + target_ref[variant_offset + 1 :]
    guide_alt = target_alt.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    first_sequence = f"{'G' * 35}{target_ref}{'C' * 130}"
    sibling_seed = guide_alt[1:8].translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sibling = list(f"{'A' * 35}{sibling_seed}{'C' * 130}")
    variant_offset_in_transcript = 35 + variant_offset
    sibling[variant_offset_in_transcript] = ref
    second_sequence = "".join(sibling)

    fasta = tmp_path / "transcripts.fa"
    fasta.write_text(f">ENST00000000001.1\n{first_sequence}\n>ENST00000000002.1\n{second_sequence}\n")
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        "".join(
            f"1\tfixture\t{feature}\t1\t{len(sequence)}\t.\t+\t.\t"
            f'gene_id "ENSG00000000001"; transcript_id "{transcript_id}"; gene_name "SYNTHETIC_TARGET_VARIANT";\n'
            for transcript_id, sequence in (
                ("ENST00000000001.1", first_sequence),
                ("ENST00000000002.1", second_sequence),
            )
            for feature in ("transcript", "exon")
        )
    )
    coordinate = f"chr1:{1 + variant_offset_in_transcript}:{ref}:{alt}"
    cache = tmp_path / "cache"
    resolver = VariantResolver(
        min_af=0.001,
        clinvar_filters=[ClinVarSignificance.PATHOGENIC, ClinVarSignificance.LIKELY_PATHOGENIC],
        cache_dir=cache / "variants",
        variant_mode="target",
    )
    query = resolver.parse_identifier(coordinate)
    resolver.cache.put(
        resolver._get_cache_key(query),
        VariantRecord(
            id="synthetic-target-variant",
            chr="chr1",
            pos=1 + variant_offset_in_transcript,
            ref=ref,
            alt=alt,
            sources=[VariantSource.LOCAL_VCF],
            af=0.05,
            clinvar_significance=ClinVarSignificance.PATHOGENIC,
        ),
    )
    intent = tmp_path / "target-intent.json"
    intent.write_text(
        json.dumps(
            {
                "selectivity": "isoform_selective",
                "target_species": ["human"],
                "offtarget_screen_species": [],
                "required_transcript_ids": ["ENST00000000001"],
                "excluded_transcript_ids": ["ENST00000000002"],
            }
        )
    )
    mirna_reference = tmp_path / "mature.fa"
    mirna_reference.write_text(">hsa-synthetic-fixture\nACGUACGUACGUACGUACGUA\n")
    mirna_manifest = tmp_path / "mirna-references.json"
    mirna_manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "references": [
                    {
                        "source": "mirbase",
                        "species": "human",
                        "release": "target-variant-fixture",
                        "sha256": hashlib.sha256(mirna_reference.read_bytes()).hexdigest(),
                        "path": mirna_reference.name,
                    }
                ],
            }
        )
    )
    output = tmp_path / "output"
    result = subprocess.run(
        [
            "sirnaforge",
            "workflow",
            "SYNTHETIC_TARGET_VARIANT",
            "--input-fasta",
            str(fasta),
            "--output-dir",
            str(output),
            "--target-intent-config",
            str(intent),
            "--run-mode",
            "qualified",
            "--species",
            "human",
            "--query-species",
            "human",
            "--transcriptome-fasta",
            f"human:{fasta}",
            "--mirna-db",
            "mirbase",
            "--mirna-reference-manifest",
            str(mirna_manifest),
            "--snp",
            coordinate,
            "--variant-mode",
            "target",
            "--min-af",
            "0.001",
            "--transcript-seed-species",
            "human",
            "--transcript-seed-required-species",
            "human",
            "--max-transcript-seed-sites",
            "1000",
            "--filter-action",
            "max_repeat_transcript_fraction=warn",
        ],
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
        env={
            **os.environ,
            "SIRNAFORGE_CACHE_DIR": str(cache),
            "SIRNAFORGE_ANNOTATION_PROVIDER": "gencode",
            "SIRNAFORGE_GENCODE_GTF": str(gtf),
            "SIRNAFORGE_GENCODE_ASSEMBLY": "GRCh38",
        },
    )
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"

    candidates = pd.read_csv(output / "sirnaforge/candidates_all.csv")
    alternate_candidates = candidates[
        (candidates["guide_sequence"] == guide_alt)
        & (candidates["transcript_id"] == "ENST00000000001.1")
        & candidates["targeted_alleles"].map(lambda value: json.loads(value) == ["alt"])
    ]
    assert len(alternate_candidates) == 1
    candidate = alternate_candidates.iloc[0]
    guide_observation = json.loads((output / "target_intent/guide_observations.json").read_text())["guides"][
        candidate["id"]
    ]
    assert guide_observation["verdict"] == "pass"
    assert guide_observation["coverage_numerator_ids"] == ["ENST00000000001"]
    assert guide_observation["excluded_hit_transcript_ids"] == []
    seed_summary = json.loads((output / "logs/workflow_summary.json").read_text())["transcript_seed"]
    assert seed_summary["per_candidate"][candidate["id"]]["by_species"]["human"]["site_count"] > 0
