"""Installed CLI coverage for target intent across isoforms."""

import hashlib
import json
import os
import subprocess

import pandas as pd
import pytest

pytestmark = [pytest.mark.integration, pytest.mark.runs_in_container]


@pytest.mark.parametrize(
    ("selectivity", "required_ids", "expected_selective_pass"),
    [
        ("pan_isoform", [], False),
        ("isoform_selective", ["ENST00000000001"], True),
    ],
)
def test_cli_target_intent_changes_candidate_eligibility(tmp_path, selectivity, required_ids, expected_selective_pass):
    """The public workflow preserves the full annotation denominator and selective eligibility."""
    guide = "AACCTGTAAGCTGATCGTACG"
    target = guide.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    first_sequence = f"{'G' * 35}{target}{'C' * 130}"
    second_sequence = "A" * len(first_sequence)

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
                "excluded_transcript_ids": [],
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
    assert result.returncode == 0, f"{result.stdout}\n{result.stderr}"

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

    workflow_summary = json.loads((output / "logs/workflow_summary.json").read_text())
    seed_summary = workflow_summary["transcript_seed"]
    assert seed_summary["status"] == "complete"
    assert seed_summary["required_species"] == ["human"]
    seed_observations = seed_summary["per_candidate"]
    assert all(
        seed_observations[candidate_id]["by_species"]["human"]["status"] == "complete"
        for candidate_id in target_candidates["id"]
    )
    assert "max_transcript_seed_sites_verdict" in candidates
    assert set(candidates.loc[target_candidates.index, "max_transcript_seed_sites_verdict"]) <= {"pass", "fail"}

    qualified = pd.read_csv(output / "sirnaforge/candidates_qualified.csv")
    qualified_fasta = output / "sirnaforge/candidates_qualified.fasta"
    if expected_selective_pass:
        assert guide in set(qualified["guide_sequence"]), target_candidates[
            ["id", "selection_state", "passes_filters", "quality_issues"]
        ].to_dict(orient="records")
        assert qualified_fasta.is_file()
        assert guide in qualified_fasta.read_text()
    else:
        assert guide not in set(qualified["guide_sequence"])
