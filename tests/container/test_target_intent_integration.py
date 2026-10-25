"""Installed CLI coverage for target intent across isoforms."""

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
            "--species",
            "human",
            "--skip-off-targets",
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
    target_verdicts = [observations[candidate_id]["verdict"] for candidate_id in target_candidates["id"]]
    assert ("pass" in target_verdicts) is expected_selective_pass
    assert ("fail" in target_verdicts) is not expected_selective_pass
