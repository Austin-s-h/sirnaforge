"""Unit tests for the ZFN Nextflow bridge helpers."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from Bio.Seq import Seq

from sirnaforge.models.zfn import MatchOrientation, Strand, ZFNOffTargetSite
from sirnaforge.zfn.nextflow_bridge import make_zfn_shard_manifest, run_zfn_shard_search
from sirnaforge.zfn.search import ExhaustiveZFNOffTargetSearcher

LEFT = "GCGTACGTA"
RIGHT = "TACGGCATA"


def _write_fasta(tmp_path: Path, sequence: str, name: str = "chr3") -> Path:
    fasta = tmp_path / "search_space.fa"
    fasta.write_text(f">{name}\n{sequence}\n", encoding="utf-8")
    return fasta


def _canonical_site(left: str, right: str, spacer: str = "AAAAA") -> str:
    return f"{left}{spacer}{str(Seq(right).reverse_complement())}"


def test_make_zfn_shard_manifest_chunks_single_contig_when_enabled(tmp_path: Path) -> None:
    """Nextflow shard manifests should chunk single-contig references like the direct engine."""
    fasta = _write_fasta(tmp_path, "A" * 200)
    output_tsv = tmp_path / "zfn_shards.tsv"

    result = make_zfn_shard_manifest(
        genome_fasta=fasta,
        left_half_site=LEFT,
        right_half_site=RIGHT,
        spacer_lengths="5",
        max_mismatches=0,
        sharding_enabled=True,
        shard_chunk_mb=0.000025,
        shard_overlap_bp=0,
        shard_chromosomes="chr3",
        output_tsv=output_tsv,
    )

    with output_tsv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    assert result["sharding_active"] is True
    assert result["shards"] == len(rows)
    assert len(rows) > 1
    assert rows[0]["chrom"] == "chr3"
    assert rows[0]["core_start_1"] == "1"
    assert rows[1]["core_start_1"] == "26"


def test_run_zfn_shard_search_uses_region_search_api(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Shard execution should call the shard-bounded direct engine API rather than whole-search."""
    fasta = _write_fasta(tmp_path, _canonical_site(LEFT, RIGHT))
    output_sites_csv = tmp_path / "sites.csv"
    output_summary_json = tmp_path / "summary.json"

    def _fail_whole_search(*args: object, **kwargs: object) -> list[ZFNOffTargetSite]:
        raise AssertionError("whole-search API should not be used for shard execution")

    def _fake_region_search(
        self: ExhaustiveZFNOffTargetSearcher,
        params: object,
        chrom: str,
        scan_start0: int,
        scan_end0: int,
        core_start0: int | None = None,
        core_end0: int | None = None,
        annotation: object | None = None,
        top_n_sites: int | None = None,
    ) -> list[ZFNOffTargetSite]:
        assert chrom == "chr3"
        assert scan_start0 == 0
        assert core_start0 == 0
        assert top_n_sites is None
        return [
            ZFNOffTargetSite(
                site_id="chr3:1-23:L...R:mm0+0",
                chrom="chr3",
                start_1based=1,
                end_1based=23,
                strand=Strand.PLUS,
                orientation=MatchOrientation.LR,
                spacer_len=5,
                sequence=_canonical_site(LEFT, RIGHT),
                left_mismatches=0,
                right_mismatches=0,
                left_seed_mismatches=0,
                right_seed_mismatches=0,
                left_mismatch_positions=[],
                right_mismatch_positions=[],
                total_mismatches=0,
                score=100.0,
                score_components={},
                dimer_compatible=True,
                region="unknown",
                nearest_gene=None,
                left_aligned=LEFT,
                right_aligned=RIGHT,
            )
        ]

    monkeypatch.setattr(ExhaustiveZFNOffTargetSearcher, "search", _fail_whole_search)
    monkeypatch.setattr(ExhaustiveZFNOffTargetSearcher, "search_region", _fake_region_search)

    result = run_zfn_shard_search(
        shard_id="chr3:1-23",
        shard_chrom="chr3",
        scan_start_1=1,
        scan_end_1=23,
        core_start_1=1,
        core_end_1=23,
        shard_max_mismatches=0,
        left_half_site=LEFT,
        right_half_site=RIGHT,
        genome_fasta=fasta,
        algorithm="zfn_v2",
        dimer_mode="heterodimer_only",
        spacer_lengths="5",
        annotation_file=None,
        output_sites_csv=output_sites_csv,
        output_summary_json=output_summary_json,
    )

    summary = json.loads(output_summary_json.read_text(encoding="utf-8"))
    csv_lines = output_sites_csv.read_text(encoding="utf-8").strip().splitlines()

    assert result == {"shard_id": "chr3:1-23", "sites": 1, "chrom": "chr3"}
    assert summary["summary"]["off_target_sites"] == 1
    assert len(summary["candidates"]) == 1
    assert len(csv_lines) == 2
