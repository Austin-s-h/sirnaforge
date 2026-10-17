"""``benchmark design``'s filter-accounting contract: one run, two verdict sets (#109, bm-design slice).

Builds a prepared artifact directory by hand -- ``manifest.json``/``observations.csv`` written
through the real :mod:`sirnaforge.benchmark.artifact` models, ``design_inputs.fasta`` written
directly -- rather than going through ``benchmark prepare`` (a sibling slice's module), so this
test exercises only what :func:`sirnaforge.benchmark.design.design_artifact` itself is responsible
for. Guide sequences are synthetic and constructed so a specific gate is the one under test:

* ``LOW_GC_OBS``: GC content below the shipped ``gc_min`` floor (30.0) but above a widened floor,
  isolating the widening-vs-narrowing distinction the whole module exists to get right.
* ``HOMOPOLYMER_OBS``: a 4-nt run, GC well inside range, isolating ``max_poly_runs`` -- the one gate
  #109 requires to stay active under every policy this surface can resolve.
* ``TOO_SHORT_OBS``: its ``design_context_id`` names a FASTA record shorter than the paired length,
  so ``SiRNADesigner`` enumerates zero windows for it -- the "no candidate was ever built" case.

No override is needed for ``max_repeat_transcript_fraction``: that gate's ``FilterSpec`` declares
stage ``design`` (always "evaluable" per the resolver), but the only code that ever records an
observed value for it (``SiRNADesigner.stamp_repeat_verdict``) is invoked from ``workflow.py``'s
post-screening step against a real cDNA reference, never from ``design_from_file`` -- native
transcript mapping is explicitly out of scope for #109. ``design_artifact`` excludes it from both
verdict sets for exactly that reason (``_UNMEASURED_BY_DESIGN_FROM_FILE``), the same way it already
excludes the nine ``post_screen`` gates a design-only run holds no evidence for.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from Bio.Seq import Seq

from sirnaforge.benchmark.artifact import (
    ACCOUNTING_FILENAME,
    DESIGN_INPUTS_FASTA_FILENAME,
    MANIFEST_FILENAME,
    OBSERVATIONS_FILENAME,
    BenchmarkArtifactCounts,
    BenchmarkArtifactManifest,
    BenchmarkArtifactOutputs,
    BenchmarkObservation,
    FilterExclusionCounts,
    GCWideningBlock,
    GCWideningEntry,
    ManifestOutputEntry,
    ManifestPanelBlock,
    PolynucleotideRunRequirementBlock,
    read_accounting,
    read_manifest,
    write_manifest,
    write_observations,
)
from sirnaforge.benchmark.design import design_artifact
from sirnaforge.config.run_policy import RunPolicyError, default_for
from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation, SettingSource

PAIRED_LENGTH = 21

# GC ~9.5%: below the shipped gc_min floor (30.0), above a widened floor of 5.0. No homopolymer run
# longer than 1, so max_poly_runs cannot confound this row's verdict under either policy.
LOW_GC_GUIDE = "ATATATATATATATATGCATA"[:PAIRED_LENGTH]

# GC ~42.9%, comfortably inside [30, 60] under either policy; a leading 4-nt run of A is the one
# thing either policy can reject it for.
HOMOPOLYMER_GUIDE = "AAAAGCATGCATGCATGCATG"[:PAIRED_LENGTH]

assert len(LOW_GC_GUIDE) == PAIRED_LENGTH
assert len(HOMOPOLYMER_GUIDE) == PAIRED_LENGTH


def _observation(
    *, observation_id: str, guide: str, design_context_id: str, paired_length: int = PAIRED_LENGTH
) -> BenchmarkObservation:
    """A minimally-varying, fully-populated ``BenchmarkObservation`` for one synthetic guide."""
    return BenchmarkObservation(
        observation_id=observation_id,
        panel_id="synthetic_fixture",
        source_row_index=0,
        architecture="fully_complementary",
        assay_label="synthetic",
        measured_endpoint="inhibition_fraction",
        measured_value=None,
        full_guide_sequence=guide,
        guide_length=len(guide),
        paired_guide_sequence=guide,
        paired_slice_start_1based=1,
        paired_length=paired_length,
        passenger_sequence=None,
        guide_3p_overhang=None,
        passenger_3p_overhang=None,
        duplex_pairing_status="unstated",
        source_citation="Synthetic fixture, no real citation (#109 bm-design slice)",
        source_redistribution=None,
        source_file="<synthetic>",
        source_sha256="0" * 64,
        target_transcript_id=None,
        target_identity_status="unavailable",
        target_start_1based=None,
        target_end_1based=None,
        target_strand=None,
        design_context_id=design_context_id,
        design_context_source="measured_target_site",
        split=None,
        compatibility_status="compatible",
        compatibility_reason="",
    )


def _placeholder_manifest(
    *, observations_path: Path, fasta_path: Path, observation_count: int
) -> BenchmarkArtifactManifest:
    """A manifest shaped like one ``benchmark prepare`` would have written, before ``design`` runs.

    Every field ``BenchmarkArtifactManifest`` declares is required at once (frozen, no partial
    construction), so a "prepare-only" manifest still needs syntactically valid ``run_policy``/
    ``gc_widening``/``polynucleotide_run_requirement`` blocks -- ``design_artifact`` overwrites every
    one of them, so their content here is a placeholder, not a claim.
    """
    return BenchmarkArtifactManifest(
        tool_version="0.0.0-fixture",
        created_utc="2026-01-01T00:00:00Z",
        invoked_command=("sirnaforge", "benchmark", "prepare", "--panel", "synthetic_fixture"),
        panel=ManifestPanelBlock(
            panel_id="synthetic_fixture",
            display_name="Synthetic fixture panel",
            architecture="fully_complementary",
            citation="Synthetic fixture, no real citation (#109 bm-design slice)",
            redistribution=None,
            data_present=True,
            descriptor_hash="sha256:" + "0" * 64,
        ),
        paired_length=PAIRED_LENGTH,
        requested_length=PAIRED_LENGTH,
        split_rule=None,
        inputs=(),
        outputs=BenchmarkArtifactOutputs(
            observations_csv=ManifestOutputEntry(
                path=str(observations_path),
                exists=True,
                size_bytes=observations_path.stat().st_size,
                sha256="0" * 64,
                rows=observation_count,
            ),
            design_inputs_fasta=ManifestOutputEntry(
                path=str(fasta_path),
                exists=True,
                size_bytes=fasta_path.stat().st_size,
                sha256="0" * 64,
                sequences=2,
            ),
        ),
        counts=BenchmarkArtifactCounts(
            observations_in_source=observation_count,
            observations_kept=observation_count,
            observations_incompatible=0,
            mapped_native=0,
            mapped_panel_local=0,
            mapping_unavailable=observation_count,
            entered_design=0,
            no_candidate=0,
            default_pass=0,
            benchmark_pass=0,
            excluded_by_filter={},
        ),
        polynucleotide_run_requirement=PolynucleotideRunRequirementBlock(
            comparator=FilterComparator.LE,
            threshold=3,
            action=FilterAction.FAIL,
            evaluated=True,
            excluded=FilterExclusionCounts(default=0, benchmark=0),
        ),
        run_policy={},
        default_run_policy={},
        gc_widening=GCWideningBlock(
            gc_min=GCWideningEntry(default=30.0, benchmark=30.0, source=SettingSource.BUILTIN_PROFILE, widened=False),
            gc_max=GCWideningEntry(default=60.0, benchmark=60.0, source=SettingSource.BUILTIN_PROFILE, widened=False),
        ),
    )


@pytest.fixture
def artifact_dir(tmp_path: Path) -> Path:
    """A prepared artifact directory with three observations: widenable-GC, homopolymer, no-context."""
    directory = tmp_path / "synthetic_fixture__len21"
    directory.mkdir()

    low_gc_context = str(Seq(LOW_GC_GUIDE).reverse_complement())
    homopolymer_context = str(Seq(HOMOPOLYMER_GUIDE).reverse_complement())
    too_short_context = "ATATATATAT"  # 10 nt: shorter than PAIRED_LENGTH, so zero windows enumerate
    assert len(low_gc_context) == PAIRED_LENGTH
    assert len(homopolymer_context) == PAIRED_LENGTH
    assert len(too_short_context) < PAIRED_LENGTH

    fasta_path = directory / DESIGN_INPUTS_FASTA_FILENAME
    fasta_path.write_text(
        f">low_gc_ctx\n{low_gc_context}\n>homopolymer_ctx\n{homopolymer_context}\n>too_short_ctx\n{too_short_context}\n"
    )

    observations = [
        _observation(observation_id="synthetic_fixture:000000", guide=LOW_GC_GUIDE, design_context_id="low_gc_ctx"),
        _observation(
            observation_id="synthetic_fixture:000001", guide=HOMOPOLYMER_GUIDE, design_context_id="homopolymer_ctx"
        ),
        _observation(observation_id="synthetic_fixture:000002", guide=LOW_GC_GUIDE, design_context_id="too_short_ctx"),
    ]
    observations_path = directory / OBSERVATIONS_FILENAME
    write_observations(observations, observations_path)

    manifest = _placeholder_manifest(
        observations_path=observations_path, fasta_path=fasta_path, observation_count=len(observations)
    )
    write_manifest(manifest, directory / MANIFEST_FILENAME)
    return directory


def test_widened_gc_min_passes_benchmark_but_fails_default(artifact_dir: Path) -> None:
    """A candidate below the shipped GC floor enters design, passes wide, fails narrow (#109 bm-design)."""
    design_artifact(artifact_dir=artifact_dir, gc_min=5.0)

    rows = {row.observation_id: row for row in read_accounting(artifact_dir / ACCOUNTING_FILENAME)}
    row = rows["synthetic_fixture:000000"]

    assert row.entered_design is True
    assert row.benchmark_filter_status is FilterEvaluation.PASS
    assert row.default_filter_status is FilterEvaluation.FAIL
    assert row.default_filter_reasons == "gc_content_min"


def test_homopolymer_fails_both_policies_identically(artifact_dir: Path) -> None:
    """max_poly_runs is not GC-widenable, so a 4-nt run fails identically under both policies."""
    design_artifact(artifact_dir=artifact_dir, gc_min=5.0)

    rows = {row.observation_id: row for row in read_accounting(artifact_dir / ACCOUNTING_FILENAME)}
    row = rows["synthetic_fixture:000001"]

    assert row.entered_design is True
    assert "max_poly_runs" in row.default_filter_reasons.split(";")
    assert "max_poly_runs" in row.benchmark_filter_reasons.split(";")
    assert row.default_filter_status is FilterEvaluation.FAIL
    assert row.benchmark_filter_status is FilterEvaluation.FAIL


def test_no_candidate_keeps_the_observation_row_not_evaluated(artifact_dir: Path) -> None:
    """A context shorter than the paired length yields zero candidates; the observation still gets a row."""
    design_artifact(artifact_dir=artifact_dir, gc_min=5.0)

    rows = {row.observation_id: row for row in read_accounting(artifact_dir / ACCOUNTING_FILENAME)}
    row = rows["synthetic_fixture:000002"]

    assert row.entered_design is False
    assert row.candidate_id is None
    assert row.default_filter_status is FilterEvaluation.NOT_EVALUATED
    assert row.benchmark_filter_status is FilterEvaluation.NOT_EVALUATED
    assert row.default_filter_reasons == ""
    assert row.benchmark_filter_reasons == ""


def test_gc_min_narrowing_is_refused_before_any_design_work(artifact_dir: Path) -> None:
    """A gc_min tighter than the shipped default raises RunPolicyError, and touches nothing on disk."""
    accounting_path = artifact_dir / ACCOUNTING_FILENAME
    manifest_bytes_before = (artifact_dir / MANIFEST_FILENAME).read_bytes()

    with pytest.raises(RunPolicyError):
        design_artifact(artifact_dir=artifact_dir, gc_min=default_for("gc_min") + 1.0)

    assert not accounting_path.exists()
    assert (artifact_dir / MANIFEST_FILENAME).read_bytes() == manifest_bytes_before


def test_gc_max_narrowing_is_refused(artifact_dir: Path) -> None:
    """The same refusal in the other direction: gc_max may not be lowered below the shipped default."""
    with pytest.raises(RunPolicyError):
        design_artifact(artifact_dir=artifact_dir, gc_max=default_for("gc_max") - 1.0)


def test_manifest_policy_blocks_share_profile_hash_and_default_reports_shipped_gc_min(artifact_dir: Path) -> None:
    """Both policy blocks resolve against the same built-in profile; only the default reports its gc_min."""
    design_artifact(artifact_dir=artifact_dir, gc_min=5.0)

    manifest = read_manifest(artifact_dir / MANIFEST_FILENAME)
    benchmark_block = manifest.run_policy
    default_block = manifest.default_run_policy

    assert benchmark_block["profile"]["content_hash"] == default_block["profile"]["content_hash"]

    default_gc_min = next(record["value"] for record in default_block["resolved_settings"] if record["key"] == "gc_min")
    assert default_gc_min == default_for("gc_min")

    benchmark_gc_min = next(
        record["value"] for record in benchmark_block["resolved_settings"] if record["key"] == "gc_min"
    )
    assert benchmark_gc_min == 5.0
    assert manifest.gc_widening.gc_min.widened is True
    assert manifest.gc_widening.gc_min.default == default_for("gc_min")
