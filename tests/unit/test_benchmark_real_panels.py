"""End-to-end ``prepare``+``design`` over the *vendored* OligoGym panels (#109, execution pass).

``tests/unit/test_benchmark_panels.py`` already re-derives every geometry claim in the registry
straight from ``tests/data/benchmarks/oligogym/records.csv`` (2,850/907/356, the revcomp relations,
the ``100 - label`` flip). Nothing here repeats that. What no file covered before is whether the two
*surfaces* -- ``sirnaforge.benchmark.prepare`` and ``sirnaforge.benchmark.design`` -- actually carry
those 4,113 real rows through to an artifact, and with which tallies. That is what this file pins:
per-panel prepare counts reconciled against the vendored ``manifest.json``, the design-stage
``entered_design``/``no_candidate``/``default_pass``/``benchmark_pass`` for a real panel, and the
invariance of the default verdict set under a widened benchmark run (#109's core accounting claim).

Two facts about *how* it reaches those surfaces, both deliberate:

1. Every panel is prepared by its **registry id**, through the shipped vendored path. This file used
   to declare a *probe* descriptor per panel -- the registry's own descriptor re-validated under a
   ``data_present=False`` id, which was the one legal way to hand ``prepare_artifact`` a
   ``--panel-csv`` -- and to pre-filter the shared table itself, because ``prepare.py`` named only
   ``huesken_subset`` as vendored and never called
   :meth:`~sirnaforge.benchmark.panels.PanelDescriptor.selects_row`. It now reads the path off
   :attr:`~sirnaforge.benchmark.panels.PanelDescriptor.vendored_csv` and honours the selector, so both
   work-arounds are gone and the tallies below are the ones ``sirnaforge benchmark prepare --panel
   ichihara`` itself reports.
2. The design stage here runs over ``design_inputs.fasta`` as ``prepare`` writes it -- one record per
   compatible observation, holding the *measured target site* (``design_context_source=
   measured_target_site``, the 21 nt reverse complement of the measured guide), not the 161 nt
   synthetic transcripts in ``tests/data/benchmarks/oligogym/inputs/``. So ``entered_design`` here
   means "the designer rebuilt this row's 19 nt core from its own measured site", and every one of
   the numbers below is evidence about design enumeration, never about native accessibility -- the
   coordinates these rows carry are positions in a fabricated context (see
   ``tests/data/benchmarks/README.md``).
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from sirnaforge.benchmark.artifact import read_accounting, read_manifest, read_observations
from sirnaforge.benchmark.design import design_artifact
from sirnaforge.benchmark.prepare import prepare_artifact

#: ``tests/unit/`` -> ``tests/`` -> repo root, the walk ``prepare.py`` itself does from ``src/``.
REPO_ROOT = Path(__file__).resolve().parents[2]
OLIGOGYM_DIR = REPO_ROOT / "tests" / "data" / "benchmarks" / "oligogym"
RECORDS_CSV = OLIGOGYM_DIR / "records.csv"
VENDORED_MANIFEST = OLIGOGYM_DIR / "manifest.json"

#: The bytes every count in this file rests on. Pinned independently of ``manifest.json`` so a
#: regenerated ``records.csv`` fails here loudly instead of quietly shifting a tally, the same
#: two-files-must-agree rule ``test_benchmark_prepare.py`` applies to the Huesken subset.
RECORDS_CSV_SHA256 = "2fb449362985c52c89ef02d8c7606da641b8cb314d86dac1d7f0cba74611a9c1"

#: Per panel: the ``manifest.json`` datasets it selects, and ``(in_source, kept, incompatible)`` as
#: ``prepare`` reports them at paired length 19. Ichihara is the sum of its two datasets
#: (2,431 + 419); Shmushkovich keeps nothing at any length, being the asymmetric panel #109 excludes.
_EXPECTED: dict[str, tuple[tuple[str, ...], tuple[int, int, int]]] = {
    "ichihara": (("ichihara_2007_1", "ichihara_2007_2"), (2850, 2850, 0)),
    "martinelli": (("martinelli_2023_1",), (907, 907, 0)),
    "shmushkovich": (("shmushkovich_2018_1",), (356, 0, 356)),
}

#: The two target-identity values #109 may honestly emit for these rows. ``prepare.py`` writes
#: ``unavailable`` today and will write ``synthetic_context_local`` once it passes
#: ``derive_observation``'s answer through (that value is already what
#: ``panels.derive_observation`` returns for all 4,113 rows). Both are asserted as acceptable so this
#: file pins the load-bearing half of the claim -- that neither is ever a *native* or *confirmed*
#: identity -- without pinning which of the two honest values the writer currently chooses.
_HONEST_IDENTITIES = frozenset({"unavailable", "synthetic_context_local"})


def _prepare_vendored(panel_id: str, tmp_path: Path, *, paired_length: int = 19) -> Path:
    """Prepare one vendored panel through the real ``prepare_artifact``, returning its artifact dir.

    No ``--panel-csv`` and no descriptor substitution: the panel is named by its registry id, the
    table comes off its own ``vendored_csv``, and ``prepare`` applies the row selector, so the artifact
    this returns is the one the CLI writes for the same invocation. ``paired_length`` is always passed
    because ``shmushkovich`` is asymmetric and declares none.
    """
    prepare_artifact(
        panel_id=panel_id,
        paired_length=paired_length,
        out_dir=tmp_path / "artifacts",
        invoked_command=("sirnaforge", "benchmark", "prepare", "--panel", panel_id),
    )
    return tmp_path / "artifacts" / f"{panel_id}__len{paired_length}"


def _vendored_dataset_counts() -> dict[str, int]:
    """``manifest.json``'s own per-dataset row counts -- the authority these tests reconcile against."""
    return dict(json.loads(VENDORED_MANIFEST.read_text())["datasets"])


@pytest.mark.unit
def test_the_vendored_table_is_the_bytes_every_count_here_was_measured_on() -> None:
    """A regenerated ``records.csv`` must break loudly here, not silently move a tally.

    Every number in this file was measured by running ``prepare``/``design`` over these exact bytes;
    a checksum mismatch means the numbers below describe a file this repository no longer holds.
    """
    assert hashlib.sha256(RECORDS_CSV.read_bytes()).hexdigest() == RECORDS_CSV_SHA256
    assert _vendored_dataset_counts() == {
        "ichihara_2007_1": 2431,
        "ichihara_2007_2": 419,
        "martinelli_2023_1": 907,
        "shmushkovich_2018_1": 356,
    }


@pytest.mark.unit
@pytest.mark.parametrize("panel_id", sorted(_EXPECTED))
def test_prepare_carries_every_selected_vendored_row_into_the_artifact(panel_id: str, tmp_path: Path) -> None:
    """Per-panel prepare counts, reconciled against the vendored manifest's own dataset totals.

    Two independent statements in one place, because a disagreement between them is the failure mode
    this pass exists to catch: ``observations_in_source`` must equal the sum of the row counts
    ``manifest.json`` declares for the datasets this panel's selector names, and the kept/incompatible
    split must be the one the declared architecture implies -- everything for the two paired-core
    panels, nothing for the asymmetric one. ``design_inputs.fasta`` carries exactly the kept rows,
    so an incompatible row can never reach the designer.
    """
    datasets, (in_source, kept, incompatible) = _EXPECTED[panel_id]
    declared = _vendored_dataset_counts()

    artifact_dir = _prepare_vendored(panel_id, tmp_path)
    manifest = read_manifest(artifact_dir / "manifest.json")
    counts = manifest.counts

    assert counts.observations_in_source == sum(declared[dataset] for dataset in datasets) == in_source
    assert (counts.observations_kept, counts.observations_incompatible) == (kept, incompatible)
    assert counts.observations_kept + counts.observations_incompatible == counts.observations_in_source
    assert manifest.outputs.design_inputs_fasta.sequences == kept
    assert len(read_observations(artifact_dir / "observations.csv")) == in_source


@pytest.mark.unit
@pytest.mark.parametrize("panel_id", sorted(_EXPECTED))
def test_no_prepared_row_claims_native_or_confirmed_target_evidence(panel_id: str, tmp_path: Path) -> None:
    """A site at 1-based 71 in a fabricated context must never be written as native evidence.

    The whole reason ``TargetIdentityStatus`` grew a third value (#109) rather than reusing
    ``panel_local``: these rows do carry an exact coordinate, and it indexes ``70 A +
    revcomp(guide) + 70 A``, not a transcript. #110 owns native identity, so no row an artifact
    written here holds may read as one: no identity value outside the two honest ones, no transcript
    id, no strand, and ``mapped_native`` zero in the manifest's own tally. The word "native" is not
    banned from the file's bytes -- the asymmetric rows' ``compatibility_reason`` says "#110 owns
    native asymmetric handling", which is a pointer to the issue that owns the work, not a claim of
    evidence -- but "confirmed" is, since #109's vocabulary has no such value to write.
    """
    artifact_dir = _prepare_vendored(panel_id, tmp_path)
    observations = read_observations(artifact_dir / "observations.csv")

    assert {observation.target_identity_status for observation in observations} <= _HONEST_IDENTITIES
    assert all(observation.target_transcript_id is None for observation in observations)
    assert all(observation.target_strand is None for observation in observations)
    assert read_manifest(artifact_dir / "manifest.json").counts.mapped_native == 0
    assert "confirmed" not in (artifact_dir / "observations.csv").read_text()


@pytest.mark.unit
@pytest.mark.parametrize(
    ("panel_id", "dataset", "column", "spot_checks"),
    [
        ("martinelli", "martinelli_2023_1", "efficacy_higher_is_better", (0, 453, 906)),
        # 74.00336740698765 is the last Ichihara row: a full-precision label, and the one most likely
        # to be silently rounded by a float round-trip through the CSV writer.
        ("ichihara", "ichihara_2007_1", "efficacy_higher_is_better", (0, 1215, 2430)),
        ("shmushkovich", "shmushkovich_2018_1", "label_processed", (0, 178, 355)),
    ],
)
def test_prepared_measured_values_are_the_rows_own_label_string_verbatim(
    panel_id: str,
    dataset: str,
    column: str,
    spot_checks: tuple[int, ...],
    tmp_path: Path,
) -> None:
    """``measured_value`` must survive prepare as the string ``records.csv`` holds, not a re-format.

    ``derive_observation`` reading the right column is already pinned in
    ``test_benchmark_panels.py``; what is pinned here is the rest of the trip -- ``float()`` and back
    out through ``write_observations`` -- because a value that arrives rounded, rescaled or
    direction-flipped is indistinguishable in the artifact from a value that was measured that way.
    Every selected row is checked against its own source row, in order, and three are additionally
    checked against ``label_processed`` (the label the panel published) so the two columns cannot
    drift apart unnoticed for Ichihara/Martinelli, where they are byte-identical, or be confused with
    the derived ``100 - label`` flip that only Shmushkovich's
    ``efficacy_higher_is_better`` carries.
    """
    with RECORDS_CSV.open(newline="") as handle:
        source = [row for row in csv.DictReader(handle) if row["dataset"] == dataset]

    artifact_dir = _prepare_vendored(panel_id, tmp_path)
    with (artifact_dir / "observations.csv").open(newline="") as handle:
        written = list(csv.DictReader(handle))
    # `write_observations` sorts by the zero-padded `observation_id`, so written order is source
    # order. Ichihara's panel spans two datasets and `records.csv` holds ichihara_2007_1's 2,431 rows
    # first, so its rows are this panel's leading block; every other panel is one whole dataset.
    written = written[: len(source)]
    assert len(written) == len(source)

    assert [row["measured_value"] for row in written] == [row[column] for row in source]
    for index in spot_checks:
        assert written[index]["measured_value"] == source[index]["label_processed"]
        assert written[index]["full_guide_sequence"] == source[index]["guide_sequence"].upper()


@pytest.mark.unit
def test_designing_the_martinelli_panel_enters_all_907_rows_and_agrees_with_itself(tmp_path: Path) -> None:
    """The design-stage tally for a real vendored panel, at the shipped default policy.

    Every row's measured site reproduces its own 19 nt core (``guide_match=paired_core_exact`` on all
    907), so ``no_candidate`` is 0 and the two verdict sets must be *identical* -- an unwidened
    benchmark run resolves the same thresholds as the default one, so any gap between
    ``default_pass`` and ``benchmark_pass`` here would mean the two policies were not resolved from
    the same profile. 683 of 907 pass; the polynucleotide gate #109 requires to stay active accounts
    for 179 of the exclusions under *both* policies, which is what "kept active and explicitly
    recorded" reduces to in a tally.
    """
    artifact_dir = _prepare_vendored("martinelli", tmp_path)
    result = design_artifact(artifact_dir=artifact_dir, invoked_command=("sirnaforge", "benchmark", "design"))

    assert (result.entered_design, result.no_candidate) == (907, 0)
    assert result.default_pass == result.benchmark_pass == 683
    rows = read_accounting(artifact_dir / "accounting.csv")
    assert len(rows) == 907
    assert {row.guide_match for row in rows} == {"paired_core_exact"}

    manifest = read_manifest(artifact_dir / "manifest.json")
    assert manifest.polynucleotide_run_requirement.evaluated is True
    assert manifest.polynucleotide_run_requirement.excluded.default == 179
    assert manifest.polynucleotide_run_requirement.excluded.benchmark == 179
    assert manifest.gc_widening.gc_min.widened is False
    assert manifest.gc_widening.gc_max.widened is False


@pytest.mark.unit
def test_all_356_shmushkovich_rows_reach_the_accounting_and_none_reaches_the_designer(tmp_path: Path) -> None:
    """The asymmetric panel is recorded and refused, never sliced -- end to end, not just in ``panels.py``.

    #109 requires every measured observation to appear in the artifact even when it cannot be
    designed, so the honest result for this panel is 356 accounting rows that all say
    ``entered_design=False``/``not_evaluated`` and an empty ``design_inputs.fasta`` -- not 356 missing
    rows, and not 356 rows quietly truncated to a 19 nt core. ``design_artifact`` must also survive a
    zero-record FASTA rather than raising out of the designer.
    """
    artifact_dir = _prepare_vendored("shmushkovich", tmp_path)
    result = design_artifact(artifact_dir=artifact_dir, invoked_command=("sirnaforge", "benchmark", "design"))

    assert (result.entered_design, result.no_candidate) == (0, 356)
    assert result.default_pass == result.benchmark_pass == 0
    rows = read_accounting(artifact_dir / "accounting.csv")
    assert len(rows) == 356
    assert all(row.entered_design is False for row in rows)
    assert {row.default_filter_status.value for row in rows} == {"not_evaluated"}
    assert {row.benchmark_filter_status.value for row in rows} == {"not_evaluated"}
    assert (artifact_dir / "design_inputs.fasta").read_text() == ""

    observations = read_observations(artifact_dir / "observations.csv")
    assert {observation.compatibility_status for observation in observations} == {"incompatible"}
    assert all(observation.paired_guide_sequence == observation.full_guide_sequence for observation in observations)


@pytest.mark.unit
def test_a_widened_gc_run_leaves_the_default_verdict_set_byte_identical(tmp_path: Path) -> None:
    """Widening the benchmark's GC bounds must move the benchmark tally and nothing else (#109).

    The claim the whole two-verdict-sets design rests on: the default set is re-derived from the
    widened run's own observed values, so widening may only *add* candidates the default policy then
    still judges by unchanged thresholds. Checked on the same real artifact, designed twice -- once
    unwidened, once at 10/90 -- comparing the serialised ``default_run_policy`` block byte for byte
    and the default tally row by row. 28 more rows pass under the widened benchmark policy (683 ->
    711) and the 40 rows the default GC ceiling excludes stay excluded under it, which is the
    asymmetry the design predicts: the benchmark's own ceiling exclusions drop to 0.
    """
    artifact_dir = _prepare_vendored("martinelli", tmp_path)
    argv = ("sirnaforge", "benchmark", "design")

    unwidened = design_artifact(artifact_dir=artifact_dir, invoked_command=argv)
    before = json.loads((artifact_dir / "manifest.json").read_text())
    widened = design_artifact(artifact_dir=artifact_dir, gc_min=10.0, gc_max=90.0, invoked_command=argv)
    after = json.loads((artifact_dir / "manifest.json").read_text())

    assert json.dumps(before["default_run_policy"], sort_keys=True) == json.dumps(
        after["default_run_policy"], sort_keys=True
    )
    assert unwidened.default_pass == widened.default_pass == 683
    assert unwidened.benchmark_pass == 683
    assert widened.benchmark_pass == 711
    assert after["gc_widening"]["gc_min"] == {"default": 30.0, "benchmark": 10.0, "source": "explicit", "widened": True}
    assert after["gc_widening"]["gc_max"] == {"default": 60.0, "benchmark": 90.0, "source": "explicit", "widened": True}

    exclusions = after["counts"]["excluded_by_filter"]
    assert exclusions["gc_content_max"] == {"default": 40, "benchmark": 0}
    assert exclusions["max_poly_runs"] == {"default": 179, "benchmark": 179}
