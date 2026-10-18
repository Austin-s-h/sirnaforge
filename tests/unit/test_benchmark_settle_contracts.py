"""Cross-slice contracts the #109 settle pass fixed, pinned against the one vendored panel.

Everything here is an *integration* assertion between slices that landed separately
(``benchmark/artifact.py``, ``panels.py``, ``prepare.py``, ``design.py``) and could not be checked
inside any one of them. Three defects the settle pass found are pinned so they cannot come back:

1. ``ManifestOutputEntry.path`` recorded ``str(<absolute path>)``, so two runs of the same panel
   bytes under two ``--out-dir`` roots wrote manifests that could not be compared -- which is the
   whole content of #109's "the manifest is sufficient to reproduce the artifact from checksummed
   inputs" criterion. It now records the fixed inner filename.
2. ``design_artifact`` *appended* a ``prepared_artifact`` input entry whose checksum was of
   ``manifest.json`` -- a file it overwrites moments later. Re-running ``benchmark design`` therefore
   grew ``inputs`` by one stale entry per run. It now records the two files it actually reads, keyed
   by role and replaced rather than appended.
3. ``FilterScope``'s ``species``/``hit_classes`` are ``frozenset[str]`` and dumped in hash order,
   which varies with ``PYTHONHASHSEED``. Two identical runs wrote manifests differing only in the
   order of ``hit_classes``, defeating the same comparison as (1).

The panel exercised is ``huesken_subset`` -- ``tests/unit/data/sirna_efficacy_subset.csv``, the only
benchmark panel whose bytes this repository vendors (see ``tests/unit/data/README.md``). Ichihara,
Martinelli, Shmushkovich and OligoGym are named in #109/#110 and are **not present**; nothing in this
file is evidence about them.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from sirnaforge.benchmark.artifact import (
    ACCOUNTING_FILENAME,
    CANDIDATES_ALL_FILENAME,
    DESIGN_INPUTS_FASTA_FILENAME,
    MANIFEST_FILENAME,
    OBSERVATIONS_FILENAME,
    read_manifest,
)
from sirnaforge.benchmark.design import design_artifact
from sirnaforge.benchmark.prepare import prepare_artifact
from sirnaforge.config.run_policy import DEFAULT_PROFILE_NAME, EntryPoint, resolve_run_policy
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.models.policy import FilterScope

#: Every file in a fully designed artifact, and whether its bytes are expected to be reproducible.
#: ``manifest.json`` is excluded because ``created_utc`` is the one field the schema deliberately
#: leaves nondeterministic (:class:`~sirnaforge.benchmark.artifact.BenchmarkArtifactManifest`).
_REPRODUCIBLE_FILES = (
    OBSERVATIONS_FILENAME,
    DESIGN_INPUTS_FASTA_FILENAME,
    CANDIDATES_ALL_FILENAME,
    ACCOUNTING_FILENAME,
)

_ARGV = ("sirnaforge", "benchmark", "prepare", "--panel", "huesken_subset")


def _prepare_and_design(out_dir: Path) -> Path:
    """One full ``prepare`` then ``design`` pass over the vendored panel, at its declared length."""
    prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=out_dir, invoked_command=_ARGV)
    artifact_dir = out_dir / "huesken_subset__len19"
    design_artifact(artifact_dir=artifact_dir, invoked_command=_ARGV)
    return artifact_dir


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.unit
def test_manifest_output_paths_are_the_fixed_inner_filenames(tmp_path: Path) -> None:
    """No output entry may carry a caller's directory, or the manifest is not comparable between runs.

    Recording ``str(artifact_dir / "observations.csv")`` made the manifest a function of ``--out-dir``
    rather than of the panel bytes and the policy, which is exactly what #109's reproducibility
    criterion needs it not to be.
    """
    artifact_dir = _prepare_and_design(tmp_path / "deeply" / "nested" / "root")
    outputs = read_manifest(artifact_dir / MANIFEST_FILENAME).outputs

    assert outputs.observations_csv.path == OBSERVATIONS_FILENAME
    assert outputs.design_inputs_fasta.path == DESIGN_INPUTS_FASTA_FILENAME
    assert outputs.candidates_all_csv is not None
    assert outputs.candidates_all_csv.path == CANDIDATES_ALL_FILENAME
    assert outputs.accounting_csv is not None
    assert outputs.accounting_csv.path == ACCOUNTING_FILENAME
    assert str(tmp_path) not in (artifact_dir / MANIFEST_FILENAME).read_text()


@pytest.mark.unit
def test_manifest_output_checksums_match_the_bytes_on_disk(tmp_path: Path) -> None:
    """A recorded checksum has to be of the file that is actually there, not of an earlier version."""
    artifact_dir = _prepare_and_design(tmp_path)
    outputs = read_manifest(artifact_dir / MANIFEST_FILENAME).outputs

    for entry in (outputs.observations_csv, outputs.design_inputs_fasta):
        assert entry.sha256 == _sha256(artifact_dir / entry.path)
    for optional in (outputs.candidates_all_csv, outputs.accounting_csv):
        assert optional is not None
        assert optional.sha256 == _sha256(artifact_dir / optional.path)


@pytest.mark.unit
def test_re_running_design_does_not_grow_the_inputs_list(tmp_path: Path) -> None:
    """``role`` is a key in ``inputs``: a second design replaces its entries, never appends them.

    The pre-fix behaviour appended a fresh ``prepared_artifact`` entry per run, so a twice-designed
    artifact carried two same-role entries with different checksums and no rule saying which one
    described the artifact as it then stood.
    """
    artifact_dir = _prepare_and_design(tmp_path)
    first = read_manifest(artifact_dir / MANIFEST_FILENAME).inputs

    design_artifact(artifact_dir=artifact_dir, invoked_command=_ARGV)
    second = read_manifest(artifact_dir / MANIFEST_FILENAME).inputs

    roles = [entry.role for entry in second]
    assert roles == [entry.role for entry in first]
    assert len(roles) == len(set(roles))
    assert second == first


@pytest.mark.unit
def test_recorded_prepared_inputs_are_the_files_design_actually_read(tmp_path: Path) -> None:
    """The design stage's provenance names ``observations.csv``/``design_inputs.fasta``, not the manifest.

    Checksumming ``manifest.json`` was worthless as provenance: ``design_artifact`` overwrites it, so
    the hash named bytes the finished artifact no longer held.
    """
    artifact_dir = _prepare_and_design(tmp_path)
    by_role = {entry.role: entry for entry in read_manifest(artifact_dir / MANIFEST_FILENAME).inputs}

    assert set(by_role) == {"panel_csv", "prepared_observations_csv", "prepared_design_inputs_fasta"}
    assert by_role["prepared_observations_csv"].sha256 == _sha256(artifact_dir / OBSERVATIONS_FILENAME)
    assert by_role["prepared_design_inputs_fasta"].sha256 == _sha256(artifact_dir / DESIGN_INPUTS_FASTA_FILENAME)


@pytest.mark.unit
def test_filter_scope_serialises_both_set_axes_sorted() -> None:
    """A ``frozenset`` dumps in hash order; a manifest field must not depend on ``PYTHONHASHSEED``."""
    scope = FilterScope(species=frozenset({"mouse", "human"}), hit_classes=frozenset({"undetermined", "off_target"}))
    payload = scope.model_dump(mode="json")

    assert payload["species"] == ["human", "mouse"]
    assert payload["hit_classes"] == ["off_target", "undetermined"]


@pytest.mark.unit
def test_a_resolved_policy_manifest_carries_sorted_hit_classes(tmp_path: Path) -> None:
    """The same guarantee where it actually bit: inside ``run_policy``'s per-filter scope blocks."""
    artifact_dir = _prepare_and_design(tmp_path)
    payload = json.loads((artifact_dir / MANIFEST_FILENAME).read_text())

    scoped = [
        entry["scope"]["hit_classes"]
        for policy in ("run_policy", "default_run_policy")
        for entry in payload[policy]["filters"]
        if entry["scope"]["hit_classes"]
    ]
    assert scoped, "no filter in either policy block declares a hit-class scope; this test proves nothing"
    assert all(classes == sorted(classes) for classes in scoped)


@pytest.mark.unit
def test_two_full_passes_differ_only_in_the_declared_nondeterministic_field(tmp_path: Path) -> None:
    """The reproducibility criterion, end to end: same bytes and same argv, same artifact.

    Byte equality for all four data files, and manifest equality once ``created_utc`` -- the single
    field the schema declares nondeterministic, and the only one that lives outside every CSV -- is
    dropped from both sides.
    """
    first = _prepare_and_design(tmp_path / "pass1")
    second = _prepare_and_design(tmp_path / "pass2")

    for filename in _REPRODUCIBLE_FILES:
        assert (first / filename).read_bytes() == (second / filename).read_bytes(), filename

    payloads = []
    for artifact_dir in (first, second):
        payload = json.loads((artifact_dir / MANIFEST_FILENAME).read_text())
        del payload["created_utc"]
        payloads.append(payload)
    assert payloads[0] == payloads[1]


@pytest.mark.unit
def test_prepared_design_inputs_fasta_is_consumable_by_the_unmodified_design_path(tmp_path: Path) -> None:
    """#109's "consumable by the current fixed-length design/workflow path", executed rather than asserted.

    ``design_inputs.fasta`` is handed to a plain ``SiRNADesigner`` resolved under
    ``EntryPoint.DESIGN_COMMAND`` -- the same object ``sirnaforge design`` builds, with no benchmark
    policy and no benchmark code in the path -- and must enumerate candidates for every record.
    ``candidates_all.csv`` is then the *superset* ``design.py``'s docstring promises: the plain path's
    survivors plus the enumeration-time rejects a fail-action gate dropped.
    """
    prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path, invoked_command=_ARGV)
    fasta = tmp_path / "huesken_subset__len19" / DESIGN_INPUTS_FASTA_FILENAME

    policy = resolve_run_policy(
        entry_point=EntryPoint.DESIGN_COMMAND, profile_name=DEFAULT_PROFILE_NAME, stated={"sirna_length": 19}
    )
    result = SiRNADesigner(policy.design_parameters).design_from_file(str(fasta))

    manifest = read_manifest(tmp_path / "huesken_subset__len19" / MANIFEST_FILENAME)
    assert result.total_sequences == manifest.counts.observations_kept
    assert result.candidates, "the prepared FASTA produced no candidates through the unmodified design path"
    assert all(len(candidate.guide_sequence) == 19 for candidate in result.candidates)
