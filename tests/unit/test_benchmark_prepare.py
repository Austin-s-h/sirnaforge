"""``benchmark prepare``: panel ingest into a checksummed artifact (#109, bm-prepare slice).

``sirnaforge.benchmark.prepare.prepare_artifact`` is exercised here against ``huesken_subset``
(``tests/unit/data/sirna_efficacy_subset.csv`` -- see ``tests/unit/data/README.md``) plus the three
synthetic ``--panel-csv`` fixtures in ``tests/unit/data/benchmark/``, each read under the
architecture-level ``user_supplied_*`` descriptor that names its geometry:
``user_supplied_fully_complementary``, ``user_supplied_asymmetric`` and
``user_supplied_paired_core_with_overhang`` (the same architecture as Huesken). Reusing those
fixtures rather than inventing new ones means this file exercises the exact compatibility branches
``tests/unit/data/benchmark/README.md`` documents (too-short guide, strand-length mismatch,
asymmetric geometry) without adding a second copy of any of them.

These tests used to hand those same fixtures to ``martinelli``/``shmushkovich``/``ichihara``, which
was legal while those ids were ``data_present=False``. ``f4beab7`` vendored their bytes, so
``prepare`` now -- correctly -- refuses a ``--panel-csv`` for them: a run must not be able to read
different bytes than the ones its own manifest names as vendored
(``test_vendored_panel_refuses_an_explicit_panel_csv_override``). The ``user_supplied_*`` ids are
where a caller's own table belongs, and they vendor nothing, so the fixtures reach the same
compatibility branches with no claim attached to a published panel's name.

The issue's cited source (``docs/prd_benchmark_artifacts_and_variable_length.md``) does not exist in
this repository; the issue body (#109) is the entire specification. A passing test on a
``user_supplied_*`` id is evidence about the artifact/prepare *contract* for one declared duplex
geometry, and about nothing else: a synthetic fixture's ``efficacy`` value is an arbitrary
placeholder, not a measurement. The vendored panels' own rows are carried through
``prepare``/``design`` in ``tests/unit/test_benchmark_real_panels.py``.
"""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pytest
from Bio.Seq import Seq

from sirnaforge.benchmark.artifact import BenchmarkArtifactError, read_observations
from sirnaforge.benchmark.panels import PANEL_REGISTRY, predeclared_split
from sirnaforge.benchmark.prepare import BenchmarkPrepareError, prepare_artifact

DATA_DIR = Path(__file__).parent / "data"
HUESKEN_CSV = DATA_DIR / "sirna_efficacy_subset.csv"
HUESKEN_TRANSCRIPTS = DATA_DIR / "sirna_efficacy_subset_transcripts.fa"
BENCHMARK_FIXTURES = DATA_DIR / "benchmark"

# Recorded in tests/unit/data/README.md's checksum table (issue #102). Pinned here, independently
# of that file, so a byte-changed fixture is caught by two files agreeing to disagree.
HUESKEN_CSV_SHA256 = "6e0a2efea505a1df7110136a8818e5bbbb28243dc12becdfbc73402fa8a5e10d"


@pytest.mark.unit
def test_prepare_huesken_subset_writes_the_three_declared_files(tmp_path: Path) -> None:
    """A prepare at the panel's declared length (19) writes exactly the three #109 files."""
    manifest = prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path)

    artifact_dir = tmp_path / "huesken_subset__len19"
    assert {path.name for path in artifact_dir.iterdir()} == {
        "manifest.json",
        "observations.csv",
        "design_inputs.fasta",
    }
    assert manifest.panel.panel_id == "huesken_subset"
    assert manifest.paired_length == 19


@pytest.mark.unit
def test_manifest_input_sha256_matches_the_readme_checksum(tmp_path: Path) -> None:
    """``inputs[0]`` is the panel csv, and its checksum must match the one pinned in the vendored-data README.

    That README (``tests/unit/data/README.md``, issue #102) is the authority for this fixture's
    bytes; a mismatch here means either this fixture or that table changed without the other.
    """
    manifest = prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path)

    assert manifest.inputs[0].role == "panel_csv"
    assert manifest.inputs[0].sha256 == HUESKEN_CSV_SHA256
    assert manifest.inputs[0].sha256 == hashlib.sha256(HUESKEN_CSV.read_bytes()).hexdigest()


@pytest.mark.unit
def test_counts_partition_the_source_into_kept_and_incompatible(tmp_path: Path) -> None:
    """``observations_in_source`` always equals ``observations_kept + observations_incompatible``.

    Checked at two lengths: 19 (every 21 nt Huesken guide is compatible, so incompatible is 0 -- a
    trivially-true equation is not a real check) and 22 (every guide is too short for a 22 nt paired
    core, since the panel's guides are uniformly 21 nt, so incompatible equals the full count).
    """
    manifest_19 = prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path / "a")
    counts_19 = manifest_19.counts
    assert counts_19.observations_in_source == counts_19.observations_kept + counts_19.observations_incompatible
    assert counts_19.observations_incompatible == 0
    assert counts_19.observations_kept == counts_19.observations_in_source == 180

    manifest_22 = prepare_artifact(panel_id="huesken_subset", paired_length=22, out_dir=tmp_path / "b")
    counts_22 = manifest_22.counts
    assert counts_22.observations_in_source == counts_22.observations_kept + counts_22.observations_incompatible
    assert counts_22.observations_kept == 0
    assert counts_22.observations_incompatible == counts_22.observations_in_source == 180


@pytest.mark.unit
def test_paired_length_24_raises_naming_19_23(tmp_path: Path) -> None:
    """#109 is fixed-length 19-23 only; a length outside that range is refused before any file is touched."""
    with pytest.raises(BenchmarkArtifactError, match="19") as excinfo:
        prepare_artifact(panel_id="huesken_subset", paired_length=24, out_dir=tmp_path)
    assert "23" in str(excinfo.value)
    assert not any(tmp_path.iterdir())


@pytest.mark.unit
def test_second_prepare_into_the_same_directory_refuses_without_overwrite(tmp_path: Path) -> None:
    """Re-running prepare into an existing artifact directory must be refused, not silently replace it."""
    prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path)

    with pytest.raises(BenchmarkPrepareError, match="already exists"):
        prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path)

    # overwrite=True is the escape hatch, and it must actually succeed.
    prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path, overwrite=True)


@pytest.mark.unit
def test_two_prepares_produce_byte_identical_observations_csv(tmp_path: Path) -> None:
    """Two independent prepares of the same checksummed inputs must write identical bytes.

    No column anywhere carries a timestamp or an absolute path, and rows are always written sorted
    by ``observation_id`` (``sirnaforge.benchmark.artifact.write_observations``) -- this is the
    regression assertion the manifest's reproducibility claim rests on, not merely a claim.
    """
    prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path / "run1")
    prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path / "run2")

    observations_1 = (tmp_path / "run1" / "huesken_subset__len19" / "observations.csv").read_bytes()
    observations_2 = (tmp_path / "run2" / "huesken_subset__len19" / "observations.csv").read_bytes()
    assert observations_1 == observations_2

    fasta_1 = (tmp_path / "run1" / "huesken_subset__len19" / "design_inputs.fasta").read_bytes()
    fasta_2 = (tmp_path / "run2" / "huesken_subset__len19" / "design_inputs.fasta").read_bytes()
    assert fasta_1 == fasta_2


@pytest.mark.unit
def test_unknown_panel_id_refuses_naming_the_registry(tmp_path: Path) -> None:
    """An unregistered panel_id is refused before anything is opened, naming valid ids."""
    with pytest.raises(ValueError, match="unknown benchmark panel_id") as excinfo:
        prepare_artifact(panel_id="not_a_real_panel", out_dir=tmp_path)
    assert "huesken_subset" in str(excinfo.value)


@pytest.mark.unit
@pytest.mark.parametrize(
    "panel_id",
    [
        "huesken_full",
        "user_supplied_paired_core_with_overhang",
        "user_supplied_fully_complementary",
        "user_supplied_asymmetric",
    ],
)
def test_registered_panel_without_vendored_bytes_requires_panel_csv(panel_id: str, tmp_path: Path) -> None:
    """Every registered-but-unvendored panel refuses with a reason, naming the option that fixes it.

    All four ids this repository declares no bytes for: ``huesken_full``, whose table
    (``work/sirna_bench.csv``) is untracked here, and the three architecture-level ``user_supplied_*``
    ids, which ship no bytes by construction. ``paired_length`` is passed for all of them because
    ``user_supplied_asymmetric`` declares none, and its missing-length refusal is a different one
    (:func:`~sirnaforge.benchmark.prepare._resolve_paired_length`) that would mask this test's.

    The message must also name that descriptor's own required columns. For a ``user_supplied_*`` id
    "which columns" is the caller's next question -- the table is theirs and this repository has never
    seen it -- and the descriptor is the only thing that can answer it.
    """
    descriptor = PANEL_REGISTRY[panel_id]
    assert descriptor.data_present is False
    with pytest.raises(BenchmarkPrepareError, match="no vendored bytes") as excinfo:
        prepare_artifact(panel_id=panel_id, paired_length=19, out_dir=tmp_path)
    assert "--panel-csv" in str(excinfo.value)
    for column in descriptor.required_source_columns():
        assert column in str(excinfo.value)
    assert not any(tmp_path.iterdir())


@pytest.mark.unit
@pytest.mark.parametrize("panel_id", ["huesken_subset", "ichihara", "martinelli", "shmushkovich"])
def test_vendored_panel_refuses_an_explicit_panel_csv_override(panel_id: str, tmp_path: Path) -> None:
    """A panel that ships its own bytes: a caller may not silently substitute different ones.

    Matches the CLI's own help text for ``--panel-csv`` (``cli.py::benchmark_prepare``): a run must
    not be able to read bytes other than the ones its manifest names as vendored. Checked on all four
    ingestible vendored ids, and with the panel's *own* table as the override for ``huesken_subset``
    -- even identical bytes are refused, because the guard is on provenance, not on content. This is
    also why the synthetic fixtures below are prepared under ``user_supplied_*`` ids: since ``f4beab7``
    vendored the OligoGym bytes, handing one of those panels a crafted table is exactly what this
    refusal exists to stop.
    """
    with pytest.raises(BenchmarkPrepareError, match="vendored bytes"):
        prepare_artifact(panel_id=panel_id, paired_length=19, panel_csv=HUESKEN_CSV, out_dir=tmp_path)
    assert not any(tmp_path.iterdir())


@pytest.mark.unit
def test_user_supplied_fully_complementary_via_panel_csv(tmp_path: Path) -> None:
    """A fully-complementary table: guide and passenger pair blunt over their full, equal length.

    ``synthetic_fully_complementary.csv`` (``tests/unit/data/benchmark/README.md``) carries three
    compatible 21 nt records and one (``FC0004``) whose passenger is 20 nt against a 21 nt guide --
    the "declared fully complementary but the strands differ in length" incompatible branch. No
    ``--paired-length`` is passed, so this also pins the descriptor's declared default of 21 (the
    canonical blunt duplex) through to the artifact directory's name.
    """
    manifest = prepare_artifact(
        panel_id="user_supplied_fully_complementary",
        panel_csv=BENCHMARK_FIXTURES / "synthetic_fully_complementary.csv",
        out_dir=tmp_path,
    )
    assert manifest.panel.data_present is False
    assert manifest.paired_length == 21
    assert manifest.counts.observations_in_source == 4
    assert manifest.counts.observations_kept == 3
    assert manifest.counts.observations_incompatible == 1

    observations = read_observations(tmp_path / "user_supplied_fully_complementary__len21" / "observations.csv")
    by_row = {observation.source_row_index: observation for observation in observations}
    compatible = by_row[0]
    assert compatible.compatibility_status == "compatible"
    assert compatible.guide_3p_overhang == ""  # blunt: the whole guide is the paired core
    assert compatible.duplex_pairing_status == "fully_complementary"
    # A caller's table locates nothing: this repository knows nothing about their target context, not
    # even whether one exists, so `unavailable` is the only honest value it can record for these rows.
    assert compatible.target_identity_status == "unavailable"
    assert compatible.target_transcript_id is None

    incompatible = by_row[3]
    assert incompatible.compatibility_status == "incompatible"
    assert "differ in length" in incompatible.compatibility_reason
    assert incompatible.paired_guide_sequence == incompatible.full_guide_sequence  # never truncated to fit


@pytest.mark.unit
def test_user_supplied_asymmetric_is_incompatible_by_construction(tmp_path: Path) -> None:
    """An asymmetric (15/20) duplex is incompatible with #109's fixed-length interface unconditionally.

    ``synthetic_asymmetric_15_20.csv`` records are kept in ``observations.csv`` -- #109 never drops
    a measured row -- but every one is incompatible, with both strand lengths named in its reason, and
    the guide is reported verbatim rather than truncated or padded to force it into a 19-23 nt paired
    core (#110 owns giving this a real path). ``paired_length`` is mandatory here, because an
    asymmetric descriptor declares none and inventing one would tag the artifact with a length nothing
    measured; it only names the directory, and cannot make a single row compatible.
    """
    manifest = prepare_artifact(
        panel_id="user_supplied_asymmetric",
        panel_csv=BENCHMARK_FIXTURES / "synthetic_asymmetric_15_20.csv",
        paired_length=19,
        out_dir=tmp_path,
    )
    assert manifest.counts.observations_kept == 0
    assert manifest.counts.observations_incompatible == manifest.counts.observations_in_source == 2
    assert manifest.outputs.design_inputs_fasta.sequences == 0

    observations = read_observations(tmp_path / "user_supplied_asymmetric__len19" / "observations.csv")
    assert len(observations) == 2
    for observation in observations:
        assert observation.compatibility_status == "incompatible"
        assert "asymmetric" in observation.compatibility_reason
        assert "20 nt" in observation.compatibility_reason and "15 nt" in observation.compatibility_reason
        assert observation.paired_guide_sequence == observation.full_guide_sequence
        assert observation.guide_3p_overhang is None
        assert observation.target_identity_status == "unavailable"


@pytest.mark.unit
def test_paired_core_with_overhang_too_short_guide_is_incompatible(tmp_path: Path) -> None:
    """A guide shorter than the requested paired core is incompatible, never stretched to fit.

    ``synthetic_paired_core_with_overhang.csv``'s ``SYN0005`` is an 18 nt guide, deliberately shorter
    than any #109 fixed length; requested at 19 it must be incompatible while its four 21 nt siblings
    (a measured 2 nt 3' overhang, the same shape as the vendored Huesken subset) are compatible. That
    fixture states no passenger at all, which is the ``user_supplied_*`` optional-column case: the
    mapped ``passenger_sequence`` column is simply absent, and every row's field is ``None`` rather
    than an error.
    """
    manifest = prepare_artifact(
        panel_id="user_supplied_paired_core_with_overhang",
        panel_csv=BENCHMARK_FIXTURES / "synthetic_paired_core_with_overhang.csv",
        paired_length=19,
        out_dir=tmp_path,
    )
    assert manifest.counts.observations_kept == 4
    assert manifest.counts.observations_incompatible == 1

    observations = read_observations(tmp_path / "user_supplied_paired_core_with_overhang__len19" / "observations.csv")
    by_row = {observation.source_row_index: observation for observation in observations}
    too_short = by_row[4]  # SYN0005, the fifth row
    assert too_short.full_guide_sequence == "ACGTACGTACGTACGTAC"
    assert too_short.compatibility_status == "incompatible"
    assert "18" in too_short.compatibility_reason and "19" in too_short.compatibility_reason

    compatible = by_row[0]
    assert compatible.compatibility_status == "compatible"
    assert compatible.guide_3p_overhang == "TT"  # last 2 nt of a 21 nt guide sliced to a 19 nt core
    assert all(observation.passenger_sequence is None for observation in observations)


@pytest.mark.unit
def test_incompatible_observation_is_kept_but_never_enters_design_inputs_fasta(tmp_path: Path) -> None:
    """#109's "never drop a measured observation" rule, checked against the two files it spans.

    The incompatible row must still have a design_context_id (the schema requires one on every row),
    but no FASTA record may exist for it, since ``design_inputs.fasta`` is one record per *compatible*
    observation only.
    """
    manifest = prepare_artifact(
        panel_id="user_supplied_paired_core_with_overhang",
        panel_csv=BENCHMARK_FIXTURES / "synthetic_paired_core_with_overhang.csv",
        paired_length=19,
        out_dir=tmp_path,
    )
    artifact_dir = tmp_path / "user_supplied_paired_core_with_overhang__len19"
    observations = read_observations(artifact_dir / "observations.csv")
    too_short = next(observation for observation in observations if observation.source_row_index == 4)
    assert too_short.design_context_id  # non-empty, per the schema

    fasta_records = (artifact_dir / "design_inputs.fasta").read_text()
    assert too_short.design_context_id not in fasta_records
    assert manifest.outputs.design_inputs_fasta.sequences == 4


@pytest.mark.unit
def test_split_matches_the_predeclared_accession_parity_rule(tmp_path: Path) -> None:
    """Every huesken_subset observation's ``split`` matches ``predeclared_split`` on its raw accession.

    Cross-checks ``observations.csv`` against the raw source table plus the shared rule
    (``sirnaforge.benchmark.panels.predeclared_split``, itself pinned against
    ``tests/unit/data/README.md`` in ``test_benchmark_panels.py``) rather than re-asserting the
    rule's own accession lists, so this test fails if -- and only if -- ``prepare`` stops applying it.
    """
    with HUESKEN_CSV.open(newline="") as handle:
        raw_rows = list(csv.DictReader(handle))

    manifest = prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path)
    assert manifest.split_rule == "sha256_accession_parity_v1"

    observations = read_observations(tmp_path / "huesken_subset__len19" / "observations.csv")
    by_row = {observation.source_row_index: observation for observation in observations}
    for index, raw_row in enumerate(raw_rows):
        expected = predeclared_split(raw_row["accession"])
        assert by_row[index].split == expected


@pytest.mark.unit
def test_design_context_defaults_to_the_measured_target_site(tmp_path: Path) -> None:
    """Without ``--panel-transcripts``, every compatible observation's design context is the measured target site.

    The reverse complement of the *full* measured guide, never trimmed to the requested paired
    length: the target site is a fact about the guide as measured, independent of which fixed length
    this particular run happens to request.
    """
    manifest = prepare_artifact(panel_id="huesken_subset", paired_length=19, out_dir=tmp_path)
    artifact_dir = tmp_path / "huesken_subset__len19"
    observations = {
        observation.source_row_index: observation
        for observation in read_observations(artifact_dir / "observations.csv")
    }
    first = observations[0]
    assert first.design_context_source == "measured_target_site"

    fasta_records = _read_fasta(artifact_dir / "design_inputs.fasta")
    expected = str(Seq(first.full_guide_sequence).reverse_complement())
    assert fasta_records[first.design_context_id] == expected
    assert manifest.outputs.design_inputs_fasta.sequences == 180


@pytest.mark.unit
def test_panel_transcripts_supplies_the_real_transcript_when_the_accession_matches(tmp_path: Path) -> None:
    """``--panel-transcripts`` flips ``design_context_source`` to ``panel_transcript`` for a matched accession.

    ``tests/unit/data/sirna_efficacy_subset_transcripts.fa`` (issue #95) carries the full-length
    transcript for every accession in this subset, headered by accession -- exactly the lookup key
    ``prepare`` uses.
    """
    manifest = prepare_artifact(
        panel_id="huesken_subset",
        paired_length=19,
        panel_transcripts=HUESKEN_TRANSCRIPTS,
        out_dir=tmp_path,
    )
    assert any(entry.role == "panel_transcripts_fasta" for entry in manifest.inputs)

    artifact_dir = tmp_path / "huesken_subset__len19"
    observations = {
        observation.source_row_index: observation
        for observation in read_observations(artifact_dir / "observations.csv")
    }
    first = observations[0]  # accession NM_003969, per tests/unit/data/sirna_efficacy_subset.csv row 0
    assert first.design_context_source == "panel_transcript"

    fasta_records = _read_fasta(artifact_dir / "design_inputs.fasta")
    transcripts = _read_fasta(HUESKEN_TRANSCRIPTS)
    assert fasta_records[first.design_context_id] == transcripts["NM_003969"]


@pytest.mark.unit
def test_a_matched_user_supplied_transcript_supplies_a_context_but_not_an_identity(tmp_path: Path) -> None:
    """``--panel-transcripts`` may give a caller's row its design context, never a target identity.

    A matched accession flips ``design_context_source`` to ``panel_transcript`` -- the caller asked for
    their own sequence to be designed over, and it is -- while ``target_identity_status`` stays
    ``unavailable`` and ``target_transcript_id`` stays unset. Matching a header in a FASTA the caller
    also supplied is not evidence about a native transcript, and keeping those two apart is the point:
    #110 owns native mapping. Both inputs are written here rather than vendored because
    ``synthetic_context.fa``'s headers are context ids, not the accessions this lookup keys on.
    """
    table = tmp_path / "caller_table.csv"
    table.write_text("accession,guide_sequence,efficacy\nCALLER_TX_1,ACGTACGTACGTACGTACGTT,0.61\n")
    fasta = tmp_path / "caller_transcripts.fa"
    fasta.write_text(">CALLER_TX_1\nTTGGCCAAGGTTCCACGTACGTACGTACGTACGTTAAGGCCTTAAGGTTCC\n")

    prepare_artifact(
        panel_id="user_supplied_paired_core_with_overhang",
        paired_length=19,
        panel_csv=table,
        panel_transcripts=fasta,
        out_dir=tmp_path / "out",
    )
    artifact_dir = tmp_path / "out" / "user_supplied_paired_core_with_overhang__len19"
    observation = read_observations(artifact_dir / "observations.csv")[0]

    assert observation.design_context_source == "panel_transcript"
    designed = _read_fasta(artifact_dir / "design_inputs.fasta")[observation.design_context_id]
    assert designed == _read_fasta(fasta)["CALLER_TX_1"]
    assert observation.target_identity_status == "unavailable"
    assert observation.target_transcript_id is None
    assert observation.target_start_1based is None and observation.target_end_1based is None


def _read_fasta(path: Path) -> dict[str, str]:
    """A tiny local FASTA parser, so this test file does not need a Biopython round trip to assert on."""
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records
