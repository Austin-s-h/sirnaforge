"""One resolved screening reference, one door, and the rename that made both possible (#99).

Every test here fails against 0.7.0's two doors: ``genome_indices_override`` wrote the pipeline
parameter directly and skipped ``_configure_transcriptome_inputs``, where the transcript index the
classifier reads was built -- so a cDNA file handed to the "genome" door aligned fine and then
classified against nothing (211,359 candidate-rows, 100% of mouse classifications unevidenced).
"""

from __future__ import annotations

import asyncio
import gzip
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from sirnaforge.config.reference_policy import (
    UNRESOLVED_SPECIES,
    ReferenceKind,
    ReferenceKindError,
    ReferencePolicyResolver,
    ReferenceState,
    SpeciesAuthority,
    WorkflowInputSpec,
    parse_index_entries,
)
from sirnaforge.data.transcriptome_manager import TranscriptomeManager
from sirnaforge.models.sirna import DesignMode, DesignParameters
from sirnaforge.workflow import INDEX_BUILD_ERROR_KEY, SiRNAWorkflow, WorkflowConfig

WORKFLOWS = Path("src/sirnaforge/pipeline/nextflow/workflows")

MOUSE_CDNA_HEADER = (
    ">ENSMUST00000108658.11 cdna chromosome:GRCm39:11:69471109:69482701:1 "
    "gene:ENSMUSG00000059552.15 gene_symbol:Trp53 transcript_biotype:protein_coding"
)


def _mouse_cdna(tmp_path: Path, name: str = "mm_cdna.fa") -> Path:
    """A one-transcript mouse cDNA file whose headers name its own assembly."""
    path = tmp_path / name
    path.write_text(f"{MOUSE_CDNA_HEADER}\nACGTACGTACGTACGTACGTA\n")
    return path


def _workflow(
    tmp_path: Path,
    name: str,
    *,
    transcriptome_indices: str | None = None,
    transcriptome_fasta: str | None = None,
    screen_species: list[str] | None = None,
) -> SiRNAWorkflow:
    """A workflow whose only screening reference is the one the test names."""
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        design_params=DesignParameters(),
        screen_species=screen_species,
        transcriptome_indices=transcriptome_indices,
        transcriptome_fasta=transcriptome_fasta,
    )
    return SiRNAWorkflow(config)


def _returning(fasta: Path) -> Callable[..., dict[str, Path]]:
    """A ``get_custom_transcriptome`` stand-in: the fetch and cache layer is not under test here."""

    def _custom(_self: TranscriptomeManager, _reference: str, **_kwargs: object) -> dict[str, Path]:
        return {"fasta": fasta}

    return _custom


def _resolve(workflow: SiRNAWorkflow) -> tuple[bool, dict[str, Any]]:
    """Run the resolver, returning ``(configured, pipeline parameters)``."""
    params: dict[str, Any] = {}
    configured = asyncio.run(workflow._resolve_screening_references(params))
    return configured, params


# ---------------------------------------------------------------------------
# One door: an override and a default differ only in provenance
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_an_explicit_index_override_builds_the_classifier_metadata_too(tmp_path: Path) -> None:
    """The override used to skip the only place a transcript index was ever built."""
    fasta = _mouse_cdna(tmp_path)
    workflow = _workflow(tmp_path, "override_out", transcriptome_indices=f"mouse:{fasta}")

    configured, params = _resolve(workflow)

    assert configured is True
    assert params["transcriptome_indices"] == f"mouse:{fasta}"
    assert params["transcriptome_species"] == "mouse"
    # The metadata defect 1 is about: without it every mouse hit publishes species_index_missing.
    species_index = workflow._transcript_index.for_species("mouse")
    assert species_index is not None and species_index.transcript_count == 1
    # And the cDNA repeat detection reuses, which the override also used to skip.
    assert workflow._species_cdna_fasta["mouse"] == fasta


@pytest.mark.unit
def test_an_override_and_a_default_differ_only_in_provenance(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Same reference through both doors: same species, kind, identity and index; different state."""
    fasta = _mouse_cdna(tmp_path)
    monkeypatch.setattr(TranscriptomeManager, "get_custom_transcriptome", _returning(fasta))

    override = _workflow(tmp_path, "prov_override", transcriptome_indices=f"mouse:{fasta}")
    default = _workflow(tmp_path, "prov_default", transcriptome_fasta=str(fasta))
    _resolve(override)
    _resolve(default)

    through_override = override._screening_references.references[0]
    through_default = default._screening_references.references[0]

    assert (through_override.species, through_override.kind, through_override.index) == (
        "mouse",
        ReferenceKind.TRANSCRIPTOME,
        str(fasta),
    )
    assert (through_override.species, through_override.kind, through_override.index) == (
        through_default.species,
        through_default.kind,
        through_default.index,
    )
    assert through_override.state is ReferenceState.EXPLICIT
    assert through_override.form.value == "prebuilt_index"
    assert through_default.form.value == "reference"
    assert override._transcript_index.for_species("mouse") is not None
    assert default._transcript_index.for_species("mouse") is not None


@pytest.mark.unit
def test_a_reference_with_no_index_yet_travels_on_the_fasta_parameter(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A FASTA named as an index prefix aligns nothing, and the pipeline calls that a completed screen.

    The pipeline indexes ``transcriptome_fastas`` and treats ``transcriptome_indices`` as prefixes, so
    a reference whose index does not exist yet has to travel on the first.
    """
    fasta = _mouse_cdna(tmp_path)
    monkeypatch.setattr(TranscriptomeManager, "get_custom_transcriptome", _returning(fasta))
    workflow = _workflow(tmp_path, "unindexed_out", transcriptome_fasta=str(fasta))

    configured, params = _resolve(workflow)

    assert configured is True
    assert params["transcriptome_fastas"] == f"mouse:{fasta}"
    assert "transcriptome_indices" not in params, "a FASTA must not be named where a prefix is read"
    assert workflow._screening_references.references[0].needs_index_build is True


@pytest.mark.unit
def test_an_index_prefix_with_no_readable_sequence_is_refused_not_screened(tmp_path: Path) -> None:
    """Aligning against an index whose sequence cannot be read can only publish undetermined hits.

    Refused on the same grounds as a failed index build: the species is unscreened, with a reason,
    rather than screened for hours against something nothing can classify.
    """
    workflow = _workflow(tmp_path, "no_companion_out", transcriptome_indices="mouse:/nonexistent/mm_cdna")

    configured, params = _resolve(workflow)

    assert configured is False
    assert "transcriptome_indices" not in params, "an unusable reference must not reach the aligner"
    reason = workflow._species_screening_shortfalls["mouse"]
    assert "resolved to" in reason and "--transcriptome-fasta" in reason, reason
    rejection = workflow._screening_references.rejections[0]
    assert rejection.species == "mouse"
    assert workflow._screening_references.references == ()


# ---------------------------------------------------------------------------
# Species comes from the resolved reference, and a declared species is intent
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_gzipped_or_binary_neighbour_of_an_index_prefix_is_refused_not_read(tmp_path: Path) -> None:
    """Header inference used to abort the whole run on such a file, and then to build an empty index.

    Both outcomes are wrong: an unreadable neighbour means the hits cannot be classified, which is a
    per-species completeness fact, not a crash and not a screen.
    """
    compressed = tmp_path / "mm_cdna.fa.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write(f"{MOUSE_CDNA_HEADER}\nACGTACGTACGTACGTACGTA\n")
    binary = tmp_path / "mm_index"
    binary.write_bytes(bytes(range(256)) * 4)

    for name, prefix in (("gz_out", compressed), ("binary_out", binary)):
        workflow = _workflow(tmp_path, name, transcriptome_indices=f"mouse:{prefix}")

        configured, params = _resolve(workflow)

        assert configured is False, prefix
        assert params == {}
        assert "gzipped or binary" in workflow._species_screening_shortfalls["mouse"]
        assert workflow._transcript_index.for_species("mouse") is None, "an empty index classifies nothing"


@pytest.mark.unit
def test_an_index_build_failure_names_the_species_of_the_reference_it_failed_for(tmp_path: Path) -> None:
    """Wave 1's per-species rejection is only useful if it names the species that went unscreened.

    The species has to be resolved before the refusal: keyed on the declaration alone, an undeclared
    reference recorded its shortfall against ``unknown`` and the real species vanished from the
    published record.
    """
    workflow = _workflow(tmp_path, "build_failed_out", transcriptome_fasta="ensembl_mouse_cdna")
    fasta = _mouse_cdna(tmp_path)

    async def _prepared(*_args: object, **_kwargs: object) -> dict[str, object]:
        return {
            "source_species": "mouse",
            "fasta": fasta,
            INDEX_BUILD_ERROR_KEY: "BWA-MEM2 index build failed for mm_cdna.fa; ran out of memory",
        }

    workflow._prepare_transcriptome_database = _prepared  # type: ignore[method-assign]

    configured, _ = _resolve(workflow)

    assert configured is False
    assert "ran out of memory" in workflow._species_screening_shortfalls["mouse"]
    assert UNRESOLVED_SPECIES not in workflow._species_screening_shortfalls
    assert workflow._screening_references.rejections[0].species == "mouse"


@pytest.mark.unit
@pytest.mark.parametrize("declared", ["transcriptome", "moose"])
def test_a_species_the_registry_does_not_recognise_is_refused_on_the_index_door(declared: str) -> None:
    """The index door has no header fallback for its label, so an unrecognised one becomes the species.

    ``transcriptome`` as a species is the defect this module exists to close, and ``--species``
    already refuses a name the registry does not know.
    """
    with pytest.raises(ValueError, match="unsupported species"):
        parse_index_entries(f"{declared}:/idx/prefix", option="--transcriptome-indices", reason="test")


@pytest.mark.unit
def test_a_declared_species_wins_over_the_headers_and_the_disagreement_is_reported(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Header inference is the fallback, not the authority: a caller who names a species means it."""
    fasta = _mouse_cdna(tmp_path)
    workflow = _workflow(tmp_path, "declared_out", transcriptome_indices=f"rat:{fasta}")

    with caplog.at_level("WARNING"):
        _resolve(workflow)

    reference = workflow._screening_references.references[0]
    assert reference.species == "rat"
    assert reference.species_authority is SpeciesAuthority.DECLARED
    assert any("disagree" in record.message for record in caplog.records), "a contradiction must be reported"


@pytest.mark.unit
def test_a_species_declared_on_a_custom_path_is_honoured(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``--transcriptome-fasta rat:/path`` states the species; the path itself stays the reference."""
    fasta = _mouse_cdna(tmp_path)
    seen: list[str] = []

    def _custom(_self: TranscriptomeManager, reference: str, **_kwargs: object) -> dict[str, Path]:
        seen.append(reference)
        return {"fasta": fasta}

    monkeypatch.setattr(TranscriptomeManager, "get_custom_transcriptome", _custom)
    workflow = _workflow(tmp_path, "declared_path_out", transcriptome_fasta=f"rat:{fasta}")

    _resolve(workflow)

    assert seen == [str(fasta)], "the species prefix must not be passed on as part of the path"
    reference = workflow._screening_references.references[0]
    assert (reference.species, reference.species_authority) == ("rat", SpeciesAuthority.DECLARED)


@pytest.mark.unit
def test_an_undeclared_custom_reference_takes_its_species_from_its_headers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With nothing declared, the reference's own assembly token decides -- never the parameter name."""
    fasta = _mouse_cdna(tmp_path)
    monkeypatch.setattr(TranscriptomeManager, "get_custom_transcriptome", _returning(fasta))
    workflow = _workflow(tmp_path, "headers_out", transcriptome_fasta=str(fasta))

    _, params = _resolve(workflow)

    reference = workflow._screening_references.references[0]
    assert (reference.species, reference.species_authority) == ("mouse", SpeciesAuthority.REFERENCE_HEADERS)
    assert params["transcriptome_species"] == "mouse"


@pytest.mark.unit
def test_a_reference_whose_species_nothing_states_is_labelled_unknown(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The placeholder is ``unknown``, not ``transcriptome``: that word is the kind, not a species."""
    fasta = tmp_path / "mystery.fa"
    fasta.write_text(">contig_1 no assembly here\nACGTACGTACGTACGTACGTA\n")
    monkeypatch.setattr(TranscriptomeManager, "get_custom_transcriptome", _returning(fasta))
    workflow = _workflow(tmp_path, "unknown_out", transcriptome_fasta=str(fasta))

    _resolve(workflow)

    reference = workflow._screening_references.references[0]
    assert reference.species == UNRESOLVED_SPECIES
    assert reference.species_authority is SpeciesAuthority.UNRESOLVED


# ---------------------------------------------------------------------------
# kind is keyed to modality, and a mismatch fails before expensive work
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_a_genomic_assembly_named_for_an_sirna_screen_fails_before_any_work(tmp_path: Path) -> None:
    """``ensembl_human_hg38_primary`` is a whole genome: refuse it before downloading one."""
    output_dir = tmp_path / "kind_out"

    with pytest.raises(ReferenceKindError, match="genome reference"):
        WorkflowConfig(
            output_dir=output_dir,
            gene_query="TP53",
            design_params=DesignParameters(),
            transcriptome_fasta="ensembl_human_hg38_primary",
        )

    assert not output_dir.exists(), "validation must precede even the output tree"


@pytest.mark.unit
def test_a_zfn_run_resolves_no_transcriptome_reference(tmp_path: Path) -> None:
    """ZFN screens genomic DNA, so the cDNA defaults are not its references -- and are not fatal.

    A default nobody asked for is disabled with a reason; naming one explicitly is a mismatch.
    """
    resolved_defaults = ReferencePolicyResolver(WorkflowInputSpec()).resolve_transcriptomes()
    config = WorkflowConfig(
        output_dir=tmp_path / "zfn_out",
        gene_query="TP53",
        design_params=DesignParameters(design_mode=DesignMode.ZFN),
        transcriptome_selection=resolved_defaults,
    )

    assert config.screening_kind is ReferenceKind.GENOME
    assert config.screening_requests == ()
    assert len(config.screening_request_rejections) == len(resolved_defaults.choices)
    assert "screens against a genome" in config.screening_request_rejections[0].reason

    with pytest.raises(ReferenceKindError, match="screens against a genome"):
        WorkflowConfig(
            output_dir=tmp_path / "zfn_explicit_out",
            gene_query="TP53",
            design_params=DesignParameters(design_mode=DesignMode.ZFN),
            transcriptome_fasta="ensembl_human_cdna",
        )


@pytest.mark.unit
def test_a_zfn_run_refuses_a_prebuilt_transcriptome_index(tmp_path: Path) -> None:
    """The kind check has to cover both doors, or the override is the one that skips validation.

    A prebuilt index arrives on ``--transcriptome-indices``, so it is a transcriptome index by
    construction; naming one on a ZFN run used to resolve a genome-kind reference onto the
    transcriptome parameters, with a plan entry contradicting its own kind.
    """
    output_dir = tmp_path / "zfn_index_out"

    with pytest.raises(ReferenceKindError, match="screens against a genome"):
        WorkflowConfig(
            output_dir=output_dir,
            gene_query="TP53",
            design_params=DesignParameters(design_mode=DesignMode.ZFN),
            transcriptome_indices=f"mouse:{_mouse_cdna(tmp_path)}",
        )

    assert not output_dir.exists(), "validation must precede even the output tree"


# ---------------------------------------------------------------------------
# The rename: nothing silently maps an old name onto a new one
# ---------------------------------------------------------------------------


@pytest.mark.unit
@pytest.mark.parametrize(
    ("argument", "replacement"),
    [("genome_species", "screen_species"), ("genome_indices_override", "transcriptome_indices")],
)
def test_a_removed_argument_names_its_replacement(tmp_path: Path, argument: str, replacement: str) -> None:
    """A hard break, but never a bare "unexpected keyword argument"."""
    with pytest.raises(TypeError, match=f"renamed to '{replacement}'"):
        WorkflowConfig(
            output_dir=tmp_path / "renamed_out",
            gene_query="TP53",
            design_params=DesignParameters(),
            **{argument: ["human"] if argument == "genome_species" else "human:/idx"},
        )


@pytest.mark.unit
def test_an_unknown_argument_is_still_reported_as_unknown(tmp_path: Path) -> None:
    """A typo must not read as a rename."""
    with pytest.raises(TypeError, match="unexpected keyword argument 'gnome_species'"):
        WorkflowConfig(
            output_dir=tmp_path / "typo_out",
            gene_query="TP53",
            design_params=DesignParameters(),
            gnome_species=["human"],
        )


@pytest.mark.unit
def test_a_stale_pipeline_parameter_is_refused_rather_than_ignored(tmp_path: Path) -> None:
    """Nextflow ignores an unknown --param, so a stale one would screen against nothing."""
    with pytest.raises(ValueError, match="'genome_indices' is now 'transcriptome_indices'"):
        WorkflowConfig(
            output_dir=tmp_path / "stale_out",
            gene_query="TP53",
            design_params=DesignParameters(),
            nextflow_config={"genome_indices": "human:/idx"},
        )


@pytest.mark.unit
def test_the_pipeline_refuses_the_removed_parameters_by_name() -> None:
    """The same refusal on the Nextflow side, where an unknown --param is otherwise accepted."""
    main_nf = (WORKFLOWS / "main.nf").read_text()

    for old_name, new_name in (
        ("genome_fastas", "transcriptome_fastas"),
        ("genome_indices", "transcriptome_indices"),
        ("genome_species", "transcriptome_species"),
    ):
        assert f"{old_name} : '{new_name}'" in main_nf or f"{old_name}: '{new_name}'" in main_nf, old_name
        assert f"params.{old_name}" not in main_nf, f"{old_name} must not still configure anything"
    # Every supplied key is folded before it is compared, because Nextflow files a hyphenated
    # --genome-indices under the camelCase key genomeIndices: matching the three snake_case
    # spellings alone let that form through, and it then screened nothing and reported success.
    assert "params.keySet()" in main_nf
    assert "replaceAll('-', '_')" in main_nf
    assert "([a-z0-9])([A-Z])" in main_nf


@pytest.mark.unit
def test_the_per_species_results_are_published_under_transcriptome() -> None:
    """The publishDir and the workflow's own fallback reader have to agree on one directory name."""
    module = (WORKFLOWS / "modules/local/offtarget_analysis.nf").read_text()

    assert 'publishDir "${params.outdir}/transcriptome"' in module
    assert "/genome" not in module


# ---------------------------------------------------------------------------
# The #104 types this issue wires: an explicit species scope and a screening plan
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_the_published_summary_carries_the_scope_and_the_plan(tmp_path: Path) -> None:
    """``FilterScope`` and ``ScreeningPlanEntry`` landed unwired; the resolver is what fills them."""
    fasta = _mouse_cdna(tmp_path)
    workflow = _workflow(tmp_path, "plan_out", transcriptome_indices=f"mouse:{fasta}")
    _resolve(workflow)
    guides = tmp_path / "candidates.fasta"
    guides.write_text(">cand_1\nACGTACGTACGTACGTACGTA\n")

    workflow._record_screening_plan(guides, {"max_hits": 100, "bwa_k": 12})
    summary = workflow._summarize_screening_references()

    assert summary["scope"]["species"] == ["mouse"]
    assert summary["screening"]["resolved_species"] == ["mouse"]
    entry = summary["screening_plan"]["entries"][0]
    assert entry["channel"] == "transcriptome"
    assert entry["species"] == "mouse"
    assert entry["reference_id"] == str(fasta)
    assert entry["search_settings"] == {"max_hits": 100, "bwa_k": 12}
    assert len(entry["guide_set_digest"]) == 16, "the digest identifies the guide set actually submitted"


@pytest.mark.unit
def test_the_species_parameter_reaches_nextflow_exactly_once(tmp_path: Path) -> None:
    """The resolver and the runner both name the species; two flags leave the last one silently winning."""
    from sirnaforge.pipeline.nextflow.config import NextflowConfig  # noqa: PLC0415

    input_file = tmp_path / "candidates.fasta"
    input_file.write_text(">cand_1\nACGTACGTACGTACGTACGTA\n")

    args = NextflowConfig(work_dir=tmp_path / "work").get_nextflow_args(
        input_file=input_file,
        output_dir=tmp_path / "out",
        screen_species=["human", "mouse"],
        additional_params={"transcriptome_species": "human"},
    )

    assert args.count("--transcriptome_species") == 1
    assert args[args.index("--transcriptome_species") + 1] == "human"
