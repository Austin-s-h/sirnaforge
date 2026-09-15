"""The manifest must name the data and the build the science rests on, or say it cannot.

``manifest.json`` already records the resolved run policy, every gate as data and a sha256 for three
of four outputs. What it does not record is what the run was *screened against* and what *built* it:
no Ensembl release, no genome build, no transcriptome checksum, no aligner index digest, no git SHA,
no ViennaRNA version -- and no digest at all for ``report.html``, the 31 MB artifact everyone reads.
Two runs therefore cannot be shown to have screened against the same references.

These tests pin the one new top-level ``provenance`` key and the single rule that makes it worth
having: **the manifest must never assert an identity it cannot verify.** ``"unknown"`` is not an
answer; an explicit, named non-availability is. This is #100's rule for screening evidence applied
to references and to the build.

Contracts this suite assumes of the implementation (it is written first, so these are the spec):

``sirnaforge.provenance``
    ``build_identity() -> dict`` returns the whole ``build`` sub-block, wrapped in ``lru_cache`` so
    it exposes ``build_identity.cache_clear()``. Two seams exist for tests, because both hazards they
    guard are unobservable otherwise: ``_git(*args, cwd=...)``, the single thin wrapper around the
    git subprocess, and ``_package_repo_root()``, repo discovery from the package source.
    ``write_report_manifest(report_html, manifest_json, out) -> Path`` writes the sidecar.

``ScreeningReference.identity_evidence``
    An optional frozen view of the reference's ``CacheMetadata``, keyed
    ``url``/``checksum``/``checksum_algorithm``/``file_size``/``size_scope``/``downloaded_at``/
    ``filters``/``local_content_hash``. It is the only route by which the bytes' identity survives
    ``_prepare_transcriptome_database``, where a fresh ``TranscriptomeManager`` is built and
    garbage-collected on return.

All fixtures are synthetic. TP53 is the documented public example gene.
"""

from __future__ import annotations

import asyncio
import hashlib
import importlib
import json
import re
import subprocess
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

import pytest

from sirnaforge import __version__
from sirnaforge.config.reference_policy import (
    ReferenceForm,
    ReferenceKind,
    ReferenceRejection,
    ReferenceState,
    ScreeningReference,
    ScreeningReferenceSet,
    SpeciesAuthority,
)
from sirnaforge.config.run_policy import EntryPoint, resolve_run_policy
from sirnaforge.core.scoring import COMPOSITE_TERMS
from sirnaforge.data.orthology import OrthologueMapping
from sirnaforge.models.scoring_profile import TERM_REGISTRY, ExperimentalAUPostScreenSiRNAWeights
from sirnaforge.models.sirna import (
    DesignParameters,
    DesignResult,
    DesignWeights,
    PostScreenMiRNAWeights,
    ScoringWeights,
    SiRNACandidate,
)
from sirnaforge.reporting import build_payload
from sirnaforge.utils.cache_utils import (
    ARTIFACT_TRANSCRIPTOME_INDEX,
    artifact_stamp_path,
    fingerprint_inputs,
    read_artifact_stamp,
    write_artifact_stamp,
)
from sirnaforge.workflow import SiRNAWorkflow, WorkflowConfig

GUIDE = "ACGTACGTACGTACGTACGTA"
REPO_ROOT = Path(__file__).resolve().parents[2]

#: Keys a resolved screening reference must carry. Present-and-null is the contract; absent is not,
#: because an absent key reads as "not applicable" where a null triple reads as "not knowable, and
#: here is why".
_REQUIRED_REFERENCE_KEYS = frozenset(
    {
        "species",
        "species_authority",
        "requiredness",
        "source_url",
        "assembly",
        "ensembl_release",
        "bytes",
        "index",
        "classification_index",
        "filters",
    }
)


# ---------------------------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------------------------


def _provenance_module() -> Any:
    """Import ``sirnaforge.provenance`` at call time, so its absence fails one test, not collection."""
    return importlib.import_module("sirnaforge.provenance")


@pytest.fixture(autouse=True)
def _clear_build_identity_cache() -> Iterator[None]:
    """``build_identity`` caches a git subprocess for the process, so env-sensitive tests must reset it."""
    try:
        module = _provenance_module()
    except ImportError:
        yield
        return
    module.build_identity.cache_clear()
    yield
    module.build_identity.cache_clear()


def _candidate(candidate_id: str = "c1") -> SiRNACandidate:
    """A passing candidate with enough component scores that step6 can re-score and write a CSV."""
    return SiRNACandidate(
        id=candidate_id,
        transcript_id="ENST00000000001",
        position=1,
        guide_sequence=GUIDE,
        passenger_sequence=GUIDE.translate(str.maketrans("ATCG", "TAGC"))[::-1],
        length=len(GUIDE),
        gc_content=47.6,
        asymmetry_score=0.5,
        mfe=-4.2,
        structure="." * len(GUIDE),
        duplex_stability=-39.0,
        paired_fraction=0.1,
        design_score=50.0,
        component_scores={
            "target_accessibility": 0.5,
            "asymmetry": 0.5,
            "gc_content": 0.5,
            "duplex_stability_score": 0.6,
            "melting_temp_c": 70.0,
            "dg_5p": -8.0,
            "dg_3p": -9.0,
            "delta_dg_end": 1.0,
        },
    )


def _qualified_workflow(tmp_path: Path, name: str) -> SiRNAWorkflow:
    """A qualified human run, resolved through the one resolver as every real entry point does."""
    policy = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        query_species="human",
        screen_species=["human"],
    )
    config = WorkflowConfig(
        output_dir=tmp_path / name,
        gene_query="TP53",
        resolved_policy=policy,
        screen_species=["human"],
        query_species="human",
    )
    workflow = SiRNAWorkflow(config)
    workflow._gene_transcript_ids = {"ENST00000000001"}
    workflow._query_gene_symbols = {"TP53"}
    return workflow


def _design_only_workflow(tmp_path: Path, name: str) -> SiRNAWorkflow:
    """A design-only run: no screen, no orthology lookup, no Nextflow."""
    policy = resolve_run_policy(entry_point=EntryPoint.DESIGN_COMMAND)
    config = WorkflowConfig(output_dir=tmp_path / name, gene_query="TP53", resolved_policy=policy)
    return SiRNAWorkflow(config)


def _cached_reference(
    tmp_path: Path,
    *,
    species: str = "human",
    with_index_stamp: bool = False,
) -> tuple[ScreeningReference, dict[str, str]]:
    """A resolved reference whose bytes exist on disk, plus the md5 of those bytes.

    The FASTA is real so a digest can actually be computed and compared; the recorded sizes are the
    container's, because that is the shape the run under review produced.
    """
    cache = tmp_path / "cache" / "transcriptomes"
    cache.mkdir(parents=True, exist_ok=True)
    fasta = cache / "152cb327e903.fa"
    fasta.write_text(">ENST00000000001.1 gene_symbol:TP53\nACGTACGTACGTACGTACGTAAA\n")
    index_prefix = cache / "152cb327e903_index"

    inputs = fingerprint_inputs([fasta])
    checksum = next(iter(inputs.values()))

    if with_index_stamp:
        members = []
        for suffix in (".0123", ".bwt.2bit.64"):
            member = index_prefix.with_name(index_prefix.name + suffix)
            member.write_bytes(b"\x00" * 2048)
            members.append(member)
        write_artifact_stamp(
            ARTIFACT_TRANSCRIPTOME_INDEX,
            index_prefix,
            inputs=inputs,
            outputs=members,
        )
        # Every bwa-mem2 member of a real human index exceeds the full-digest budget, so at least one
        # member in the wild is sampled. Force that here: a sampled digest must never be republished
        # without its mode, or it reads as a whole-file one.
        stamp_path = artifact_stamp_path(index_prefix)
        stamp = json.loads(stamp_path.read_text())
        stamp["outputs"][members[1].name]["digest_mode"] = "sampled"
        stamp_path.write_text(json.dumps(stamp, indent=2))

    identity_evidence = {
        "url": "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "checksum": checksum,
        "checksum_algorithm": "md5",
        "file_size": fasta.stat().st_size,
        "size_scope": "decompressed_fasta",
        "downloaded_at": "2026-09-09T16:55:17.703607",
        "filters": None,
        "local_content_hash": checksum,
    }
    kwargs: dict[str, Any] = {
        "species": species,
        "kind": ReferenceKind.TRANSCRIPTOME,
        "identity": f"ensembl_{species}_cdna",
        "index": str(index_prefix),
        "species_authority": SpeciesAuthority.SOURCE_REGISTRY,
        "form": ReferenceForm.REFERENCE,
        "state": ReferenceState.DEFAULT,
        "reason": "auto-selected for the query species",
        "fasta": str(fasta),
        "needs_index_build": False,
        "identity_evidence": identity_evidence,
    }
    try:
        reference = ScreeningReference(**kwargs)
    except TypeError as exc:  # pragma: no cover - the pre-implementation state
        raise AssertionError(
            "ScreeningReference must accept an optional `identity_evidence` mapping so the bytes' "
            "url/checksum/size/downloaded_at survive `_prepare_transcriptome_database`; without it "
            f"the manifest can never name which transcriptome was screened ({exc})"
        ) from exc

    return reference, {"checksum": checksum, "fasta": str(fasta), "index": str(index_prefix)}


def _with_resolved_human(workflow: SiRNAWorkflow, reference: ScreeningReference) -> SiRNAWorkflow:
    workflow._screening_references = ScreeningReferenceSet(
        kind=ReferenceKind.TRANSCRIPTOME,
        references=(reference,),
        requested_species=("human",),
    )
    workflow._species_cdna_fasta = {reference.species: Path(str(reference.fasta))}
    return workflow


def _with_rejected_human(workflow: SiRNAWorkflow) -> SiRNAWorkflow:
    workflow._screening_references = ScreeningReferenceSet(
        kind=ReferenceKind.TRANSCRIPTOME,
        references=(),
        rejections=(
            ReferenceRejection(
                species="human",
                identity="ensembl_human_cdna",
                reason="BWA-MEM2 index build failed: MemoryError",
            ),
        ),
        requested_species=("human",),
    )
    return workflow


def _manifest(workflow: SiRNAWorkflow) -> dict[str, Any]:
    """Build the manifest and round-trip it through JSON, as the written artifact is."""
    base = workflow.config.output_dir / "sirnaforge"
    base.mkdir(parents=True, exist_ok=True)
    manifest = workflow._build_fair_manifest(
        all_csv=base / "candidates_all.csv",
        pass_csv=base / "candidates_pass.csv",
        pass_fasta=base / "candidates_pass.fasta",
    )
    written: dict[str, Any] = json.loads(json.dumps(manifest))
    return written


def _provenance(workflow: SiRNAWorkflow) -> dict[str, Any]:
    """The provenance block of a JSON round-tripped manifest, with its schema version asserted."""
    manifest = _manifest(workflow)
    assert "provenance" in manifest, "the manifest must carry a top-level `provenance` block"
    block: dict[str, Any] = manifest["provenance"]
    assert block.get("schema_version"), "the provenance block must name its own schema version"
    return block


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _leaves(node: Any, path: str = "provenance") -> Iterator[tuple[str, str | None, Any, Any]]:
    """Yield ``(path, key, value, parent)`` for every scalar in the block."""
    if isinstance(node, Mapping):
        for key, value in node.items():
            child = f"{path}.{key}"
            if isinstance(value, (Mapping, list)):
                yield from _leaves(value, child)
            else:
                yield child, key, value, node
    elif isinstance(node, list):
        for index, value in enumerate(node):
            child = f"{path}[{index}]"
            if isinstance(value, (Mapping, list)):
                yield from _leaves(value, child)
            else:
                yield child, None, value, node


def _mappings(node: Any, path: str = "provenance") -> Iterator[tuple[str, Mapping[str, Any]]]:
    if isinstance(node, Mapping):
        yield path, node
        for key, value in node.items():
            yield from _mappings(value, f"{path}.{key}")
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from _mappings(value, f"{path}[{index}]")


def _run_step6(workflow: SiRNAWorkflow, candidates: list[SiRNACandidate] | None = None) -> Path:
    """Run design ranking then step6, exactly as a real run does, and return ``sirnaforge/``."""
    guides = candidates or [_candidate()]
    result = DesignResult(
        input_file="<test>",
        parameters=workflow.config.design_params,
        candidates=list(guides),
        top_candidates=list(guides),
        total_sequences=1,
        total_candidates=len(guides),
        filtered_candidates=len(guides),
        processing_time=0.0,
    )
    workflow._apply_post_screen_ranking(result)
    asyncio.run(workflow.step6_generate_reports(result))
    return workflow.config.output_dir / "sirnaforge"


# ---------------------------------------------------------------------------------------------
# Invariant 1: a qualified run cannot emit an unnamed required reference
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_a_qualified_run_names_every_reference_it_required(tmp_path):
    """Every required channel/species pair must resolve to a named, identified reference.

    This is the headline gap: the report scopes its off-target liability rows to "human, nm<=2" and
    never says *which* human transcriptome, so two runs cannot be shown to have screened the same
    bytes.
    """
    workflow = _qualified_workflow(tmp_path, "qualified")
    reference, expected = _cached_reference(tmp_path)
    _with_resolved_human(workflow, reference)

    required = workflow.config.resolved_policy.evidence_requirements.required_pairs
    assert ("transcriptome", "human") in required, "fixture guard: this run must require human transcriptome"

    provenance = _provenance(workflow)
    screening = provenance["references"]["screening"]
    resolved = {entry["species"]: entry for entry in screening["resolved"]}

    for channel, species in sorted(required):
        assert channel == "transcriptome"
        assert species in resolved, f"required pair ({channel}, {species}) has no resolved reference in the manifest"
        entry = resolved[species]
        missing = _REQUIRED_REFERENCE_KEYS - set(entry)
        assert not missing, f"resolved reference for {species} omits {sorted(missing)}"
        assert entry["requiredness"] == "required"
        assert entry["species_authority"] == "source_registry"

    human = resolved["human"]
    assert human["source_url"]["value"] == (
        "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    )
    assert human["assembly"]["value"] == "GRCh38"
    # The release genuinely floated: /pub/current_fasta pins none, and the manifest must say that
    # rather than guess or omit.
    assert human["ensembl_release"]["value"] is None
    assert human["ensembl_release"]["pinned"] is False
    # `algorithm` is mandatory: this digest is md5 while `files.*` is sha256. `size_scope` is
    # mandatory: the recorded size is decompressed and the remote Content-Length is not.
    assert human["bytes"]["digest"]["algorithm"] == "md5"
    assert human["bytes"]["digest"]["value"] == expected["checksum"]
    assert human["bytes"]["size_scope"] == "decompressed_fasta"
    # Explicit null, never an absent key: null means the whole reference was screened.
    assert human["filters"] is None

    assert provenance["coverage"]["unnamed_required_references"] == []
    assert provenance["coverage"]["required_pairs"] == [["transcriptome", "human"]]
    assert provenance["coverage"]["unidentified_required_references"] == []


@pytest.mark.unit
def test_a_required_reference_that_did_not_resolve_is_declared_not_hidden(tmp_path):
    """A shortfall is a fact about the run; erasing it is how "unscreened" came to read as "clean"."""
    workflow = _with_rejected_human(_qualified_workflow(tmp_path, "rejected"))

    provenance = _provenance(workflow)
    screening = provenance["references"]["screening"]

    # The block is never omitted: absent reads "a screen happened and we forgot".
    assert screening["resolved"] == []
    assert screening["requested_species"] == ["human"]
    assert screening["rejected"] == [
        {
            "species": "human",
            "identity": "ensembl_human_cdna",
            "reason": "BWA-MEM2 index build failed: MemoryError",
        }
    ]

    coverage = provenance["coverage"]
    assert "unnamed_required_references" in coverage, "the shortfall key must be present, not omitted"
    assert coverage["unnamed_required_references"] == [["transcriptome", "human"]]


@pytest.mark.unit
@pytest.mark.parametrize("shape", ["resolved", "rejected", "design_only"])
def test_the_provenance_block_fabricates_nothing(tmp_path, shape):
    """The recursive no-fabrication walk: absent evidence must not read like present evidence.

    One test, applied to every leaf, is what stops the *next* field from regressing to silence --
    which is how ``design_parameters.genome_index: null`` and a bare ``tool_version: "0.7.1"``
    survived this long.
    """
    workflow = _qualified_workflow(tmp_path, f"walk_{shape}")
    if shape == "resolved":
        reference, _ = _cached_reference(tmp_path)
        _with_resolved_human(workflow, reference)
    elif shape == "rejected":
        _with_rejected_human(workflow)
    else:
        workflow = _design_only_workflow(tmp_path, f"walk_{shape}")

    provenance = _provenance(workflow)
    leaves = list(_leaves(provenance))
    assert len(leaves) > 30, "a provenance block this small cannot be recording build and reference identity"

    for path, key, value, parent in leaves:
        if not isinstance(value, str):
            continue
        # 1. The bare string "unknown" is forbidden, except the pre-existing UNRESOLVED_SPECIES
        #    literal, which may only appear paired with species_authority: "unresolved".
        if value == "unknown":
            assert key == "species" and parent.get("species_authority") == "unresolved", (
                f'{path} is the bare string "unknown"; name the reason-class instead'
            )
        # 2. A synthesised fingerprint is never a commit.
        assert not value.startswith("nogit-"), f"{path} republishes runner.py's synthesised nogit fingerprint"
        # 3. An mtime-derived fingerprint may live only under a key that says so.
        if value.startswith("mtime-"):
            assert key == "workflow_dir_fingerprint", (
                f"{path} carries an mtime fingerprint under a key that reads like a real revision"
            )

    for path, mapping in _mappings(provenance):
        if "value" not in mapping or "state" not in mapping:
            continue
        if mapping["value"] is not None:
            continue
        state = mapping["state"]
        reason = mapping.get("reason")
        assert isinstance(state, str) and state and state != "unknown", (
            f"{path} is a null value with no reason-class in `state`"
        )
        assert isinstance(reason, str) and len(reason) >= 20, (
            f"{path} is a null value whose `reason` does not explain the non-availability: {reason!r}"
        )


# ---------------------------------------------------------------------------------------------
# Invariant 2: build identity, and the release gate that fails closed
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_a_dev_build_declares_itself(tmp_path):
    """``tool_version: "0.7.1"`` is a release string, not a build identity.

    The run under review was a dev build and said nothing about it. A dev build must SAY it is one,
    naming which counts it fails the release gate on.
    """
    head = subprocess.run(  # noqa: S603
        ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
        timeout=10,
    ).stdout.strip()

    build = _provenance(_qualified_workflow(tmp_path, "dev_build"))["build"]

    assert build["release"] is False
    assert build["kind"] == "dev_worktree"
    assert build["tool_version"] == __version__
    vcs = build["vcs"]
    assert vcs["repo_verified"] is True
    assert vcs["commit"]["value"] == head
    assert re.fullmatch(r"[0-9a-f]{40}", vcs["commit"]["value"])
    # describe never falls back to __version__: disagreement between them is the signal.
    describe = vcs["describe"]["value"]
    assert describe and describe != __version__
    assert describe in build["detail"], "detail must quote the describe string it judged"
    assert __version__ in build["detail"], "detail must name the __version__ that matched no tag"


@pytest.mark.unit
def test_the_release_gate_fails_closed_when_git_cannot_be_run(tmp_path, monkeypatch):
    """An undecidable build is not a release build, and a failed `git status` is not cleanliness."""
    module = _provenance_module()

    def _explode(*_args: Any, **_kwargs: Any) -> str:
        raise subprocess.TimeoutExpired(cmd="git", timeout=5)

    monkeypatch.setattr(module, "_git", _explode)
    module.build_identity.cache_clear()

    build = _provenance(_qualified_workflow(tmp_path, "no_git"))["build"]

    assert build["release"] is False
    assert build["kind"] in {"installed_no_vcs", "unknown_environment"}
    assert build["detail"]
    vcs = build["vcs"]
    assert vcs["commit"]["value"] is None
    assert vcs["commit"]["state"] == "git_binary_missing"
    # NOT False: False is a positive claim that the tree was clean.
    assert vcs["dirty"]["value"] is None
    assert vcs["dirty"]["state"] == "git_binary_missing"


@pytest.mark.unit
def test_a_foreign_repository_is_never_attributed_to_this_build(tmp_path, monkeypatch):
    """Three ancestor directories of the installed package carry ``.git``, including the sandbox root.

    Walking up until a ``.git`` appears therefore attributes somebody else's commit to this build.
    The guard is ``git ls-files --error-unmatch src/sirnaforge/__init__.py``: a repo that does not
    track the package cannot vouch for it.
    """
    module = _provenance_module()
    foreign = tmp_path / "not_sirnaforge"
    foreign.mkdir()
    monkeypatch.setattr(module, "_package_repo_root", lambda: foreign)
    module.build_identity.cache_clear()

    vcs = _provenance(_qualified_workflow(tmp_path, "foreign_repo"))["build"]["vcs"]

    assert vcs["repo_verified"] is False
    assert vcs["commit"]["value"] is None, "an unverified repo must not supply a commit"
    assert vcs["commit"]["state"] == "repo_does_not_track_package"


@pytest.mark.unit
def test_a_container_run_says_so_and_claims_no_image_digest(tmp_path, monkeypatch):
    """A tag is mutable and is never substituted for a digest.

    ``0.0.0+local`` is the honest recorded value for the docker-dev image, which passes no VERSION
    build-arg -- and it must never be promoted into ``image_digest`` or ``image_ref``.
    """
    module = _provenance_module()
    monkeypatch.setenv("SIRNAFORGE_IN_CONTAINER", "1")
    monkeypatch.setenv("BUILD_VERSION", "0.0.0+local")
    module.build_identity.cache_clear()

    container = _provenance(_qualified_workflow(tmp_path, "container"))["build"]["container"]

    assert container["in_container"] is True
    assert container["detection"] == "env:SIRNAFORGE_IN_CONTAINER"
    assert container["image_digest"]["value"] is None
    assert container["image_digest"]["state"] == "not_injected"
    assert container["image_ref"]["value"] is None
    assert container["build_version"]["value"] == "0.0.0+local"
    assert container["build_version"]["state"] == "env:BUILD_VERSION"


@pytest.mark.unit
def test_the_build_block_records_the_runtime_the_report_drew_with(tmp_path):
    """The report draws secondary structures, so the ViennaRNA that drew them is part of the result.

    The aligner stays an explicit non-availability: no code path captures a bwa-mem2 version, and a
    conda pin on PATH is a constraint, not an observation.
    """
    build = _provenance(_qualified_workflow(tmp_path, "runtime"))["build"]

    runtime = build["runtime"]
    assert runtime["python"]["version"] and runtime["python"]["implementation"]
    assert runtime["viennarna"]["binary_probed"] is False, "no code execs RNAfold; do not quote its banner"
    assert runtime["aligner"]["value"] is None
    assert runtime["aligner"]["state"] == "not_captured"

    dependencies = build["dependencies"]
    assert re.fullmatch(r"[0-9a-f]{64}", dependencies["distributions_sha256"])
    assert dependencies["distribution_count"] > 1
    # Derived from importlib.metadata.requires("sirnaforge") -- the same drift class as
    # reported_not_scored. A declared-but-uninstalled dependency maps to null, never to omission.
    declared = dependencies["declared"]
    assert isinstance(declared, dict) and declared
    assert "biopython" in declared, "the declared set must come from the package metadata, not a hand-list"


# ---------------------------------------------------------------------------------------------
# Invariant 3: reported_not_scored is derived, not a literal
# ---------------------------------------------------------------------------------------------


def _expected_split(weights: ScoringWeights) -> tuple[list[str], list[str]]:
    scored = {term for vector in weights.all_vectors() for term in vector.terms}
    return (
        [term for term in TERM_REGISTRY if term in scored],
        [term for term in TERM_REGISTRY if term not in scored],
    )


def _experimental_au_workflow(tmp_path: Path, name: str) -> SiRNAWorkflow:
    """A run whose post-screen vector pays weight to ``au_1_5``.

    ``model_construct`` is deliberate: ``ScoringWeights.postscreen_sirna`` is typed to the shipped
    vector, so the experimental AU vector is the only way to exercise a run whose scored set differs
    from the module-level ``COMPOSITE_TERMS`` the hard-coded literal was copied from. Only
    ``design_params`` is passed -- ``WorkflowConfig`` refuses a resolved policy alongside different
    design parameters, and its adapter resolves the policy from these.
    """
    weights = ScoringWeights.model_construct(
        design=DesignWeights(),
        postscreen_sirna=ExperimentalAUPostScreenSiRNAWeights(),
        postscreen_mirna=PostScreenMiRNAWeights(),
    )
    params = DesignParameters()
    params.scoring = weights
    config = WorkflowConfig(output_dir=tmp_path / name, gene_query="TP53", design_params=params)
    return SiRNAWorkflow(config)


@pytest.mark.unit
def test_reported_not_scored_is_derived_from_this_runs_vectors(tmp_path):
    """The literal was wrong in three ways, and each way is a real finding.

    ``au_1_5`` and ``pos1_mismatch`` were missing -- registered terms the shipped profile does not
    score. ``paired_fraction`` was present but has no ``TermRecord`` at all: it is a gate input, and
    its record lives in ``design_parameters.filters.max_paired_fraction``.
    """
    scoring = _manifest(_qualified_workflow(tmp_path, "derived"))["scoring"]
    expected_scored, expected_unscored = _expected_split(DesignParameters().scoring)

    assert scoring["reported_not_scored"] == expected_unscored
    assert expected_unscored == ["pos1_mismatch", "au_1_5", "empirical", "isoform_coverage", "conservation"]
    # Keeps tests/unit/test_manifest_parameters.py:82-83 green: both lists stay flat.
    assert "empirical" in scoring["reported_not_scored"]
    assert "conservation" in scoring["reported_not_scored"]
    assert scoring["scored_terms"] == expected_scored
    # A reader must be able to tell a derived list from an older run's literal without diffing releases.
    assert "TERM_REGISTRY" in scoring["reported_not_scored_source"]


@pytest.mark.unit
def test_promoting_a_term_moves_it_with_no_source_edit(tmp_path):
    """The #96 drift test: composite scoring once applied two hidden normalisations nobody declared.

    Promote a term to scored and the hard-coded literal does not notice. This test is what the
    literal would have failed silently.
    """
    scoring = _manifest(_experimental_au_workflow(tmp_path, "au"))["scoring"]

    assert "au_1_5" in scoring["scored_terms"], "a term this run pays weight to is scored, not merely reported"
    assert "au_1_5" not in scoring["reported_not_scored"]
    # The unmoved terms stay put: promotion moves one term, not the vocabulary.
    for term in ("pos1_mismatch", "empirical", "isoform_coverage", "conservation"):
        assert term in scoring["reported_not_scored"]
        assert term not in scoring["scored_terms"]


@pytest.mark.unit
def test_no_scoring_term_is_invented_or_lost(tmp_path):
    """The two lists must partition the registry, for the shipped and the experimental weight sets."""
    default_scoring = _manifest(_qualified_workflow(tmp_path, "partition_default"))["scoring"]

    experimental_manifest = _manifest(_experimental_au_workflow(tmp_path, "partition_au"))

    for scoring in (default_scoring, experimental_manifest["scoring"]):
        scored = set(scoring["scored_terms"])
        unscored = set(scoring["reported_not_scored"])
        assert scored | unscored == set(TERM_REGISTRY)
        assert not scored & unscored

    # Pins the existing assertion at test_manifest_parameters.py:70.
    assert default_scoring["scored_terms"] == list(COMPOSITE_TERMS)
    # The correction: paired_fraction has no TermRecord, so it was never a candidate for weight. Its
    # record lives where it is actually applied.
    assert "paired_fraction" not in default_scoring["reported_not_scored"]
    assert "max_paired_fraction" in experimental_manifest["design_parameters"]["filters"]


# ---------------------------------------------------------------------------------------------
# Invariant 4: the report reader is unaffected, and the block is JSON-clean and stable
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_the_added_key_is_inert_for_the_report_reader(tmp_path):
    """``payload.py`` reads only ``gene_query``, ``run_policy``, ``tool_version`` and ``run_timestamp``.

    Adding one top-level key must change nothing it sees. This test never writes to payload.py; it
    asserts the boundary holds.
    """
    workflow = _qualified_workflow(tmp_path, "payload")
    reference, _ = _cached_reference(tmp_path)
    _with_resolved_human(workflow, reference)
    base = _run_step6(workflow)

    manifest = json.loads((base / "manifest.json").read_text())
    assert "provenance" in manifest

    payload = build_payload(workflow.config.output_dir)

    assert payload.provenance["tool_version"] == manifest["tool_version"] != "unknown"
    assert payload.provenance["run_timestamp"] == manifest["run_timestamp"] != "unknown"
    # The gate panel must still come from the run's own policy, not the library defaults.
    assert payload.provenance["policy_source"] == "the run's own manifest"
    assert payload.run["gene_query"] == "TP53"
    assert not [caveat for caveat in payload.caveats if "manifest" in caveat.lower()]


@pytest.mark.unit
def test_the_provenance_block_is_json_clean_and_stable_within_a_process(tmp_path):
    """No Path, set or Enum may leak in, and one process must run the git subprocess once."""
    module = _provenance_module()
    module.build_identity.cache_clear()

    calls: list[Any] = []
    real_git = module._git

    def _counted(*args: Any, **kwargs: Any) -> Any:
        calls.append(args)
        return real_git(*args, **kwargs)

    module._git = _counted
    try:
        workflow = _qualified_workflow(tmp_path, "stable")
        reference, _ = _cached_reference(tmp_path)
        _with_resolved_human(workflow, reference)

        first = workflow._build_fair_manifest(
            all_csv=tmp_path / "all.csv", pass_csv=tmp_path / "pass.csv", pass_fasta=tmp_path / "pass.fasta"
        )
        json.dumps(first)  # no Path, set, Enum or numpy scalar anywhere
        after_first = len(calls)
        second = workflow._build_fair_manifest(
            all_csv=tmp_path / "all.csv", pass_csv=tmp_path / "pass.csv", pass_fasta=tmp_path / "pass.fasta"
        )
    finally:
        module._git = real_git

    assert first["provenance"] == second["provenance"]
    assert after_first > 0, "the build block must actually consult git"
    assert len(calls) == after_first, "build_identity must be cached; the git subprocess runs once per process"


# ---------------------------------------------------------------------------------------------
# Invariant 5: the report attestation chain
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_the_report_is_attested_by_a_sidecar_not_by_rewriting_the_manifest(tmp_path):
    """``report.html`` is the artifact everyone reads and it carries no digest at all.

    A second write pass is not available: the report embeds the manifest's content, so rewriting
    manifest.json would leave the report quoting v1 while v2 sits on disk -- the artifact everyone
    reads disagreeing with the artifact that attests it. The sidecar attests both, so the pair is
    mutually verifiable with no rewrite and no circularity.
    """
    workflow = _qualified_workflow(tmp_path, "attested")
    reference, _ = _cached_reference(tmp_path)
    _with_resolved_human(workflow, reference)
    base = _run_step6(workflow)

    report = base / "report.html"
    manifest_path = base / "manifest.json"
    sidecar_path = base / "report_manifest.json"
    assert report.exists()
    assert sidecar_path.exists(), "write_report_manifest must run immediately after write_report succeeds"

    sidecar = json.loads(sidecar_path.read_text())
    assert sidecar["report_html"]["sha256"] == _sha256(report)
    assert sidecar["report_html"]["size_bytes"] == report.stat().st_size
    # manifest.json AS WRITTEN: a mismatch here means the manifest was rewritten after rendering.
    assert sidecar["manifest_json"]["sha256"] == _sha256(manifest_path)

    manifest = json.loads(manifest_path.read_text())
    artifact = manifest["provenance"]["artifacts"]["report_html"]
    assert artifact["attested"] is False
    assert artifact["state"] == "attested_by_sidecar"
    assert artifact["attestation_artifact"] == "report_manifest.json"
    assert len(artifact["reason"]) >= 20
    # `files` means "digested". The report is not one of those, and must not be listed there.
    assert "report_html" not in manifest["files"]


@pytest.mark.unit
def test_a_failed_report_is_not_logged_as_a_failed_sidecar(tmp_path, monkeypatch, caplog):
    """#103's two-handler rule, first half: a failed report is reported as exactly that.

    One shared handler once logged "Failed to write self-contained HTML report" for a report that had
    in fact succeeded, because the registration beside it was what failed.
    """
    monkeypatch.setattr(
        "sirnaforge.workflow.write_report",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("render exploded")),
    )
    workflow = _qualified_workflow(tmp_path, "no_report")

    with caplog.at_level("WARNING"):
        base = _run_step6(workflow)

    assert (base / "manifest.json").exists(), "a failed report must not cost the manifest"
    assert not (base / "report_manifest.json").exists()
    messages = caplog.text
    assert "Failed to write self-contained HTML report" in messages
    assert "report manifest" not in messages.lower(), "no sidecar failure happened; none may be reported"


@pytest.mark.unit
def test_a_failed_sidecar_is_not_logged_as_a_failed_report(tmp_path, monkeypatch, caplog):
    """The other half of the two-handler rule: the sidecar has its own handler."""
    monkeypatch.setattr(
        "sirnaforge.workflow.write_report_manifest",
        lambda *_a, **_k: (_ for _ in ()).throw(OSError("read-only run root")),
    )
    workflow = _qualified_workflow(tmp_path, "no_sidecar")

    with caplog.at_level("WARNING"):
        base = _run_step6(workflow)

    assert (base / "report.html").exists(), "the report succeeded and must be reported as such"
    assert not (base / "report_manifest.json").exists()
    assert "Failed to write self-contained HTML report" not in caplog.text
    assert "report manifest" in caplog.text.lower()


# ---------------------------------------------------------------------------------------------
# Invariant 6: the second manifest caller
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_the_offtarget_only_entry_point_is_provenanced_and_honest_about_what_it_lacks(tmp_path):
    """``write_offtarget_only_exports`` bypasses step5 and step6 and renders no report at all."""
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW, query_species="human")
    config = WorkflowConfig(output_dir=tmp_path / "offtarget_only", gene_query="TP53", resolved_policy=policy)
    workflow = SiRNAWorkflow(config)

    workflow.write_offtarget_only_exports([_candidate()])
    manifest = json.loads((config.output_dir / "sirnaforge" / "manifest.json").read_text())
    provenance = manifest["provenance"]

    # The build block needs no DesignResult, so it is fully populated here too.
    build = provenance["build"]
    assert build["release"] is False
    assert build["kind"]
    assert build["tool_version"] == __version__
    assert build["runtime"]["python"]["version"]
    # No Nextflow screen ran, so there is no pipeline revision to record -- omitted, not nulled.
    assert "pipeline_revision" not in build

    screening = provenance["references"]["screening"]
    assert screening["resolved"] == []
    assert isinstance(screening["disabled_reason"], str) and screening["disabled_reason"]

    orthology = provenance["databases"]["orthology"]
    assert orthology["attempted"] is False
    assert orthology["source"] is None, "empty() must stop defaulting to compara"
    assert len(orthology["reason"]) >= 20

    # Nothing was expected here, so "not attested" would be the wrong word.
    assert provenance["artifacts"]["report_html"]["state"] == "not_rendered_by_this_entry_point"


# ---------------------------------------------------------------------------------------------
# Invariant 7: orthology never claims a call it did not make
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_orthology_never_claims_a_compara_call_it_did_not_make(tmp_path):
    """A pre-existing defect, fixed before anything is published.

    Today ``empty()`` inherits the ``SOURCE_COMPARA`` field default, so every single-species screen,
    every design-only run and every run whose gene yielded neither ID nor symbol already publishes
    ``{"source": "ensembl_compara", "orthologue_types": [all three]}`` into workflow_summary.json --
    a REST call that never happened. Copying ``summary()`` into manifest.json without this fix
    promotes fabricated provenance into the attested artifact.
    """
    summary = OrthologueMapping.empty().summary()

    assert summary["source"] is None
    assert "orthologue_types" not in summary, "types accepted by a call that never ran are not provenance"

    orthology = _provenance(_design_only_workflow(tmp_path, "design_only_orthology"))["databases"]["orthology"]
    assert orthology["attempted"] is False
    assert orthology["source"] is None
    assert len(orthology["reason"]) >= 20


# ---------------------------------------------------------------------------------------------
# Invariant 8: index attestation is three-state, and digest_mode survives
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_an_unstamped_index_is_named_unattested_rather_than_assumed_good(tmp_path):
    """bwa-mem2 aligns against the index files alone; an index nobody stamped is unattested bytes."""
    workflow = _qualified_workflow(tmp_path, "no_stamp")
    reference, paths = _cached_reference(tmp_path, with_index_stamp=False)
    _with_resolved_human(workflow, reference)
    assert read_artifact_stamp(Path(paths["index"])) is None, "fixture guard: no stamp beside this index"

    index = _provenance(workflow)["references"]["screening"]["resolved"][0]["index"]

    assert index["attested"] is False
    assert index["state"] == "no_stamp_for_this_index"
    assert index["members"] == []
    assert index["built_from_digest"] is None
    assert index["aligner"] == "bwa-mem2"
    assert index["aligner_version"]["value"] is None
    assert len(index["reason"]) >= 20


@pytest.mark.unit
def test_a_stamped_index_proves_it_agrees_with_its_reference_and_keeps_digest_mode(tmp_path):
    """``built_from_digest`` is what proves the index and the reference are the same release.

    ``digest_mode`` is mandatory and never dropped: every bwa-mem2 member of a human index exceeds
    the full-digest budget, so a sampled head/tail digest must not read as a whole-file one. The
    members are never collapsed into one "index digest" -- a four-file artifact proves its own
    completeness only member by member.
    """
    workflow = _qualified_workflow(tmp_path, "stamped")
    reference, paths = _cached_reference(tmp_path, with_index_stamp=True)
    _with_resolved_human(workflow, reference)

    entry = _provenance(workflow)["references"]["screening"]["resolved"][0]
    index = entry["index"]

    assert index["attested"] is True
    assert index["state"] == "artifact_stamp"
    assert index["built_from_digest"] == entry["bytes"]["digest"]["value"] == paths["checksum"]
    assert "digest" not in index, "members are never collapsed into one index digest"

    members = {member["name"]: member for member in index["members"]}
    assert len(members) == 2
    for member in members.values():
        assert member["digest"]
        assert member["digest_algorithm"] == "md5"
        assert member["digest_mode"] in {"full", "sampled"}
        assert member["size"] == 2048
    assert {member["digest_mode"] for member in members.values()} == {"full", "sampled"}


# ---------------------------------------------------------------------------------------------
# Invariant 9: bookkeeping never fails a run
# ---------------------------------------------------------------------------------------------


@pytest.mark.unit
def test_a_provenance_failure_is_named_rather_than_silent(tmp_path, monkeypatch):
    """Provenance is bookkeeping: it must never fail a run, and it must never fail quietly either."""
    monkeypatch.setattr(
        "sirnaforge.workflow.build_identity",
        lambda *_a, **_k: (_ for _ in ()).throw(RuntimeError("git went missing mid-run")),
    )
    workflow = _qualified_workflow(tmp_path, "assembly_failed")
    base = _run_step6(workflow)

    manifest = json.loads((base / "manifest.json").read_text())

    assert manifest["files"], "the rest of the manifest must survive a provenance failure"
    assert manifest["provenance"] == {
        "schema_version": "1.0",
        "state": "provenance_assembly_failed",
        "reason": "git went missing mid-run",
    }
