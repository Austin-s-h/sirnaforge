"""siRNAforge Workflow Orchestrator.

Coordinates the complete siRNA design pipeline:
1. Transcript retrieval and validation
2. ORF validation and reporting
3. siRNA candidate generation and scoring
4. Top-N candidate selection and reporting
5. Off-target analysis with Nextflow pipeline
"""

from __future__ import annotations

import asyncio
import csv
import json
import math
import os
import re
import shutil
import tempfile
import time
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, NamedTuple, cast

import pandas as pd
from Bio.Seq import Seq
from pandera.typing import DataFrame
from rich.console import Console
from rich.progress import Progress

from sirnaforge import __version__
from sirnaforge.config import (
    DEFAULT_TRANSCRIPTOME_SOURCES,
    UNRESOLVED_SPECIES,
    ReferenceChoice,
    ReferenceForm,
    ReferencePolicyResolver,
    ReferenceRejection,
    ReferenceRequest,
    ReferenceSelection,
    ScreeningReference,
    ScreeningReferenceSet,
    SpeciesAuthority,
    WorkflowInputSpec,
    build_screening_requests,
    screening_kind_for_design_mode,
)
from sirnaforge.config.reference_policy import parse_index_entries, resolve_reference_species
from sirnaforge.config.run_policy import (
    FILTER_SPEC_BY_ID,
    EntryPoint,
    ResolvedRunPolicy,
    RunPolicyError,
    describe_parameters,
    resolve_run_policy,
)
from sirnaforge.core.design import (
    MiRNADesigner,
    SiRNADesigner,
    biogenesis_features,
    record_gate_outcome,
)
from sirnaforge.core.filtering import (
    Comparator,
    GateSpec,
    evaluate_gate,
    evaluate_gates,
    first_rejection,
)
from sirnaforge.core.hit_annotation import (
    CLASSIFICATION_COLUMNS,
    LIABILITY_CLASSES,
    HitAnnotator,
    accumulate_hit_class,
    annotate_hit_row,
    count_persisted_classes,
    is_annotated,
    liabilities_counted,
    write_classified_hits,
)
from sirnaforge.core.hit_classification import (
    ClassificationContext,
    HitClass,
    HitClassCounts,
    classify_hit,
)
from sirnaforge.core.off_target import OffTargetAnalysisManager
from sirnaforge.core.repeat_detection import (
    RepeatDetector,
    normalize_guide_sequence,
)
from sirnaforge.core.scoring import (
    SCORING_WEIGHT_SET_VERSION,
    ScoringError,
    compute_composite,
    conservation_sub_score,
    isoform_coverage_sub_score,
    off_target_sub_score,
)
from sirnaforge.core.screening_evidence import (
    RECONCILIATION_FILENAME,
    EvidenceEnvelope,
    EvidenceProducer,
    EvidenceSource,
    Reconciliation,
    build_plan,
    collect_evidence,
    completed_pairs,
    guide_set_digest,
    not_requested_entry,
    parse_reconciliation_payload,
    reconcile,
    reconciliation_payload,
)
from sirnaforge.core.selection import (
    CandidateView,
    SelectionInputs,
    select,
)
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.data.base import DatabaseType, FastaUtils, TranscriptInfo
from sirnaforge.data.ensembl_references import ENSEMBL_ASSEMBLIES, infer_species_from_cdna_headers
from sirnaforge.data.gene_search import GeneSearcher
from sirnaforge.data.orf_analysis import ORFAnalyzer
from sirnaforge.data.orthology import (
    ENSEMBL_BASE_URL,
    SOURCE_COMPARA,
    OrthologTable,
    OrthologueMapping,
    load_ortholog_table,
    mapping_from_table,
    resolve_orthologues,
)
from sirnaforge.data.species_registry import (
    CANONICAL_SPECIES_REGISTRY,
    MIRGENEDB_SPECIES_TABLE,
    normalize_species_name,
)
from sirnaforge.data.transcript_annotation import EnsemblTranscriptModelClient
from sirnaforge.data.transcript_index import TranscriptGeneIndex
from sirnaforge.data.transcriptome_manager import INDEX_BUILD_ERROR_KEY, TranscriptomeManager
from sirnaforge.models.evidence import (
    EvidenceStatus,
    ScreeningEvidenceEntry,
    ScreeningPlan,
)
from sirnaforge.models.policy import (
    EvidenceRequirements,
    FilterAction,
    FilterEvaluation,
    Requiredness,
    RunMode,
    RunStatus,
    ScreeningChannel,
)
from sirnaforge.models.schemas import ORFValidationSchema, SiRNACandidateSchema
from sirnaforge.models.scoring_profile import TERM_REGISTRY
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    DesignResult,
    OffTargetFilterCriteria,
    SelectionState,
    SiRNACandidate,
    build_candidate_row,
    ranking_score,
)
from sirnaforge.models.sirna import SiRNACandidate as _ModelSiRNACandidate
from sirnaforge.models.variant import VariantRecord
from sirnaforge.models.zfn import (
    GenomicAnnotationConfig,
    ZFNDesignParameters,
    ZFNDesignResult,
    ZFNShardingConfig,
)
from sirnaforge.pipeline import NextflowConfig, NextflowRunner
from sirnaforge.provenance import (
    PROVENANCE_SCHEMA_VERSION,
    absent,
    assembly_identity,
    build_identity,
    digest_scope,
    index_attestation,
    pipeline_revision_identity,
    present,
    report_html_artifact,
    reported_term_partition,
    write_report_manifest,
)
from sirnaforge.reporting import ReportPayload, build_payload, write_quilt_summarize, write_report
from sirnaforge.utils.cache_utils import resolve_cache_subdir, stable_cache_key
from sirnaforge.utils.control_candidates import DIRTY_CONTROL_LABEL, inject_dirty_controls
from sirnaforge.utils.ensembl_ids import strip_version
from sirnaforge.utils.hashing import file_sha256
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.modification_patterns import apply_modifications_to_candidate
from sirnaforge.utils.parsing import parse_csv
from sirnaforge.utils.resource_resolver import InputSource, resolve_input_source
from sirnaforge.utils.species import is_human_species
from sirnaforge.validation import ValidationConfig, ValidationMiddleware
from sirnaforge.workflow_variant import (
    VariantWorkflowConfig,
    normalize_variant_mode,
    parse_clinvar_filter_string,
    resolve_workflow_variants,
)
from sirnaforge.zfn import emit_zfn_experimental_warning
from sirnaforge.zfn.design import ZFNDesigner

logger = get_logger(__name__)
console = Console(record=True, force_terminal=False, legacy_windows=True)

#: Keys of one ``filtering_stats.per_species`` bucket: the five classes plus the two shortfall flags
#: ``accumulate_hit_class`` tallies. Named once so the pre-seeded bucket and the tally cannot drift.
_PER_SPECIES_COUNTERS: tuple[str, ...] = (
    *(member.value for member in HitClass),
    "symbol_lookup_missing",
    "species_index_missing",
)


#: Guide-set digest stamped on evidence a run synthesized without ever recording a plan (a direct
#: call into result processing). Never a real 16-character sha256 prefix, so it cannot collide with
#: one and can only ever explain an absence.
_UNRECORDED_GUIDE_SET_DIGEST = "unrecorded-guide-set"


def describe_shortfall_reasons(reasons: Mapping[str, int]) -> str:
    """One human sentence naming what was actually missing, most costly first.

    ``no_evidence:<channel>:<species>`` and ``unknown:<filter_id>`` are different problems with
    different fixes -- a screen to re-run against, versus an annotation the gate could not read --
    so the console says which rather than attributing both to screening.
    """
    parts: list[str] = []
    for reason, count in reasons.items():
        kind, _, detail = reason.partition(":")
        if kind == "no_evidence":
            channel, _, species = detail.partition(":")
            parts.append(f"{count}x no {channel} evidence for {species}")
        elif kind == "unknown":
            parts.append(f"{count}x undecided {detail}")
        else:  # pragma: no cover - defensive: an unrecognised reason is still reported verbatim
            parts.append(f"{count}x {reason}")
    return "; ".join(parts) if parts else "no reason recorded"


def _describe_design_failure(error: BaseException) -> str:
    """Why a transcript was dropped, in one line that is never blank.

    ``str(exc)`` is empty for a bare ``RuntimeError()`` and for most cancellations, and a shortfall
    recorded with no reason is indistinguishable from one recorded by mistake -- the same rule
    ``models/evidence.py``'s ``failure_carries_a_reason`` validator enforces on screening evidence.
    """
    detail = str(error).strip()
    return f"{type(error).__name__}: {detail}" if detail else type(error).__name__


#: Reason the route through result processing reports when the pipeline exited 0 and staged no output
#: directory at all. Named here because that route publishes ``partial`` like any thin result, and
#: writing nothing is a failure to execute rather than a partial answer.
_EXECUTION_ERROR_REASON = "missing_output"

#: What each ``reason`` a screening exit publishes means in :class:`RunStatus` terms. Keyed on the
#: reason rather than on ``status`` because the exits that spell themselves ``skipped`` mean opposite
#: things: a channel the user switched off was never requested, while an engine that could not be
#: executed is an error the caller must act on. A reason absent from here does not fall back to a
#: guess -- see :func:`screening_run_status`.
_RUN_STATUS_BY_REASON: Mapping[str, RunStatus] = {
    "user_disabled": RunStatus.NOT_REQUESTED,
    "no_candidates": RunStatus.NO_ELIGIBLE,
    "nextflow_failed": RunStatus.EXECUTION_ERROR,
    "nextflow_unavailable": RunStatus.EXECUTION_ERROR,
    _EXECUTION_ERROR_REASON: RunStatus.EXECUTION_ERROR,
}


def screening_run_status(*, status: str, reason: str = "", required_evidence_missing: bool = False) -> RunStatus:
    """Translate one screening exit's ``{status, reason}`` pair into the one run-status vocabulary.

    Pure and total, because it is the single definition of what those seven string pairs mean (#100)
    and every exit is stamped through it. The result is published beside ``status``/``reason``, never
    instead of them.

    ``required_evidence_missing`` is what lets a step whose own word is ``completed`` report
    ``INCOMPLETE``: ``status`` is decided from the aggregate's self-report, which says nothing was
    missing whenever aggregation itself is what failed, while the reconciled evidence knows the unit
    was planned and published nothing. An unrecognised status word claims no completeness rather
    than inheriting one, so a new exit cannot arrive silently reading ``COMPLETED``.
    """
    known = _RUN_STATUS_BY_REASON.get(reason)
    if known is not None:
        return known
    if status == "missing":
        return RunStatus.EXECUTION_ERROR
    if status == "completed":
        return RunStatus.INCOMPLETE if required_evidence_missing else RunStatus.COMPLETED
    return RunStatus.INCOMPLETE


class CandidateTablePaths(NamedTuple):
    """Where a run's candidate tables and manifest live, named once for both entry points.

    ``step6_generate_reports`` and :meth:`SiRNAWorkflow.write_offtarget_only_exports` publish the same
    six files into the same ``sirnaforge/`` sub-directory, which is what lets ``sirnaforge report
    <output_dir>`` read either kind of run with no special case.
    """

    all_csv: Path
    pass_csv: Path
    pass_fasta: Path
    qualified_csv: Path
    qualified_fasta: Path
    provisional_csv: Path
    manifest_json: Path

    @classmethod
    def under(cls, base: Path) -> CandidateTablePaths:
        """The table paths under one run's ``sirnaforge/`` directory."""
        return cls(
            all_csv=base / "candidates_all.csv",
            pass_csv=base / "candidates_pass.csv",
            pass_fasta=base / "candidates_pass.fasta",
            qualified_csv=base / "candidates_qualified.csv",
            qualified_fasta=base / "candidates_qualified.fasta",
            provisional_csv=base / "candidates_provisional.csv",
            manifest_json=base / "manifest.json",
        )


class OffTargetGateCounts(NamedTuple):
    """The counts the off-target gates compare, in the order ``_check_offtarget_filters`` reads them.

    Named so the no-hit and with-hit paths can hand the same eight numbers to the same gate call.
    Every field defaults to zero and :data:`_ZERO_OFFTARGET_COUNTS` is that all-zero instance: a
    completed screen that found nothing is a *measurement* of zero, and giving it a name is what
    stopped the clean-screen path from skipping the gates entirely (issue #106).
    """

    transcriptome_0mm: int = 0
    transcriptome_1mm: int = 0
    transcriptome_2mm: int = 0
    transcriptome_seed_0mm: int = 0
    mirna_0mm_seed: int = 0
    mirna_high_risk: int = 0
    total_hits: int = 0
    genuine_off_target_count: int = 0


_ZERO_OFFTARGET_COUNTS = OffTargetGateCounts()

_TRANSCRIPTOME_CHANNEL = frozenset({ScreeningChannel.TRANSCRIPTOME})
_MIRNA_SEED_CHANNEL = frozenset({ScreeningChannel.MIRNA_SEED})

#: Which screening channel supplies each post-screen gate's count. One mapping, because
#: ``_check_offtarget_filters`` and ``_apply_post_screen_ranking`` both need this answer and must not
#: diverge. A declared filter absent from it reads design-time or annotation evidence instead
#: (``min_isoform_coverage``, every design-stage gate), so an incomplete channel never excuses it.
POST_SCREEN_FILTER_CHANNELS: Mapping[str, frozenset[ScreeningChannel]] = {
    "max_transcriptome_hits_0mm": _TRANSCRIPTOME_CHANNEL,
    "max_transcriptome_hits_1mm": _TRANSCRIPTOME_CHANNEL,
    "max_transcriptome_hits_2mm": _TRANSCRIPTOME_CHANNEL,
    "max_transcriptome_seed_perfect": _TRANSCRIPTOME_CHANNEL,
    "max_off_target_count": _TRANSCRIPTOME_CHANNEL,
    "max_mirna_perfect_seed": _MIRNA_SEED_CHANNEL,
    "max_mirna_1mm_seed": _MIRNA_SEED_CHANNEL,
    "fail_on_high_risk_mirna": _MIRNA_SEED_CHANNEL,
    # Counts two channels at once, so it is decided by neither alone.
    "max_total_offtarget_hits": _TRANSCRIPTOME_CHANNEL | _MIRNA_SEED_CHANNEL,
}

#: Which run statistic each off-target rejection increments. Shared by the no-hit and with-hit paths.
_OFFTARGET_REJECTION_STATS: Mapping[SiRNACandidate.FilterStatus, str] = {
    SiRNACandidate.FilterStatus.TRANSCRIPTOME_PERFECT_MATCH: "failed_perfect_match",
    SiRNACandidate.FilterStatus.TRANSCRIPTOME_1MM: "failed_transcriptome_1mm",
    SiRNACandidate.FilterStatus.TRANSCRIPTOME_2MM: "failed_transcriptome_2mm",
    SiRNACandidate.FilterStatus.TRANSCRIPTOME_SEED_PERFECT: "failed_transcriptome_seed_perfect",
    SiRNACandidate.FilterStatus.MIRNA_PERFECT_SEED: "failed_mirna_seed",
    SiRNACandidate.FilterStatus.HIGH_RISK_MIRNA: "failed_high_risk_mirna",
    SiRNACandidate.FilterStatus.EXCESS_OFF_TARGETS: "failed_excess_off_targets",
}


#: 0.7.0 keyword arguments the #99 rename removed, and what to pass instead. A hard break: nothing
#: silently maps an old name onto the new one, because the two were never quite the same thing --
#: ``genome_indices_override`` bypassed reference resolution, and ``transcriptome_indices`` does not.
#: They are still named here so the failure quotes the replacement instead of "unexpected keyword".
RENAMED_ARGUMENTS: Mapping[str, str] = {
    "genome_species": "screen_species",
    "genome_indices_override": "transcriptome_indices",
}


#: Pipeline parameters the rename removed. Refused rather than ignored: Nextflow accepts an unknown
#: ``--param`` silently, so a stale ``genome_indices`` in a raw ``nextflow_config`` would configure
#: nothing and the run would screen against nothing and report success.
RENAMED_NEXTFLOW_PARAMS: Mapping[str, str] = {
    "genome_indices": "transcriptome_indices",
    "genome_fastas": "transcriptome_fastas",
    "genome_species": "transcriptome_species",
}


def _policy_filter_actions(policy: ResolvedRunPolicy | None) -> dict[str, FilterAction] | None:
    """The per-filter actions a resolved policy decided, or None when there is no policy to ask.

    None rather than an empty dict on purpose: empty would read as "every filter has no action" and
    the designer would fall back to declared defaults for a run that did resolve a policy.
    """
    if policy is None:
        return None
    return {resolved.filter_id: resolved.descriptor.action for resolved in policy.filters}


def refuse_renamed_arguments(supplied: Mapping[str, Any]) -> None:
    """Refuse a removed ``genome_*`` keyword, naming its replacement.

    Raises:
        TypeError: Any keyword was supplied. Unknown keywords are reported the same way Python
            would, so a typo does not read as a rename.
    """
    if not supplied:
        return
    renamed = sorted(name for name in supplied if name in RENAMED_ARGUMENTS)
    unknown = sorted(name for name in supplied if name not in RENAMED_ARGUMENTS)
    parts = [f"{name!r} was renamed to {RENAMED_ARGUMENTS[name]!r} in 0.7.1" for name in renamed]
    parts += [f"unexpected keyword argument {name!r}" for name in unknown]
    detail = "; ".join(parts)
    if renamed:
        detail += (
            ". Screening references are transcriptomes: 'genome' now means genomic DNA and belongs to ZFN only (#99)"
        )
    raise TypeError(detail)


class WorkflowConfig:
    """Configuration for the complete siRNA design workflow."""

    def __init__(
        self,
        output_dir: Path,
        gene_query: str,
        input_fasta: Path | None = None,
        database: DatabaseType = DatabaseType.ENSEMBL,
        design_params: DesignParameters | None = None,
        # off-target selection now always equals design_params.top_n
        nextflow_config: Mapping[str, Any] | None = None,
        transcriptome_indices: str | None = None,
        screen_species: list[str] | None = None,
        query_species: str | None = None,
        mirna_database: str = "mirgenedb",
        mirna_species: Sequence[str] | None = None,
        transcriptome_fasta: str | None = None,
        transcriptome_filter: str | None = None,
        transcriptome_selection: ReferenceSelection | None = None,
        ortholog_mapping_file: Path | str | None = None,
        validation_config: ValidationConfig | None = None,
        log_file: str | None = None,
        write_json_summary: bool = True,
        num_threads: int | None = None,
        input_source: InputSource | None = None,
        keep_nextflow_work: bool = False,
        variant_config: VariantWorkflowConfig | None = None,
        zfn_config: ZFNWorkflowConfig | None = None,
        resolved_policy: ResolvedRunPolicy | None = None,
        **renamed: Any,
    ):
        """Initialize workflow configuration.

        ``resolved_policy`` is the run policy resolved once by ``config/run_policy.py``. A caller
        that supplies its own ``design_params`` instead gets the same object built by the adapter,
        so every path -- CLI, Python API, off-target-only, direct WorkflowConfig -- carries a policy
        and none of them re-derives a threshold. Resolution happens before the output directories
        below are created, so an invalid configuration costs nothing.

        ``screen_species`` are the species this run asks to screen; ``transcriptome_indices`` are
        ``species:index_prefix`` references the caller has already built. Both feed one resolver
        (#99), so an override and a default differ only in provenance. ``**renamed`` exists solely
        to refuse the 0.7.0 ``genome_*`` spellings with the new name rather than a bare TypeError.
        """
        refuse_renamed_arguments(renamed)
        self.output_dir = Path(output_dir)
        self.input_source = input_source
        self.zfn_config = zfn_config

        resolved_input = input_source.local_path if input_source else (Path(input_fasta) if input_fasta else None)
        self.input_fasta = resolved_input
        # Preserve the user-supplied gene_query as the logical label even when using an input FASTA
        self.gene_query = gene_query
        self.database = database
        if (
            resolved_policy is not None
            and design_params is not None
            and design_params != resolved_policy.design_parameters
        ):
            raise RunPolicyError(
                "WorkflowConfig received both a resolved policy and different design_params; "
                "pass one or the other so there is a single resolved configuration"
            )
        self.design_params = (
            resolved_policy.design_parameters if resolved_policy else (design_params or DesignParameters())
        )
        self.resolved_policy = resolved_policy or describe_parameters(
            self.design_params, entry_point=EntryPoint.SCREENING_WORKFLOW, query_species=query_species
        )
        # single source of truth: number of candidates selected everywhere
        self.top_n = self.design_params.top_n
        self.nextflow_config: dict[str, Any] = dict(nextflow_config) if nextflow_config else {}
        stale_params = sorted(set(self.nextflow_config) & set(RENAMED_NEXTFLOW_PARAMS))
        if stale_params:
            raise ValueError(
                "nextflow_config carries pipeline parameters the #99 rename removed: "
                + ", ".join(f"{name!r} is now {RENAMED_NEXTFLOW_PARAMS[name]!r}" for name in stale_params)
                + ". Nextflow ignores an unknown parameter, so leaving them would screen against nothing."
            )

        # Explicit index entries are references like any other: they are parsed here, before any
        # directory is created, and resolved by the same resolver as the defaults rather than being
        # written straight into the Nextflow parameters. Writing them straight through is what let a
        # cDNA file align successfully and then classify against nothing (#99 defect 1).
        index_requests = parse_index_entries(
            transcriptome_indices,
            option="--transcriptome-indices",
            reason="explicit index override (--transcriptome-indices)",
        )
        override_species = [request.declared_species or "" for request in index_requests] or None

        # Species this run asks to screen against. Track whether they were requested or defaulted.
        species_explicitly_provided = screen_species is not None or override_species is not None
        requested_species = screen_species or ["human", "rat", "rhesus"]
        if override_species:
            requested_species = override_species

        # Normalize all species names to canonical form for consistent comparisons
        self.screen_species: list[str] = list(dict.fromkeys(normalize_species_name(s) for s in requested_species))
        self.species_explicitly_requested = species_explicitly_provided
        # Organism of the TARGET transcripts, when the caller states it outright. None means
        # "derive it from where the transcripts actually came from" -- see SiRNAWorkflow.__init__.
        # Deliberately NOT defaulted from screen_species: that list is an unordered set of
        # references to screen against, so no position in it identifies the target.
        stated_query_species = (query_species or "").strip()
        self.query_species: str | None = normalize_species_name(stated_query_species) if stated_query_species else None
        self.mirna_database = mirna_database
        # One species vocabulary for the miRNA channel too, in the caller's own order (#100). The CLI
        # resolves --species into miRNA *database codes* -- `--species human,rhesus` becomes
        # ['hsa', 'mml'] -- and unnormalised codes make the plan/evidence join key on 'human' never
        # match. Normalizing here is enough because this attribute is the single source of the plan
        # entry, the pipeline parameter and the published summary, and MiRNADatabaseManager resolves a
        # canonical name to the same source as its code (identical cache_key for mirgenedb; the mirbase
        # sources are keyed on the canonical name outright), so the screen reads the same database.
        if mirna_species:
            normalized_mirna = [normalize_species_name(value) for value in mirna_species if value]
            self.mirna_species = list(dict.fromkeys(normalized_mirna))
        else:
            self.mirna_species = []
        # Store transcriptome filter for later use
        self.transcriptome_filter = transcriptome_filter
        if transcriptome_selection is None and transcriptome_fasta:
            transcriptome_selection = ReferenceSelection(
                choices=(ReferenceChoice.explicit(transcriptome_fasta, reason="legacy transcriptome argument"),)
            )
        self.transcriptome_selection = transcriptome_selection or ReferenceSelection.disabled(
            "no transcriptome configured"
        )
        self.transcriptome_references = [
            choice.value for choice in self.transcriptome_selection.choices if choice.value
        ]
        # The one door. Every screening reference -- explicit index override or resolved default --
        # becomes a request here, and a modality mismatch (a genomic assembly handed to an siRNA run,
        # a transcriptome handed to ZFN) raises now, before a multi-gigabyte download.
        self.screening_kind = screening_kind_for_design_mode(self.design_params.design_mode)
        self.screening_requests, self.screening_request_rejections = build_screening_requests(
            kind=self.screening_kind,
            selection=self.transcriptome_selection,
            index_requests=index_requests,
        )
        # Offline orthologue evidence (#101): when set, cross-species classification reads this file
        # instead of calling Ensembl Compara, so an air-gapped run -- and every fixture -- is
        # deterministic and never waits on REST.
        self.ortholog_mapping_file = Path(ortholog_mapping_file) if ortholog_mapping_file else None
        self.validation_config = validation_config or ValidationConfig()
        self.log_file = log_file
        self.write_json_summary = write_json_summary
        self.keep_nextflow_work = keep_nextflow_work
        # Parallelism for design stage (cap at 4 CPUs for better efficiency)
        requested_threads = num_threads if num_threads is not None else (os.cpu_count() or 4)
        self.num_threads = max(1, min(4, requested_threads))
        # Variant targeting configuration
        self.variant_config = variant_config

        if self.mirna_database and self.mirna_species:
            self.nextflow_config.setdefault("mirna_db", self.mirna_database)
            self.nextflow_config.setdefault("mirna_species", ",".join(self.mirna_species))

        # Create output structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "sirnaforge").mkdir(exist_ok=True)
        (self.output_dir / "off_target").mkdir(exist_ok=True)
        (self.output_dir / "logs").mkdir(exist_ok=True)
        if self.design_params.design_mode != DesignMode.ZFN:
            (self.output_dir / "transcripts").mkdir(exist_ok=True)
            (self.output_dir / "orf_reports").mkdir(exist_ok=True)


class ZFNWorkflowConfig:
    """Configuration for ZFN pair evaluation and off-target search workflow.

    This carries the scientifically distinct ZFN parameters:
    - Left/right half-site sequences (IUPAC-validated, 9-18 bp)
    - Genomic search space (whole-genome FASTA)
    - Algorithm choice (homology / conserved_g / zfn_v2)
    - Spacer/dimer/mismatch constraints
    - Optional genomic annotation for region classification
    """

    def __init__(
        self,
        zfn_params: ZFNDesignParameters,
        annotation: GenomicAnnotationConfig | None = None,
    ):
        """Initialize ZFN workflow configuration."""
        self.zfn_params = zfn_params
        self.annotation = annotation


class SiRNAWorkflow:
    """Main workflow orchestrator for siRNA/miRNA/ZFN design pipeline."""

    def __init__(self, config: WorkflowConfig):
        """Initialize the workflow orchestrator."""
        self.config = config
        self.gene_searcher = GeneSearcher()
        self.orf_analyzer = ORFAnalyzer()
        self.validation = ValidationMiddleware(config.validation_config)

        # Select designer based on design mode
        self.zfn_designer: ZFNDesigner | None = None
        self.sirnaforgeer: SiRNADesigner | MiRNADesigner | None = None
        # The design gates need the ACTIONS as well as the thresholds. DesignParameters carries only
        # the numbers, so without this a `--filter-action` override would reach the off-target gates
        # and be silently ignored by the design ones.
        design_actions = _policy_filter_actions(getattr(config, "resolved_policy", None))
        if config.design_params.design_mode == DesignMode.ZFN:
            self.zfn_designer = ZFNDesigner()
        elif config.design_params.design_mode == DesignMode.MIRNA:
            self.sirnaforgeer = MiRNADesigner(config.design_params, filter_actions=design_actions)
        else:
            self.sirnaforgeer = SiRNADesigner(config.design_params, filter_actions=design_actions)

        self.results: dict[str, Any] = {}
        self._nextflow_cache_info: dict[str, Any] | None = None
        # What step5's orthologue lookup produced. None is not "nothing resolved": write_offtarget_only
        # and design-only runs never reach step5, so the manifest reports orthology as not attempted
        # rather than publishing a Compara call that never happened.
        self._orthology_mapping: OrthologueMapping | None = None
        self._annotation_summary: dict[str, Any] = {}
        self._gene_transcript_ids: set[str] = set()
        self._query_gene_ids: set[str] = set()
        self._query_gene_symbols: set[str] = set()
        self._protein_coding_transcript_ids: set[str] = set()
        self._protein_coding_transcript_count: int = 0
        self._transcript_index = TranscriptGeneIndex()
        self._species_explicitly_requested: bool = False
        # Species actually handed to Nextflow, which is a superset of config.screen_species when
        # extra species arrive on transcriptome_indices/transcriptome_fastas. Conservation is scored
        # against this list, so its denominator can never be smaller than the set screened.
        self._active_screen_species: list[str] = []
        # What the one resolver produced: {species, kind, identity, index} per reference, plus the
        # requests it could not use. Empty until _resolve_screening_references runs.
        self._screening_references = ScreeningReferenceSet(kind=config.screening_kind)
        # What this run intends to screen, digest-keyed to the guide set actually submitted. Recorded
        # from the REQUESTED references, before resolution can drop one: a plan derived from what
        # survived resolution cannot represent a missing reference at all (#100). None until the
        # screening stage records one; a run that asked for nothing records an empty plan, which
        # claims nothing.
        self._screening_plan: ScreeningPlan | None = None
        # The plan reconciled against what the run actually published, and the authority for
        # ``_completed_evidence_pairs``. None until screening reports, which is not the same as
        # reporting that nothing completed.
        self._screening_evidence: Reconciliation | None = None
        # Digest of the guide set the plan and every evidence entry are keyed on. Held so a
        # synthesized entry joins the plan entry it answers for rather than opening a second key.
        self._guide_set_digest: str | None = None
        # Distinct guides handed to the aligner, as counted when the FASTA was written. None means
        # nothing was submitted through this workflow, never "all of them".
        self._submitted_guide_count: int | None = None
        # Species requested for screening that never reached Nextflow, and why.
        self._species_screening_shortfalls: dict[str, str] = {}
        # Transcripts a design failure dropped, and why -- the design-stage twin of the record above
        # (#100). The parallel design path catches a failed batch and carries on, so the batch's
        # transcripts leave the candidate pool while DesignResult.total_sequences keeps counting them.
        self._design_input_shortfalls: dict[str, str] = {}
        self._representative_to_candidates: dict[str, list[SiRNACandidate]] = {}
        self._candidate_id_to_representative: dict[str, str] = {}
        # Channel/species pairs with completed evidence, spelled as ``ChannelRequirement.key`` so they
        # compare directly against ``EvidenceRequirements``. None means screening has not reported yet,
        # which is not the same as reporting that nothing completed.
        self._completed_evidence_pairs: frozenset[tuple[str, str]] | None = None
        # Species a gate with an unrestricted scope counts, i.e. what this run screened. Set alongside
        # the pair record so the two cannot describe different runs.
        self._screened_species_scope: frozenset[str] = frozenset()
        # What the last selection decided and why, so an empty shortlist names its cause.
        self._selection_summary: dict[str, Any] = {}
        # Every eligible candidate, not the top_n slice. None until a selection has run.
        self._eligible_candidate_ids: frozenset[str] | None = None
        # Single authoritative query species, set once (not re-inferred per call site), and never
        # screen_species[0] -- that list's order carries no meaning.
        # The real answer is a property of where the target transcripts came from: the gene-query
        # database (GeneSearcher.query_species). An input FASTA states no organism, so it takes the
        # same answer, and WorkflowConfig(query_species=...) states it outright when it differs.
        self._query_species: str = self.config.query_species or self.gene_searcher.query_species(self.config.database)
        # Parsed here, not at first use: a malformed mapping file must be reported before a screen
        # runs, not after it. None means "no file supplied, resolve orthologues over REST".
        self._ortholog_table: OrthologTable | None = (
            load_ortholog_table(self.config.ortholog_mapping_file) if self.config.ortholog_mapping_file else None
        )
        self._species_cdna_fasta: dict[str, Path] = {}
        self._guide_to_transcripts: dict[str, frozenset[str]] = {}
        self._repeat_summary: dict[str, Any] = {"status": "not_run"}

        # Optional: Initialize transcript annotation client (not used by default yet)
        # This can be enabled via environment variable or config flag in the future
        try:
            self._annotation_client: EnsemblTranscriptModelClient | None = EnsemblTranscriptModelClient()
        except Exception:
            self._annotation_client = None
        self._dirty_controls_added: int = 0
        # The payload and report path step 6 registered with Quilt, kept so the registration can be
        # re-issued once logs/workflow_summary.json exists. None until a report has been written, which
        # is what stops the re-issue from naming a summary no run produced.
        self._report_registration: tuple[ReportPayload, Path] | None = None

    async def run_complete_workflow(self) -> dict[str, Any]:
        """Run the complete design workflow (siRNA/miRNA or ZFN)."""
        # ── ZFN mode: fundamentally different execution path ──
        if self.config.design_params.design_mode == DesignMode.ZFN:
            return await self._run_zfn_workflow()

        console.print("\n🧬 [bold cyan]Starting siRNAforge Workflow[/bold cyan]")
        console.print(f"Gene Query: [yellow]{self.config.gene_query}[/yellow]")
        console.print(f"Output Directory: [blue]{self.config.output_dir}[/blue]")

        start_time = time.perf_counter()

        # Validate input parameters (quiet: avoid verbose warnings in console)
        _ = self.validation.validate_input_parameters(self.config.design_params)

        with Progress(console=console) as progress:
            main_task = progress.add_task("[cyan]Overall Progress", total=6)

            # Step 1: Transcript Retrieval
            progress.update(main_task, description="[cyan]Retrieving transcripts...")
            transcripts = await self.step1_retrieve_transcripts(progress)
            progress.advance(main_task)

            # All isoforms of the queried gene are "on-target" - a guide hitting a sibling
            # isoform (shared exon) perfectly is not an off-target, it's the intended gene.
            self._gene_transcript_ids = {strip_version(t.transcript_id) for t in transcripts if t.transcript_id}
            self._query_gene_ids = {strip_version(t.gene_id) for t in transcripts if t.gene_id}
            self._query_gene_symbols = {t.gene_name.strip().upper() for t in transcripts if t.gene_name}
            if not self._query_gene_symbols and self.config.gene_query:
                self._query_gene_symbols = {self.config.gene_query.strip().upper()}

            # Record protein-coding transcript set for isoform coverage scoring
            protein_coding_transcripts = {
                strip_version(t.transcript_id)
                for t in transcripts
                if t.transcript_id and t.transcript_type == "protein_coding"
            }
            self._protein_coding_transcript_ids = protein_coding_transcripts
            self._protein_coding_transcript_count = len(protein_coding_transcripts)

            # Track species request: if config shows explicit request, use it; otherwise default was applied
            self._species_explicitly_requested = self.config.species_explicitly_requested

            # Variant Resolution (optional, after transcript retrieval)
            # Save resolved variants on the workflow instance for later use
            if self.config.variant_config and self.config.variant_config.has_variants:
                progress.update(main_task, description="[cyan]Resolving variants...")
                self.resolved_variants = await self.resolve_variants_step(progress)
                progress.advance(main_task)
            else:
                # Skip variant resolution step
                self.resolved_variants = []
                progress.advance(main_task)

            # Step 2: ORF Validation
            progress.update(main_task, description="[cyan]Validating ORFs...")
            orf_results = await self.step2_validate_orfs(transcripts, progress)
            progress.advance(main_task)

            # Step 3: siRNA Design
            progress.update(main_task, description="[cyan]Designing siRNAs...")
            design_results = await self.step3_design_sirnas(transcripts, progress)
            progress.advance(main_task)

            # Step 4: Off-target Analysis and Scoring. Repeat detection now runs inside this
            # step (see step5_offtarget_analysis) so it can reuse the transcriptome reference
            # already materialized for screening, instead of fetching it a second time.
            progress.update(main_task, description="[cyan]Running off-target analysis...")
            offtarget_results = await self.step5_offtarget_analysis(design_results)
            progress.advance(main_task)

            # Step 5: Generate Reports (after off-target analysis completes)
            progress.update(main_task, description="[cyan]Generating reports...")
            await self.step6_generate_reports(design_results)
            progress.advance(main_task)

        total_time = max(0.0, time.perf_counter() - start_time)

        # Compile final results
        # Serialize authoritative design parameters into the workflow summary.
        # Dumped wholesale so a newly added threshold cannot go unrecorded.
        design_parameters: dict[str, Any] = self.config.design_params.model_dump(mode="json")

        final_results: dict[str, Any] = {
            "workflow_config": {
                "gene_query": self.config.gene_query,
                "database": self.config.database.value,
                "output_dir": str(self.config.output_dir),
                "processing_time": total_time,
                "mirna_reference": {
                    "database": self.config.mirna_database,
                    "species": self.config.mirna_species,
                },
            },
            "transcript_summary": self._summarize_transcripts(transcripts),
            "transcript_annotation_summary": self._annotation_summary or {"enabled": False},
            "orf_summary": self._summarize_orf_results(orf_results),
            "design_summary": self._summarize_design_results(design_results),
            # Why the shortlist is the size it is. Beside design_summary, not inside it: this is a
            # statement about evidence and run mode. An empty shortlist from an incomplete screen and
            # one from a complete screen that found nothing eligible must not read alike.
            "selection_summary": self._selection_summary,
            "design_parameters": design_parameters,
            "repeat_summary": self._repeat_summary,
            "offtarget_summary": offtarget_results,
            "reference_summary": self._summarize_screening_references(),
        }

        # Optionally save workflow summary JSON (store in logs/)
        if self.config.write_json_summary:
            summary_file = self.config.output_dir / "logs" / "workflow_summary.json"
            with summary_file.open("w") as f:
                json.dump(final_results, f, indent=2, default=str)
            # Re-issue the Quilt registration now the summary is on disk. write_quilt_summarize omits
            # an artifact that does not exist, and step 6 runs before this write, so without this the
            # "Run manifest / Workflow summary" row collapses to manifest.json alone (#103). Re-issuing
            # rather than reordering keeps the summary's processing_time measuring the whole run, and
            # gating on the write means a write_json_summary=False or failed run never registers a path
            # it does not have.
            if summary_file.exists() and self._report_registration is not None:
                self._write_quilt_summarize(*self._report_registration)

        console.print(f"\n✅ [bold green]Workflow completed in {total_time:.2f}s[/bold green]")
        console.print(f"📊 Results saved to: [blue]{self.config.output_dir}[/blue]")

        # Persist the Rich console stream to a log file for auditing
        try:
            stream_log = self.config.output_dir / "logs" / "workflow_stream.log"
            # Append the captured console output
            with stream_log.open("a", encoding="utf-8") as lf:
                lf.write(console.export_text(clear=False))
        except Exception:
            # Do not fail the workflow if log export fails
            logger.warning("Failed to export console stream to workflow_stream.log")

        return final_results

    # ──────────────────────────────────────────────────────
    #  ZFN workflow: pair evaluation + exhaustive off-target
    # ──────────────────────────────────────────────────────

    async def _run_zfn_workflow(self) -> dict[str, Any]:
        """Execute the ZFN pair evaluation and off-target search workflow.

        ZFN is scientifically distinct from siRNA/miRNA:
        - Input is a user-provided half-site pair (not transcript FASTA)
        - Off-target is exhaustive sliding-window on whole-genome FASTA
          with FokI seed-region penalties and paired hit assembly
        - Steps 1–2 (transcript fetch, ORF validation) are skipped

                Sharding/chunking behavior is configured via ``zfn_params.sharding``.
                The search layer applies generic contig-aware planning and will avoid
                chunk sharding when the selected FASTA resolves to a single contig.
        """
        if self.zfn_designer is None:
            raise RuntimeError("ZFN designer not initialized for design_mode=zfn")

        zfn_cfg = self.config.zfn_config
        if zfn_cfg is None:
            raise RuntimeError(
                "ZFNWorkflowConfig not provided — ensure --zfn-left-half-site and --zfn-right-half-site are set"
            )

        zfn_params = zfn_cfg.zfn_params
        annotation = zfn_cfg.annotation
        annotation_source: str | None = None
        if annotation is not None:
            if annotation.annotation_path:
                annotation_source = annotation.annotation_path
            elif annotation.annotation_reference:
                annotation_source = annotation.annotation_reference

        # Covers callers who drive SiRNAWorkflow directly. Under the CLI the entry point
        # already emitted it, and the once-per-process latch keeps this from repeating.
        emit_zfn_experimental_warning(console)

        console.print("\n🧬 [bold cyan]Starting ZFN Pair Evaluation Workflow[/bold cyan]")
        console.print(f"Left half-site:  [yellow]{zfn_params.left_half_site}[/yellow]")
        console.print(f"Right half-site: [yellow]{zfn_params.right_half_site}[/yellow]")
        console.print(f"Algorithm: [blue]{zfn_params.algorithm.value}[/blue]")
        console.print(f"Dimer mode: [blue]{zfn_params.dimer_mode.value}[/blue]")
        console.print(f"Spacer lengths: [blue]{zfn_params.spacer_constraints.allowed_spacer_lengths}[/blue]")
        console.print(f"Output: [blue]{self.config.output_dir}[/blue]")

        start_time = time.perf_counter()

        with Progress(console=console) as progress:
            main_task = progress.add_task("[cyan]ZFN Workflow Progress", total=3)

            # Step 1: Evaluate pair (design + exhaustive off-target search)
            progress.update(main_task, description="[cyan]Running ZFN pair evaluation & off-target search...")
            zfn_result: ZFNDesignResult = self.zfn_designer.evaluate_pair(
                params=zfn_params,
                annotation=annotation,
            )
            progress.advance(main_task)

            # Step 2: Generate reports
            progress.update(main_task, description="[cyan]Generating ZFN reports...")
            zfn_output = self.config.output_dir / "sirnaforge"
            zfn_output.mkdir(parents=True, exist_ok=True)

            # Off-target sites CSV
            offtarget_csv = zfn_output / "offtarget_sites.csv"
            zfn_result.save_offtargets_csv(str(offtarget_csv))
            console.print(f"  Off-target sites: [green]{offtarget_csv}[/green]")

            # Candidate summary JSON
            candidate_json = zfn_output / "candidate_summary.json"
            candidate_payloads: list[dict[str, Any]] = []
            for cand in zfn_result.candidates:
                payload = cand.model_dump(mode="json")
                payload["on_target_result"] = cand.component_scores.get("on_target_quality")
                candidate_payloads.append(payload)

            candidate_data: dict[str, Any] = {
                "schema_version": "zfn_candidate_summary.v1",
                "search_contract": zfn_params.canonical_search_contract().model_dump(mode="json"),
                "candidates": candidate_payloads,
                "summary": zfn_result.get_summary(),
            }
            candidate_json.write_text(json.dumps(candidate_data, indent=2, default=str))
            console.print(f"  Candidate summary: [green]{candidate_json}[/green]")
            progress.advance(main_task)

            # Step 3: Write workflow summary
            progress.update(main_task, description="[cyan]Writing workflow summary...")
            total_time = max(0.0, time.perf_counter() - start_time)

            summary = zfn_result.get_summary()
            summary.update(
                {
                    "workflow_mode": "zfn",
                    "left_half_site": zfn_params.left_half_site,
                    "right_half_site": zfn_params.right_half_site,
                    "dimer_mode": zfn_params.dimer_mode.value,
                    "algorithm": zfn_params.algorithm.value,
                    "spacer_lengths": zfn_params.spacer_constraints.allowed_spacer_lengths,
                    "max_mismatches": zfn_params.half_site_constraints.max_mismatches,
                    "seed_len_from_foki": zfn_params.half_site_constraints.seed_len_from_fokI,
                    "seed_max_mismatches": zfn_params.half_site_constraints.seed_max_mismatches,
                    "search_backend": zfn_params.search_backend.value,
                    "search_space_reference": zfn_params.search_space_reference,
                    "search_space_fasta": zfn_params.search_space_fasta,
                    "search_space_index": zfn_params.search_space_index,
                    "annotation_source": annotation_source,
                    "total_workflow_time_s": round(total_time, 3),
                    "output_dir": str(self.config.output_dir),
                    "offtarget_csv": str(offtarget_csv),
                    "candidate_json": str(candidate_json),
                    "sharding_enabled": zfn_params.sharding.enabled,
                    "shard_chunk_size_bp": zfn_params.sharding.chunk_size_bp,
                    "shard_overlap_bp": zfn_params.sharding.overlap_bp,
                    "shard_chromosomes": zfn_params.sharding.chromosomes,
                    "shard_max_workers": zfn_params.sharding.max_workers,
                }
            )

            if zfn_result.candidates:
                cand = zfn_result.candidates[0]
                summary["composite_score"] = cand.composite_score
                summary["predicted_sites_total"] = cand.predicted_sites_total
                summary["predicted_sites_exonic"] = cand.predicted_sites_exonic
                summary["passes_filters"] = cand.passes_offtarget_filters
                summary["on_target_result"] = cand.component_scores.get("on_target_quality")

            if self.config.write_json_summary:
                log_dir = self.config.output_dir / "logs"
                log_dir.mkdir(parents=True, exist_ok=True)
                summary_path = log_dir / "workflow_summary.json"
                summary_path.write_text(json.dumps(summary, indent=2, default=str))
                console.print(f"  Workflow summary: [green]{summary_path}[/green]")

            progress.advance(main_task)

        console.print(f"\n✅ [bold green]ZFN workflow completed in {total_time:.1f}s[/bold green]")

        return summary

    async def step1_retrieve_transcripts(self, progress: Progress) -> list[TranscriptInfo]:
        """Step 1: Retrieve and validate transcript sequences."""
        task = progress.add_task("[yellow]Fetching transcripts...", total=3)
        # If an input FASTA was provided, read sequences directly and create TranscriptInfo objects
        if self.config.input_fasta:
            if self.config.input_source:
                origin = self.config.input_source
                prefix = "🌐 Downloaded" if origin.downloaded else "📂 Local"
                console.print(f"{prefix} input FASTA: [blue]{origin.original}[/blue]")
            sequences = FastaUtils.read_fasta(self.config.input_fasta)
            progress.advance(task)

            transcripts: list[TranscriptInfo] = []
            for header, seq in sequences:
                # header may contain transcript id and metadata; use first token as id
                tid = header.split()[0]
                transcripts.append(
                    TranscriptInfo(
                        transcript_id=tid,
                        transcript_name=None,
                        transcript_type="unknown",
                        gene_id=self.config.gene_query,
                        gene_name=self.config.gene_query,
                        sequence=seq,
                        length=len(seq),
                        database=self.config.database,
                    )
                )

            # Save a normalized transcripts FASTA in the output directory
            transcript_file = self.config.output_dir / "transcripts" / f"{self.config.gene_query}_transcripts.fasta"
            sequences_out = [(f"{t.transcript_id} {t.gene_name}", t.sequence or "") for t in transcripts]
            FastaUtils.save_sequences_fasta(sequences_out, transcript_file)
            progress.advance(task)

            console.print(f"📄 Loaded {len(transcripts)} sequences from FASTA: {self.config.input_fasta}")

            # Quiet transcript validation (no verbose console warnings)
            _ = self.validation.validate_transcripts(transcripts)

            return transcripts

        # Otherwise perform a gene search
        gene_result = await self.gene_searcher.search_gene(
            self.config.gene_query, self.config.database, include_sequence=True
        )
        progress.advance(task)

        if not gene_result.success:
            raise ValueError(f"No results found for gene '{self.config.gene_query}' in {self.config.database}")

        # Get transcripts
        transcripts = gene_result.transcripts
        progress.advance(task)

        # Filter for protein-coding transcripts
        protein_transcripts = [t for t in transcripts if t.transcript_type == "protein_coding" and t.sequence]

        if not protein_transcripts:
            raise ValueError("No protein-coding transcripts found with sequences")

        # Save transcripts to file
        transcript_file = self.config.output_dir / "transcripts" / f"{self.config.gene_query}_transcripts.fasta"
        sequences = [
            (f"{t.transcript_id} {t.gene_name} type:{t.transcript_type} length:{t.length}", t.sequence or "")
            for t in protein_transcripts
            if t.sequence is not None
        ]

        FastaUtils.save_sequences_fasta(sequences, transcript_file)
        progress.advance(task)

        console.print(f"📄 Retrieved {len(protein_transcripts)} protein-coding transcripts")

        # If canonical transcripts are present, save them separately
        canonical_transcripts = [t for t in transcripts if getattr(t, "is_canonical", False) and t.sequence]
        if canonical_transcripts:
            canonical_file = self.config.output_dir / "transcripts" / f"{self.config.gene_query}_canonical.fasta"
            canonical_sequences = [
                (
                    f"{t.transcript_id} {t.gene_name} type:{t.transcript_type} length:{t.length} canonical:true",
                    t.sequence or "",
                )
                for t in canonical_transcripts
                if t.sequence is not None
            ]
            FastaUtils.save_sequences_fasta(canonical_sequences, canonical_file)
            console.print(f"⭐ Canonical transcripts saved: {canonical_file.name}")

        # Quiet transcript validation (no verbose console warnings)
        _ = self.validation.validate_transcripts(protein_transcripts)

        # Optional: Enrich with genomic annotations if client is available
        await self._enrich_transcript_annotations(protein_transcripts)

        return protein_transcripts

    async def _enrich_transcript_annotations(self, transcripts: list[TranscriptInfo]) -> None:
        """Optionally enrich transcripts with genomic annotations.

        This is a non-breaking enhancement that fetches additional genomic metadata
        for transcripts when the annotation client is available.
        Results are logged to workflow summary but do not modify transcript objects.
        """
        if not self._annotation_client:
            return

        # Only try to annotate if we have Ensembl transcript IDs
        transcript_ids = [t.transcript_id for t in transcripts if t.transcript_id.startswith("ENST")]
        if not transcript_ids:
            return

        try:
            # Use default reference for annotation
            reference = ReferenceChoice.default("GRCh38", reason="auto-selected for annotation")

            bundle = await self._annotation_client.fetch_by_ids(
                ids=transcript_ids[:10],  # Limit to first 10 to avoid excessive API calls
                species="human",
                reference=reference,
            )

            # Store summary for workflow output
            self._annotation_summary = {
                "enabled": True,
                "provider": "ensembl_rest",
                "transcripts_queried": len(transcript_ids[:10]),
                "transcripts_resolved": bundle.resolved_count,
                "transcripts_unresolved": bundle.unresolved_count,
                "reference": reference.to_metadata(),
            }

            if bundle.resolved_count > 0:
                console.print(
                    f"📊 Genomic annotations: {bundle.resolved_count}/{len(transcript_ids[:10])} transcripts enriched"
                )
        except Exception as e:
            logger.debug(f"Transcript annotation enrichment failed (non-critical): {e}")
            self._annotation_summary = {"enabled": False, "error": str(e)}

    async def resolve_variants_step(self, progress: Progress) -> list[VariantRecord]:
        """Resolve variants for targeting or avoidance (optional workflow step).

        This step runs after transcript retrieval and before siRNA design,
        resolving and filtering variants based on the workflow configuration.

        This step is run after transcript retrieval and before ORF validation and siRNA design.
        Variants are resolved using ClinVar, Ensembl Variation, and/or VCF files.

        Args:
            progress: Rich progress tracker

        Returns:
            List of resolved VariantRecords that passed filters
        """
        if not self.config.variant_config or not self.config.variant_config.has_variants:
            return []

        task = progress.add_task("[yellow]Resolving variants...", total=2)

        # Resolve variants using the workflow variant module
        variants = await resolve_workflow_variants(
            config=self.config.variant_config,
            gene_name=self.config.gene_query,
            output_dir=self.config.output_dir,
        )
        progress.advance(task)

        if variants:
            console.print(
                f"🧬 Resolved {len(variants)} variant(s) for {self.config.variant_config.variant_mode.value} mode"
            )
            for variant in variants[:5]:  # Show first 5
                console.print(f"  • {variant.id or variant.to_vcf_style()}")
            if len(variants) > 5:
                console.print(f"  ... and {len(variants) - 5} more")
        else:
            console.print("⚠️  No variants passed filters")

        progress.advance(task)
        return variants

    async def step2_validate_orfs(self, transcripts: list[TranscriptInfo], progress: Progress) -> dict[str, Any]:
        """Step 2: Validate ORFs and generate validation report."""
        task = progress.add_task("[yellow]Analyzing ORFs...", total=len(transcripts) + 1)

        orf_results: dict[str, Any] = {}
        valid_transcripts: list[TranscriptInfo] = []

        for transcript in transcripts:
            try:
                analysis = await self.orf_analyzer.analyze_transcript(transcript)
                orf_results[transcript.transcript_id] = analysis

                if analysis.has_valid_orf:
                    valid_transcripts.append(transcript)

                progress.advance(task)

            except Exception as e:
                logger.warning(f"ORF analysis failed for {transcript.transcript_id}: {e}")
                progress.advance(task)

        # Generate ORF validation report
        report_file = self.config.output_dir / "orf_reports" / "orf_validation.txt"
        self._generate_orf_report(orf_results, report_file)
        progress.advance(task)

        console.print(f"🔍 ORF validation: {len(valid_transcripts)}/{len(transcripts)} transcripts have valid ORFs")
        return {"results": orf_results, "valid_transcripts": valid_transcripts}

    async def step3_design_sirnas(self, transcripts: list[TranscriptInfo], progress: Progress) -> DesignResult:
        """Step 3: Design siRNA candidates for valid transcripts.

        Parallelizes per-transcript design when not running from a user-provided input FASTA,
        to preserve backward-compatibility with tests and monkeypatching of design_from_file.
        Set env SIRNAFORGE_PARALLEL_DESIGN=1 to force parallel mode.

        The two branches have inverse failure modes, and only the parallel one needs a record
        (#100). The single-call branch has no try/except at all: a design failure propagates and the
        run stops, so no result can claim a target set the designer never covered. The parallel
        branch catches a failed batch (and, inside the batch, a failed transcript) and carries on,
        which is the right call for a 40-transcript gene -- but the loss has to be written down, so
        every dropped transcript id goes to ``_design_input_shortfalls`` and from there into
        ``design_summary`` and ``SelectionInputs.design_input_shortfalls``.
        """
        # Create temporary FASTA file for siRNA design (preserves original behavior)
        temp_fasta = self.config.output_dir / "transcripts" / "temp_for_design.fasta"
        sequences = [(f"{t.transcript_id}", t.sequence) for t in transcripts if t.sequence]
        FastaUtils.save_sequences_fasta(sequences, temp_fasta)

        use_parallel = (self.config.input_fasta is None) or (os.getenv("SIRNAFORGE_PARALLEL_DESIGN", "0") == "1")

        if not use_parallel:
            # Original single-call path (compatible with tests that patch design_from_file)
            task = progress.add_task("[yellow]Designing siRNAs...", total=2)
            progress.advance(task)
            assert self.sirnaforgeer is not None, "designer not initialised for siRNA/miRNA mode"
            design_result = self.sirnaforgeer.design_from_file(str(temp_fasta))
            self._store_guide_to_transcripts(self.sirnaforgeer.last_guide_to_transcripts)
            added_controls = inject_dirty_controls(design_result)
            self._dirty_controls_added = len(added_controls)
            if added_controls:
                console.print(
                    f"🧪 Added {len(added_controls)} {DIRTY_CONTROL_LABEL} candidates for signal verification"
                )
            progress.advance(task)
            _ = self.validation.validate_design_results(design_result)
            temp_fasta.unlink(missing_ok=True)
            console.print(f"🎯 Generated {len(design_result.candidates)} siRNA candidates")
            console.print(f"   Top {len(design_result.top_candidates)} candidates selected for further analysis")
            return design_result

        # Parallel per-transcript path
        start = time.perf_counter()
        total = len(sequences)
        task = progress.add_task("[yellow]Designing siRNAs...", total=total if total > 0 else 1)

        results: list[DesignResult] = []
        guide_to_transcripts: dict[str, set[str]] = {}

        # Batch transcripts for more efficient threading
        transcript_batches = self._batch_transcripts(transcripts)

        with ThreadPoolExecutor(max_workers=self.config.num_threads) as executor:
            futures = {executor.submit(self._process_transcript_batch, batch): batch for batch in transcript_batches}

            for fut in as_completed(futures):
                try:
                    batch_results, batch_guide_mapping, batch_shortfalls = fut.result()
                    results.extend(batch_results)
                    # Merge guide-to-transcript mappings
                    for guide_seq, transcript_set in batch_guide_mapping.items():
                        guide_to_transcripts.setdefault(guide_seq, set()).update(transcript_set)
                    # Transcripts the batch itself lost, recorded here rather than inside the worker so
                    # only this thread ever writes the run's shortfall record.
                    self._record_design_input_shortfalls(batch_shortfalls)
                    # Advance progress by the number of transcripts in this batch
                    batch = futures[fut]
                    progress.advance(task, len(batch))
                except Exception as e:
                    batch = futures[fut]
                    batch_transcript_ids = [t.transcript_id for t in batch]
                    logger.exception(f"Design failed for transcript batch {batch_transcript_ids}: {e}")
                    # The whole batch's results are discarded with the exception, so every transcript
                    # in it is unrepresented -- including any that had already been designed.
                    self._record_design_input_shortfalls(
                        dict.fromkeys(batch_transcript_ids, _describe_design_failure(e))
                    )
                    progress.advance(task, len(batch))

        self._store_guide_to_transcripts(guide_to_transcripts)

        if self._design_input_shortfalls:
            # Said out loud at the step that lost them, not only in the JSON: the candidate counts
            # printed two lines below are over the transcripts that survived, and nothing else on the
            # console distinguishes a gene designed whole from one designed in part.
            dropped = sorted(self._design_input_shortfalls)
            console.print(
                f"⚠️  {len(dropped)} of {total} transcript(s) were dropped by a design failure and are "
                f"covered by no candidate: {', '.join(dropped)}"
            )
            logger.error(f"Design input shortfalls: {self._design_input_shortfalls}")

        # Merge candidates
        all_candidates: list[SiRNACandidate] = [c for dr in results for c in dr.candidates]
        rejected_pool: list[SiRNACandidate] = [c for dr in results for c in getattr(dr, "rejected_candidates", [])]

        # Recompute transcript hit metrics across all inputs
        total_seqs = total
        for c in all_candidates:
            hits = len(guide_to_transcripts.get(c.guide_sequence, {c.transcript_id}))
            c.transcript_hit_count = hits
            c.transcript_hit_fraction = (hits / total_seqs) if total_seqs > 0 else 0.0

        # Sort, compute top-N (prefer passing candidates)
        all_candidates.sort(key=ranking_score, reverse=True)
        passing = [
            c
            for c in all_candidates
            if (c.passes_filters is True)
            or (
                hasattr(_ModelSiRNACandidate, "FilterStatus")
                and c.passes_filters == _ModelSiRNACandidate.FilterStatus.PASS
            )
        ]
        top_candidates = (passing or all_candidates)[: self.config.top_n]

        processing_time = max(0.0, time.perf_counter() - start)
        filtered_count = len(passing)
        tool_versions = results[0].tool_versions if results else {}

        combined = DesignResult(
            input_file="<parallel_transcripts>",
            parameters=self.config.design_params,
            candidates=all_candidates,
            top_candidates=top_candidates,
            total_sequences=total_seqs,
            total_candidates=len(all_candidates),
            filtered_candidates=filtered_count,
            processing_time=processing_time,
            tool_versions=tool_versions,
            rejected_candidates=rejected_pool,
        )

        added_controls = inject_dirty_controls(combined)
        self._dirty_controls_added = len(added_controls)
        if added_controls:
            console.print(f"🧪 Added {len(added_controls)} {DIRTY_CONTROL_LABEL} candidates for signal verification")

        _ = self.validation.validate_design_results(combined)

        temp_fasta.unlink(missing_ok=True)
        console.print(f"🎯 Generated {len(combined.candidates)} siRNA candidates (threads={self.config.num_threads})")
        console.print(f"   Top {len(combined.top_candidates)} candidates selected for further analysis")
        return combined

    def _store_guide_to_transcripts(self, guide_to_transcripts: dict[str, set[str]] | None) -> None:
        """Cache a guide -> source-transcripts mapping for post-screen isoform coverage scoring.

        `design_from_sequence` and the miRNA designer never populate this, so absence here
        (an empty guide_to_transcripts, or one that is None) leaves the coverage term inactive
        for those candidates rather than computing it from the wrong numerator (see D4/D8).

        Args:
            guide_to_transcripts: Raw mapping from guide sequence to source transcript IDs, as
                produced by SiRNADesigner.design_from_file or the parallel per-transcript path.
        """
        if not guide_to_transcripts:
            return
        self._guide_to_transcripts = {
            normalize_guide_sequence(guide): frozenset(strip_version(tid) for tid in tids)
            for guide, tids in guide_to_transcripts.items()
        }

    def _batch_transcripts(
        self, transcripts: list[TranscriptInfo], batch_size: int | None = None
    ) -> list[list[TranscriptInfo]]:
        """Group transcripts into batches for more efficient threading.

        Args:
            transcripts: List of transcripts to batch
            batch_size: Number of transcripts per batch. If None, automatically calculated
                       based on transcript lengths to aim for ~2 seconds of work per batch.

        Returns:
            List of transcript batches
        """
        if batch_size is None:
            # Estimate batch size based on transcript lengths
            total_length = sum(len(t.sequence or "") for t in transcripts)
            if total_length == 0 or len(transcripts) == 0:
                batch_size = 1
            else:
                avg_length = total_length / len(transcripts)
                # Rough estimate: aim for batches with ~2000 candidates each
                # (1000bp transcript ≈ 980 candidates for 21nt siRNAs)
                target_candidates_per_batch = 2000
                # Subtract siRNA length from transcript length when estimating candidates
                # (default siRNA length is 21nt, but use configured value)
                sirna_length = self.config.design_params.sirna_length
                batch_size = max(1, int(target_candidates_per_batch / max(1, avg_length - sirna_length + 1)))
                # Cap batch size to avoid memory issues
                batch_size = min(batch_size, 20)

        batches: list[list[TranscriptInfo]] = []
        for i in range(0, len(transcripts), batch_size):
            batch = transcripts[i : i + batch_size]
            if batch:  # Only add non-empty batches
                batches.append(batch)

        return batches

    def _process_transcript_batch(
        self, batch: list[TranscriptInfo]
    ) -> tuple[list[DesignResult], dict[str, set[str]], dict[str, str]]:
        """Process a batch of transcripts and return results, guide mapping and what it lost.

        Runs in a worker thread, so it *returns* its shortfalls rather than writing them onto the
        workflow: the caller merges them, and only the collecting thread mutates run state. A
        transcript whose design raised here is dropped exactly as a whole failed batch is, so it is
        the same loss and gets the same record (#100).

        Args:
            batch: List of transcripts to process

        Returns:
            Tuple of (design_results, guide_to_transcripts_mapping, dropped_transcript_id -> reason)
        """
        results: list[DesignResult] = []
        guide_to_transcripts: dict[str, set[str]] = {}
        shortfalls: dict[str, str] = {}

        for transcript in batch:
            if not transcript.sequence:
                continue

            try:
                assert self.sirnaforgeer is not None, "designer not initialised"
                dr = self.sirnaforgeer.design_from_sequence(transcript.sequence, transcript.transcript_id)
                results.append(dr)

                # Build guide-to-transcript mapping for this batch
                for c in dr.candidates:
                    guide_to_transcripts.setdefault(c.guide_sequence, set()).add(c.transcript_id)

            except Exception as e:
                logger.exception(f"Design failed for transcript {transcript.transcript_id}: {e}")
                shortfalls[transcript.transcript_id] = _describe_design_failure(e)
                continue

        return results, guide_to_transcripts, shortfalls

    def _record_design_input_shortfalls(self, shortfalls: Mapping[str, str]) -> None:
        """Record transcripts a design failure dropped, so a lost target cannot read as a designed one.

        First reason wins per transcript, matching ``record_filter_verdict``'s rule: the per-transcript
        reason a surviving batch reports is more specific than any later re-statement, and a transcript
        belongs to exactly one batch, so the two records here cannot describe the same loss twice.
        """
        for transcript_id, reason in shortfalls.items():
            self._design_input_shortfalls.setdefault(transcript_id, reason)

    def _apply_modifications_to_results(self, design_results: DesignResult) -> None:
        """Apply chemical modification patterns to all candidates in design results.

        Args:
            design_results: DesignResult containing candidates to modify
        """
        pattern = self.config.design_params.modification_pattern
        overhang = self.config.design_params.default_overhang

        # Apply modifications to all candidates
        for candidate in design_results.candidates:
            apply_modifications_to_candidate(
                candidate,
                pattern_name=pattern,
                overhang=overhang,
                target_gene=self.config.gene_query,
            )

        console.print(f"✨ Applied {pattern} modification pattern with {overhang} overhangs to all candidates")

    def _run_repeat_detection(self, candidates: list[SiRNACandidate]) -> dict[str, Any]:
        """Detect repeat elements in candidate guide sequences before scoring.

        Scans all distinct guide sequences against the query species' cDNA reference to flag
        guides overlapping repeat elements. Costs ~47s for a ~1GB human reference. Runs from
        within step5_offtarget_analysis (after screening but before scoring) so it can reuse
        the transcriptome reference already materialized for off-target screening rather than
        fetching a second, redundant copy just to find the query species' FASTA. Its caller
        does not reach it at all when the user disabled off-target analysis, because there is
        then no reference to reuse.
        """
        distinct_guides = {normalize_guide_sequence(c.guide_sequence) for c in candidates}
        # The RESOLVED threshold, not the module constant. They agree by default, so a user who set
        # `max_repeat_transcript_fraction` got the constant silently -- and the run summary published
        # the value they asked for beside a scan that never used it.
        threshold = self.config.design_params.filters.max_repeat_transcript_fraction
        action = self._filter_action_for("max_repeat_transcript_fraction", threshold)

        query_cdna_fasta = self._species_cdna_fasta.get(self._query_species)
        if not query_cdna_fasta:
            console.print("⚠️  Query species cDNA reference not available; skipping repeat detection")
            return {
                "status": "skipped",
                "reason": "reference_unavailable",
                "query_species": self._query_species,
                "repeat_flagged_count": 0,
                "threshold_fraction": threshold,
            }

        console.print(
            f"🔍 Scanning {len(distinct_guides)} distinct guides against {self._query_species} cDNA (~47s)..."
        )
        detector = RepeatDetector(threshold_fraction=threshold)
        scan_result = detector.scan(distinct_guides, query_cdna_fasta)

        # Stamp repeat verdicts on candidates
        observations = scan_result.observations
        for candidate in candidates:
            SiRNADesigner.stamp_repeat_verdict(candidate, observations, action)

        repeat_flagged_count = sum(1 for c in candidates if c.repeat_flagged)
        console.print(
            f"🔍 Repeat detection: {repeat_flagged_count}/{len(candidates)} candidates flagged "
            f"(threshold: {scan_result.threshold_fraction:.1%}, reference: {scan_result.reference_transcript_count} transcripts)"
        )

        return {
            "status": "completed",
            "query_species": self._query_species,
            "distinct_guides_scanned": len(distinct_guides),
            "reference_transcript_count": scan_result.reference_transcript_count,
            "threshold_fraction": scan_result.threshold_fraction,
            "repeat_flagged_count": repeat_flagged_count,
            "repeat_sequences": list(scan_result.repeat_sequences),
        }

    async def step6_generate_reports(self, design_results: DesignResult) -> None:  # noqa: C901, PLR0912
        """Step 6: Generate comprehensive reports.

        Four candidate tables, and the difference between them is the point (#100). ``candidates_all``
        is every row; ``candidates_pass`` is the gate verdict and nothing else, kept exactly as it was
        because it is the compat deliverable (narrowing it once deleted the whole output of an
        unscreened run); ``candidates_qualified`` and ``candidates_provisional`` are the *selection*,
        named as files so "which of these can I order?" does not depend on the reader knowing to filter
        a column.
        """
        # No user-facing top-candidates FASTA or text/json summaries are produced anymore.
        # We only keep canonical CSV outputs (ALL + PASS) for candidates. Off-target analysis
        # prepares its own internal FASTA input under off_target/.

        if self.config.design_params.apply_modifications:
            self._apply_modifications_to_results(design_results)

        base = self.config.output_dir / "sirnaforge"
        tables = CandidateTablePaths.under(base)
        report_file = self.config.output_dir / "orf_reports" / "orf_validation.txt"

        try:
            self._write_candidate_tables(design_results.candidates, tables)
        except Exception as e:  # Do not fail workflow for reporting extras
            logger.warning(f"Failed to write all/pass CSVs: {e}")

        try:
            variant_links_path = self.config.output_dir / "logs" / "candidate_variants.json"
            self._write_candidate_variant_links(design_results.candidates, variant_links_path)
        except Exception as e:
            logger.warning(f"Failed to write candidate variant links: {e}")

        manifest_path = tables.manifest_json
        report_html_path = base / "report.html"
        try:
            manifest = self._build_fair_manifest(
                all_csv=tables.all_csv,
                pass_csv=tables.pass_csv,
                pass_fasta=tables.pass_fasta,
                orf_report=report_file,
                qualified_csv=tables.qualified_csv,
                qualified_fasta=tables.qualified_fasta,
                provisional_csv=tables.provisional_csv,
                report_html=report_html_path,
            )
            with manifest_path.open("w") as mf:
                json.dump(manifest, mf, indent=2)
        except Exception as e:
            logger.warning(f"Failed to write FAIR manifest: {e}")

        # Rendering the report and registering it are two writes with two failure modes, so they get
        # two handlers: a read-only run root left report.html written and quilt_summarize.json not,
        # and one shared handler then logged "Failed to write self-contained HTML report" for a report
        # that had in fact succeeded (#103).
        payload: ReportPayload | None = None
        report_path: Path | None = None
        render_error: str | None = None
        try:
            # Pass the resolved policy: rediscovering it from the manifest works, but this run holds
            # the gates it actually applied, and a default panel would report other thresholds.
            payload = build_payload(self.config.output_dir, policy=self.config.resolved_policy)
            report_path = write_report(payload, report_html_path)
        except Exception as e:
            render_error = f"{type(e).__name__}: {e}"
            logger.warning(f"Failed to write self-contained HTML report: {e}")

        summarize_path: Path | None = None
        if payload is not None and report_path is not None:
            self._report_registration = (payload, report_path)
            summarize_path = self._write_quilt_summarize(payload, report_path)

        # Unconditional, because manifest.json above claims report.html only as PENDING and the sidecar
        # is the one place the outcome can be recorded without rewriting the bytes the report quoted. A
        # sidecar written only on success left the manifest attesting a report.html a failed render
        # never wrote. Its own handler, per #103's two-handler rule: a failed sidecar must not be logged
        # as a failed report.
        try:
            write_report_manifest(
                report_path or report_html_path,
                manifest_path,
                base / "report_manifest.json",
                render_error=render_error,
            )
        except Exception as e:
            logger.warning(f"Failed to write the report manifest sidecar: {e}")

        console.print("📋 Generated comprehensive reports and FAIR metadata")
        # Each line only for an artifact that is really there. An operator told a file exists when it
        # does not has to discover the gap from the catalog instead (#103), and candidates_pass.fasta
        # is deliberately absent when no candidate passed.
        if report_file.exists():
            console.print("   - ORF validation report: orf_reports/")
        if tables.all_csv.exists() or tables.pass_csv.exists():
            console.print("   - siRNA candidate CSVs: sirnaforge/ (candidates_all.csv, candidates_pass.csv)")
        if tables.pass_fasta.exists():
            console.print("   - siRNA candidate FASTA: sirnaforge/ (candidates_pass.fasta)")
        # The selection exports are #100's, and they get the same treatment: a qualified run with nothing
        # eligible writes no qualified file, and saying otherwise would advertise an empty shortlist.
        if tables.qualified_csv.exists() or tables.qualified_fasta.exists() or tables.provisional_csv.exists():
            console.print(
                "   - Selection exports: sirnaforge/ (candidates_qualified.csv/.fasta, candidates_provisional.csv)"
            )
        if report_path is not None and report_path.exists():
            console.print("   - Self-contained HTML report: sirnaforge/report.html")
        if summarize_path is not None and summarize_path.exists():
            console.print("   - Quilt package summary: quilt_summarize.json")

    def _write_quilt_summarize(self, payload: ReportPayload, report_path: Path) -> Path | None:
        """Register the run's artifacts where Quilt reads them; return the path, or ``None`` if it failed.

        At the run root, not beside the report: Quilt reads a summarize file only at the package root,
        and the run directory is what gets published, so a report written without this shows nothing at
        all in a package view (#103). Called again after ``logs/workflow_summary.json`` lands, because
        the writer registers only artifacts that exist by the time it runs.
        """
        out_path = Path(self.config.output_dir) / "quilt_summarize.json"
        try:
            return write_quilt_summarize(payload, report_path, self.config.output_dir, out_path=out_path)
        except Exception as e:
            logger.warning(f"Failed to write Quilt package summary: {e}")
            return None

    def _build_candidate_frame(self, candidates: Sequence[SiRNACandidate]) -> pd.DataFrame:
        """The validated candidate table every export slices, from the one shared row builder.

        Shared by step6 and by #100's off-target-only exports, so the two entry points cannot drift on
        which columns they publish or on how ``passes_filters`` is spelled -- the reason the off-target
        entry point could not simply grow its own writer.

        Two dtype repairs before validation, both cases of a column that is *entirely* null:
        ``seed_7mer_hits``/``seed_8mer_hits`` become ``Int64``, and a float-typed column arrives as
        ``object`` because pandas has no float to infer from. The schema's class-level ``coerce`` never
        reaches those fields, so an all-null ``mfe`` fails validation outright. A column with even one
        value present already infers ``float64``, so this repair can only ever fire on a column the run
        measured nowhere.

        Raises:
            pandera.errors.SchemaError: The rows do not satisfy ``SiRNACandidateSchema``. Callers that
                treat a missing table as a reporting extra catch this themselves.
        """
        # Single shared row-builder (models/sirna.py) so this CSV and DesignResult.save_csv
        # can never drift on which columns they emit (issue #80 F2).
        rows: list[dict[str, Any]] = [build_candidate_row(candidate) for candidate in candidates]

        if rows:
            all_df = pd.DataFrame(rows)
        else:
            template_cols = list(SiRNACandidateSchema.to_schema().columns.keys())
            all_df = pd.DataFrame(columns=template_cols)

        for col in ("seed_7mer_hits", "seed_8mer_hits"):
            if col in all_df.columns:
                all_df[col] = all_df[col].astype("Int64")

        for name, column in SiRNACandidateSchema.to_schema().columns.items():
            dtype = str(column.dtype)
            if not dtype.startswith("float") or name not in all_df.columns:
                continue
            if all_df[name].dtype == object and all_df[name].isna().all():
                all_df[name] = all_df[name].astype(dtype)

        if "passes_filters" not in all_df.columns:
            all_df["passes_filters"] = pd.Series(dtype="object")

        validated_all = SiRNACandidateSchema.validate(all_df)

        if self.config.design_params.design_mode != DesignMode.MIRNA:
            mirna_cols = [
                "guide_pos1_base",
                "pos1_pairing_state",
                "seed_class",
                "supp_13_16_score",
                "seed_7mer_hits",
                "seed_8mer_hits",
                "seed_hits_weighted",
                "off_target_seed_risk_class",
            ]
            existing = [col for col in mirna_cols if col in validated_all.columns]
            if existing:
                validated_all = validated_all.drop(columns=existing)

        def _normalize_pass(value: Any) -> str:
            normalized = "FAIL"
            try:
                if value is True or (isinstance(value, int | float) and value == 1):
                    normalized = "PASS"
                elif value is False or (isinstance(value, int | float) and value == 0):
                    normalized = "FAIL"
                elif isinstance(value, str):
                    cleaned = value.strip().upper()
                    if cleaned in {"PASS", "TRUE", "YES"}:
                        normalized = "PASS"
                    elif cleaned in {"FAIL", "FALSE", "NO"}:
                        normalized = "FAIL"
                    else:
                        normalized = cleaned
                else:
                    normalized = "PASS" if bool(value) else "FAIL"
            except Exception:
                normalized = "FAIL"
            return normalized

        validated_all["passes_filters"] = [_normalize_pass(value) for value in validated_all["passes_filters"]]
        return cast(pd.DataFrame, validated_all)

    def _write_selection_exports(
        self,
        validated_all: pd.DataFrame,
        *,
        qualified_csv: Path,
        qualified_fasta: Path,
        provisional_csv: Path,
    ) -> None:
        """Write the two selection-named candidate tables from the resolved selection (#100).

        Both are slices of the already-validated all-candidates frame, in its order -- which is the
        selection's own ranking, because :meth:`_apply_post_screen_ranking` re-sorted the candidate list
        step6 built these rows from. Re-filtering on something else, as ``candidates_pass.csv`` does on
        ``passes_filters``, is what let a withheld candidate lead a "pass" list.

        The CSVs are written even when empty, header and all: a qualified run that could qualify nobody
        must publish that answer, and an absent file reads as "this version does not report it".
        The FASTA is removed instead, mirroring ``candidates_pass.fasta`` -- an empty FASTA is not a
        valid order list, and a stale one from a previous run in the same directory would be worse.
        """
        if "selection_state" not in validated_all.columns:
            # An older row builder, or a frame that lost the column: say nothing rather than guess.
            logger.warning("No selection_state column: candidates_qualified/provisional not written")
            return

        qualified_df = validated_all[validated_all["selection_state"] == SelectionState.ELIGIBLE.value].copy()
        provisional_df = validated_all[validated_all["selection_state"] == SelectionState.PROVISIONAL.value].copy()
        qualified_df.to_csv(qualified_csv, index=False)
        provisional_df.to_csv(provisional_csv, index=False)

        if qualified_df.empty:
            qualified_fasta.unlink(missing_ok=True)
        else:
            try:
                self._write_pass_candidates_fasta(qualified_df, qualified_fasta)
            except Exception as e:
                logger.warning(f"Failed to write qualified candidates FASTA: {e}")

    def _write_pass_candidates_fasta(self, pass_df: pd.DataFrame, output_path: Path) -> None:
        """Write passing candidates to FASTA, each header naming its selection state.

        This file leaves the tool as a list of sequences to order, so a guide the run could not qualify
        must not look like one it could. On the header rather than filtered out: filtering deleted the
        whole deliverable of an unscreened run (152 guides to 0 on a measured ``--input-fasta`` run), and
        a reader who splits on the state gets the narrower list whenever they want it.

        Also writes ``candidates_qualified.fasta`` (#100), whose rows are already narrowed to the
        eligible selection: one header format for both files, so the two order lists cannot disagree on
        what a header means.

        Args:
            pass_df: DataFrame of the candidates to write, in the order to write them
            output_path: Path to write the FASTA file
        """
        try:
            sequences: list[tuple[str, str]] = []
            for _, row in pass_df.iterrows():
                # Header carries whichever score the row actually has: composite_score is null
                # until screening, and a design-only run has design_score instead.
                score = row.get("composite_score")
                if score is None or pd.isna(score):
                    score = row.get("design_score")
                header = (
                    f"{row['id']} score={score:.1f}" if score is not None and not pd.isna(score) else str(row["id"])
                )
                state = row.get("selection_state")
                if isinstance(state, str) and state and state != SelectionState.ELIGIBLE.value:
                    header = f"{header} selection={state}"
                # A gate the run was told to only *warn* about still failed, and this file is an order
                # list -- so the header says which. Without it a guide outside a GC window the user
                # widened to `warn` sat here indistinguishable from one inside it.
                warned = sorted(
                    str(key)[: -len("_verdict")]
                    for key, value in row.items()
                    if str(key).endswith("_verdict") and value == FilterEvaluation.FAIL.value
                )
                if warned and row.get("passes_filters") == "PASS":
                    header = f"{header} warned={','.join(warned)}"
                sequence = str(row["guide_sequence"])
                sequences.append((header, sequence))

            # Use FastaUtils to write the sequences
            FastaUtils.save_sequences_fasta(sequences, output_path)
            logger.info(f"Saved {len(sequences)} candidates to FASTA: {output_path}")

        except Exception as e:
            logger.error(f"Failed to write PASS candidates FASTA: {e}")
            raise

    def _write_candidate_variant_links(self, candidates: Sequence[Any], output_path: Path) -> None:
        """Persist mapping between candidates and overlapped variants for observability."""
        entries: list[dict[str, Any]] = []
        for candidate in candidates:
            overlapped: list[Any] = list(cast(Sequence[Any], getattr(candidate, "overlapped_variants", None) or []))
            if not overlapped:
                continue
            entry: dict[str, Any] = {
                "id": getattr(candidate, "id", None),
                "transcript_id": getattr(candidate, "transcript_id", None),
                "variant_mode": getattr(candidate, "variant_mode", None),
                "allele_specific": bool(getattr(candidate, "allele_specific", False)),
                "targeted_alleles": list(getattr(candidate, "targeted_alleles", [])),
                "overlapped_variants": overlapped,
            }
            entries.append(entry)

        payload: dict[str, Any] = {
            "gene": self.config.gene_query,
            "total_candidates": len(candidates),
            "variant_annotated_candidates": len(entries),
            "candidates": entries,
        }

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w") as fh:
            json.dump(payload, fh, indent=2)
        logger.info(f"Wrote candidate variant links to {output_path}")

    def _count_fasta_sequences(self, path: Path) -> int:
        try:
            # Simple FASTA count: lines starting with '>'
            with path.open("r") as fh:
                return sum(1 for line in fh if line.startswith(">"))
        except Exception:
            return 0

    def _build_fair_manifest(
        self,
        *,
        all_csv: Path,
        pass_csv: Path,
        pass_fasta: Path,
        orf_report: Path | None = None,
        qualified_csv: Path | None = None,
        qualified_fasta: Path | None = None,
        provisional_csv: Path | None = None,
        report_html: Path | None = None,
    ) -> dict[str, Any]:
        """Create a manifest JSON describing generated outputs (checksums, sizes, counts).

        The three selection exports (#100) default to ``None`` and are then absent from ``files``
        entirely, which is the honest record for a caller that published no selection: an
        ``exists: false`` entry would claim the run tried to write one. ``orf_report`` follows the same
        rule for the off-target-only entry point, which validates no ORF and so never asks for one.

        ``report_html`` is the path the caller is *about* to render, and is deliberately not added to
        ``files``, where every key means "digested": the report is rendered after this manifest is
        written and embeds its content, so it is attested by a sidecar instead. A caller that renders
        no report passes ``None``, and the artifact record says so rather than reading "not attested".
        """
        now = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}"
        files: dict[str, dict[str, Any]] = {}

        def add_file(key: str, p: Path, ftype: str, extra: dict[str, Any] | None = None) -> None:
            if not p.exists():
                files[key] = {"path": str(p), "type": ftype, "exists": False}
                return
            entry: dict[str, Any] = {
                "path": str(p),
                "type": ftype,
                "exists": True,
                "size_bytes": p.stat().st_size,
                "sha256": file_sha256(p),
            }
            if extra:
                entry.update(extra)
            files[key] = entry

        # Row counts for CSVs
        def csv_rows(p: Path) -> int:
            try:
                # subtract header if file has at least one line
                with p.open("r") as fh:
                    lines = sum(1 for _ in fh)
                return max(0, lines - 1)
            except Exception:
                return 0

        add_file("candidates_all_csv", all_csv, "csv", {"rows": csv_rows(all_csv)})
        add_file("candidates_pass_csv", pass_csv, "csv", {"rows": csv_rows(pass_csv)})
        add_file("candidates_pass_fasta", pass_fasta, "fasta", {"sequences": self._count_fasta_sequences(pass_fasta)})
        if qualified_csv is not None:
            add_file("candidates_qualified_csv", qualified_csv, "csv", {"rows": csv_rows(qualified_csv)})
        if qualified_fasta is not None:
            add_file(
                "candidates_qualified_fasta",
                qualified_fasta,
                "fasta",
                {"sequences": self._count_fasta_sequences(qualified_fasta)},
            )
        if provisional_csv is not None:
            add_file("candidates_provisional_csv", provisional_csv, "csv", {"rows": csv_rows(provisional_csv)})
        if orf_report is not None:
            add_file("orf_validation_report", orf_report, "tsv")

        # Scoring metadata: every named vector with its weights, so a row's weight_vector column
        # resolves to the exact numbers that produced it. Weights are never altered at runtime, so
        # what is recorded here is what applied.
        scoring_weights = self.config.design_params.scoring
        # Derived from THIS run's vectors, not hand-listed: whether a term is scored is a per-run
        # property, and a literal cannot notice a term being promoted -- #96's drift class. The literal
        # missed pos1_mismatch and au_1_5, registered terms the shipped profile does not score. Its one
        # right instinct was naming paired_fraction, which has no TermRecord and is still reported and
        # unscored, so the universe is the registry PLUS those terms, never the registry alone.
        scored_terms = {term for vector in scoring_weights.all_vectors() for term in vector.terms}

        return {
            "tool": "sirnaforge",
            "tool_version": __version__,
            "gene_query": self.config.gene_query,
            "run_timestamp": now,
            # Dumped wholesale: hand-listing fields silently dropped thresholds
            # (min_asymmetry_score, max_poly_runs, ...) from the run record.
            "design_parameters": self.config.design_params.model_dump(mode="json"),
            # What was requested, what it resolved to, and which authority won for each setting --
            # a threshold in design_parameters above says what applied but not why.
            "run_policy": self.config.resolved_policy.as_manifest(),
            "scoring": {
                "weight_set_version": SCORING_WEIGHT_SET_VERSION,
                "vectors": scoring_weights.as_manifest(),
                "vector_terms": {vector.name: list(vector.terms) for vector in scoring_weights.all_vectors()},
                # scored_terms / reported_not_scored / their source and definitions: one partition of one
                # universe, so a reported term cannot be absent from both lists.
                **reported_term_partition(TERM_REGISTRY, scored_terms),
            },
            "files": files,
            "provenance": self._build_provenance_block(report_html=report_html),
        }

    def _build_provenance_block(self, *, report_html: Path | None) -> dict[str, Any]:
        """Name the build and the references this run's science rests on, or say why it cannot.

        The manifest recorded every gate as data and no reference identity at all: no Ensembl release,
        no transcriptome checksum, no index digest, no git SHA. Two runs could not be shown to have
        screened the same bytes. Everything here is read from state this process already holds -- no
        manager is instantiated, because a manager built at manifest time would attest bytes this run
        may never have touched.

        Wrapped whole: provenance is bookkeeping and must never fail a run. A named failure is not
        silence, so the degraded block says what went wrong instead of dropping the key.
        """
        try:
            build = dict(build_identity())
            # Omitted rather than nulled when no Nextflow screen ran: there is no pipeline to revise.
            if self._nextflow_cache_info is not None:
                build["pipeline_revision"] = pipeline_revision_identity(
                    self._nextflow_cache_info.get("pipeline_revision")
                )
            references = self._references_provenance()
            return {
                "schema_version": PROVENANCE_SCHEMA_VERSION,
                "build": build,
                "references": references,
                "databases": self._databases_provenance(),
                # The two design_parameters inputs that are bare nulls with no consumer, named.
                "design_inputs": self._design_input_provenance(),
                # `expected`, not `rendered`: this manifest is written before the render runs, so the
                # only honest claim here is that a report was asked for. report_manifest.json states
                # the outcome.
                "artifacts": {"report_html": report_html_artifact(expected=report_html is not None)},
                # Fed the reference entries this block just published, so coverage reports on what the
                # manifest actually says rather than re-deriving it from the same state twice.
                "coverage": self._coverage_provenance(references["screening"]["resolved"]),
            }
        except Exception as exc:
            logger.warning(f"Provenance block not assembled: {exc}")
            return {
                "schema_version": PROVENANCE_SCHEMA_VERSION,
                "state": "provenance_assembly_failed",
                "reason": str(exc),
            }

    def _screening_disabled_reason(self) -> str | None:
        """Why no reference resolved, or ``None`` when the screen was attempted.

        Never omitted and never silently empty: an absent screening block reads "a screen happened and
        we forgot", while an empty ``resolved`` beside a named reason reads "no screen was attempted".
        A rejection is not a disabled screen -- the attempt is recorded in ``rejected`` instead.
        """
        references = self._screening_references
        if references.references or references.rejections:
            return None
        if self.config.resolved_policy.run_mode is RunMode.DESIGN_ONLY:
            return "design-only mode requested, so no screening reference was resolved and none was expected"
        if not getattr(self.config.design_params, "check_off_targets", True):
            return "off-target screening was switched off by request, so no screening reference was resolved"
        return (
            "no screening reference resolved and none was rejected, so no screen was attempted through this entry point"
        )

    def _references_provenance(self) -> dict[str, Any]:
        """Both independently floating Ensembl surfaces: the screened bytes, and the REST annotation.

        Deliberately not merged. The screening references are cached FASTAs fetched from
        ``/pub/current_fasta``; the annotation client talks to ``rest.ensembl.org``, which serves
        whatever release is current. Two surfaces, two releases, and neither pins one.
        """
        references = self._screening_references
        requirements = self.config.resolved_policy.evidence_requirements
        return {
            "reference_identity_schema_version": PROVENANCE_SCHEMA_VERSION,
            "screening": {
                "kind": references.kind.value,
                "requested_species": list(references.requested_species),
                "resolved": [
                    self._reference_provenance(reference, requirements=requirements)
                    for reference in references.references
                ],
                "rejected": [
                    {"species": rejection.species, "identity": rejection.identity, "reason": rejection.reason}
                    for rejection in references.rejections
                ],
                "disabled_reason": self._screening_disabled_reason(),
            },
            "annotation": self._annotation_provenance(),
        }

    def _annotation_provenance(self) -> dict[str, Any]:
        """The REST annotation surface, but only as far as a call was actually made.

        ``EnsemblTranscriptModelClient`` is constructed unconditionally in ``__init__``, so keying this
        block on the client merely *existing* published ``provider: ensembl_rest``, an endpoint and an
        ``assembly_requested`` into every design-only run -- runs that annotate nothing. ``attempted``
        comes from ``_annotation_summary``, which only ``_enrich_transcript_annotations`` writes, and
        the assembly is read back out of that record rather than restated: the literal ``"GRCh38"``
        here was a second copy of the one at the call site, the drift class #96 named.
        """
        client = self._annotation_client
        if client is None:
            return {
                "attempted": False,
                "provider": absent(
                    "annotation_client_not_configured",
                    "no Ensembl REST annotation client was constructed, so nothing was annotated",
                ),
                "endpoint": absent(
                    "annotation_client_not_configured",
                    "no annotation client exists in this run, so there is no endpoint it could have called",
                ),
                "rest_release": absent(
                    "annotation_client_not_configured",
                    "no annotation client exists in this run, so no release could be served to it",
                ),
                "assembly_requested": absent(
                    "annotation_client_not_configured",
                    "no annotation call was possible, so no reference assembly was ever requested",
                ),
            }

        summary = self._annotation_summary
        endpoint = present(
            client.base_url,
            "client_configuration",
            "the base URL of the client this run constructed; a URL is not evidence that it was called",
        )
        rest_release = absent(
            "not_probed",
            "this run made no /info/data call, and the endpoint serves whatever release is current",
        )
        if not summary:
            unattempted = (
                "an Ensembl REST client was constructed, but this run reached no annotation call: nothing was "
                "annotated and no provider supplied anything"
            )
            return {
                "attempted": False,
                "provider": absent("no_annotation_call_made", unattempted),
                "endpoint": endpoint,
                "rest_release": rest_release,
                "assembly_requested": absent("no_annotation_call_made", unattempted),
            }

        if not summary.get("enabled"):
            failed = (
                f"the annotation call did not complete ({summary.get('error') or 'no error recorded'}), so no "
                "annotation reached this run"
            )
            return {
                "attempted": True,
                "provider": absent("annotation_call_failed", failed),
                "endpoint": endpoint,
                "rest_release": rest_release,
                "assembly_requested": absent("annotation_call_failed", failed),
            }

        # Read back from the call's own record, so the assembly is this run's request rather than a
        # second hard-coded copy of the literal at the call site.
        reference: Mapping[str, Any] = summary.get("reference") or {}
        requested = reference.get("value")
        return {
            "attempted": True,
            "provider": present(
                summary.get("provider"),
                "annotation_call_record",
                f"recorded by the annotation call this run made for {summary.get('transcripts_queried')} "
                f"transcripts, of which {summary.get('transcripts_resolved')} resolved",
            ),
            "endpoint": endpoint,
            "rest_release": rest_release,
            "assembly_requested": present(
                requested,
                "annotation_call_record",
                f"the ReferenceChoice the annotation call recorded ({reference.get('state')}: "
                f"{reference.get('reason')}) -- requested, never checked against what came back",
            )
            if requested
            else absent(
                "reference_not_recorded_by_call",
                "the annotation call completed but recorded no ReferenceChoice, so the assembly it asked for "
                "cannot be named",
            ),
        }

    def _reference_provenance(
        self, reference: ScreeningReference, *, requirements: EvidenceRequirements
    ) -> dict[str, Any]:
        """One resolved reference: the provenance of the choice, and the identity of the bytes.

        The choice half is free from the resolver. The bytes half comes only from
        ``identity_evidence``, so a reference that resolved without a cache entry (a caller's prebuilt
        index) is named unidentified rather than described with borrowed numbers.
        """
        species = reference.species
        evidence: Mapping[str, Any] = reference.identity_evidence or {}
        # Two hops, because the taxid table is keyed by MirGeneDB slug while everything upstream of it
        # speaks canonical species names.
        registry_entry = CANONICAL_SPECIES_REGISTRY.get(species, {})
        slug = registry_entry.get("mirgenedb_slug")
        taxonomy_id = MIRGENEDB_SPECIES_TABLE.get(str(slug), {}).get("taxonomy_id") if slug else None
        requiredness = requirements.requiredness_of(ScreeningChannel.TRANSCRIPTOME, species)

        entry: dict[str, Any] = {
            "species": species,
            # The channel this reference answers for, so `coverage` can match required pairs against
            # these published entries instead of re-deriving them from the reference set.
            "kind": reference.kind.value,
            "species_authority": reference.species_authority.value,
            "taxonomy_id": present(str(taxonomy_id), "species_registry", "NCBI taxid from MIRGENEDB_SPECIES_TABLE")
            if taxonomy_id
            else absent(
                "species_not_in_registry",
                f"'{species}' has no entry in the canonical species registry, so no taxid can be named",
            ),
            "registry_member": species in CANONICAL_SPECIES_REGISTRY,
            # Joined from the run policy, which is the single authority on requiredness.
            "requiredness": requiredness.value if isinstance(requiredness, Requiredness) else "not_declared",
            # A source NAME, not a version: the identity says which source, never which release.
            "source_id": reference.identity,
            "form": reference.form.value,
            "state": reference.state.value,
            "reason": reference.reason,
            "needs_index_build": reference.needs_index_build,
            "local_path": {
                "fasta": reference.fasta,
                "index_prefix": reference.index,
                # The 12-hex cache stem is an md5 of a release-floating URL, so it identifies the
                # request and not the bytes; it must never be read as a digest.
                "path_scope": "container" if build_identity()["container"]["in_container"] else "host",
            },
        }
        entry.update(self._reference_bytes_provenance(reference, evidence))
        entry["index"] = index_attestation(
            Path(reference.index),
            reference_fasta_name=Path(reference.fasta).name if reference.fasta else None,
        )
        entry["classification_index"] = self._classification_index_provenance(species)
        return entry

    def _reference_bytes_provenance(self, reference: ScreeningReference, evidence: Mapping[str, Any]) -> dict[str, Any]:
        """Which bytes were screened: url, assembly, release, digest and the filters applied.

        Every field is present whether or not it is knowable. Everything here is an identity of the
        *bytes*, and a reference that never carried one is named unidentified rather than described
        with numbers borrowed from the source it was named after -- including the assembly, which the
        bundled table can only associate with a source *name* and which
        :func:`~sirnaforge.provenance.assembly_identity` therefore reads out of the FASTA's own headers
        first.
        """
        tabled_assembly = next(
            (entry.assembly for entry in ENSEMBL_ASSEMBLIES if entry.transcriptome_source_key == reference.identity),
            None,
        )
        url = evidence.get("url")
        assembly_record = assembly_identity(
            tabled=tabled_assembly,
            url=url,
            fasta=Path(reference.fasta) if reference.fasta else None,
        )
        if not evidence:
            unidentified = absent(
                "cache_metadata_not_carried",
                "this reference resolved without a cache entry, so the bytes' own identity never reached the run",
            )
            return {
                "source_url": unidentified,
                "assembly": assembly_record,
                "ensembl_release": {**unidentified, "pinned": False},
                "bytes": {
                    "size_bytes": None,
                    "size_scope": "not_recorded",
                    # Both labels say "not_recorded" rather than one saying null and the other "full".
                    "digest": {"algorithm": None, "scope": digest_scope(None, filtered=False), **unidentified},
                    "cached_at": None,
                },
                "remote_observed": dict(unidentified),
                "filters": None,
            }

        filtered = bool(evidence.get("filters"))
        size_scope = evidence.get("size_scope") or "not_recorded"
        return {
            # A filtered entry's recorded URI is the base download plus a synthetic `#filters=`
            # fragment, so it must not be published as the URI these bytes were fetched from.
            "source_url": present(
                url,
                "cache_metadata_derived_uri",
                "the base download URL plus the synthetic #filters= fragment this subset is cached under: these "
                "bytes were cut from that download, not fetched from this URI",
            )
            if filtered
            else present(url, "cache_metadata", "ReferenceSource.url recorded when these bytes were fetched"),
            "assembly": assembly_record,
            "ensembl_release": self._ensembl_release_provenance(url),
            "bytes": {
                "size_bytes": evidence.get("file_size"),
                # Mandatory: the recorded size is decompressed and the remote Content-Length is not,
                # so an unlabelled "size" invites a false mismatch.
                "size_scope": size_scope,
                "digest": {
                    # Mandatory: this is md5 while every files.* digest is sha256.
                    "algorithm": evidence.get("checksum_algorithm"),
                    # Derived from the same evidence as size_scope, never a literal: the two labels
                    # contradicted each other for an uncompressed reference, and "full" was asserted
                    # over a filtered subset.
                    "scope": digest_scope(size_scope, filtered=filtered),
                    **present(
                        evidence.get("checksum"),
                        "cache_metadata",
                        "whole-file digest of the cached artifact, computed when these bytes were cached and "
                        "re-verified on reuse; bytes.digest.scope says which artifact that is",
                    ),
                },
                "cached_at": evidence.get("downloaded_at"),
            },
            "remote_observed": present(
                evidence.get("remote_observed"),
                "http_response_headers",
                "ETag/Last-Modified/Content-Length of the response these bytes arrived in",
            )
            if evidence.get("remote_observed")
            else absent(
                "cache_hit_not_refetched",
                "a valid cache entry was served with no HTTP request, and a speculative HEAD would describe "
                "upstream now rather than the bytes screened",
            ),
            # Explicit null, never an absent key: null means the whole reference was screened.
            "filters": evidence.get("filters"),
        }

    @staticmethod
    def _ensembl_release_provenance(url: str | None) -> dict[str, Any]:
        """Whether the release these bytes came from is pinned anywhere, and where.

        ``/pub/current_fasta`` pins nothing: the release is absent from the URL, from the cache
        metadata and from the FASTA headers alike. Recording that the pointer floated is the honest
        answer; stopping it floating is a behaviour change and a separate issue.
        """
        release = re.search(r"/release-(\d+)/", url or "")
        if release:
            return {
                **present(release.group(1), "release_pinned_url", "the download URL names the Ensembl release"),
                "pinned": True,
            }
        if url and "/current_fasta" in url:
            return {
                **absent(
                    "release_floating_url",
                    "fetched via /pub/current_fasta, which pins no release: it is absent from the URL, the "
                    "cache metadata and the FASTA headers",
                ),
                "pinned": False,
            }
        return {
            **absent(
                "release_not_recorded",
                "no Ensembl release is recoverable from this reference's URL, cache metadata or headers",
            ),
            "pinned": False,
        }

    def _classification_index_provenance(self, species: str) -> dict[str, Any]:
        """The confidence bound on every symbol-tier verdict for one species.

        ``TranscriptGeneIndex.build`` never raises, so a bare ``0`` here would read as an empty
        reference rather than as an index that was never built.
        """
        index = self._transcript_index.for_species(species)
        if index is None:
            return {
                "built": False,
                "transcript_count": None,
                "state": "not_built_by_this_run",
                "reason": "no transcript index was built for this species, so no hit could be resolved to a gene",
            }
        if index.transcript_count == 0:
            return {
                "built": False,
                "transcript_count": 0,
                "state": "index_build_returned_empty",
                "reason": "the index build parsed no transcript header, so the reference classified nothing",
            }
        fasta = self._species_cdna_fasta.get(species)
        return {
            "built": True,
            "transcript_count": index.transcript_count,
            "records_without_gene_symbol": index.missing_symbol_count,
            "symbol_coverage": round(1 - index.missing_symbol_count / index.transcript_count, 4),
            "built_from_fasta": str(fasta) if fasta else None,
        }

    def _databases_provenance(self) -> dict[str, Any]:
        """The miRNA corpus and the orthology surface, each named only as far as it is knowable."""
        databases: dict[str, Any] = {}
        # Omitted entirely unless both were set: the miRNA channel did not run at all otherwise, and an
        # empty block would read as a corpus we failed to identify.
        if self.config.mirna_database and self.config.mirna_species:
            databases["mirna"] = self._mirna_provenance()
        databases["orthology"] = self._orthology_provenance()
        return databases

    def _mirna_channel_evidence(self) -> dict[str, Any]:
        """What this run's reconciled evidence says about the miRNA seed channel, per species.

        The configuration says what was asked for; only the reconciliation says whether the channel
        published anything. Both are needed, because a configured database and a screened corpus are
        different claims and the block was emitting the first as if it were the second.
        """
        channel = ScreeningChannel.MIRNA_SEED.value
        if self._screening_evidence is None:
            return {
                "completed": False,
                "per_unit": absent(
                    "no_screening_reconciliation",
                    "this run reconciled no screening evidence at all, so whether the miRNA channel ran cannot "
                    "be read from it",
                ),
            }
        statuses = {
            species: status.value
            for (entry_channel, species, _digest), status in self._screening_evidence.statuses().items()
            if entry_channel == channel
        }
        if not statuses:
            return {
                "completed": False,
                "per_unit": absent(
                    "channel_absent_from_reconciliation",
                    "the reconciled evidence names no miRNA seed unit, so this run neither planned nor observed "
                    "one for any species",
                ),
            }
        return {
            "completed": EvidenceStatus.COMPLETE.value in statuses.values(),
            "per_unit": present(
                [{"species": species, "status": status} for species, status in sorted(statuses.items())],
                "screening_evidence_reconciliation",
                "the reconciled status of each miRNA seed unit this run planned or observed",
            ),
        }

    def _mirna_provenance(self) -> dict[str, Any]:
        """What the miRNA channel screened against, which the producer does not publish.

        The corpus digest, url, download date and mature-sequence count exist only inside the
        MIRNA_SEED_ANALYSIS task, whose cache is the container's. Instantiating a manager here would
        read the *host* cache and attest bytes this run never touched, so it is not done (#111 tier 3).

        ``source_name`` is therefore gated on the channel having completed: it was emitted from
        configuration alone, so a run that merely *named* a database published a corpus identity for a
        scan that never happened. ``per_species`` was a bare empty list, which reads as "no species"
        where the truth is that no per-species record exists to read.
        """
        from sirnaforge.data.mirna_manager import MiRNADatabaseManager  # noqa: PLC0415

        database = self.config.mirna_database
        known = database in MiRNADatabaseManager.SOURCES
        evidence = self._mirna_channel_evidence()
        if not known:
            source_name = absent(
                "not_a_known_source",
                f"'{database}' names no MiRNADatabaseManager.SOURCES entry, so no corpus can be identified",
            )
        elif not evidence["completed"]:
            source_name = absent(
                "channel_published_no_completed_evidence",
                f"'{database}' names a MiRNADatabaseManager.SOURCES entry, but this run reconciled no completed "
                "miRNA seed evidence, so no corpus was screened for that name to identify",
            )
        else:
            source_name = present(
                database,
                "sources_table",
                "MiRNADatabaseManager.SOURCES entry for the database this run's completed miRNA seed scan was "
                "configured with",
            )
        return {
            "requested": {"database": database, "species": list(self.config.mirna_species)},
            "channel_evidence": evidence["per_unit"],
            "source_name": source_name,
            "release": absent(
                "unversioned_endpoint",
                "the miRNA sources are served from release-less endpoints, and the FASTA headers carry no "
                "version either",
            ),
            "resolved": absent(
                "not_published_by_producer",
                "the corpus digest, url, download date and mature-sequence count exist only inside the "
                "MIRNA_SEED_ANALYSIS task, which does not publish them",
            ),
            "per_species": absent(
                "not_published_by_producer",
                "the MIRNA_SEED_ANALYSIS task publishes no per-species corpus record, so the mature-sequence "
                "count and digest behind each requested species cannot be named",
            ),
        }

    def _orthology_provenance(self) -> dict[str, Any]:
        """Whether an orthology lookup happened, and what it is entitled to claim.

        ``attempted: false`` with a named reason, never a defaulted source: a design-only run, a
        single-species screen and a gene that yielded neither ID nor symbol all reach step5's early
        exit, and publishing ``ensembl_compara`` for them named a REST call that never happened.
        """
        mapping = self._orthology_mapping
        if mapping is None or mapping.source is None:
            if self.config.resolved_policy.run_mode is RunMode.DESIGN_ONLY:
                reason = "design-only mode requested, so no orthology lookup was attempted and none was expected"
            elif mapping is None:
                reason = (
                    "no orthologue mapping was recorded: this entry point never reaches the step that resolves "
                    "orthologues"
                )
            else:
                reason = (
                    "the lookup returned before any request was built -- a single-species screen, or a gene that "
                    "yielded neither a gene ID nor a symbol"
                )
            return {
                "attempted": False,
                "source": None,
                "endpoint": None,
                "reason": reason,
                # Computed even here: conservation is divided by the species that were *screened*, so a
                # run with no orthology lookup can still publish fractions.
                "conservation_scope": self._conservation_scope(mapping),
            }

        summary = mapping.summary()
        from_compara = mapping.source == SOURCE_COMPARA
        lookup_species = sorted(mapping.resolved_species | mapping.unresolved_species)
        return {
            "attempted": True,
            "source": mapping.source,
            "endpoint": ENSEMBL_BASE_URL if from_compara else str(self.config.ortholog_mapping_file),
            # Compara's own vocabulary is published only for a Compara call. A user-supplied mapping
            # file was being described with a release reason that names rest.ensembl.org and with all
            # three ORTHOLOGUE_TYPES, as if a REST call had filtered them -- #101's warning exactly.
            "compara_release": absent(
                "unversioned_endpoint",
                "resolved against rest.ensembl.org with no release pin; naming the endpoint is the honest "
                "substitute for a version",
            )
            if from_compara
            else absent(
                "not_a_compara_lookup",
                "these orthologues came from the supplied mapping file named under `endpoint`, so no Compara "
                "release is involved in them at all",
            ),
            "relationship_types_accepted": present(
                summary.get("orthologue_types"),
                "compara_response_filter",
                "the homology types the resolver accepted from the Compara response; anything else was dropped",
            )
            if from_compara
            else absent(
                "not_a_compara_lookup",
                "the supplied mapping file declares no relationship type, and ORTHOLOGUE_TYPES filters Compara "
                "responses only, so nothing here was type-filtered",
            ),
            "relationship_type_observed": absent(
                "not_retained",
                "the resolver tests homology['type'] for membership and then discards it",
            )
            if from_compara
            else absent(
                "not_carried_by_mapping_file",
                "the supplied mapping file carries species and gene IDs only, so no relationship type was ever "
                "read for these orthologues",
            ),
            # A symbol must never be reprinted under a gene-ID key.
            "queried_identifiers": {
                "as_gene_ids": summary["queried_gene_ids"],
                "as_symbols": summary["queried_symbols"],
            },
            # Said out loud, so a per-species record is not mistaken for a per-species test.
            "classification_input": {
                "granularity": "union_across_species",
                "gene_id_count": len(mapping.all_gene_ids),
            },
            "per_species": [
                {
                    "species": species,
                    "requiredness": self._orthology_requiredness(species),
                    "status": self._orthology_status(mapping, species),
                    "gene_ids": sorted(mapping.gene_ids_by_species.get(species, frozenset())),
                    "route": absent(
                        "not_retained",
                        "_lookup_route logs which route answered and returns only the ids it found",
                    )
                    if from_compara
                    else absent(
                        "not_applicable_to_a_mapping_file",
                        "the supplied mapping file is read in one pass, so there is no ID-then-symbol route for "
                        "a species to have been answered by",
                    ),
                }
                # The species the *lookup* covered, which is not the conservation denominator.
                for species in lookup_species
            ],
            "confidence": absent(
                "not_returned_by_endpoint",
                "the resolver requests format=condensed, which returns no score; the per-hit evidence tier is "
                "the honest categorical confidence",
            )
            if from_compara
            else absent(
                "not_carried_by_mapping_file",
                "a supplied mapping file asserts membership with no score, so no confidence accompanies these "
                "orthologues",
            ),
            "conservation_scope": self._conservation_scope(mapping),
        }

    def _completed_transcriptome_species(self) -> frozenset[str] | None:
        """Species whose transcriptome evidence reconciled complete, or None when nothing reconciled.

        None is not the empty set: a run that reconciled nothing cannot say a species was unscreened,
        and #100's rule is that a unit with no evidence must never read as a clean one.
        """
        if self._screening_evidence is None:
            return None
        return frozenset(
            species
            for channel, species in self._completed_pairs_from_evidence(self._screening_evidence)
            if channel == ScreeningChannel.TRANSCRIPTOME.value
        )

    def _conservation_scope(self, mapping: OrthologueMapping | None) -> dict[str, Any]:
        """The denominator conservation was actually divided by, and whether its fractions are floors.

        Both fields were wrong. ``denominator_species`` published the *orthology lookup's* species set,
        which is not what ``_score_candidate_post_screen`` divides by -- it divides by the species
        handed to the aligner minus the query species, so the two disagree whenever a lookup covered a
        species the screen did not (or the reverse), and the manifest named a denominator nothing was
        scored against. ``is_lower_bound`` was ``bool(unresolved_species)``, which is true in exactly
        the case where the scorer publishes *no* conservation at all: an unresolved species in the
        denominator nulls ``conservation_score`` rather than lowering it. A floor is instead what a
        species that stayed in the denominator with no alignment evidence produces.
        """
        # The scorer's own expression, replicated verbatim so the two cannot drift.
        requested = frozenset(
            normalize_species_name(s) for s in (self._active_screen_species or self.config.screen_species)
        )
        denominator = sorted(requested - {self._query_species})
        unresolved = sorted(
            set(denominator) & {normalize_species_name(s) for s in (mapping.unresolved_species if mapping else ())}
        )
        screened = self._completed_transcriptome_species()
        source = (
            "the species handed to the aligner (_active_screen_species, else config.screen_species) minus the "
            "query species -- the same expression _score_candidate_post_screen divides by"
        )

        if screened is None:
            unscreened_record = absent(
                "no_screening_reconciliation",
                "this run reconciled no screening evidence, so which of the denominator's species produced "
                "alignments is not recoverable here",
            )
            missing: list[str] = []
        else:
            missing = sorted(set(denominator) - screened)
            unscreened_record = present(
                missing,
                "screening_evidence_reconciliation",
                "denominator species with no completed transcriptome evidence; they stay in the denominator, so "
                "they can only lower a conservation fraction",
            )

        if not denominator:
            is_lower_bound = absent(
                "conservation_term_inactive",
                "no non-query species was screened, so conservation_sub_score returns None for every candidate "
                "and there is no published fraction for a bound to apply to",
            )
        elif unresolved:
            is_lower_bound = absent(
                "conservation_not_published",
                f"orthology is unresolved for {unresolved} inside the denominator, so the scorer publishes no "
                "conservation at all rather than a fraction that could be a bound",
            )
        elif screened is None:
            is_lower_bound = absent(
                "screening_coverage_unknown",
                "no reconciled evidence says which denominator species produced alignments, so whether every "
                "published fraction is a floor cannot be decided",
            )
        elif missing:
            is_lower_bound = present(
                True,
                "denominator_species_without_alignment_evidence",
                f"{missing} stayed in the conservation denominator with no completed transcriptome evidence, so "
                "every published conservation fraction is a floor",
            )
        else:
            is_lower_bound = present(
                False,
                "every_denominator_species_screened",
                "every species in the conservation denominator produced completed transcriptome evidence, so a "
                "published fraction is exact rather than a floor",
            )

        return {
            "denominator_species": denominator,
            "denominator_source": source,
            "unresolved_in_denominator": unresolved,
            "unscreened_in_denominator": unscreened_record,
            "is_lower_bound": is_lower_bound,
        }

    def _orthology_requiredness(self, species: str) -> str:
        """Requiredness of one species' transcriptome evidence, read from the run policy only."""
        requiredness = self.config.resolved_policy.evidence_requirements.requiredness_of(
            ScreeningChannel.TRANSCRIPTOME, species
        )
        return requiredness.value if isinstance(requiredness, Requiredness) else "not_declared"

    @staticmethod
    def _orthology_status(mapping: OrthologueMapping, species: str) -> str:
        """Three-way, never a boolean: resolved-empty is evidence of absence, unresolved is not."""
        if species in mapping.unresolved_species:
            return "unresolved"
        return "resolved" if mapping.gene_ids_by_species.get(species) else "resolved_empty"

    def _coverage_provenance(self, resolved_entries: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
        """The headline invariant, made machine-readable rather than left a convention.

        A non-empty list is not a crash: it is the manifest declaring, in one place, that a reference
        this run depended on cannot be shown to be the reference another run used.

        Identification is read off the entries this block just published, and it takes *both* halves:
        bwa-mem2 aligns against the index files alone, so a reference whose FASTA digest is recorded
        while its index bytes are unattested was being reported as fully identified even though the
        off-target claims rest on the index. Which half is missing is named beside the pair.
        """
        required = sorted(self.config.resolved_policy.evidence_requirements.required_pairs)
        entries = {
            str(entry.get("species")): entry
            for entry in resolved_entries
            if entry.get("kind") == ScreeningChannel.TRANSCRIPTOME.value
        }
        unnamed = [
            [channel, species]
            for channel, species in required
            if channel != ScreeningChannel.TRANSCRIPTOME.value or species not in entries
        ]
        unidentified: list[list[str]] = []
        detail: list[dict[str, Any]] = []
        for channel, species in required:
            entry = entries.get(species)
            if entry is None:
                continue
            missing: list[str] = []
            if not ((entry.get("bytes") or {}).get("digest") or {}).get("value"):
                missing.append("fasta_digest_not_recorded")
            if not (entry.get("index") or {}).get("attested"):
                missing.append("index_bytes_unattested")
            if missing:
                unidentified.append([channel, species])
                detail.append({"channel": channel, "species": species, "missing": missing})
        return {
            "required_pairs": [[channel, species] for channel, species in required],
            "unnamed_required_references": unnamed,
            "unidentified_required_references": unidentified,
            # Which evidence each unidentified pair lacks, so the list is actionable rather than a flag.
            "unidentified_required_reference_detail": detail,
        }

    def _design_input_provenance(self) -> dict[str, Any]:
        """The two ``design_parameters`` inputs that are bare nulls, and what a null there means.

        ``genome_index`` and ``snp_file`` are CLI options forwarded onto ``DesignParameters`` and read
        by nothing in 0.7.1, so ``model_dump`` emits them as bare nulls that read exactly like "not
        applicable". They were named in the gap this block was written to close, so each is stated here
        as a named absence -- with the manifest path where the bare null lives, the way
        ``reported_not_scored_definitions`` names where a term is defined.
        """
        parameters = self.config.design_params
        return {
            "genome_index": {
                "design_parameters_path": "design_parameters.genome_index",
                **(
                    present(
                        parameters.genome_index,
                        "supplied_but_unread",
                        "--genome-index was recorded on design_parameters, but no 0.7.1 code path reads it: the "
                        "screen aligns against provenance.references.screening, so this constrained nothing",
                    )
                    if parameters.genome_index
                    else absent(
                        "not_supplied",
                        "no --genome-index was given, and no 0.7.1 code path reads the field either; the indexes "
                        "actually aligned against are named per species in provenance.references.screening",
                    )
                ),
            },
            "snp_file": {
                "design_parameters_path": "design_parameters.snp_file",
                **(
                    present(
                        parameters.snp_file,
                        "supplied_but_unread",
                        "--snp-file was recorded on design_parameters, but avoid_snps is unimplemented in 0.7.1 "
                        "(no code reads the flag), so no variant was avoided on account of this file",
                    )
                    if parameters.snp_file
                    else absent(
                        "not_supplied",
                        "no --snp-file was given, and avoid_snps is unimplemented in 0.7.1 (no code reads the "
                        "flag), so no SNP avoidance happened whether or not one had been",
                    )
                ),
            },
        }

    def _write_candidate_tables(self, candidates: Sequence[SiRNACandidate], tables: CandidateTablePaths) -> None:
        """Write the four candidate tables and the PASS/qualified FASTAs from one candidate frame.

        Error-neutral on purpose: ``step6_generate_reports`` must never fail a run over a reporting
        extra and wraps this call, while ``write_offtarget_only_exports`` lets it raise, because for
        that entry point these tables ARE the deliverable.
        """
        validated_all = self._build_candidate_frame(candidates)
        pass_df = validated_all[validated_all["passes_filters"] == "PASS"].copy()

        validated_all.to_csv(tables.all_csv, index=False)
        pass_df.to_csv(tables.pass_csv, index=False)

        if pass_df.empty:
            tables.pass_fasta.unlink(missing_ok=True)
        else:
            try:
                self._write_pass_candidates_fasta(pass_df, tables.pass_fasta)
            except Exception as e:
                logger.warning(f"Failed to write PASS candidates FASTA: {e}")

        self._write_selection_exports(
            validated_all,
            qualified_csv=tables.qualified_csv,
            qualified_fasta=tables.qualified_fasta,
            provisional_csv=tables.provisional_csv,
        )

    def write_offtarget_only_exports(self, candidates: Sequence[SiRNACandidate]) -> dict[str, str]:
        """Publish the candidate tables and the manifest for a run of pre-designed guides (#100).

        ``run_offtarget_only_workflow`` bypasses step5 and step6 entirely. This writes the same four
        tables step6 writes, from the same :meth:`_build_candidate_frame` into the same ``sirnaforge/``
        sub-directory, so one column set and one directory shape serve both entry points and
        ``sirnaforge report <output_dir>`` needs no special case.

        What is *not* published is as deliberate. A pre-designed guide arrives with no design score and
        no transcript context, so ``design_score``, ``composite_score`` and every accessibility column
        stay null and ``scored_after_screening`` stays False: score-unavailable is the honest cell.
        Deriving accessibility from the guide alone would be a fabricated target measurement, and
        falling back to a design score would invent a potency claim for a guide this tool did not
        design. Nothing here writes an ORF report either, because no ORF was validated.

        Returns:
            Written file name -> path, for the caller's own result dictionary. Only files that exist:
            the qualified FASTA is absent when nothing qualified.
        """
        base = self.config.output_dir / "sirnaforge"
        base.mkdir(parents=True, exist_ok=True)
        tables = CandidateTablePaths.under(base)

        self._write_candidate_tables(candidates, tables)

        manifest = self._build_fair_manifest(
            all_csv=tables.all_csv,
            pass_csv=tables.pass_csv,
            pass_fasta=tables.pass_fasta,
            qualified_csv=tables.qualified_csv,
            qualified_fasta=tables.qualified_fasta,
            provisional_csv=tables.provisional_csv,
        )
        with tables.manifest_json.open("w") as mf:
            json.dump(manifest, mf, indent=2)

        written = {
            "candidates_all_csv": tables.all_csv,
            "candidates_pass_csv": tables.pass_csv,
            "candidates_pass_fasta": tables.pass_fasta,
            "candidates_qualified_csv": tables.qualified_csv,
            "candidates_qualified_fasta": tables.qualified_fasta,
            "candidates_provisional_csv": tables.provisional_csv,
            "manifest_json": tables.manifest_json,
        }
        return {name: str(path) for name, path in written.items() if path.exists()}

    async def step5_offtarget_analysis(self, design_results: DesignResult) -> dict[str, Any]:
        """Step 5: Detect repeat elements, then run off-target analysis via the Nextflow pipeline.

        Repeat detection runs here rather than as its own workflow step so it can reuse the
        transcriptome reference this step already materializes for screening, instead of
        fetching it a second time just to locate the query species' FASTA. The corollary is
        that ``check_off_targets=False`` (``--skip-off-targets``) skips repeat detection as
        well: both are reference-based scans, and the reference is what the flag exists to
        avoid paying for.
        """
        candidates_for_offtarget = self._select_candidates_for_offtarget(design_results)

        if not candidates_for_offtarget:
            console.print("⚠️  No candidates available for off-target analysis")
            self._repeat_summary = {
                "status": "skipped",
                "reason": "no_candidates",
                "repeat_flagged_count": 0,
                "threshold_fraction": self.config.design_params.filters.max_repeat_transcript_fraction,
            }
            # Nothing to rank, but every step5 exit still publishes a selection outcome.
            self._apply_post_screen_ranking(design_results)
            return self._publish_run_status({"status": "skipped", "reason": "no_candidates"})

        # Honour the skip request BEFORE touching any reference: materializing the default
        # transcriptomes downloads and indexes multi-gigabyte cDNA files, and repeat detection
        # scans against that same reference (~47s). Doing either first made
        # --skip-off-targets/check_off_targets=False cost nearly as much as a real screen, and
        # made "skipped by user request" a lie about work already done.
        if not getattr(self.config.design_params, "check_off_targets", True):
            console.print("⚠️  Off-target analysis skipped by user request")
            console.print("   ↳ repeat-element detection skipped too: it needs the same cDNA reference")
            self._repeat_summary = {
                "status": "skipped",
                "reason": "user_disabled",
                "repeat_flagged_count": 0,
                "threshold_fraction": self.config.design_params.filters.max_repeat_transcript_fraction,
            }
            # Still rebuild top_candidates: repeat flags may have been stamped elsewhere, and
            # downstream reporting expects a ranked list on every path (issue #80 F4/F5).
            self._apply_post_screen_ranking(design_results)
            return self._publish_run_status({"status": "skipped", "reason": "user_disabled"})

        # One `finally`, so selection runs on every exit: the returns, the caught Nextflow failure, and
        # an exception that propagates. Not a try/except, because `_prepare_offtarget_input`'s
        # duplicate-id refusal must stay a refusal rather than becoming a skip. Without this a qualified
        # run whose screen produced no evidence still published the pre-screen shortlist (#100).
        try:
            # Prepare input files
            input_fasta = await self._prepare_offtarget_input(candidates_for_offtarget)

            # Resolve the screening references once; repeat detection and Nextflow both reuse them
            # instead of each fetching/indexing their own copy.
            additional_params: dict[str, Any] = dict(self.config.nextflow_config)
            # Before resolution, so a reference that cannot be fetched is still a planned unit that
            # reconciles as failed rather than a species that was never mentioned again (#100).
            self._record_screening_plan(input_fasta, additional_params)
            has_transcriptome = await self._resolve_screening_references(additional_params)
            self._repeat_summary = self._run_repeat_detection(candidates_for_offtarget)

            # Try Nextflow pipeline first. We do NOT run the simplistic sequence-based fallback
            # (it produces low-value results) when Nextflow is unavailable. Instead mark as skipped
            # so downstream steps/users can see the explicit reason.
            try:
                return await self._run_nextflow_offtarget_analysis(
                    candidates_for_offtarget, input_fasta, additional_params, has_transcriptome
                )
            except Exception as e:
                console.print(f"⚠️  Nextflow execution failed: {e}")
                logger.exception("Nextflow pipeline execution error")
                # An aborted screen leaves the gates undecided, not unrun (#106). Only the candidates
                # that reached no gate: the failure can be raised after integration already gated some.
                self._gate_without_screening_evidence(candidates_for_offtarget)
                return self._publish_run_status({"status": "skipped", "reason": "nextflow_failed", "error": str(e)})
        finally:
            self._apply_post_screen_ranking(design_results)

    def _select_candidates_for_offtarget(self, design_results: DesignResult) -> list[SiRNACandidate]:
        """Return ALL candidates plus any dirty controls for off-target analysis.

        Screening now happens for every distinct candidate (user story 10), not just top_n.
        Dirty controls remain included as sentinels to verify that downstream aligners,
        Nextflow modules, and reports are actually running.
        """
        # Select all candidates (every distinct candidate must be screened)
        selected: list[SiRNACandidate] = list(design_results.candidates)

        # Ensure dirty controls are included (they should already be in candidates, but double-check)
        dirty_controls = [c for c in design_results.candidates if self._is_dirty_control_candidate(c)]
        for control in dirty_controls:
            if control not in selected:
                selected.append(control)

        return selected

    @staticmethod
    def _is_dirty_control_candidate(candidate: SiRNACandidate) -> bool:
        """Identify dirty control sequences injected for observability."""
        status = getattr(candidate, "passes_filters", True)
        issues = getattr(candidate, "quality_issues", []) or []

        status_is_dirty = False
        if isinstance(status, bool):
            status_is_dirty = False
        elif isinstance(status, SiRNACandidate.FilterStatus):
            status_is_dirty = status == SiRNACandidate.FilterStatus.DIRTY_CONTROL
        else:
            status_is_dirty = str(status) == DIRTY_CONTROL_LABEL

        return status_is_dirty or (DIRTY_CONTROL_LABEL in issues)

    @staticmethod
    def _passes_filters(candidate: SiRNACandidate) -> bool:
        """True when a candidate is currently passing, under either passing representation.

        Mirrors the two-way test used elsewhere in the repo (see core/design.py
        design_from_file): passes_filters can be the bare bool True or the PASS enum member.
        """
        return candidate.passes_filters is True or candidate.passes_filters == SiRNACandidate.FilterStatus.PASS

    def _apply_post_screen_ranking(self, design_results: DesignResult) -> None:
        """Re-rank candidates and rebuild ``top_candidates`` from the ones the run can stand behind.

        Screening activates the off_target, isoform_coverage and conservation terms, so the design-time
        order no longer reflects the final ranking. Re-sorts ``design_results.candidates`` in place
        (step6 writes the CSVs in this order). Called after repeat detection and again on every path out
        of screening, failures included; idempotent.

        Two independent reasons a gate-passing candidate is still held out, kept separate because
        neither implies the other:

        1. **Comparability.** A design-time score is on a different vector and is systematically the
           more optimistic number, so an unscored candidate cannot compete against post-screen
           neighbours. Applies only to a batch holding both; a wholly unscored batch is internally
           consistent on its own terms.
        2. **Required evidence**, decided by run mode and independently of score availability (#100).

        - ``DESIGN_ONLY``: nothing was screened, so nothing is missing; the design shortlist stands.
        - ``QUALIFIED``: requires the policy's declared required evidence, and no ``UNKNOWN`` verdict on
          a gate reading it.
        - ``EXPLORATORY``: incomplete evidence is retained, sorted below complete evidence.

        Every exclusion is counted by reason in :attr:`_selection_summary`, because a qualified run that
        found nothing eligible and one that never produced the evidence both leave the shortlist empty.

        The decision itself is :func:`sirnaforge.core.selection.select`, reached through
        :meth:`_apply_selection` (#100), which knows nothing about ``DesignResult`` -- the
        off-target-only entry point has none and must reach the same answer through the same code.
        ``EXPLORATORY``'s retained-but-incomplete candidates carry ``SelectionState.PROVISIONAL``: they
        still rank, below complete evidence, and still enter ``top_candidates``.
        """
        design_results.top_candidates = self._apply_selection(design_results.candidates)

    def _apply_selection(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Decide the selection over ``candidates``, apply it to them in place, and return the shortlist.

        Everything :meth:`_apply_post_screen_ranking` does except naming the object that holds the
        shortlist, so ``run_offtarget_only_workflow`` -- which bypasses step5 and therefore holds no
        ``DesignResult`` -- publishes a selection under the same rules, the same states and the same
        console messages (#100).

        Re-sorts ``candidates`` in place (the CSV writers emit them in this order) and returns the
        first ``top_n`` eligible candidates.
        """
        # Before the views are built, because eligibility reads this gate's verdict off the candidate
        # and a gate that never ran has none to read (#105).
        self._complete_isoform_coverage_verdicts(candidates)
        views = tuple(self._candidate_view(candidate, ordinal) for ordinal, candidate in enumerate(candidates))
        result = select(views, self._selection_inputs())

        # Ordinals, not ids or identity: duplicate ids are refused before screening, but selection still
        # runs on such a list via step5's `finally`, and a dict keyed by id would collapse two answers.
        by_ordinal = list(candidates)
        # Re-sorted in place, because step6 writes the CSVs in this list's order.
        candidates[:] = [by_ordinal[ordinal] for ordinal in result.order]
        top_candidates = [by_ordinal[ordinal] for ordinal in result.top_ordinals]
        # On the candidate, not bolted onto a dataframe later, so both CSV writers emit one column set
        # from build_candidate_row and the FASTA header can read it.
        for selected in result.per_candidate:
            by_ordinal[selected.ordinal].selection_state = selected.state.value
        self._eligible_candidate_ids = frozenset(by_ordinal[ordinal].id for ordinal in result.eligible_ordinals)
        self._selection_summary = result.summary

        summary = result.summary
        withheld = int(summary["evidence_excluded"])
        scored_count = int(summary["scored_after_screening"])
        provisional = int(summary["provisional_candidates"])
        logger.info(
            f"Re-ranked {len(result.eligible_ordinals)} eligible candidates after screening (excluded "
            f"{summary['repeat_excluded']} repeat-flagged, {summary['filter_excluded']} failing a gate, "
            f"{withheld} for incomplete required evidence, {summary['unscored_excluded']} for an "
            "incomparable score)"
        )
        if withheld:
            reasons = cast(dict[str, int], summary["evidence_shortfall_reasons"])
            logger.error(
                f"{withheld} of {len(candidates)} candidates are held out of the qualified shortlist for "
                f"incomplete required evidence: {reasons}. Re-run with the missing evidence, or use "
                "run_mode=exploratory to retain them with that scope stated."
            )
            # The reasons, not a generic sentence about screening: an unknown caused by a coverage
            # annotation gap is not a screening failure and does not have the same fix.
            console.print(
                f"⚠️  {withheld} candidate(s) excluded from the qualified shortlist "
                f"({describe_shortfall_reasons(reasons)}). This run cannot qualify them; "
                "re-run with the missing evidence, or use run_mode=exploratory to keep them labelled."
            )
            console.print(
                "   ↳ they remain in candidates_all.csv and candidates_pass.csv with "
                f"selection_state={SelectionState.WITHHELD.value}"
            )
        if provisional:
            # Retained, ranked and exported -- but named, because an exploratory result that reads like a
            # qualified one is the confusion this state exists to remove.
            console.print(
                f"ℹ️  {provisional} candidate(s) are provisional: retained on incomplete evidence and "
                "ranked below every fully evidenced candidate (candidates_provisional.csv)."
            )
        if 0 < scored_count < len(candidates):
            logger.error(
                f"{len(candidates) - scored_count} of {len(candidates)} candidates could not be scored after "
                "screening; they keep design-time scores and are excluded from top_candidates because the two "
                "scores are not comparable."
            )
        return top_candidates

    def _complete_isoform_coverage_verdicts(self, candidates: Sequence[SiRNACandidate]) -> None:
        """Apply the coverage gate to any candidate that never reached it, so its verdict exists (#105).

        ``_apply_isoform_coverage_gate`` is reached from :meth:`_score_and_gate`, which only the
        integration path calls -- and only for a candidate whose query species was actually screened. So
        on every other route out of screening (``nextflow_unavailable``, ``nextflow_failed``, the basic
        sequence-only fallback, a missing output directory, and any unscreened candidate on the
        integration path itself) a configured floor reached no candidate at all: ``filter_verdicts``
        held no ``min_isoform_coverage`` key and ``build_candidate_row`` exported ``not_evaluated``,
        which reads as "no floor was configured" rather than "the floor could not be checked" (#105).

        Selection is the one point every one of those paths passes through, so the gate is completed
        here. Only where no verdict was recorded yet: the call that held the measurement owns the
        answer, and re-deciding would log the same rejection twice and re-offer the rejection label.

        Coverage is only computed during post-screen scoring, so on those paths there is nothing to
        compare and the verdict is ``UNKNOWN`` -- which never rejects, but does withhold a candidate
        from a QUALIFIED shortlist, because a run that cannot check a floor it was given has not met
        it. With no floor configured the gate still makes no claim and the verdict stays
        ``NOT_EVALUATED``, so a default run is unaffected.
        """
        for candidate in candidates:
            if "min_isoform_coverage" in (candidate.filter_verdicts or {}):
                continue
            self._apply_isoform_coverage_gate(candidate)

    def _candidate_view(self, candidate: SiRNACandidate, ordinal: int) -> CandidateView:
        """One candidate as the facts selection reads, and nothing else.

        ``unknown_filter_ids`` is computed here rather than inside selection so that module never has
        to know how a verdict was spelled -- only that it was undecided.
        """
        verdicts = candidate.filter_verdicts or {}
        return CandidateView(
            candidate_id=candidate.id,
            ordinal=ordinal,
            passes_filters=self._passes_filters(candidate),
            repeat_flagged=bool(candidate.repeat_flagged),
            scored_after_screening=bool(candidate.scored_after_screening),
            off_target_screened=bool(candidate.off_target_screened),
            ranking_score=ranking_score(candidate),
            unknown_filter_ids=frozenset(
                filter_id for filter_id, verdict in verdicts.items() if verdict == FilterEvaluation.UNKNOWN.value
            ),
        )

    def _selection_inputs(self) -> SelectionInputs:
        """Everything selection needs that is a property of the run rather than of one candidate.

        One builder, so every entry point that selects (step5's three exits, and #100's
        offtarget-only path) asks the same question of the same policy and the same evidence record.

        ``repeat_rejects`` consults the repeat gate's own action, so
        ``max_repeat_transcript_fraction=warn`` means something here.

        ``design_input_shortfalls`` is snapshotted, not aliased: :class:`SelectionInputs` is frozen and
        must describe the run as it was when the decision was made. It reaches the summary as
        run-level ``required_evidence_missing`` entries only -- a batch failure disqualifies no
        individual candidate, because the guides that *were* designed are as well evidenced as before;
        what it denies is the claim that the target set was covered (#100).
        """
        policy = getattr(self.config, "resolved_policy", None)
        repeat_action = self._filter_action_for(
            "max_repeat_transcript_fraction", self.config.design_params.filters.max_repeat_transcript_fraction
        )
        return SelectionInputs(
            run_mode=policy.run_mode if policy is not None else None,
            requirements=policy.evidence_requirements if policy is not None else None,
            completed_pairs=self._completed_evidence_pairs,
            query_species=self._query_species,
            filter_channels=POST_SCREEN_FILTER_CHANNELS,
            repeat_rejects=repeat_action is FilterAction.FAIL,
            top_n=self.config.top_n,
            design_input_shortfalls=dict(self._design_input_shortfalls),
        )

    async def _prepare_offtarget_input(self, candidates: list[SiRNACandidate]) -> Path:
        """Prepare FASTA input file for off-target analysis with deduplication.

        Deduplicates by normalized guide sequence: writes one FASTA record per distinct
        sequence using a representative candidate ID, and keeps a mapping from
        representative ID to all candidates sharing that sequence for later fan-out.
        """
        input_fasta = self.config.output_dir / "off_target" / "input_candidates.fasta"

        # Candidate ids must be globally unique: a collision would silently fan one
        # candidate's screening results onto a different candidate's sequence below.
        id_counts: dict[str, int] = {}
        for candidate in candidates:
            id_counts[candidate.id] = id_counts.get(candidate.id, 0) + 1
        duplicate_ids = sorted(cid for cid, count in id_counts.items() if count > 1)
        if duplicate_ids:
            raise ValueError(
                f"Duplicate candidate id(s) before off-target screening: {duplicate_ids[:5]}"
                f"{' (+more)' if len(duplicate_ids) > 5 else ''}. Refusing to screen: results "
                "cannot be safely attributed to the correct sequence."
            )

        # Deduplicate by normalized guide sequence
        sequence_to_candidates: dict[str, list[SiRNACandidate]] = {}
        for candidate in candidates:
            norm_guide = normalize_guide_sequence(candidate.guide_sequence)
            sequence_to_candidates.setdefault(norm_guide, []).append(candidate)

        # Pick one representative per distinct sequence and build FASTA
        sequences: list[tuple[str, str]] = []
        representative_to_candidates: dict[str, list[SiRNACandidate]] = {}
        candidate_id_to_representative: dict[str, str] = {}
        for _norm_guide, cand_list in sequence_to_candidates.items():
            # Use the first candidate as the representative
            representative = cand_list[0]
            sequences.append((representative.id, representative.guide_sequence))
            representative_to_candidates[representative.id] = cand_list
            for cand in cand_list:
                candidate_id_to_representative[cand.id] = representative.id
                # The FASTA header the aligner sees, hence qname on every hit row. Written onto the
                # candidate so the join is an id match rather than a byte comparison of sequences.
                cand.screen_query_id = representative.id

        FastaUtils.save_sequences_fasta(sequences, input_fasta)

        # Store the mappings for fan-out during integration. Duplicate ids are rejected
        # above, so this id->id lookup is O(1) and never needs to compare candidates by value.
        self._representative_to_candidates = representative_to_candidates
        self._candidate_id_to_representative = candidate_id_to_representative
        # What "submitted" means in every evidence entry this run synthesizes: the FASTA records, not
        # the candidate count. Recorded here because this is the only place that knows it, and an
        # entry claiming a submitted count nothing observed is the defect the counts model exists for.
        self._submitted_guide_count = len(sequences)

        console.print(
            f"📝 Prepared off-target input: {len(candidates)} candidates → "
            f"{len(sequences)} distinct sequences (deduplication ratio: {len(candidates) / max(1, len(sequences)):.2f}x)"
        )
        return input_fasta

    @staticmethod
    def _carry_index_build_error(
        manager: TranscriptomeManager, prepared: Mapping[str, Any], payload: dict[str, Any]
    ) -> None:
        """Copy a failed index build onto the prepared reference, so the caller can refuse it.

        The manager records the reason against the FASTA rather than in its result dict, whose type
        is shared with the genome and annotation managers.
        """
        error = manager.index_build_errors.get(str(prepared.get("fasta")))
        if error:
            payload[INDEX_BUILD_ERROR_KEY] = error

    async def _prepare_transcriptome_database(
        self, transcriptome_ref: str, filter_spec: list[str] | None = None
    ) -> dict[str, Any] | None:
        """Fetch and index one transcriptome reference, reporting what each species authority saw.

        Deliberately does **not** decide the species: it reports the bundled source's own label and
        what the reference's headers say, and :meth:`_resolve_screening_reference` applies the
        authority order. A method that both fetched and labelled is how a parameter name came to
        stand in for a species.

        Args:
            transcriptome_ref: Bundled source name (e.g. ``ensembl_human_cdna``), local path, or
                HTTP(S)/FTP URL.
            filter_spec: Optional list of filter names (e.g. ``['protein_coding']``).

        Returns:
            ``fasta``/``index`` paths plus ``source_species``, ``header_species`` and the manager's
            ``identity`` view of the cached bytes, or None when preparation failed. The identity is
            forwarded rather than re-derived: the manager built here is garbage-collected on return,
            and it is the only thing that ever knew the url, size, md5 and download date.
        """
        try:
            manager = TranscriptomeManager()

            # Check if it's a pre-configured source
            if transcriptome_ref in manager.SOURCES:
                logger.info(f"Using pre-configured transcriptome source: {transcriptome_ref}")

                # Apply filters if specified
                if filter_spec and len(filter_spec) > 0:
                    logger.info(f"Applying filters: {', '.join(filter_spec)}")
                    raw_result = manager.get_filtered_transcriptome(
                        transcriptome_ref, filters=filter_spec, build_index=True
                    )
                else:
                    raw_result = manager.get_transcriptome(transcriptome_ref, build_index=True)

                if raw_result is None:
                    return None

                enriched_result: dict[str, Any] = {"source_species": manager.SOURCES[transcriptome_ref].species}
                enriched_result.update(raw_result)
                self._carry_index_build_error(manager, raw_result, enriched_result)
                return enriched_result

            # Otherwise treat as custom path/URL
            logger.info(f"Processing custom transcriptome reference: {transcriptome_ref}")
            if filter_spec and len(filter_spec) > 0:
                logger.warning("Filtering is not supported for custom transcriptome paths; ignoring filters")
            raw_custom = manager.get_custom_transcriptome(transcriptome_ref, build_index=True)
            if raw_custom is None:
                return None

            enriched_custom: dict[str, Any] = {
                "header_species": infer_species_from_cdna_headers(Path(raw_custom["fasta"]))
            }
            enriched_custom.update(raw_custom)
            self._carry_index_build_error(manager, raw_custom, enriched_custom)
            return enriched_custom

        except Exception as e:
            logger.exception(f"Failed to prepare transcriptome database from {transcriptome_ref}")
            console.print(f"⚠️  Transcriptome preparation error: {e}")
            return None

    @staticmethod
    def _readable_cdna_fasta(path: Path) -> bool:
        """Whether a file can be read as plain-text FASTA, which is what the classifier needs.

        A compressed or binary neighbour of an index prefix parses as an empty transcript index, and
        an empty index classifies nothing while reporting a completed screen.
        """
        try:
            with path.open() as handle:
                for line in handle:
                    if line.strip():
                        return line.startswith(">")
        except (OSError, UnicodeDecodeError):
            return False
        return False

    def _adopt_prebuilt_index(self, request: ReferenceRequest) -> dict[str, Any]:
        """Adopt a caller-built index, locating the sequence file its hits must be classified against.

        bwa-mem2 writes its index files beside the FASTA they were built from, so an index prefix is
        usually that FASTA's path; the plain suffixes are tried after it. Only a file that reads as
        plain-text FASTA counts, because that is what the transcript->gene index is built from. The
        index files at the prefix are never checked -- they may exist only inside the container -- but
        a host-readable sequence file is required, so an index bundle mounted into the container alone
        is refused rather than screened against something nothing can classify.
        """
        prefix = Path(request.value)
        candidates = [prefix, *(prefix.with_name(prefix.name + suffix) for suffix in (".fa", ".fasta", ".fna"))]
        companion = next((path for path in candidates if path.is_file() and self._readable_cdna_fasta(path)), None)
        return {
            "index": prefix,
            "fasta": companion,
            "header_species": infer_species_from_cdna_headers(companion) if companion else None,
        }

    def _reject_reference(self, request: ReferenceRequest, species: str, reason: str) -> ReferenceRejection:
        """Record a reference the run cannot use, so the species is unscreened rather than assumed clean."""
        self._species_screening_shortfalls[species] = reason
        console.print(f"❌ No usable reference for {species}: {reason}")
        logger.error(f"Refusing to screen {species} with {request.value}: {reason}")
        return ReferenceRejection(species=species, identity=request.value, reason=reason)

    async def _resolve_screening_reference(
        self, request: ReferenceRequest
    ) -> tuple[ScreeningReference | None, ReferenceRejection | None]:
        """Resolve one request into a typed ``{species, kind, identity, index}`` reference.

        The one door: an explicit index override and a resolved default run through this method
        alike, so both build the species label, the transcript->gene index the classifier reads and
        the cDNA FASTA repeat detection reuses. An override that bypassed all three would leave a cDNA
        file aligning successfully and then classifying against nothing.
        """
        console.print(f"📚 Transcriptome reference: {request.value} ({request.state.value}, {request.form.value})")
        declared = request.declared_species
        fallback_species = normalize_species_name(declared) if declared else UNRESOLVED_SPECIES

        if request.form is ReferenceForm.PREBUILT_INDEX:
            payload: dict[str, Any] | None = self._adopt_prebuilt_index(request)
        else:
            try:
                payload = await self._prepare_transcriptome_database(request.value, self._resolve_filter_spec())
            except Exception as exc:  # pragma: no cover - defensive logging path
                logger.exception("Failed to prepare transcriptome database")
                return None, self._reject_reference(request, fallback_species, f"preparation failed: {exc}")
            if payload and not payload.get("fasta"):
                payload = None

        if not payload:
            return None, self._reject_reference(
                request, fallback_species, "the reference could not be fetched, cached or read"
            )

        # Resolved before any refusal below, so a rejection names the species that went unscreened
        # rather than 'unknown': the shortfall map is what reports the missing species, and blaming a
        # phantom one is how "unscreened" started reading as "clean" in the first place.
        species, authority, conflict = resolve_reference_species(
            declared=declared,
            source_species=payload.get("source_species"),
            header_species=payload.get("header_species"),
        )
        if conflict:
            logger.warning("Species labels disagree for %s: %s", request.value, conflict)
            console.print(f"⚠️  Species labels disagree for {request.value}: {conflict}")
        if authority is SpeciesAuthority.UNRESOLVED:
            logger.warning(
                "Could not establish the species of %s: neither the caller nor its headers name one, so it is "
                "labelled '%s'. Orthology cannot be resolved for it, so cross-species hits would be reported as "
                "unqualified off-targets.",
                request.value,
                UNRESOLVED_SPECIES,
            )
        else:
            logger.info("Species '%s' for %s, from its %s", species, request.value, authority.value)

        # An index build that was attempted and failed leaves no index, and the FASTA is not a
        # substitute: handed to Nextflow as an index prefix it aligns nothing, which the pipeline
        # reports as success. Refuse the reference and record why, so the species is unscreened
        # rather than silently screened against nothing.
        index_build_error = payload.get(INDEX_BUILD_ERROR_KEY)
        if index_build_error and not payload.get("index"):
            return None, self._reject_reference(request, species, str(index_build_error))

        fasta = payload.get("fasta")
        if not fasta:
            return None, self._reject_reference(
                request,
                species,
                "no readable plain-text cDNA FASTA was found beside the index prefix (a gzipped or binary "
                "neighbour does not count), so its hits could not be resolved to genes: put the cDNA FASTA "
                f"at '{request.value}' (or '{request.value}.fa'), or pass it to --transcriptome-fasta and "
                "let the run index it",
            )

        # A transcript->gene index for EVERY screened species (not just human), so orthologs can be
        # recognized across species, and the FASTA is stashed for repeat detection to reuse.
        if not self._transcript_index.for_species(species):
            self._transcript_index.build(species, Path(fasta))
        self._species_cdna_fasta.setdefault(species, Path(fasta))

        # An index when there is one, otherwise the FASTA -- and then the reference travels on the
        # FASTA parameter, so the pipeline indexes it. Naming a FASTA as an index prefix instead
        # aligns nothing, which the pipeline reports as a completed screen.
        prebuilt = payload.get("index")
        index = prebuilt or fasta
        console.print(
            f"✨ Screening reference resolved: {species} ({authority.value}) → {Path(index).name}"
            f"{'' if prebuilt else ' (the pipeline will build the index)'}"
        )
        return (
            ScreeningReference(
                species=species,
                kind=self.config.screening_kind,
                identity=request.value,
                index=str(index),
                species_authority=authority,
                form=request.form,
                state=request.state,
                reason=request.reason,
                fasta=str(fasta),
                needs_index_build=prebuilt is None,
                # `identity` names the source; this names the bytes. A prebuilt index adopted from the
                # caller has no cache entry, so it stays None rather than being invented here.
                identity_evidence=payload.get("identity"),
            ),
            None,
        )

    def _resolve_filter_spec(self) -> list[str] | None:
        """Parse the configured transcriptome filter names, continuing unfiltered when invalid."""
        if not self.config.transcriptome_filter:
            return None

        from sirnaforge.data.transcriptome_filter import get_filter_spec  # noqa: PLC0415

        try:
            filter_spec = get_filter_spec(self.config.transcriptome_filter)
        except ValueError as exc:
            logger.error(f"Invalid transcriptome filter specification: {exc}")
            console.print(f"⚠️  Invalid filter specification: {exc}")
            return None
        if filter_spec:
            console.print(f"🔍 Applying transcriptome filters: {', '.join(filter_spec)}")
        return filter_spec

    def _resolve_active_screen_species(self, params: Mapping[str, Any]) -> list[str]:
        """Filter the requested screen species down to those with a reference the pipeline received."""
        requested = [species.strip() for species in self.config.screen_species if species.strip()]
        available: set[str] = set()
        # Both keys are read because main.nf mixes them into ONE alignment channel: transcriptome
        # FASTAs are indexed in the container, prebuilt indices are used as they are, and either way
        # the species really is screened. Omitting one dropped species that had in fact been aligned
        # from the conservation denominator.
        for key in ("transcriptome_indices", "transcriptome_fastas"):
            raw_value = params.get(key) or self.config.nextflow_config.get(key)
            available.update(self._parse_species_entries(raw_value))

        if available:
            filtered = [species for species in requested if species in available]
            # Dropping a requested species here is a completeness fact about the run, not a detail
            # of list construction: recorded so it reaches the warnings and the published summary.
            for species in requested:
                if species not in available:
                    self._species_screening_shortfalls[species] = (
                        "requested for screening but no transcriptome index or FASTA was resolved for it, "
                        "so it was never submitted to the aligner"
                    )
                    logger.warning(
                        f"Species '{species}' was requested for off-target screening but no reference was "
                        "resolved for it; it will not be screened and its hit counts are unknown."
                    )
            for species in sorted(available):
                if species not in filtered:
                    filtered.append(species)
            # A species that reached the pipeline parameters without going through the resolver -- a
            # raw nextflow_config passthrough -- has no transcript index here, so its alignments can
            # only ever publish as undetermined. Reported rather than refused: the caller wrote a
            # pipeline parameter by hand, which the resolver is not asked to police.
            for species in filtered:
                if not self._transcript_index.for_species(species):
                    logger.warning(
                        f"Species '{species}' reaches the aligner without a reference resolved on this host, so no "
                        "transcript index was built for it: its hits cannot be resolved to genes and will publish "
                        "as undetermined with species_index_missing set."
                    )
        else:
            # No reference resolved for ANYTHING, so no transcriptome screen was configured at all
            # (miRNA-only, or --skip-off-targets). Recording a per-species shortfall here would
            # blame each species for a screen nobody asked for; the aggregate reports that case as
            # "TRANSCRIPTOME ANALYSIS STATUS: NOT PERFORMED" and the run still reports partial via
            # the "no alignment evidence for any species" guard in _process_nextflow_results.
            filtered = requested

        # Remembered because this list, not config.screen_species, is what gets screened: scoring
        # conservation against the shorter config list yielded a denominator smaller than its
        # numerator, which aborted post-screen scoring for the whole candidate.
        self._active_screen_species = list(filtered)
        return filtered

    @staticmethod
    def _parse_species_entries(raw_value: Any) -> set[str]:
        """Extract species identifiers from 'species:path' style strings."""
        species: set[str] = set()
        if not raw_value:
            return species

        values: list[str]
        if isinstance(raw_value, str):
            values = [raw_value]
        elif isinstance(raw_value, list | tuple | set):
            iterable = cast(Iterable[Any], raw_value)
            values = [str(entry) for entry in iterable]
        else:
            values = [str(raw_value)]

        for value in values:
            for token in value.split(","):
                entry = token.strip()
                if not entry:
                    continue
                if ":" in entry:
                    species.add(entry.split(":", 1)[0].strip())
                else:
                    species.add(entry)
        return species

    def _prepare_nextflow_cache(
        self,
        nf_config: NextflowConfig,
        screen_species: Sequence[str],
        additional_params: Mapping[str, Any],
        pipeline_revision: str,
    ) -> dict[str, Any]:
        """Configure cached work and home directories for Nextflow runs.

        Args:
            nf_config: Nextflow configuration
            screen_species: Species the screen covers (used in cache key)
            additional_params: Additional pipeline parameters
            pipeline_revision: Git revision of pipeline

        Returns:
            Cache metadata dictionary
        """
        cache_root = resolve_cache_subdir("nextflow")
        home_dir = cache_root / "home"
        work_root = cache_root / "work"
        home_dir.mkdir(parents=True, exist_ok=True)
        work_root.mkdir(parents=True, exist_ok=True)

        payload: dict[str, Any] = {
            "pipeline_revision": pipeline_revision,
            "profile": nf_config.profile,
            "max_cpus": nf_config.max_cpus,
            "max_memory": nf_config.max_memory,
            "max_time": nf_config.max_time,
            "screen_species": sorted(screen_species),
            # ``evidence_plan`` is deliberately excluded: it is a path inside this run's own output
            # directory, so keying the shared work dir on it would give every run a fresh one and
            # rebuild every index. The plan still reaches the pipeline as a staged path input, so
            # Nextflow re-runs the aggregation task by itself when its content changes (#100).
            "additional_params": self._normalize_param_dict(
                {key: value for key, value in additional_params.items() if key != "evidence_plan"}
            ),
            "extra_params": self._normalize_param_dict(nf_config.extra_params),
        }
        cache_key = stable_cache_key(payload)
        work_dir = work_root / cache_key
        work_dir.mkdir(parents=True, exist_ok=True)

        metadata: dict[str, Any] = {
            "payload": payload,
            "work_dir": str(work_dir),
            "nxf_home": str(home_dir),
            "created_at": time.time(),
        }
        metadata_file = work_dir / "cache_metadata.json"
        try:
            metadata_file.write_text(json.dumps(metadata, indent=2))
        except OSError as exc:
            logger.debug(f"Unable to write Nextflow cache metadata: {exc}")

        nf_config.work_dir = work_dir
        nf_config.nxf_home = home_dir

        cache_info = {
            "cache_key": cache_key,
            "work_dir": str(work_dir),
            "nxf_home": str(home_dir),
            "pipeline_revision": pipeline_revision,
        }
        self._nextflow_cache_info = cache_info
        return cache_info

    @staticmethod
    def _normalize_param_dict(params: Mapping[str, Any]) -> dict[str, Any]:
        """Convert values to JSON-friendly primitives for hashing."""
        normalized: dict[str, Any] = {}
        for key in sorted(params):
            value = params[key]
            if isinstance(value, Path):
                normalized[key] = str(value)
            elif isinstance(value, list | tuple | set):
                iterable = cast(Iterable[Any], value)
                normalized[key] = [SiRNAWorkflow._stringify_param(entry) for entry in iterable]
            else:
                normalized[key] = SiRNAWorkflow._stringify_param(value)
        return normalized

    @staticmethod
    def _stringify_param(value: Any) -> Any:
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, list | tuple | set):
            iterable = cast(Iterable[Any], value)
            return [SiRNAWorkflow._stringify_param(entry) for entry in iterable]
        return value

    def _publish_nextflow_work_reference(self) -> None:
        """Write workdir pointer and optionally expose a symlink inside output."""
        if not self._nextflow_cache_info:
            return

        info = self._nextflow_cache_info
        work_dir = Path(info["work_dir"])
        off_target_dir = self.config.output_dir / "off_target"
        off_target_dir.mkdir(parents=True, exist_ok=True)

        ref_file = off_target_dir / "NEXTFLOW_WORKDIR.txt"
        lines = [
            "Nextflow intermediate cache",
            f"Work directory: {work_dir}",
            f"Cache key: {info.get('cache_key')}",
            f"Pipeline revision: {info.get('pipeline_revision', 'unknown')}",
        ]
        try:
            ref_file.write_text("\n".join(lines) + "\n")
        except OSError as exc:
            logger.debug(f"Unable to write Nextflow work reference: {exc}")

        link_path = off_target_dir / "nextflow_work"
        if self.config.keep_nextflow_work:
            self._ensure_symlink(link_path, work_dir)
        elif link_path.exists() or link_path.is_symlink():
            try:
                if link_path.is_dir() and not link_path.is_symlink():
                    shutil.rmtree(link_path)
                else:
                    link_path.unlink()
            except OSError:
                pass

    @staticmethod
    def _ensure_symlink(link_path: Path, target: Path) -> None:
        """Ensure link_path points at target, replacing existing artifacts."""
        try:
            if link_path.exists() or link_path.is_symlink():
                try:
                    if link_path.resolve() == target.resolve():
                        return
                except OSError:
                    pass
                if link_path.is_dir() and not link_path.is_symlink():
                    shutil.rmtree(link_path)
                else:
                    link_path.unlink()
            link_path.symlink_to(target, target_is_directory=True)
        except OSError as exc:
            logger.debug(f"Unable to create Nextflow workdir symlink: {exc}")

    @staticmethod
    def _read_results_json(results_dir: Path, filename: str) -> dict[str, Any] | None:
        """One JSON object published by the pipeline, or None when there is nothing readable.

        The aggregated subdirectory wins over the results root, because that is where publishDir puts
        the aggregate; an unreadable or non-object payload is None, never a partially-trusted dict.
        """
        results_path = Path(results_dir)
        search_roots = [results_path / "aggregated", results_path]
        for root in search_roots:
            candidate = root / filename
            if not candidate.exists():
                continue
            try:
                with candidate.open() as fh:
                    payload = json.load(fh)
            except Exception as exc:  # pragma: no cover - defensive logging path
                logger.warning(f"Failed to read aggregated summary {candidate}: {exc}")
                return None
            if isinstance(payload, dict):
                return cast(dict[str, Any], payload)
            logger.warning(f"Aggregated summary {candidate} is not a JSON object; skipping")
            return None
        return None

    def _load_offtarget_aggregates(self, results_dir: Path) -> dict[str, Any]:
        """Load aggregated Nextflow summary JSON files when available."""
        aggregated: dict[str, Any] = {}

        transcriptome_summary = self._read_results_json(results_dir, "combined_summary.json")
        if transcriptome_summary:
            aggregated["transcriptome"] = transcriptome_summary

        mirna_summary = self._read_results_json(results_dir, "combined_mirna_summary.json")
        if mirna_summary:
            aggregated["mirna"] = mirna_summary

        return aggregated

    async def _resolve_screening_references(self, additional_params: dict[str, Any]) -> bool:
        """Resolve every screening reference request and write the resolved set into the parameters.

        One resolver for both doors, so the pipeline receives one species list and one index list
        that cannot disagree with the metadata the classifier holds. Requests that resolved to
        nothing are kept as rejections, not dropped.
        """
        requests = self.config.screening_requests
        rejections = list(self.config.screening_request_rejections)
        for rejection in self.config.screening_request_rejections:
            logger.info("Screening reference not resolved: %s (%s)", rejection.identity, rejection.reason)
        if not requests:
            reason = self.config.transcriptome_selection.disabled_reason or "no screening reference requested"
            console.print(f"ℹ️  Transcriptome off-target disabled ({reason})")
            self._screening_references = ScreeningReferenceSet(
                kind=self.config.screening_kind,
                rejections=tuple(rejections),
                requested_species=tuple(self.config.screen_species),
            )
            return False

        references: list[ScreeningReference] = []
        for request in requests:
            resolved, refusal = await self._resolve_screening_reference(request)
            if refusal is not None:
                rejections.append(refusal)
            if resolved is None:
                continue
            if resolved.species not in self.config.screen_species:
                self.config.screen_species.append(resolved.species)
            references.append(resolved)

        self._screening_references = ScreeningReferenceSet(
            kind=self.config.screening_kind,
            references=tuple(references),
            rejections=tuple(rejections),
            requested_species=tuple(self.config.screen_species),
        )
        if not references:
            return False

        # An already-built index and a FASTA still to be indexed travel on different parameters: the
        # pipeline reads one as a prefix and indexes the other, and naming a FASTA as a prefix aligns
        # nothing while reporting success.
        for key, resolved_value in (
            ("transcriptome_indices", self._screening_references.index_parameter),
            ("transcriptome_fastas", self._screening_references.fasta_parameter),
        ):
            existing = additional_params.get(key)
            merged = [token.strip() for token in str(existing).split(",")] if existing else []
            merged = [entry for entry in merged if entry]
            for entry in resolved_value.split(","):
                if entry and entry not in merged:
                    merged.append(entry)
            if merged:
                additional_params[key] = ",".join(merged)
        additional_params["transcriptome_species"] = ",".join(self._screening_references.species)
        return True

    @staticmethod
    def _plan_search_settings(additional_params: Mapping[str, Any]) -> dict[str, str | int | float | bool | None]:
        """The alignment settings a plan entry needs to be reproducible, as stated for this run."""
        settings: dict[str, str | int | float | bool | None] = {}
        for key in ("max_hits", "bwa_k", "bwa_T", "seed_start", "seed_end"):
            if key not in additional_params:
                continue
            value = additional_params[key]
            settings[key] = value if isinstance(value, str | int | float | bool) or value is None else str(value)
        return settings

    def _requested_screening_units(self) -> tuple[tuple[tuple[str, str | None], ...], tuple[str, ...]]:
        """What this run asked to screen: transcriptome ``(species, reference_id)`` pairs, then miRNA species.

        Requested, not resolved. Every species the run named is here whether or not a reference for it
        can be fetched, which is the whole point: a species dropped by ``_reject_reference`` has to
        reconcile as a failure, and it can only do that if the plan already holds an entry for it.

        A channel nobody asked for gets no entry at all -- that is what keeps ``NOT_REQUESTED``
        distinguishable from ``FAILED``. So transcriptome units exist only when a screening reference
        was requested, and miRNA units only when both a database and a species list were given, which
        is the same condition that populates the pipeline's ``mirna_db``/``mirna_species`` parameters.
        """
        transcriptome: dict[str, str | None] = {}
        if self.config.screening_requests:
            # The request's own identity is the reference_id: the reconciler adopts the plan's value,
            # and an in-container emitter never knows it.
            for request in self.config.screening_requests:
                if request.declared_species:
                    transcriptome.setdefault(normalize_species_name(request.declared_species), request.value)
            for species in self.config.screen_species:
                transcriptome.setdefault(normalize_species_name(species), None)

        mirna_species = (
            tuple(dict.fromkeys(normalize_species_name(species) for species in self.config.mirna_species))
            if self.config.mirna_database and self.config.mirna_species
            else ()
        )
        return tuple(transcriptome.items()), mirna_species

    def _record_screening_plan(self, input_fasta: Path, additional_params: dict[str, Any]) -> None:
        """Record what this run intends to screen, and hand the plan to the pipeline.

        Called BEFORE ``_resolve_screening_references``, because that is where a reference can be
        dropped: a plan built from what survived resolution cannot represent a missing one, so the
        species whose index failed to build simply vanished from every artifact (#100). The plan is
        also serialized to ``screening_plan.json`` and passed on as ``evidence_plan``, so the
        pipeline's own aggregation reconciles against the same expectation this workflow holds.

        The guide-set digest is the submitted FASTA's own hash: two screens of one reference with
        different guide sets are different evidence, and joining them would attribute one screen's
        counts to the other's guides.
        """
        transcriptome, mirna_species = self._requested_screening_units()
        self._guide_set_digest = guide_set_digest(input_fasta)
        self._screening_plan = build_plan(
            guide_set_digest=self._guide_set_digest,
            transcriptome=transcriptome,
            mirna_species=mirna_species,
            search_settings=self._plan_search_settings(additional_params),
        )

        plan_file = self.config.output_dir / "screening_plan.json"
        try:
            plan_file.write_text(self._screening_plan.model_dump_json(indent=2))
        except OSError as exc:  # pragma: no cover - defensive logging path
            logger.warning(f"Could not write the screening plan to {plan_file}: {exc}")
            return
        additional_params["evidence_plan"] = str(plan_file)

    def _summarize_screening_references(self) -> dict[str, Any]:
        """The published reference record: what resolved, over which species, and what did not.

        ``scope`` is the resolved screen's species set as a :class:`FilterScope`: what the run
        planned to screen, fixed at resolution. It is not a coverage report -- a species whose
        alignment published nothing is subtracted post-run in ``filtering_stats.unscreened_species``,
        not from here. No gate reads it in 0.7.1; applying a scope is #101's.

        ``screening_evidence`` is the plan reconciled against what the run published: what was asked
        for, what each unit actually did, and where each answer came from. It is additive behind the
        same ``is not None`` guard as ``screening_plan``, so it reaches every ``workflow_summary.json``
        write site for free (#100).
        """
        summary: dict[str, Any] = {
            "transcriptome": self.config.transcriptome_selection.to_metadata(),
            "screening": self._screening_references.to_metadata(),
            "scope": self._screening_references.scope().model_dump(mode="json"),
        }
        if self._screening_plan is not None:
            summary["screening_plan"] = self._screening_plan.model_dump(mode="json")
        if self._screening_evidence is not None:
            summary["screening_evidence"] = reconciliation_payload(self._screening_evidence)
        return summary

    def _log_nextflow_targets(
        self,
        active_species: Sequence[str],
        has_transcriptome: bool,
        additional_params: Mapping[str, Any],
    ) -> None:
        """Emit console updates about the screening references handed to the pipeline."""
        if active_species:
            console.print(f"🔭 Nextflow transcriptome species: {', '.join(active_species)}")
        else:
            console.print("🔭 Nextflow transcriptome species: (none)")

        if has_transcriptome:
            transcriptome_species = str(additional_params.get("transcriptome_species", ""))
            pretty = transcriptome_species or "unspecified"
            console.print(f"🗂️  Transcriptome indices resolved for: {pretty}")

    async def _run_nextflow_offtarget_analysis(
        self,
        candidates: list[SiRNACandidate],
        input_fasta: Path,
        additional_params: dict[str, Any] | None = None,
        has_transcriptome: bool | None = None,
    ) -> dict[str, Any]:
        """Run Nextflow-based off-target analysis.

        Args:
            candidates: Candidates to screen.
            input_fasta: Deduplicated FASTA input for the Nextflow pipeline.
            additional_params: Pre-configured Nextflow parameters (screening references already
                resolved, and the plan already recorded), or None to do both here.
            has_transcriptome: Paired with additional_params; ignored when that is None.
        """
        if additional_params is None:
            additional_params = dict(self.config.nextflow_config)
            # Recorded before resolution, never after: see _record_screening_plan (#100).
            self._record_screening_plan(input_fasta, additional_params)
            has_transcriptome = await self._resolve_screening_references(additional_params)
        has_transcriptome = bool(has_transcriptome)
        active_species = self._resolve_active_screen_species(additional_params)
        has_transcriptome = has_transcriptome or bool(additional_params.get("transcriptome_indices"))
        self._log_nextflow_targets(active_species, has_transcriptome, additional_params)

        if not active_species and not has_transcriptome:
            console.print("ℹ️  No transcriptome indices configured; skipping Nextflow run")
            return await self._basic_offtarget_analysis(candidates)

        runner, _ = self._setup_nextflow_runner(active_species, additional_params)

        if not self._validate_nextflow_environment(runner):
            # No screen ran, so the gates report unknown rather than never having been reached (#106).
            self._gate_without_screening_evidence(candidates)
            return self._publish_run_status({"status": "skipped", "reason": "nextflow_unavailable"})

        # Execute pipeline
        console.print("🚀 Running embedded Nextflow off-target analysis...")
        nf_output_dir = self.config.output_dir / "off_target" / "results"

        results = await runner.run_offtarget_analysis(
            input_file=input_fasta,
            output_dir=nf_output_dir,
            screen_species=active_species,
            additional_params=additional_params,
            show_progress=True,
        )

        self._publish_nextflow_work_reference()

        if results["status"] == "completed":
            return await self._process_nextflow_results(candidates, nf_output_dir, results)

        console.print(f"❌ Nextflow pipeline failed: {results}")
        return await self._basic_offtarget_analysis(candidates)

    async def run_nextflow_offtarget_analysis(
        self,
        candidates: list[SiRNACandidate],
        input_fasta: Path,
    ) -> dict[str, Any]:
        """Public wrapper for Nextflow off-target analysis execution."""
        return await self._run_nextflow_offtarget_analysis(candidates=candidates, input_fasta=input_fasta)

    def _setup_nextflow_runner(
        self,
        screen_species: Sequence[str],
        additional_params: Mapping[str, Any],
    ) -> tuple[NextflowRunner, dict[str, Any]]:
        """Configure Nextflow runner with user settings and cached workdirs.

        Args:
            screen_species: Species the screen covers
            additional_params: Additional pipeline parameters

        Returns:
            Configured NextflowRunner and cache metadata
        """
        # Auto-detect environment to use appropriate profile
        # This will automatically switch to 'local' profile when running inside a container
        nf_config = NextflowConfig.auto_configure()

        # Apply user overrides from workflow config, BUT preserve auto-detected profile
        # unless explicitly overridden AND we're not in a container
        if self.config.nextflow_config:
            for key, value in self.config.nextflow_config.items():
                # Don't allow profile override when running in container
                # (container detection takes precedence for safety)
                if key == "profile" and nf_config.is_running_in_docker():
                    logger.warning(
                        f"Ignoring user profile override '{value}' - running in container, using 'local' profile"
                    )
                    continue
                setattr(nf_config, key, value)

        # Log the execution environment for debugging
        env_info = nf_config.get_environment_info()
        logger.info(f"Nextflow execution: {env_info.get_execution_summary()}")

        runner = NextflowRunner(nf_config)
        cache_info = self._prepare_nextflow_cache(
            nf_config=nf_config,
            screen_species=screen_species,
            additional_params=additional_params,
            pipeline_revision=runner.get_pipeline_revision(),
        )
        return runner, cache_info

    def _validate_nextflow_environment(self, runner: NextflowRunner) -> bool:
        """Validate Nextflow installation and workflow files."""
        validation = runner.validate_installation()
        if not validation["nextflow"]:
            console.print("⚠️  Nextflow not available; off-target analysis will be skipped")
            return False
        if not validation["workflow_files"]:
            console.print("⚠️  Nextflow workflows not found; off-target analysis will be skipped")
            return False
        return True

    async def _process_nextflow_results(
        self, candidates: list[SiRNACandidate], output_dir: Path, results: dict[str, Any]
    ) -> dict[str, Any]:
        """Process and map Nextflow pipeline results to candidates."""
        console.print("✅ Nextflow pipeline completed successfully")
        parsed = await self._parse_nextflow_results(output_dir)
        aggregated_views = self._load_offtarget_aggregates(output_dir)
        published_status = "completed"
        workflow_warnings: list[str] = []

        tx_summary = aggregated_views.get("transcriptome") if aggregated_views else None
        if tx_summary:
            for warning_msg in self._transcriptome_shortfall_warnings(tx_summary):
                published_status = "partial"
                console.print(warning_msg)
                workflow_warnings.append(warning_msg)

        # Shortfalls decided before Nextflow ran are reported on the same footing as ones the
        # aggregate found: a species dropped for want of a reference appeared in no artifact at all.
        for species, reason in sorted(self._species_screening_shortfalls.items()):
            published_status = "partial"
            warning_msg = f"⚠️  '{species}' was not screened: {reason}"
            console.print(warning_msg)
            workflow_warnings.append(warning_msg)

        # POSITIVE evidence, deliberately not the aggregate's self-reported missing_species: this
        # method reports "completed" whenever the output directory merely exists, and
        # _load_offtarget_aggregates returns {} when combined_summary.json is absent — which is
        # exactly the shape of a run where aggregation itself failed. missing_species is then empty
        # and the run looks complete. Only a species with a published alignment file has earned the
        # reading "no hits here means clean".
        screened_species = self._species_with_alignment_evidence(tx_summary)
        if not screened_species:
            published_status = "partial"
            warning_msg = (
                "⚠️  No transcriptome alignment evidence for any species (no aggregated summary, or "
                "miRNA-only mode): off-target counts are unknown, so candidates keep their design-time scores."
            )
            console.print(warning_msg)
            workflow_warnings.append(warning_msg)

        # Integrate off-target results into candidates with filtering. The species that produced
        # alignments are passed through because a candidate with no hits from an alignment that
        # never ran must not be scored as if it had come back clean.
        filter_criteria = getattr(self.config.design_params, "offtarget_filters", None) or OffTargetFilterCriteria()
        # Orthologue evidence is resolved here, in the async layer, and handed to the synchronous
        # classifier as a plain gene-ID set: gene symbols cannot decide orthology (mouse TP53 is
        # Trp53), and classify_hit is pure by contract. Degrades to the symbol heuristic on failure.
        ortholog_mapping = await self._resolve_ortholog_mapping(screened_species, parsed)
        mirna_screened = self._mirna_channel_completed(aggregated_views.get("mirna") if aggregated_views else None)
        if not mirna_screened:
            # Warned, but not a status downgrade: ``published_status`` speaks for the evidence the run
            # *requires*, and every miRNA pair is exploratory in 0.7.1. The consequence is carried as
            # UNKNOWN on the miRNA gates and as a flag in filtering_stats.
            warning_msg = (
                "⚠️  No miRNA seed evidence was published for this run: the miRNA gates report unknown rather "
                "than a clean scan, and their observed values stay empty."
            )
            console.print(warning_msg)
            workflow_warnings.append(warning_msg)

        # Reconciled BEFORE integration, because the pairs the gates read come from it: what the run
        # planned, against what it published. Replaces _species_with_alignment_evidence as the
        # authority for completeness -- that heuristic survives inside it as the fallback for a
        # result directory written before evidence existed (#100).
        self._screening_evidence = self._reconcile_screening_evidence(
            output_dir, screened_species=screened_species, mirna_screened=mirna_screened
        )
        for key in self._screening_evidence.keys_with_status(EvidenceStatus.FAILED):
            logger.error(f"No {key[0]} evidence for '{key[1]}': that unit was planned and published nothing.")

        updated_candidates, stats = self._integrate_offtarget_results(
            candidates,
            parsed,
            filter_criteria,
            screened_species=screened_species,
            ortholog_gene_ids=ortholog_mapping.all_gene_ids,
            mirna_screened=mirna_screened,
            unresolved_orthology_species=ortholog_mapping.unresolved_species,
        )
        stats["mirna_channel_screened"] = mirna_screened
        stats["orthology"] = ortholog_mapping.summary()
        # Held for the manifest: this is the only place that knows whether a lookup happened.
        self._orthology_mapping = ortholog_mapping
        # A species dropped before Nextflow is unscreened in exactly the sense this field names, so
        # it belongs in it; the reasons are published beside it rather than only logged.
        stats["species_screening_shortfalls"] = dict(self._species_screening_shortfalls)
        stats["unscreened_species"] = sorted(
            set(cast(list[str], stats.get("unscreened_species") or [])) | set(self._species_screening_shortfalls)
        )
        workflow_warnings.extend(self._persist_hit_classifications(parsed))
        self._log_offtarget_statistics(stats, aggregated_views, output_dir)

        # Map parsed results for return structure, keyed by candidate id. The LOOKUP is by
        # screen_query_id: the aligner saw one record per distinct guide, so a deduplicated
        # candidate's evidence lives under its representative's id and looking it up by `id`
        # published a fabricated zero for every non-representative (32,463 of 34,863 ids on the
        # frozen baseline) while candidates_all.csv carried the real count for the same ids.
        # Two scalars per candidate, never the rows behind them: carrying `hits` produced a 2.9 GiB
        # workflow_summary.json whose other keys totalled 14 KB. Nothing reads the rows from this map --
        # every consumer of the alignments works off `parsed` -- and the detail stays on disk in the
        # tables named by `detail_files`.
        mapped = {}
        for c in updated_candidates:
            entry = parsed.get("results", {}).get(c.screen_query_id or c.id)
            mapped[c.id] = {
                # The candidate's own liability count, not the ingest tally: the raw entry counts
                # on-target, ortholog, repeat and miRNA rows too, so taking it here gave one field
                # name two different definitions across two artifacts of the same run.
                "off_target_count": c.off_target_count,
                "off_target_score": entry.get("off_target_score", 0.0) if entry else 0.0,
            }

        return self._publish_run_status(
            {
                "status": published_status,
                "method": "embedded_nextflow",
                "output_dir": str(output_dir),
                "results": mapped,
                "detail_files": self._offtarget_detail_files(parsed, output_dir),
                "execution_metadata": results,
                "filtering_stats": stats,
                "aggregated": aggregated_views,
                "warnings": workflow_warnings,
            },
            # An output directory that does not exist is a failure to execute, not a thin result: the
            # runner reported exit 0 and staged nothing, so there is no partial answer to read. The
            # published `status` stays whatever this method decided, per #100's beside-not-instead rule.
            reason=_EXECUTION_ERROR_REASON if parsed.get("status") == "missing" else None,
        )

    @staticmethod
    def _offtarget_detail_files(parsed: Mapping[str, Any], output_dir: Path) -> dict[str, list[str]]:
        """Where the per-alignment rows actually live, so the summary can point instead of copy.

        Paths are the ones the parser read, not a glob a reader has to re-derive: a guessed pattern
        goes stale the moment the pipeline's layout changes, and the summary would then name files
        that do not exist.
        """
        transcriptome = [
            str(table["path"])
            for table in cast(list[dict[str, Any]], parsed.get("transcriptome_hit_tables") or [])
            if table.get("path")
        ]
        mirna = [str(path) for path in cast(list[Path], parsed.get("mirna_hit_files") or [])]
        return {
            "transcriptome": sorted(transcriptome),
            "mirna": sorted(mirna),
            "root": [str(output_dir)],
        }

    @staticmethod
    def _persist_hit_classifications(parsed: Mapping[str, Any]) -> list[str]:
        """Write the class, symbols, shortfall flags and ortholog evidence tier onto every table read.

        Every transcriptome file the parser read is rewritten with the classification columns appended, so
        a reviewer can tell an on-target isoform alignment from a liability without re-deriving
        anything. A header-only table is rewritten too: the published schema must not depend on
        whether the run produced any hits, or a consumer selecting ``hit_class`` fails on exactly
        the clean runs.

        Returns the reconciliation warnings, which the caller surfaces alongside the run status.
        Rows are annotated by ``_classify_orphan_hit_rows`` before this runs, so there is no
        unannotated-row branch here: an unannotated row is a reconciliation failure, reported as
        one, rather than a silent skip.
        """
        tables = cast(list[dict[str, Any]], parsed.get("transcriptome_hit_tables") or [])
        persisted = 0
        tally: dict[str, int] = {}
        for table in tables:
            tsv_path = table.get("path")
            rows = cast(list[dict[str, Any]], table.get("rows") or [])
            fieldnames = cast(list[str], table.get("fieldnames") or [])
            if not tsv_path or not fieldnames:
                continue
            try:
                written = write_classified_hits(Path(tsv_path), rows, fieldnames)
            except OSError as exc:
                logger.warning(f"Could not persist hit classification to {tsv_path}: {exc}")
                continue
            persisted += written
            for name, count in count_persisted_classes(rows).items():
                tally[name] = tally.get(name, 0) + count
            logger.info(f"Wrote {'/'.join(CLASSIFICATION_COLUMNS)} onto {tsv_path} for {written} rows")

        if tables:
            console.print(
                f"🏷️  Persisted hit classification for {persisted} alignment(s) across "
                f"{len(tables)} file(s) ({', '.join(f'{name}={count}' for name, count in sorted(tally.items()))})"
            )
        return SiRNAWorkflow._reconcile_persisted_hits(parsed, persisted, sum(tally.values()))

    @staticmethod
    def _reconcile_persisted_hits(parsed: Mapping[str, Any], persisted: int, classified: int) -> list[str]:
        """Check the published rows against the hits the candidate counters were built from.

        The two are counted from different structures — ``results`` by the row ingest,
        ``persisted`` by the tables that were written — so a path that feeds candidates hits it
        never publishes shows up here. That is the failure this check exists for: an aggregated
        table rejected upstream left the fallback ingesting real per-species rows into candidates
        while the published table stayed header-only, so the hit table reported no liabilities
        beside candidates carrying dozens each.
        """
        counted = 0
        for entry in cast(dict[str, dict[str, Any]], parsed.get("results") or {}).values():
            for hit in cast(list[Mapping[str, Any]], entry.get("hits") or []):
                if "mirna_id" not in hit and "database" not in hit:
                    counted += 1

        warnings: list[str] = []
        if persisted < counted:
            warnings.append(
                f"⚠️  Hit table/candidate mismatch: {counted} transcriptome hit(s) reached the candidate counters "
                f"but only {persisted} row(s) were published. The published hit table under-reports liabilities; "
                "treat off_target_count, not the table, as the count for this run."
            )
        if classified < persisted:
            warnings.append(
                f"⚠️  {persisted - classified} published hit row(s) carry no usable classification; "
                "their hit_class cell is not one of the taxonomy's values."
            )
        for warning in warnings:
            logger.error(warning)
            console.print(warning)
        return warnings

    @staticmethod
    def _species_with_alignment_evidence(tx_summary: Mapping[str, Any] | None) -> list[str]:
        """Species with a published transcriptome alignment file, per the aggregated summary.

        Returns the species for which the aggregator actually saw analysis output — the only
        positive evidence available that an alignment ran. Returns an empty list when there is no
        such evidence for any species, which covers three failures that would otherwise read as a
        clean complete run:

        - the aggregate is missing entirely (aggregation never ran, so nothing reported anything);
        - the aggregate exists but names no species (a hand-written or legacy summary);
        - miRNA-only mode, where sirna_offtarget_analysis.nf derives the species list from
          the resolved index channel and falls back to '' — no transcriptome alignment happened at all, so
          the off-target term has nothing to stand on.

        ``species_screened`` is preferred over any file tally because a file can be discovered and
        still be unusable: a 0-byte ``*_analysis.tsv`` (what offtarget_analysis.nf's stub emits) or
        one the schema rejects counts in ``species_file_counts`` while contributing no alignment.
        An explicitly empty ``species_screened`` is evidence, not a missing key, so it is honoured.
        """
        if not tx_summary:
            return []
        screened = tx_summary.get("species_screened")
        if screened is not None:
            return [str(species) for species in cast(list[Any], screened)]
        # Older summaries carry no screened list; their per-species file tally is the same evidence,
        # coarser, and cannot see a file that was discovered and then rejected.
        usable_counts = cast(dict[str, int], tx_summary.get("usable_species_file_counts") or {})
        file_counts = usable_counts or cast(dict[str, int], tx_summary.get("species_file_counts") or {})
        if file_counts:
            return [species for species, count in file_counts.items() if count]
        # Older summaries carry no per-species counts; species_analyzed minus the species the
        # aggregator itself flagged as producing no files is the same evidence, coarser.
        analyzed = [str(species) for species in cast(list[Any], tx_summary.get("species_analyzed") or [])]
        missing = {str(species) for species in cast(list[Any], tx_summary.get("missing_species") or [])}
        return [species for species in analyzed if species not in missing]

    @staticmethod
    def _mirna_channel_completed(mirna_summary: Mapping[str, Any] | None) -> bool:
        """Whether the miRNA seed scan actually ran, per its own aggregated summary.

        Channel-level, not per species: the scan is one batch over every submitted guide, so that is the
        honest granularity.

        ``total_candidates`` is the evidence, not ``analysis_files_processed`` -- the latter counts files
        the aggregator *globbed*, including a 0-byte one it then failed to read, the same weak signal
        that let an empty ``*_analysis.tsv`` publish a completed transcriptome screen. ``total_candidates``
        counts files that parsed and validated, so a header-only table (a scan that found nothing) counts
        and an unreadable one does not. An absent summary is unknown.
        """
        if not mirna_summary:
            return False
        for key in ("total_candidates", "analysis_files_processed"):
            value = mirna_summary.get(key)
            if isinstance(value, int):
                return value > 0
        return False

    def _evidence_species_scope(self) -> frozenset[str]:
        """Species a channel-level completion speaks for: what this run screened, plus the query species.

        The miRNA scan is one batch over every submitted guide, so it has no per-species outcome of
        its own; this is the scope its completion covers, and it is deliberately the same set the
        hand-rolled pair record used before evidence existed.
        """
        return frozenset(
            normalize_species_name(species) for species in (self._active_screen_species or self.config.screen_species)
        ) | {self._query_species}

    def _legacy_evidence_envelopes(
        self, digest: str, *, screened_species: Sequence[str], mirna_screened: bool
    ) -> tuple[EvidenceEnvelope, ...]:
        """Evidence inferred from an aggregate that carries no per-unit envelope of its own.

        A 0.7.0/0.7.1 result directory has no ``*_evidence.json`` anywhere, so the three-tier
        :meth:`_species_with_alignment_evidence` heuristic is the only positive evidence available.
        Those entries are marked :attr:`EvidenceSource.LEGACY_SUMMARY` rather than ``ENVELOPE``:
        nothing published them, so a caller that wants only first-hand records can drop them, while
        every existing result directory keeps the answer it has today.

        Counts stay unobserved. The heuristic establishes that an alignment ran, not what it counted,
        and writing a zero here would be exactly the fabricated clean screen #100 exists to stop.
        """
        entries = [
            ScreeningEvidenceEntry(
                channel=ScreeningChannel.TRANSCRIPTOME,
                species=species,
                guide_set_digest=digest,
                status=EvidenceStatus.COMPLETE,
                submitted_guide_digest=digest,
                submitted_guides=self._submitted_guide_count,
            )
            for species in dict.fromkeys(normalize_species_name(species) for species in screened_species)
        ]
        if mirna_screened:
            entries.extend(
                ScreeningEvidenceEntry(
                    channel=ScreeningChannel.MIRNA_SEED,
                    species=species,
                    guide_set_digest=digest,
                    status=EvidenceStatus.COMPLETE,
                    submitted_guide_digest=digest,
                    submitted_guides=self._submitted_guide_count,
                )
                for species in sorted(self._evidence_species_scope())
            )
        return tuple(
            EvidenceEnvelope(
                producer=EvidenceProducer.WORKFLOW_SYNTHESIS, source=EvidenceSource.LEGACY_SUMMARY, entry=entry
            )
            for entry in entries
        )

    def _not_requested_envelopes(
        self, digest: str, *, plan: ScreeningPlan, observed: Sequence[EvidenceEnvelope]
    ) -> tuple[EvidenceEnvelope, ...]:
        """``NOT_REQUESTED`` entries for declared pairs the plan deliberately holds no entry for.

        The plan, not a producer's silence, is what says whether a channel was asked for: a miRNA
        channel nobody requested is not the same as one that was requested and published nothing.
        Only pairs the policy takes a position on are synthesized, because an
        entry nobody declared answers no question, and only ``workflow_synthesis`` may emit this
        status -- a task that ran was requested by definition.
        """
        policy = getattr(self.config, "resolved_policy", None)
        if policy is None:
            return ()
        known = {(entry.channel.value, entry.species) for entry in plan.entries}
        known.update((envelope.entry.channel.value, envelope.entry.species) for envelope in observed)
        return tuple(
            EvidenceEnvelope(
                producer=EvidenceProducer.WORKFLOW_SYNTHESIS,
                source=EvidenceSource.SYNTHESIZED,
                entry=not_requested_entry(
                    requirement.channel,
                    requirement.species,
                    digest,
                    "declared by the run policy but never requested for screening, so nothing was searched",
                ),
            )
            for requirement in policy.evidence_requirements.channel_requirements
            if requirement.key not in known
        )

    def _reconcile_screening_evidence(
        self, output_dir: Path, *, screened_species: Sequence[str], mirna_screened: bool
    ) -> Reconciliation:
        """Reconcile the run's plan against what screening actually published.

        Three sources, in order of how directly they speak: the per-unit envelopes the tasks wrote,
        the reconciliation the pipeline's own aggregation published, and -- for a directory written
        before #100 existed -- the legacy summary heuristic. Whichever answers, every *planned* unit
        with nothing to show for it becomes a FAILED entry, so a species whose index build crashed
        and a pipeline that aborted before aggregation both leave a record instead of vanishing.
        """
        observed = collect_evidence(output_dir)
        if not observed:
            published = parse_reconciliation_payload(self._read_results_json(output_dir, RECONCILIATION_FILENAME))
            if published is not None:
                return published

        digest = self._guide_set_digest or (
            observed[0].entry.guide_set_digest if observed else _UNRECORDED_GUIDE_SET_DIGEST
        )
        plan = (
            self._screening_plan
            if self._screening_plan is not None
            # No plan recorded: a caller that reached result processing directly. The resolved set
            # restates one, rejections included, so the reconciliation is still keyed on the units
            # this run asked for rather than on whatever happened to publish.
            else self._screening_references.requested_plan(guide_set_digest=digest)
        )
        envelopes = list(observed) or list(
            self._legacy_evidence_envelopes(digest, screened_species=screened_species, mirna_screened=mirna_screened)
        )
        envelopes.extend(self._not_requested_envelopes(digest, plan=plan, observed=envelopes))
        return reconcile(plan, envelopes)

    def _completed_pairs_from_evidence(self, reconciliation: Reconciliation) -> frozenset[tuple[str, str]]:
        """The (channel, species) pairs reconciled evidence says completed: #100's integration seam.

        Spelled as ``ChannelRequirement.key`` is, so eligibility compares rather than re-derives, and
        derived from evidence rather than from a species list so a censored or failed unit cannot be
        read as clean. ``strict`` drops entries inferred from a pre-#100 summary, and applies only
        when this run's own tasks published envelopes -- where by construction there is nothing
        legacy to drop. A run is never downgraded by its own fallback, and no existing result
        directory changes the answer it already gives (#100).

        The whole reconciliation goes in, not its ``evidence``: the plan has to be in scope where the
        projection to (channel, species) drops the guide-set digest, or an envelope for a foreign
        guide set satisfies the requirement whose plan entry just failed over that mismatch.
        """
        policy = getattr(self.config, "resolved_policy", None)
        qualified = policy is not None and policy.run_mode is RunMode.QUALIFIED
        envelope_backed = any(source is EvidenceSource.ENVELOPE for source in reconciliation.sources.values())
        return completed_pairs(reconciliation, strict=qualified and envelope_backed)

    def _missing_required_evidence_units(self) -> tuple[tuple[str, str], ...]:
        """Required channel/species pairs this run reconciled no complete evidence for (#100).

        Asked of the reconciliation, not of the aggregate's own word: only a run that recorded a plan
        can know a unit was expected, and only ``EvidenceRequirements`` may say a unit was required.
        A run that reconciled nothing (a direct call into result processing, or an exit that returns
        before any producer runs) returns nothing here rather than claiming a shortfall it cannot
        substantiate -- that exit's ``reason`` already says what happened.

        Completeness comes from :meth:`_completed_pairs_from_evidence`, the same rule that fills
        :attr:`_completed_evidence_pairs` for eligibility, so the run status and the shortlist can
        never disagree about which unit was complete -- and a censored unit counts as missing in
        both, because a lower bound cannot show a ceiling was respected. Recomputed rather than read
        off that attribute because integration returns before assigning it on exactly the exit this
        question matters most for: the one where the pipeline published no results at all.
        """
        policy = getattr(self.config, "resolved_policy", None)
        if self._screening_evidence is None or policy is None:
            return ()
        completed = self._completed_pairs_from_evidence(self._screening_evidence)
        return tuple(sorted(policy.evidence_requirements.required_pairs - completed))

    def _publish_run_status(self, summary: dict[str, Any], *, reason: str | None = None) -> dict[str, Any]:
        """Stamp one screening summary with its :class:`RunStatus`, beside the keys it already has.

        Every exit out of screening goes through here, which is what makes the vocabulary one
        vocabulary rather than a seventh string. ``status``, ``reason`` and ``method`` are left
        exactly as they were: consumers compare ``status == "completed"`` (the CLI table, #103's
        report, existing result directories) and this slice publishes an additional key rather than
        redefining theirs.

        ``reason`` overrides the summary's own for an exit that carries no reason key but is still a
        named outcome -- the missing-output route reports ``partial`` because that is what its caller
        published, while what actually happened is that the pipeline wrote no output directory.
        """
        summary["run_status"] = screening_run_status(
            status=str(summary.get("status") or ""),
            reason=reason if reason is not None else str(summary.get("reason") or ""),
            required_evidence_missing=bool(self._missing_required_evidence_units()),
        ).value
        return summary

    @staticmethod
    def _transcriptome_shortfall_warnings(tx_summary: Mapping[str, Any]) -> list[str]:
        """Every way the aggregate says a requested species was not screened, as run warnings.

        A rejected file and an absent file are different facts with the same consequence, so both
        are reported with their reason.
        """
        missing = [str(species) for species in cast(list[Any], tx_summary.get("missing_species") or [])]
        unscreened = [str(species) for species in cast(list[Any], tx_summary.get("unscreened_species") or [])]
        rejected = cast(dict[str, Any], tx_summary.get("rejected_species_files") or {})

        warnings: list[str] = []
        if missing:
            warnings.append(
                "⚠️  No transcriptome alignment files were generated for: "
                f"{', '.join(missing)}. This usually means the BWA-MEM2 indexing stage ran out of memory. "
                "Increase Nextflow --max_memory (32GB+ recommended for human transcriptomes) or pre-build indices."
            )
        for species, reasons in sorted(rejected.items()):
            detail = "; ".join(str(reason) for reason in cast(list[Any], reasons))
            warnings.append(
                f"⚠️  Alignment output for '{species}' was rejected and contributed nothing to this screen: {detail}"
            )
        rejected_only = [species for species in unscreened if species not in missing]
        if rejected_only:
            warnings.append(
                f"⚠️  Not screened: {', '.join(rejected_only)}. Zero hits for these species means unknown, "
                "not clean; they keep no off-target evidence in this run."
            )
        return warnings

    def _log_offtarget_statistics(
        self,
        stats: Mapping[str, Any],
        aggregated_views: Mapping[str, Any],
        output_dir: Path,
    ) -> None:
        """Emit structured console logs for off-target statistics."""
        candidates_with_hits = stats.get("candidates_with_offtargets", 0)
        if candidates_with_hits:
            console.print(f"📊 Off-target analysis: {candidates_with_hits} candidates with hits")
            summaries = (
                ("failed_perfect_match", "❌ {} failed: perfect transcriptome matches"),
                ("failed_transcriptome_1mm", "❌ {} failed: 1mm transcriptome threshold"),
                ("failed_transcriptome_2mm", "❌ {} failed: 2mm transcriptome threshold"),
                ("failed_transcriptome_seed_perfect", "❌ {} failed: perfect transcriptome seed matches"),
                ("failed_mirna_seed", "❌ {} failed: miRNA perfect seed matches"),
                ("failed_high_risk_mirna", "❌ {} failed: high-risk miRNA hits"),
                ("failed_isoform_coverage", "❌ {} failed: protein-coding isoform coverage floor"),
            )
            for key, template in summaries:
                count = stats.get(key, 0)
                if count:
                    console.print(f"   {template.format(count)}")

            human_tx = stats.get("human_transcriptome_hits", 0)
            other_tx = stats.get("other_transcriptome_hits", 0)
            if human_tx or other_tx:
                console.print(f"   🧬 Transcriptome hits — human: {human_tx}, other: {other_tx}")

            human_mirna = stats.get("human_mirna_hits", 0)
            other_mirna = stats.get("other_mirna_hits", 0)
            if human_mirna or other_mirna:
                console.print(f"   🌱 miRNA hits — human: {human_mirna}, other: {other_mirna}")

        self._log_unscored_after_screening(stats)

        missing_on_target = stats.get("missing_on_target_hit", 0)
        self._log_missing_on_target(missing_on_target)

        if aggregated_views:
            tx_summary = aggregated_views.get("transcriptome")
            if tx_summary:
                self._log_aggregated_transcriptome(tx_summary)

            mirna_summary = aggregated_views.get("mirna")
            if mirna_summary:
                human_hits = mirna_summary.get("human_hits", 0)
                other_hits = mirna_summary.get("other_species_hits", 0)
                console.print(f"   🌱 Aggregated miRNA hits — human: {human_hits}, other: {other_hits}")

        trace_file = Path(output_dir) / "pipeline_info" / "execution_trace.txt"
        if trace_file.exists():
            console.print(f"   📘 Nextflow execution trace: {trace_file}")

    def _log_aggregated_transcriptome(self, tx_summary: Mapping[str, Any]) -> None:
        """Report the aggregate's per-species view, separating "clean" from "not screened"."""
        species_counts = cast(dict[str, int], tx_summary.get("hits_per_species", {}) or {})
        human_hits = tx_summary.get("human_hits", 0)
        other_hits = tx_summary.get("other_species_hits", 0)
        # These roll-ups are derived from the screened species' counts, so on a run that screened
        # nothing they are 0 for want of evidence. Printing that as the first line of the block was
        # the same "absent reads as clean" the per-species reporting below exists to prevent.
        if self._species_with_alignment_evidence(tx_summary):
            console.print(f"   🧾 Aggregated transcriptome hits — human: {human_hits}, other: {other_hits}")
        else:
            console.print("   🧾 Aggregated transcriptome hits — none screened, so counts are unknown (not zero)")
        if species_counts:
            formatted = ", ".join(f"{k}: {v}" for k, v in sorted(species_counts.items()))
            console.print(f"      per species: {formatted}")

        missing_species = [str(species) for species in cast(list[Any], tx_summary.get("missing_species") or [])]
        if missing_species:
            console.print(
                "      ⚠️ Transcriptome alignment files were missing for: "
                f"{', '.join(missing_species)} (likely insufficient memory during BWA indexing)."
            )

        # "No hits detected" is a clean-screen claim, so it is made only about species that were
        # screened. Reported against species_analyzed it read a rejected alignment as good news.
        screened_species = self._species_with_alignment_evidence(tx_summary)
        zero_hit_species = [species for species in screened_species if species_counts.get(species, 0) == 0]
        if zero_hit_species:
            console.print(f"      ℹ️ No transcriptome hits detected for: {', '.join(zero_hit_species)}")

        unscreened = [
            str(species)
            for species in cast(list[Any], tx_summary.get("unscreened_species") or [])
            if str(species) not in missing_species
        ]
        if unscreened:
            console.print(
                f"      ⚠️ No usable alignment evidence for: {', '.join(unscreened)} (counts unknown, not zero)."
            )

    @staticmethod
    def _log_unscored_after_screening(stats: Mapping[str, Any]) -> None:
        """Warn when screening ran but some candidates still carry their design-time score."""
        not_scored = stats.get("candidates_not_scored_after_screening", 0)
        if not_scored:
            console.print(
                f"   ⚠️ {not_scored} of {stats.get('candidates_analyzed', 0)} candidates could not be scored after "
                "screening and kept their design-time score (see scored_after_screening in the candidate CSV)"
            )

    @staticmethod
    def _log_missing_on_target(missing_on_target: int) -> None:
        """Warn when no candidates confirmed a 0-mismatch hit against their own source transcript."""
        if missing_on_target:
            console.print(
                f"   ⚠️ {missing_on_target} candidate(s) had no confirmed 0-mismatch hit against their own "
                "source transcript (on-target self-match not found in the transcriptome index)"
            )

    async def _basic_offtarget_analysis(self, candidates: list[SiRNACandidate]) -> dict[str, Any]:
        """Fallback basic off-target analysis: sequence-only, and no evidence for any planned unit.

        This path aligns nothing against a reference, so every unit the run planned FAILED, and the
        published evidence says so rather than leaving the plan unanswered. Without that record the
        fallback was indistinguishable from a screen that ran and found nothing (#100).
        """
        # Use simplified analysis when external tools are not available
        analyzer = OffTargetAnalysisManager(species="human")  # Default to human for basic analysis
        results = {}

        for candidate in candidates:
            analysis_result = analyzer.analyze_sirna_candidate(candidate)

            # Extract relevant metrics for backward compatibility
            mirna_hits = analysis_result.get("mirna_hits", [])
            transcriptome_hits = analysis_result.get("transcriptome_hits", [])

            # Calculate basic scores
            off_target_count = len(mirna_hits) + len(transcriptome_hits)
            penalty = off_target_count * 10  # Simple penalty calculation
            score = math.exp(-penalty / 50)  # Score calculation

            results[candidate.id] = {
                "off_target_count": off_target_count,
                "off_target_penalty": penalty,
                "off_target_score": score,
                "method": "sequence_analysis",
            }

        # Save results
        results_file = self.config.output_dir / "off_target" / "basic_analysis.json"
        with results_file.open("w") as f:
            json.dump(results, f, indent=2)

        # Every planned unit failed here, stated from the plan rather than left unanswered: this
        # method is also the landing place after a Nextflow failure, so "nothing was published" is
        # exactly what happened to each one.
        self._screening_evidence = reconcile(
            self._screening_plan or ScreeningPlan(),
            (),
            missing_detail=(
                "no reference alignment ran: the basic sequence-only fallback replaced the screen, "
                "so this unit produced no screening evidence"
            ),
        )
        self._completed_evidence_pairs = self._completed_pairs_from_evidence(self._screening_evidence)
        # Nothing was aligned, so every gate reading a screening channel is undecidable here rather
        # than unrun (#106): this path reaches _integrate_offtarget_results on no route at all.
        self._gate_without_screening_evidence(candidates)

        # `partial`, never `completed`: this path aligns nothing and never reaches
        # _integrate_offtarget_results, so no candidate gets a verdict at all. It is also reached after a
        # Nextflow failure, where `completed` printed as "Off-target Analysis: Complete".
        console.print(
            f"📊 Basic sequence-only off-target analysis for {len(candidates)} candidates "
            "(no reference alignment: this is not a completed screen)"
        )
        return self._publish_run_status(
            {
                "status": "partial",
                "method": "basic",
                "results": results,
                "aggregated": {},
                "screening_evidence": reconciliation_payload(self._screening_evidence),
            }
        )

    async def _parse_nextflow_results(self, output_dir: Path) -> dict[str, Any]:  # noqa: PLR0912
        """Parse results from Nextflow off-target analysis.

        Parses BOTH transcriptome AND miRNA results from their respective
        output directories and combines them into a single results structure for
        candidate filtering.
        """
        results: dict[str, dict[str, Any]] = {}

        if not output_dir.exists():
            return {"status": "missing", "method": "nextflow", "output_dir": str(output_dir), "results": results}

        # Check for combined transcriptome results in aggregated subdirectory
        aggregated_dir = output_dir / "aggregated"

        def _aggregate_path(filename: str) -> Path:
            return (aggregated_dir / filename) if aggregated_dir.exists() else (output_dir / filename)

        def _ingest_row(row: dict[str, Any]) -> None:
            qname = row.get("qname") or row.get("query") or row.get("id")
            if not qname:
                return
            try:
                score = float(row.get("offtarget_score") or row.get("score") or 0)
            except Exception:
                score = 0.0
            entry = results.setdefault(qname, {"off_target_count": 0, "off_target_score": 0.0, "hits": []})
            entry["off_target_count"] += 1
            entry["off_target_score"] = max(entry["off_target_score"], score)
            entry["hits"].append(row)

        # One table per transcriptome file read, each keeping its own ordered row list, so the
        # classification columns are written back onto exactly the rows that were read from that
        # file. Every transcriptome row reaching a candidate counter is in one of these tables.
        transcriptome_tables: list[dict[str, Any]] = []
        # Recorded, not just logged: the run summary names these files rather than copying their rows,
        # and a path only reachable from a log line is not a pointer a consumer can follow.
        mirna_hit_files: list[Path] = []

        def _record_mirna_file(path: Path) -> None:
            """Publish one ingested miRNA file as a detail pointer, at most once.

            Every path whose rows reached the miRNA counters, the aggregate included (#108). The summary
            points at these files rather than carrying their rows, so this list is the only route to
            them.

            Deduplicated by path: the aggregate and a per-file fallback can name the same file, and
            naming one twice reads as two sources of evidence.
            """
            if path not in mirna_hit_files:
                mirna_hit_files.append(path)

        def _ingest_tsv(path: Path, table: dict[str, Any] | None = None) -> bool:
            if not path.exists() or path.stat().st_size == 0:
                return False
            found = False
            with path.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                if table is not None and reader.fieldnames:
                    table["path"] = path
                    table["fieldnames"] = [name for name in reader.fieldnames if name]
                for row in reader:
                    _ingest_row(row)
                    if table is not None:
                        cast(list[dict[str, Any]], table["rows"]).append(row)
                    found = True
            return found

        def _ingest_transcriptome_tsv(path: Path) -> bool:
            """Ingest a transcriptome hit file into its own persistable table."""
            table: dict[str, Any] = {"path": None, "fieldnames": [], "rows": []}
            found = _ingest_tsv(path, table)
            if table["path"] is not None:
                transcriptome_tables.append(table)
            return found

        def _ingest_json(path: Path, table: dict[str, Any] | None = None) -> bool:
            if not path.exists() or path.stat().st_size == 0:
                return False
            raw_data: list[Any] | dict[str, Any] | str | int | float | bool | None
            try:
                with path.open() as fh:
                    raw_data = json.load(fh)
            except Exception:
                raw_data = []
            found = False
            data: list[dict[str, Any]] = []
            if isinstance(raw_data, list):
                raw_entries: list[Any] = raw_data
            else:
                raw_entries = []
            for entry in raw_entries:
                if isinstance(entry, dict):
                    data.append(cast(dict[str, Any], entry))
            for item in data:
                _ingest_row(item)
                if table is not None:
                    if not table["fieldnames"]:
                        table["fieldnames"] = [str(key) for key in item if key]
                    cast(list[dict[str, Any]], table["rows"]).append(item)
                found = True
            return found

        def _ingest_transcriptome_json(path: Path, tsv_path: Path) -> bool:
            """Ingest the JSON aggregate into a table that will be published as ``tsv_path``.

            The JSON aggregate is the one path whose rows reached the candidate counters with no
            table to be republished into, so the published hit table under-reported liabilities the
            candidates had already been charged for. Giving those rows the TSV the producer would
            have written puts them back under the same guarantee as every other path.
            """
            table: dict[str, Any] = {"path": None, "fieldnames": [], "rows": []}
            found = _ingest_json(path, table)
            if found and table["fieldnames"]:
                table["path"] = tsv_path
                transcriptome_tables.append(table)
            return found

        transcriptome_hits_found = _ingest_transcriptome_tsv(_aggregate_path("combined_offtargets.tsv"))
        if not transcriptome_hits_found:
            transcriptome_hits_found = _ingest_transcriptome_json(
                _aggregate_path("combined_offtargets.json"), _aggregate_path("combined_offtargets.tsv")
            )

        mirna_aggregate_tsv = _aggregate_path("combined_mirna_hits.tsv")
        mirna_hits_found = _ingest_tsv(mirna_aggregate_tsv)
        if mirna_hits_found:
            _record_mirna_file(mirna_aggregate_tsv)
        else:
            mirna_aggregate_json = _aggregate_path("combined_mirna_hits.json")
            mirna_hits_found = _ingest_json(mirna_aggregate_json)
            if mirna_hits_found:
                _record_mirna_file(mirna_aggregate_json)

        if not transcriptome_hits_found or not mirna_hits_found:
            transcriptome_files: list[Path] = []
            mirna_files: list[Path] = []

            if not transcriptome_hits_found:
                transcriptome_dir = output_dir / "transcriptome"
                if transcriptome_dir.exists():
                    transcriptome_files = list(transcriptome_dir.glob("*_analysis.tsv"))

            if not mirna_hits_found:
                mirna_dir = output_dir / "mirna"
                if mirna_dir.exists():
                    mirna_files = list(mirna_dir.glob("*_analysis.tsv"))

            files: list[Path] = []
            if not transcriptome_hits_found:
                files.extend(transcriptome_files)
            if not mirna_hits_found:
                files.extend(mirna_files)

            if not files and not transcriptome_hits_found and not mirna_hits_found:
                # Last resort: scan for any TSV files. Treated as transcriptome hits, which is what the
                # *_offtargets.tsv name means, so they are persisted rather than counted and lost.
                transcriptome_files = list(output_dir.glob("**/*_offtargets.tsv"))
                files = list(transcriptome_files)

            transcriptome_file_set = {Path(path) for path in transcriptome_files}
            for fpath in files:
                path = Path(fpath)
                # Genome rows carry a class and must reach the published table; miRNA rows have no
                # class of their own and feed the miRNA counters only.
                if path in transcriptome_file_set:
                    _ingest_transcriptome_tsv(path)
                elif _ingest_tsv(path):
                    # Recording the ingest here is what stops it happening twice. The glob above
                    # already matches mirna/mirna_analysis.tsv, so a trailing "if not
                    # mirna_hits_found" retry of that exact path re-read the same file and doubled
                    # every miRNA counter; there is no layout in which the retry reached a file the
                    # glob did not.
                    logger.info(f"Parsed miRNA analysis results from {path}")
                    _record_mirna_file(path)
                    mirna_hits_found = True

        return {
            "status": "completed",
            "method": "nextflow",
            "output_dir": str(output_dir),
            "results": results,
            "transcriptome_hit_tables": transcriptome_tables,
            "mirna_hit_files": mirna_hit_files,
        }

    def _gate_offtarget_counts(
        self,
        candidate: SiRNACandidate,
        *,
        counts: OffTargetGateCounts,
        filter_criteria: OffTargetFilterCriteria,
        complete_pairs: frozenset[tuple[str, str]],
        stats: dict[str, Any],
    ) -> None:
        """Apply every off-target gate to one candidate's counts, and count what rejected it.

        One entry point for both the no-hit and with-hit paths, so a gate cannot be wired into only one
        of them -- how the clean-screen case ended up ungated (#106).
        :meth:`_gate_without_screening_evidence` is the third caller: the paths that publish no
        evidence at all come through here too, with ``complete_pairs`` empty.

        ``passes_filters`` is set by the rejecting gate through ``record_filter_verdict``, which keeps
        the first failure. Not assigned here: overwriting it unconditionally let an off-target label mask
        a design verdict already on the row.
        """
        should_fail, fail_status = self._check_offtarget_filters(
            *counts,
            filter_criteria,
            candidate,
            complete_pairs=complete_pairs,
        )
        if not should_fail or fail_status is None:
            return
        logger.info(f"Candidate {candidate.id} failed off-target filter: {fail_status.value}")
        stat_key = _OFFTARGET_REJECTION_STATS.get(fail_status)
        if stat_key is not None:
            stats[stat_key] += 1

    def _requested_species_scope(self) -> frozenset[str]:
        """The species every post-screen gate's evidence is scoped over: what was asked for, plus the query.

        ``_active_screen_species`` rather than ``config.screen_species`` because extra species can
        arrive on ``transcriptome_indices``/``transcriptome_fastas``, and the query species is always in
        scope even when it was never listed. Shared by the integration path and by
        :meth:`_gate_without_screening_evidence` so a gate's evidence requirement is the same set
        whether or not a screen ran -- a scope that shrank when the screen failed would leave the gate
        with no pair to be incomplete about.
        """
        requested = frozenset(
            normalize_species_name(species) for species in (self._active_screen_species or self.config.screen_species)
        )
        return requested | {self._query_species}

    def _gate_without_screening_evidence(
        self,
        candidates: Sequence[SiRNACandidate],
        *,
        filter_criteria: OffTargetFilterCriteria | None = None,
    ) -> None:
        """Reach every off-target gate on a path that published no screening evidence at all (#106).

        Four paths -- the basic sequence-only fallback, ``nextflow_unavailable``, ``nextflow_failed``
        and an output directory that never appeared -- return without reaching
        ``_integrate_offtarget_results``, so no candidate reaches ``_gate_offtarget_counts`` and every
        channel-reading gate would export ``not_evaluated``, the same cell a run with no threshold
        configured writes. A screen that never happened and a screen that came back clean must not
        export the same nine cells.

        ``complete_pairs=frozenset()``: nothing completed, so no gate may report a pass, and each one
        records ``UNKNOWN`` with an empty observed value rather than a fabricated zero. An undecidable
        gate never rejects, so ``passes_filters`` is untouched and no rejection counter can move --
        which is why no statistics are published from here. What it does do is withhold the candidate
        from a QUALIFIED shortlist under the named reason ``unknown:<filter_id>``.

        Only candidates that reached no off-target gate: the call that held the counts owns the
        answer, and ``record_filter_verdict`` overwrites, so re-gating one would replace a decided
        PASS/FAIL with ``UNKNOWN``. That matters on the ``nextflow_failed`` path, which can be reached
        by an exception raised after integration already ran.

        Deliberately not called from the ``user_disabled`` or ``no_candidates`` exits: a channel the
        user switched off was never requested, and ``not_evaluated`` is the honest cell for it.
        """
        criteria = filter_criteria or getattr(self.config.design_params, "offtarget_filters", None)
        criteria = criteria or OffTargetFilterCriteria()
        # Set here as well as in _integrate_offtarget_results, which never ran on these paths.
        self._screened_species_scope = self._screened_species_scope or self._requested_species_scope()
        # Nothing here can reject, so this dict only exists to satisfy the shared entry point.
        discarded_stats: dict[str, Any] = dict.fromkeys(_OFFTARGET_REJECTION_STATS.values(), 0)
        gated = 0
        for candidate in candidates:
            verdicts = candidate.filter_verdicts or {}
            if any(filter_id in verdicts for filter_id in POST_SCREEN_FILTER_CHANNELS):
                continue
            self._gate_offtarget_counts(
                candidate,
                counts=_ZERO_OFFTARGET_COUNTS,
                filter_criteria=criteria,
                complete_pairs=frozenset(),
                stats=discarded_stats,
            )
            gated += 1
        if gated:
            logger.warning(
                f"No screening evidence was published for this run: the off-target gates on {gated} candidate(s) "
                "report unknown rather than a clean screen, and their observed values stay empty."
            )

    def _gate_evidence_pairs(self, filter_id: str, channels: frozenset[ScreeningChannel]) -> frozenset[tuple[str, str]]:
        """The channel x species pairs one gate needs before it may report a pass.

        The species come from the gate's own declared ``FilterScope``, so the evidence a gate requires
        is the evidence it counts. An unrestricted scope means every species this run screened.
        """
        policy = getattr(self.config, "resolved_policy", None)
        scoped: frozenset[str] = frozenset()
        if policy is not None:
            try:
                scoped = frozenset(policy.descriptor(filter_id).scope.species)
            except (KeyError, RunPolicyError):
                scoped = frozenset()
        species = scoped or self._screened_species_scope
        return frozenset((channel.value, name) for channel in channels for name in species)

    def _check_offtarget_filters(
        self,
        transcriptome_0mm: int,
        transcriptome_1mm: int,
        transcriptome_2mm: int,
        transcriptome_seed_0mm: int,
        mirna_0mm_seed: int,
        mirna_high_risk: int,
        total_hits: int,
        genuine_off_target_count: int,
        filter_criteria: OffTargetFilterCriteria,
        candidate: SiRNACandidate,
        *,
        complete_pairs: frozenset[tuple[str, str]] | None = None,
    ) -> tuple[bool, SiRNACandidate.FilterStatus | None]:
        """Record every off-target gate's verdict, and report the first that rejects.

        Each gate writes its own outcome and the value it compared onto the candidate. That matters
        more here than at the design stage: six of these gates read human-stratified counters that
        exist only as locals in the caller, and the identically named exported columns are all-species
        totals and disagree.

        The rules themselves live in :mod:`sirnaforge.core.filtering` (#100), which decides every gate
        from (threshold, action, observed, evidence completeness); this method's job is to say what
        each gate compares, which evidence its pass claim rests on, and which label it rejects under.

        ``complete_pairs`` is which channel x species pairs this candidate holds evidence for; ``None``
        means the caller has no per-species record and every gate is treated as evidenced. A gate whose channels did
        not all complete cannot report a pass -- its count is a lower bound, so a zero means "nothing was
        looked for" -- and records ``UNKNOWN`` with no observed value. Not the lower bound: #103's report
        re-derives verdicts from that column, and a 0 there would produce a confident pass the run never
        made. The count survives on the candidate's own hit counters.

        A gate whose incomplete count already exceeds its ceiling still FAILs: a lower bound above the
        threshold is proof enough. Only the pass direction needs completeness.

        ``transcriptome_seed_0mm`` is the only input that sees a *partial* hit whose seed paired
        perfectly. ``nm`` is a guide-level distance, so a clipped or gapped hit carries nm > 2 and
        lands in none of the ``transcriptome_{0,1,2}mm`` strata even when its seed is intact; the
        seed counter and ``genuine_off_target_count`` are the only signals left. Its threshold
        (``max_transcriptome_seed_perfect``) defaults to ``None``, so this check is inert until a
        user opts in. Unlike the three mismatch counts it is not species-split -- it is the
        reported ``transcriptome_hits_seed_0mm`` column verbatim, so the gate fires on exactly the
        number the user sees.

        Returns:
            Tuple of (should_fail, fail_status enum or None)
        """
        # (filter_id, threshold, observed, label). The filter_id is what ties each check to the
        # declared filter whose action decides whether exceeding the threshold rejects the candidate
        # or is merely recorded, the key the verdict is published under, and the lookup into
        # POST_SCREEN_FILTER_CHANNELS for the channels its count depends on.
        checks: list[tuple[str, int | None, int, SiRNACandidate.FilterStatus]] = [
            (
                "max_transcriptome_hits_0mm",
                filter_criteria.max_transcriptome_hits_0mm,
                transcriptome_0mm,
                SiRNACandidate.FilterStatus.TRANSCRIPTOME_PERFECT_MATCH,
            ),
            (
                "max_transcriptome_hits_1mm",
                filter_criteria.max_transcriptome_hits_1mm,
                transcriptome_1mm,
                SiRNACandidate.FilterStatus.TRANSCRIPTOME_1MM,
            ),
            (
                "max_transcriptome_hits_2mm",
                filter_criteria.max_transcriptome_hits_2mm,
                transcriptome_2mm,
                SiRNACandidate.FilterStatus.TRANSCRIPTOME_2MM,
            ),
            (
                "max_transcriptome_seed_perfect",
                filter_criteria.max_transcriptome_seed_perfect,
                transcriptome_seed_0mm,
                SiRNACandidate.FilterStatus.TRANSCRIPTOME_SEED_PERFECT,
            ),
            (
                "max_mirna_perfect_seed",
                filter_criteria.max_mirna_perfect_seed,
                mirna_0mm_seed,
                SiRNACandidate.FilterStatus.MIRNA_PERFECT_SEED,
            ),
            (
                "max_total_offtarget_hits",
                filter_criteria.max_total_offtarget_hits,
                total_hits,
                SiRNACandidate.FilterStatus.TOTAL_OFFTARGETS,
            ),
            (
                "max_off_target_count",
                filter_criteria.max_off_target_count,
                genuine_off_target_count,
                SiRNACandidate.FilterStatus.EXCESS_OFF_TARGETS,
            ),
            # The boolean flag, as a ceiling of zero, so it goes through the same path as every other
            # gate instead of a trailing special case. A high-risk hit is by definition a perfect seed
            # hit, so it must not sit behind max_mirna_perfect_seed's own ceiling of 0.
            (
                "fail_on_high_risk_mirna",
                0 if filter_criteria.fail_on_high_risk_mirna else None,
                mirna_high_risk,
                SiRNACandidate.FilterStatus.HIGH_RISK_MIRNA,
            ),
        ]

        # Every check is evaluated and recorded, by the one shared evaluator. Returning on the first
        # rejection left the remaining gates unmeasured, so their reported counts were a function of
        # this list's order; deriving the rules here left five call sites to keep in step.
        specs: list[GateSpec] = []
        observed: dict[str, float | None] = {}
        statuses: dict[str, SiRNACandidate.FilterStatus] = {}
        for filter_id, threshold, value, status in checks:
            channels = POST_SCREEN_FILTER_CHANNELS[filter_id]
            specs.append(
                GateSpec(
                    filter_id=filter_id,
                    threshold=threshold,
                    action=self._filter_action_for(filter_id, threshold),
                    # Every gate here is a ceiling, including the boolean flag as a ceiling of zero.
                    comparator=Comparator.AT_MOST,
                    channels=channels,
                    evidence_pairs=self._gate_evidence_pairs(filter_id, channels),
                )
            )
            observed[filter_id] = value
            statuses[filter_id] = status

        outcomes = evaluate_gates(specs, observed, complete_pairs=complete_pairs)
        for outcome in outcomes:
            record_gate_outcome(candidate, outcome, status=statuses[outcome.filter_id])

        rejection = first_rejection(outcomes)
        return rejection is not None, statuses[rejection.filter_id] if rejection is not None else None

    def _filter_action_for(self, filter_id: str, threshold: float | None) -> FilterAction:
        """The action in force for one post-screen filter: the run's policy, else the declared default.

        Not off-target-specific despite where it sits: ``min_isoform_coverage`` resolves through it
        too, which is what makes that gate configurable rather than hard-coded to reject.

        ``threshold`` decides the OFF case. A registry default of ``off`` means "no threshold is
        configured", not "never gate on this" -- ``max_transcriptome_seed_perfect`` ships that way
        precisely so a caller can opt in by setting one. Treating the default as authoritative would
        silently ignore a threshold the caller asked for, which is the defect that gate already had
        once. ``warn`` is the only action chosen independently of whether a threshold exists.
        """
        policy = getattr(self.config, "resolved_policy", None)
        if policy is not None:
            try:
                descriptor = policy.descriptor(filter_id)
            except (KeyError, RunPolicyError):
                descriptor = None
            if descriptor is not None:
                # OFF because the policy has no threshold is not the same as OFF because a caller
                # asked for it. Only the first is overridden by a threshold arriving on the criteria.
                if descriptor.action is FilterAction.OFF and descriptor.threshold is None and threshold is not None:
                    return FilterAction.FAIL
                return cast(FilterAction, descriptor.action)
        spec = FILTER_SPEC_BY_ID.get(filter_id)
        if spec is None:
            return FilterAction.FAIL
        if spec.default_action is FilterAction.OFF and threshold is not None:
            return FilterAction.FAIL
        return spec.default_action

    @staticmethod
    def _species_present_on_hits(offtarget_data: dict[str, Any]) -> frozenset[str]:
        """Canonical species that actually appear on a transcriptome hit row.

        A species that was screened but aligned nothing needs no orthologue lookup, because there is
        no cross-species hit to classify. Keeps a single-species screen (and any unit test) off the
        network entirely.
        """
        seen: set[str] = set()
        for entry in (offtarget_data.get("results") or {}).values():
            for hit in (entry or {}).get("hits") or []:
                if "mirna_id" in hit or "database" in hit:
                    continue
                label = hit.get("species")
                if label and str(label).strip():
                    seen.add(normalize_species_name(str(label)))
        return frozenset(seen)

    async def _resolve_ortholog_mapping(
        self, screened_species: Sequence[str] | None, offtarget_data: dict[str, Any]
    ) -> OrthologueMapping:
        """Resolve orthologue gene IDs for non-query species that actually produced hits.

        Skipped entirely for a single-species screen, which is the common case and needs no network
        call. A configured ``ortholog_mapping_file`` replaces the Compara call outright, which is the
        offline path #101 requires and the only path any test may take. A failure is not fatal: the
        mapping comes back with the species in ``unresolved_species`` and its hits fall back to the
        labelled symbol heuristic.
        """
        query_species = self._query_species
        candidates_for_lookup = (
            frozenset(normalize_species_name(s) for s in (screened_species or self._active_screen_species or ()))
            & self._species_present_on_hits(offtarget_data)
        ) - {query_species}
        if not candidates_for_lookup or not (self._query_gene_ids or self._query_gene_symbols):
            return OrthologueMapping.empty()

        if self._ortholog_table is not None:
            logger.info(
                f"Reading orthologue evidence for {sorted(candidates_for_lookup)} from "
                f"{self.config.ortholog_mapping_file} instead of Ensembl Compara"
            )
            return mapping_from_table(
                self._ortholog_table,
                frozenset(self._query_gene_ids),
                query_species,
                candidates_for_lookup,
                query_gene_symbols=frozenset(self._query_gene_symbols),
            )

        mapping = await resolve_orthologues(
            frozenset(self._query_gene_ids),
            query_species,
            candidates_for_lookup,
            # An input FASTA supplies transcript IDs, not gene IDs, and Compara answers those with an
            # empty 200. The symbol is the fallback identifier; the answer is still a stable gene ID.
            query_gene_symbols=frozenset(self._query_gene_symbols),
            # No wall-clock ceiling on a real run. Abandoning a slow-but-working Compara would report
            # a resolvable species as unresolved, and unresolved conservation still publishes 0.0
            # rather than null (#101) -- so a budget here buys speed by fabricating a measurement.
            budget=None,
        )
        if mapping.unresolved_species:
            logger.warning(
                f"Orthologue mapping unresolved for {sorted(mapping.unresolved_species)}; hits in those species "
                "fall back to gene-symbol equality, which is a heuristic and is labelled as such in the hit table."
            )
        return mapping

    def _integrate_offtarget_results(  # noqa: PLR0912, C901
        self,
        candidates: list[SiRNACandidate],
        offtarget_data: dict[str, Any],
        filter_criteria: OffTargetFilterCriteria | None = None,
        screened_species: Sequence[str] | None = None,
        ortholog_gene_ids: frozenset[str] = frozenset(),
        mirna_screened: bool | None = None,
        unresolved_orthology_species: frozenset[str] = frozenset(),
    ) -> tuple[list[SiRNACandidate], dict[str, Any]]:
        """Integrate off-target analysis results, classify hits, and score candidates.

        Workflow reordering: screening now happens before final scoring. This method:
        1. Fans out deduplicated screening results to all candidates sharing a sequence
        2. Classifies each transcriptome hit four ways (on-target, ortholog, repeat, off-target)
        3. Computes post-screen sub-scores (off-target, isoform coverage, conservation)
        4. Computes the final composite score with the full term set

        Re-ranking candidates by the new scores is the caller's job (step5_offtarget_analysis),
        which holds the DesignResult whose candidates/top_candidates need reordering.

        Args:
            candidates: List of siRNA candidates to update
            offtarget_data: Off-target results from Nextflow pipeline
            filter_criteria: Optional filtering criteria for off-targets
            screened_species: Species with positive evidence of a published alignment. A run that
                produced no alignment for the query species cannot support an off-target term, so
                its candidates keep their design-time scores instead of being awarded perfect
                specificity. None means the caller has no per-species evidence to offer (direct
                callers, the basic analysis fallback) and the requested set is assumed screened.
            ortholog_gene_ids: Version-stripped gene IDs resolved as orthologues of the query gene,
                from :func:`sirnaforge.data.orthology.resolve_orthologues`. Empty is the honest
                default: orthology then falls back to the labelled gene-symbol heuristic, which
                cannot resolve pairs like human TP53 / mouse Trp53.
            unresolved_orthology_species: Species whose orthologue lookup did not complete, from
                :attr:`OrthologueMapping.unresolved_species`. Conservation is null for a candidate
                whose denominator includes one, rather than reporting a fraction over a numerator that
                could not be populated. Empty is the honest default: a caller with no orthology outcome
                to report is not asserting that every lookup succeeded, but it has nothing better, and
                the alternative -- assuming failure -- would null out conservation on every direct call.
            mirna_screened: Whether the miRNA seed channel completed. The miRNA scan is one batch
                over every submitted guide, so completion is a property of the run and not of a
                candidate. ``False`` makes the miRNA gates ``UNKNOWN`` rather than passing on a
                fabricated zero. ``None`` means the caller has no channel evidence to offer -- the
                same convention as ``screened_species`` -- and the channel is assumed to have run.

        Returns:
            Tuple of (updated candidates, statistics dict with hit class decomposition)
        """
        if not offtarget_data or offtarget_data.get("status") != "completed":
            logger.warning("No completed off-target data available; candidates keep design-time scores")
            # Keeping the design-time score is not the whole answer: the gates still have to say they
            # could not decide, or this exit publishes `not_evaluated` for every one of them (#106).
            self._gate_without_screening_evidence(candidates, filter_criteria=filter_criteria)
            return candidates, {}

        if filter_criteria is None:
            filter_criteria = OffTargetFilterCriteria()

        results = offtarget_data.get("results", {})

        # Single authoritative query species, set once in __init__ (see comment there).
        query_species = self._query_species

        # Conservation is keyed on the set handed to the aligner, not on whether species were typed
        # on the CLI: the term goes inactive exactly when that set is query-species-only.
        # _active_screen_species, not config.screen_species, is what Nextflow was handed: extra
        # species can arrive on transcriptome_indices/transcriptome_fastas and
        # can return ortholog hits, so scoring them against the shorter config list made the
        # conservation numerator exceed its denominator.
        requested_species = frozenset(
            normalize_species_name(s) for s in (self._active_screen_species or self.config.screen_species)
        )
        # A species whose alignment never ran STAYS in this denominator: subtracting it would let a
        # degraded run out-report the complete run it degraded from. A species that was screened and
        # produced nothing can only lower conservation, never raise it. Species with no resolvable
        # index never enter this set at all, so conservation is a statement about what was compared,
        # not about the CLI species list.
        conservation_denominator = requested_species - {query_species}

        if screened_species is None:
            # No per-species evidence offered; assume what was requested was screened, which is
            # what every caller outside the Nextflow path can honestly claim.
            screened = requested_species | {query_species}
        else:
            screened = frozenset(normalize_species_name(s) for s in screened_species)
        unscreened = requested_species - screened
        query_species_unscreened = query_species not in screened
        if unscreened - {query_species}:
            logger.warning(
                f"No alignment evidence for {sorted(unscreened - {query_species})}; those species stay in the "
                "conservation denominator, so conservation is a lower bound for every candidate in this run."
            )
        if query_species_unscreened:
            logger.error(
                f"No alignment evidence for query species '{query_species}' (its alignment produced no files, or "
                "aggregation reported nothing at all); off-target counts are unknown. Refusing to compute "
                "post-screen scores: candidates keep their design-time scores."
            )
            console.print(
                f"❌ No {query_species} alignment evidence: candidates cannot be scored after screening and are "
                "reported with design-time scores. Re-run the screen before trusting the ranking."
            )

        # Channel completion, not hit presence, is what lets a zero be read as clean. The miRNA scan is
        # one batch, so its completion is run-level; the transcriptome channel is per candidate and is
        # decided inside the loop.
        unresolved_orthology = frozenset(normalize_species_name(s) for s in unresolved_orthology_species)
        mirna_channel_complete = mirna_screened is not False
        if not mirna_channel_complete:
            logger.warning(
                "No miRNA seed evidence in this run: the miRNA gates report unknown rather than passing on a "
                "count of zero, and their observed values stay empty."
            )
        # Spelled as the declared requirements are, so eligibility compares rather than re-derives.
        # Reconciled evidence is the authority whenever screening produced any (#100): a unit that
        # failed, was censored or was never requested cannot contribute a pair, however many species
        # a summary lists. The species-list derivation below remains for a caller that integrates
        # hits with no evidence at all -- a direct call, or a path that never reached a producer.
        self._completed_evidence_pairs = (
            self._completed_pairs_from_evidence(self._screening_evidence)
            if self._screening_evidence is not None
            else frozenset((ScreeningChannel.TRANSCRIPTOME.value, species) for species in screened)
            | (
                frozenset(
                    (ScreeningChannel.MIRNA_SEED.value, species) for species in requested_species | {query_species}
                )
                if mirna_channel_complete
                else frozenset()
            )
        )

        # The same pairs the gates read, so eligibility and gate evaluation cannot disagree about what
        # completed. `unscreened_pairs` is what a candidate the query-species alignment missed loses.
        self._screened_species_scope = self._requested_species_scope()
        run_complete_pairs = self._completed_evidence_pairs or frozenset()
        unscreened_pairs = frozenset(
            {(ScreeningChannel.TRANSCRIPTOME.value, query_species)} if query_species_unscreened else ()
        )

        classification_context = ClassificationContext(
            query_gene_ids=frozenset(self._query_gene_ids),
            query_gene_symbols=frozenset(self._query_gene_symbols),
            on_target_transcript_ids=frozenset(self._gene_transcript_ids),
            query_species=query_species,
            index=self._transcript_index,
            repeat_flagged_guides=frozenset(
                normalize_guide_sequence(c.guide_sequence) for c in candidates if c.repeat_flagged
            ),
            requested_species=requested_species,
            ortholog_gene_ids=ortholog_gene_ids,
        )
        # Per-row reference lookups, kept separate from the classification context: they answer
        # "what gene did this land on, and did any reference exist" for every row, whatever its class.
        annotator = HitAnnotator(index=self._transcript_index, query_species=query_species)

        # Fan out deduplicated results: each representative's results apply to all candidates sharing that sequence
        representative_results: dict[str, dict[str, Any]] = {}
        for repr_id, entry in results.items():
            representative_results[repr_id] = entry

        # Pre-seeded with every requested species (zero-filled) so a species that was screened
        # but produced no hits is distinguishable from one never requested (absent key). Species
        # seen on a hit but not requested (unexpected) still get a bucket via setdefault below.
        per_species: dict[str, dict[str, int]] = {
            species: dict.fromkeys(_PER_SPECIES_COUNTERS, 0) for species in requested_species
        }

        stats: dict[str, Any] = {
            "candidates_analyzed": len(candidates),
            "candidates_with_offtargets": 0,
            # Two different quantities, named apart. hit_classes counts ALIGNMENTS, once each, and
            # is filled in after the loop so it agrees with the published hit table and the console
            # line. The candidate-weighted variant sums each candidate's counters over the fanned-out
            # guide, which is what the gates act on: on the frozen baseline the two are 43,536 and
            # ~632,000, and publishing the second as "hit_classes" made them look like one number.
            "hit_classes": dict.fromkeys((member.value for member in HitClass), 0),
            "hit_classes_candidate_weighted": dict.fromkeys((member.value for member in HitClass), 0),
            "query_gene_transcripts_recognised": len(self._gene_transcript_ids),
            "ortholog_symbol_lookup_misses": 0,
            "species_index_misses": 0,
            "per_species": per_species,
            "failed_perfect_match": 0,
            "failed_transcriptome_1mm": 0,
            "failed_transcriptome_2mm": 0,
            "failed_transcriptome_seed_perfect": 0,
            "failed_mirna_seed": 0,
            "failed_high_risk_mirna": 0,
            "failed_excess_off_targets": 0,
            "failed_isoform_coverage": 0,
            "human_transcriptome_hits": 0,
            "other_transcriptome_hits": 0,
            "human_mirna_hits": 0,
            "other_mirna_hits": 0,
            # User intent vs the default species list, recorded (not used to gate scoring).
            "species_explicitly_requested": self._species_explicitly_requested,
            # Candidates whose final score is NOT a post-screen score, so a degraded run is
            # countable rather than merely visible in the log.
            "candidates_not_scored_after_screening": 0,
            "unscreened_species": sorted(unscreened),
        }

        for candidate in candidates:
            candidate_id = candidate.id

            # O(1) lookup: candidate id -> its representative's id -> that representative's
            # results. Falls back to the candidate's own id when it was never deduplicated.
            repr_id = self._candidate_id_to_representative.get(candidate_id, candidate_id)
            offtarget_entry = representative_results.get(repr_id)

            # Zero hits means "clean" only for a candidate that actually reached the aligner. The
            # dedup map holds every submitted candidate and is empty only when screening bypassed
            # _prepare_offtarget_input, so an id missing from a populated map was never submitted
            # and its counts are unknown, not zero (issue #78 §3).
            never_submitted = bool(self._candidate_id_to_representative) and (
                candidate_id not in self._candidate_id_to_representative
            )
            if candidate.screen_query_id is None and not never_submitted:
                # A caller that bypassed _prepare_offtarget_input still gets the join key it
                # screened under -- with an empty dedup map the candidate's own id IS the qname.
                # Guarded on never_submitted so a candidate the aligner never saw is given no join key.
                candidate.screen_query_id = repr_id
            if never_submitted:
                logger.error(
                    f"Candidate {candidate_id} was never submitted to off-target screening; "
                    "its hit counts are unknown and it keeps its design-time score."
                )

            # Mark as screened even if no hits -- but only when the alignments it would have
            # appeared in were actually produced. off_target_screened=False therefore means "this
            # candidate's screen was incomplete": any hit counts written below are a LOWER BOUND,
            # not a total. They are still written, because hits that were found are real evidence
            # and are still applied as filters -- a row showing a filter failure with zeroed counts
            # would be the more confusing contradiction (see the field docs in models/sirna.py).
            unscreened_candidate = query_species_unscreened or never_submitted
            candidate.off_target_screened = not unscreened_candidate

            # Which channel x SPECIES pairs THIS candidate holds evidence for. Pairs, not channels: a
            # gate whose scope is every screened species -- max_off_target_count is the live one --
            # passed on a lower bound whenever a secondary species failed, because completeness was
            # decided query-species-only. A candidate the aligner never saw holds nothing.
            complete_pairs: frozenset[tuple[str, str]] = (
                frozenset() if never_submitted else run_complete_pairs - unscreened_pairs
            )

            if not offtarget_entry or not offtarget_entry.get("hits"):
                # No hits: compute post-screen score with zero off-targets. When nothing was
                # aligned, that zero is an absence of evidence, and awarding
                # off_target_sub_score(0) = 1.0 would float every candidate to the top of the
                # ranking on a run that failed. Keep the design-time score instead.
                if unscreened_candidate:
                    candidate.scored_after_screening = False
                    stats["candidates_not_scored_after_screening"] += 1
                else:
                    self._score_and_gate(
                        candidate, HitClassCounts(), conservation_denominator, stats, unresolved_orthology
                    )

                # Then gate it, on the same footing as a candidate that had hits. This branch used to
                # `continue`, so a completed screen that found nothing reached no gate at all and
                # exported every one as `not_evaluated` (#106). The counts are zero and
                # ``complete_pairs`` says whether that zero was measured.
                self._gate_offtarget_counts(
                    candidate,
                    counts=_ZERO_OFFTARGET_COUNTS,
                    filter_criteria=filter_criteria,
                    complete_pairs=complete_pairs,
                    stats=stats,
                )
                continue

            stats["candidates_with_offtargets"] += 1

            # Classify each transcriptome hit and aggregate miRNA hits (unchanged)
            hit_counts = HitClassCounts()
            # transcriptome_totals/_human are stratified SUBSETS of transcriptome_off_target_total
            # / transcriptome_human_total below (nm>=3 hits count toward the totals but land in
            # no bucket), so the two families never contradict each other.
            transcriptome_totals = {0: 0, 1: 0, 2: 0}
            transcriptome_query = {0: 0, 1: 0, 2: 0}
            transcriptome_off_target_total = 0
            transcriptome_query_total = 0
            transcriptome_human_total = 0
            transcriptome_seed_0mm = 0
            mirna_total = 0
            mirna_query_total = 0
            mirna_human_total = 0
            mirna_0mm_seed_total = 0
            mirna_query_0mm_seed = 0
            mirna_1mm_seed = 0
            mirna_high_risk_total = 0
            mirna_high_risk_query = 0

            for hit in offtarget_entry.get("hits", []):
                nm = int(hit.get("nm", 0))
                seed_mismatches = int(hit.get("seed_mismatches", 0))
                offtarget_score = float(hit.get("offtarget_score", 0.0))
                species_label = hit.get("species")
                # The gates count the query species, not literally human. Testing `is_human_species`
                # here meant that on any non-human query run the mismatch-stratified gates measured
                # zero on real hits: four perfect-match mouse off-targets on a mouse-query run gave
                # observed 0 / verdict pass while the exported column said 4.
                species_is_query = (normalize_species_name(species_label) == query_species) if species_label else False
                species_is_human = is_human_species(species_label)

                is_mirna = "mirna_id" in hit or "database" in hit

                if is_mirna:
                    mirna_total += 1
                    if species_is_human:
                        mirna_human_total += 1
                    if species_is_query:
                        mirna_query_total += 1
                    if seed_mismatches == 0:
                        mirna_0mm_seed_total += 1
                        if species_is_query:
                            mirna_query_0mm_seed += 1
                        if offtarget_score < 5.0:
                            mirna_high_risk_total += 1
                            if species_is_query:
                                mirna_high_risk_query += 1
                    elif seed_mismatches == 1:
                        mirna_1mm_seed += 1
                else:
                    # Classify transcriptome hit using the four-way classifier
                    classification = classify_hit(hit, candidate.guide_sequence, classification_context)
                    # A blank/missing species label belongs to the query species (see classifier).
                    hit_species = normalize_species_name(species_label) if species_label else query_species

                    # Persist the verdict on the hit row, then count from the row that was
                    # written. The class, both symbols and both shortfall flags reach the hit
                    # table and the candidate counters from one place, so a counted hit is always
                    # a published hit. _reconcile_persisted_hits checks that as a row total; it
                    # does not check which candidate a row was attributed to.
                    # No species bucket here: per_species counts alignments, and this loop visits a
                    # deduplicated guide's rows once per candidate carrying it. See below.
                    annotate_hit_row(hit, classification, annotator)
                    hit_class = accumulate_hit_class(hit, hit_counts, None, hit_species)

                    # Only liabilities feed the mismatch-stratified counters. Letting on-target
                    # isoform hits through here would fail every guide on a multi-isoform gene
                    # against max_transcriptome_hits_0mm, which is the near-total kill switch
                    # fixed in 0.5.2. UNDETERMINED is in LIABILITY_CLASSES, so a run with no
                    # transcript index gates exactly as it did before the class existed.
                    if hit_class not in LIABILITY_CLASSES:
                        continue

                    # Every liability counts toward the totals regardless of nm, so
                    # transcriptome_hits_total agrees with off_target_count even when nm>=3 hits
                    # occur (the default exhaustive search no longer caps hits at nm<=2).
                    transcriptome_off_target_total += 1
                    treated_as_query = species_is_query or not species_label
                    if treated_as_query:
                        transcriptome_query_total += 1
                    # Kept alongside the query-species counter: the published human/other
                    # decomposition is a statement about human specifically, and stays true on a run
                    # whose query species is something else.
                    if species_is_human or not species_label:
                        transcriptome_human_total += 1

                    if nm == 0:
                        transcriptome_totals[0] += 1
                        if treated_as_query:
                            transcriptome_query[0] += 1
                    elif nm == 1:
                        transcriptome_totals[1] += 1
                        if treated_as_query:
                            transcriptome_query[1] += 1
                    elif nm == 2:
                        transcriptome_totals[2] += 1
                        if treated_as_query:
                            transcriptome_query[2] += 1

                    if seed_mismatches == 0:
                        transcriptome_seed_0mm += 1

            # Totals count every genuine off-target hit (any nm); the _totals/_human dicts above
            # are stratified nm<=2 subsets, not addends -- do not replace these with a sum().
            transcriptome_total_hits = transcriptome_off_target_total
            human_transcriptome_hits = transcriptome_human_total

            # Write per-candidate hit class fields
            candidate.on_target_hits = hit_counts.on_target
            candidate.ortholog_hits = hit_counts.ortholog
            candidate.repeat_hits = hit_counts.repeat
            # On-target, ortholog and repeat hits excluded; undecidable ones included and reported
            # separately, so a missing reference cannot loosen the screen (see LIABILITY_CLASSES).
            candidate.off_target_count = liabilities_counted(hit_counts)
            candidate.undetermined_hits = hit_counts.undetermined
            candidate.ortholog_species = ",".join(sorted(hit_counts.ortholog_species))

            # Legacy fields (still needed for reporting and miRNA filters)
            candidate.transcriptome_hits_total = transcriptome_total_hits
            candidate.transcriptome_hits_0mm = transcriptome_totals[0]
            candidate.transcriptome_hits_1mm = transcriptome_totals[1]
            candidate.transcriptome_hits_2mm = transcriptome_totals[2]
            candidate.transcriptome_hits_seed_0mm = transcriptome_seed_0mm
            # The gate inputs, written onto the row rather than left as locals: these are the numbers
            # the six query-species-stratified gates compare, and without them a client re-applying a
            # descriptor read the all-species column and disagreed with the run (#101).
            candidate.transcriptome_hits_0mm_query = transcriptome_query[0]
            candidate.transcriptome_hits_1mm_query = transcriptome_query[1]
            candidate.transcriptome_hits_2mm_query = transcriptome_query[2]
            candidate.mirna_hits_0mm_seed_query = mirna_query_0mm_seed
            candidate.mirna_hits_high_risk_query = mirna_high_risk_query
            candidate.total_offtarget_hits_query = transcriptome_query_total + mirna_query_total
            candidate.on_target_confirmed = hit_counts.on_target > 0
            candidate.mirna_hits_total = mirna_total
            candidate.mirna_hits_0mm_seed = mirna_0mm_seed_total
            candidate.mirna_hits_1mm_seed = mirna_1mm_seed
            candidate.mirna_hits_high_risk = mirna_high_risk_total
            candidate.off_target_penalty = offtarget_entry.get("off_target_score", 0.0)

            # Update global stats
            candidate_weighted = stats["hit_classes_candidate_weighted"]
            candidate_weighted["on_target"] += hit_counts.on_target
            candidate_weighted["ortholog"] += hit_counts.ortholog
            candidate_weighted["repeat"] += hit_counts.repeat
            candidate_weighted["off_target"] += hit_counts.off_target
            candidate_weighted["undetermined"] += hit_counts.undetermined
            stats["ortholog_symbol_lookup_misses"] += hit_counts.symbol_lookup_missing
            stats["species_index_misses"] += hit_counts.no_species_index
            stats["human_transcriptome_hits"] += human_transcriptome_hits
            stats["other_transcriptome_hits"] += transcriptome_total_hits - human_transcriptome_hits
            stats["human_mirna_hits"] += mirna_human_total
            stats["other_mirna_hits"] += mirna_total - mirna_human_total

            # Score candidate with post-screen terms. The hits above are real evidence even on a
            # partial run, but a run missing the query species cannot produce a trustworthy
            # off-target term, so those candidates keep their design-time score. The filters below
            # still apply: real hits can only fail a candidate, never wrongly pass one.
            if unscreened_candidate:
                candidate.scored_after_screening = False
                stats["candidates_not_scored_after_screening"] += 1
            else:
                self._score_and_gate(candidate, hit_counts, conservation_denominator, stats, unresolved_orthology)

            # Apply filtering criteria
            self._gate_offtarget_counts(
                candidate,
                counts=OffTargetGateCounts(
                    transcriptome_0mm=transcriptome_query[0],
                    transcriptome_1mm=transcriptome_query[1],
                    transcriptome_2mm=transcriptome_query[2],
                    transcriptome_seed_0mm=transcriptome_seed_0mm,
                    mirna_0mm_seed=mirna_query_0mm_seed,
                    mirna_high_risk=mirna_high_risk_query,
                    total_hits=transcriptome_query_total + mirna_query_total,
                    genuine_off_target_count=liabilities_counted(hit_counts),
                ),
                filter_criteria=filter_criteria,
                complete_pairs=complete_pairs,
                stats=stats,
            )

        stats["hits_classified_without_candidate"] = self._classify_orphan_hit_rows(
            offtarget_data, classification_context, annotator
        )

        # Counted here, after every row has a class, and over each alignment exactly once: the
        # per-candidate loop above visits a deduplicated guide's rows once per candidate carrying it.
        alignment_rows = self._alignment_rows(results)
        stats["hit_classes"] = count_persisted_classes(alignment_rows)
        # per_species is the same quantity decomposed by species, so it is counted the same way.
        # Filled from the loop it read 2x on the frozen baseline's deduplicated guides, which put
        # two contradicting decompositions of one number side by side in workflow_summary.json.
        self._tally_per_species(alignment_rows, per_species, query_species)

        # Re-ranking (excluding repeat-flagged candidates) happens in step5_offtarget_analysis,
        # where design_results is in scope to receive the reordered candidates/top_candidates.
        return candidates, stats

    @staticmethod
    def _tally_per_species(
        alignment_rows: Sequence[Mapping[str, Any]],
        per_species: dict[str, dict[str, int]],
        query_species: str,
    ) -> None:
        """Decompose the alignment-level class tally by species, in place.

        Requested species keep their zero-filled bucket, so "screened and clean" stays
        distinguishable from "never requested" (absent key). A blank species label belongs to the
        query species, exactly as the classifier reads it.
        """
        discarded = HitClassCounts()
        for hit in alignment_rows:
            if not is_annotated(hit):
                continue
            species_label = hit.get("species")
            hit_species = normalize_species_name(species_label) if species_label else query_species
            bucket = per_species.setdefault(hit_species, dict.fromkeys(_PER_SPECIES_COUNTERS, 0))
            accumulate_hit_class(hit, discarded, bucket, hit_species)

    @staticmethod
    def _alignment_rows(results: Mapping[str, Any]) -> list[Mapping[str, Any]]:
        """Every transcriptome alignment row the ingest read, once each.

        Keyed by qname, so a row shared by many candidates through their representative appears
        here once -- unlike the per-candidate loop, which sees it once per candidate. miRNA rows are
        excluded: they carry no class.
        """
        rows: list[Mapping[str, Any]] = []
        for entry in cast(dict[str, dict[str, Any]], results).values():
            for hit in cast(list[Mapping[str, Any]], entry.get("hits") or []):
                if "mirna_id" not in hit and "database" not in hit:
                    rows.append(hit)
        return rows

    @staticmethod
    def _classify_orphan_hit_rows(
        offtarget_data: Mapping[str, Any],
        classification_context: ClassificationContext,
        annotator: HitAnnotator,
    ) -> int:
        """Classify hit rows whose qname matched no candidate in this run.

        Normally there are none. A reused or stale results directory produces some, and a table
        read back from disk can carry a blank ``hit_class`` cell, so both are classified from their
        own query sequence and counted. Every table the parser read is covered, which is what lets
        the writer publish unconditionally instead of guarding against blank verdicts.
        """
        tables = cast(list[dict[str, Any]], offtarget_data.get("transcriptome_hit_tables") or [])
        orphans = [
            row
            for table in tables
            for row in cast(list[dict[str, Any]], table.get("rows") or [])
            if not is_annotated(row)
        ]
        for row in orphans:
            annotate_hit_row(row, classify_hit(row, str(row.get("qseq") or ""), classification_context), annotator)
        if orphans:
            logger.warning(
                f"{len(orphans)} aggregated off-target row(s) belong to no candidate in this run; "
                "classified from their own query sequence and excluded from candidate counters."
            )
        return len(orphans)

    def _score_and_gate(
        self,
        candidate: SiRNACandidate,
        hit_counts: HitClassCounts,
        conservation_denominator: frozenset[str],
        stats: dict[str, Any],
        unresolved_orthology_species: frozenset[str] = frozenset(),
    ) -> None:
        """Score one screened candidate, then apply the gates that need post-screen evidence.

        Both the no-hits and the with-hits branches above go through here, so the isoform-coverage
        gate cannot be wired into only one of them -- which it was, leaving the clean-screen case
        (the common one) ungated.
        """
        if not self._score_candidate_post_screen(
            candidate, hit_counts, conservation_denominator, unresolved_orthology_species
        ):
            stats["candidates_not_scored_after_screening"] += 1

        # isoform_coverage is not a scoring term; when a floor is configured it is a gate, and this
        # is the first point at which the coverage fraction exists.
        if self._apply_isoform_coverage_gate(candidate):
            stats["failed_isoform_coverage"] += 1

    def _apply_isoform_coverage_gate(self, candidate: SiRNACandidate) -> bool:
        """Record the protein-coding isoform coverage verdict, rejecting only if the action says to.

        A gate, not a scoring term: isoform coverage is reported on every candidate and read by no
        weight. The floor lives on FilterCriteria and defaults to None (off), so default behaviour
        is unchanged.

        Goes through the shared verdict recorder rather than assigning the rejection label by hand, so
        a resolved ``warn`` is honoured (#105). The recorder also gives the
        two non-rejecting outcomes somewhere to live -- PASS at or above the floor, and UNKNOWN when
        coverage could not be computed, which never rejects because an annotation gap is not evidence of
        poor coverage. Whether such a candidate can still *qualify* is decided in
        :meth:`_apply_post_screen_ranking`, not here -- which is also where
        :meth:`_complete_isoform_coverage_verdicts` calls this gate for a candidate that reached no
        scoring at all, so the verdict exists on every path and not only the one that integrated hits.

        A floor, so it declares ``Comparator.AT_LEAST`` and no evidence pairs: coverage comes from the
        transcript annotation rather than from a screening channel, which is why an incomplete channel
        never excuses it (see ``POST_SCREEN_FILTER_CHANNELS``).

        Returns True when this gate rejected the candidate, so the caller can count it.
        """
        floor = self.config.design_params.filters.min_isoform_coverage
        coverage = candidate.isoform_coverage
        spec = GateSpec(
            filter_id="min_isoform_coverage",
            threshold=floor,
            action=self._filter_action_for("min_isoform_coverage", floor),
            comparator=Comparator.AT_LEAST,
        )
        # In force nowhere still writes a word: the exported column is fixed, and an empty cell there
        # reads as a verdict rather than as "this gate was not applied". record_gate_outcome does that.
        outcome = evaluate_gate(spec, coverage)
        record_gate_outcome(candidate, outcome, status=SiRNACandidate.FilterStatus.LOW_ISOFORM_COVERAGE)
        if outcome.evaluation is not FilterEvaluation.FAIL:
            return False
        # A decided FAIL means both numbers exist, so the message can quote the comparison it made.
        logger.info(
            f"Candidate {candidate.id} failed isoform coverage ({spec.action.value}): {coverage:.3f} < {floor:.3f}"
        )
        return outcome.rejects

    def _score_candidate_post_screen(
        self,
        candidate: SiRNACandidate,
        hit_counts: HitClassCounts,
        conservation_denominator: frozenset[str],
        unresolved_orthology_species: frozenset[str] = frozenset(),
    ) -> bool:
        """Compute post-screen composite score with the full term set.

        Args:
            candidate: Candidate to score
            hit_counts: Aggregated hit class counts
            unresolved_orthology_species: Species whose orthologue lookup did not complete. Their
                conservation is unknown rather than zero, so any overlap with the denominator leaves
                ``conservation_score`` null. Empty is the honest default for a caller with no
                orthology outcome to report.
            conservation_denominator: Non-query species handed to the aligner. A species whose
                alignment failed stays in here, so a degraded run cannot outscore a complete one.

        Returns:
            True when the candidate now carries a post-screen score; False when scoring failed
            and the design-time score was kept (the caller counts these so a degraded run is
            visible, and ranking demotes them). A failure after an earlier *success* on the same
            candidate additionally clears that earlier score -- see :meth:`_clear_post_screen_score`.
        """
        # Build features from design-time component scores (reuse existing sub-scores)
        features: dict[str, float] = {}
        cs = candidate.component_scores or {}

        # Design-stage terms: reuse from component_scores, dropping any NaN. A term missing here
        # makes the vector unsatisfiable and compute_composite raises -- weights are never
        # renormalised, so the candidate keeps no post-screen score rather than a rescaled one.
        for term in ("asymmetry", "gc_content", "target_accessibility"):
            value = cs.get(term)
            if value is not None and not math.isnan(value):
                features[term] = value

        # miRNA mode scores three more declared terms. They are re-derived from the sequences
        # rather than read back off the row, so a candidate that never passed through
        # MiRNADesigner (a dirty control is a deep copy of a *rejected* candidate) is scored on
        # the same vector as everything else instead of landing on its own scale.
        vector = self.config.design_params.scoring.vector_for(
            post_screen=True, design_mode=self.config.design_params.design_mode
        )
        if self.config.design_params.design_mode == DesignMode.MIRNA:
            features.update(biogenesis_features(candidate.guide_sequence, candidate.passenger_sequence))

        # The sub-score helpers raise on a numerator exceeding its denominator. Both numerators are
        # constructed as subsets of their denominators here, so neither is reachable -- but the
        # guarantee rests on invariants several call sites away, and the cost of one being broken
        # later must not be an aborted run after screening has already been paid for. It must be
        # loud instead: see the ERROR log below and the caller's degraded-run counter.
        try:
            features["off_target"] = off_target_sub_score(liabilities_counted(hit_counts))

            # Numerator is how many of the query gene's protein-coding transcripts contain THIS
            # guide (from step3's guide->source-transcripts map), not the single transcript the
            # candidate happened to be enumerated from. Absent for design_from_sequence/miRNA
            # paths, in which case the term stays inactive rather than computing a wrong number.
            guide_transcripts = self._guide_to_transcripts.get(normalize_guide_sequence(candidate.guide_sequence))
            if guide_transcripts is not None:
                isoform_cov = isoform_coverage_sub_score(
                    len(self._protein_coding_transcript_ids & guide_transcripts),
                    self._protein_coding_transcript_count,
                )
                if isoform_cov is not None:
                    features["isoform_coverage"] = isoform_cov
                    candidate.isoform_coverage = isoform_cov

            # Numerator is intersected with the denominator on purpose: an ortholog hit in a species
            # this run never asked the aligner for (a stale cached result, a hand-edited hit table)
            # is not part of this ratio, and counting it would make the fraction exceed 1.
            conserved_species = hit_counts.ortholog_species & conservation_denominator
            unexpected_species = hit_counts.ortholog_species - conservation_denominator
            if unexpected_species:
                logger.warning(
                    f"Ortholog hits for {candidate.id} in unrequested species {sorted(unexpected_species)}; "
                    "excluded from the conservation term."
                )
            # An unresolved lookup cannot populate the numerator, so the ratio would report
            # "conserved in 0 of N" for a question nobody asked -- a 0.0 that reads as a measurement.
            # Null instead. The denominator is NOT shrunk to the resolvable species: that let a degraded
            # run outscore the complete run it degraded from (#80 F10).
            unresolved_in_denominator = conservation_denominator & unresolved_orthology_species
            if unresolved_in_denominator:
                candidate.conservation_score = None
                logger.debug(
                    f"Conservation unavailable for {candidate.id}: orthology unresolved for "
                    f"{sorted(unresolved_in_denominator)}."
                )
            else:
                conservation = conservation_sub_score(len(conserved_species), len(conservation_denominator))
                if conservation is not None:
                    features["conservation"] = conservation
                    candidate.conservation_score = conservation

            # One vector, applied exactly as declared. The miRNA biogenesis terms are inside it,
            # so there is no bonus to fold in afterwards and nothing to divide the result by.
            result = compute_composite(features, vector)
            candidate.composite_score = result.score
            candidate.weight_set_version = result.weight_set_version
            candidate.weight_vector = result.vector_name
            candidate.scored_after_screening = True
            # Write per-term contributions
            candidate.score_off_target = result.contributions.get("off_target")
            candidate.score_target_accessibility = result.contributions.get("target_accessibility")
            candidate.score_asymmetry = result.contributions.get("asymmetry")
            candidate.score_gc_content = result.contributions.get("gc_content")
            candidate.score_ago_start = result.contributions.get("ago_start")
            candidate.score_pos1_mismatch = result.contributions.get("pos1_mismatch")
            candidate.score_supp_13_16 = result.contributions.get("supp_13_16")
            return True
        except (ScoringError, ValueError) as exc:
            # ERROR, not WARNING: the candidate now carries a design-time score that is not
            # comparable with its neighbours' post-screen scores, so the run is degraded.
            logger.error(f"Post-screen scoring failed for {candidate.id}: {exc}. Keeping design-time score.")
            if candidate.scored_after_screening:
                self._clear_post_screen_score(candidate)
            candidate.scored_after_screening = False
            return False

    @staticmethod
    def _clear_post_screen_score(candidate: SiRNACandidate) -> None:
        """Drop a post-screen score this method wrote earlier, once a later attempt has failed (#100).

        Only reached when ``scored_after_screening`` was already True on entry, i.e. an earlier call
        succeeded and this one did not. Without it the row published ``scored_after_screening=False``
        beside a full set of post-screen fields from the first attempt -- a composite the run no longer
        stands behind, and the number ``ranking_score`` would still rank on.

        Deliberately narrow in two ways. It clears exactly the ten fields the success path writes, all
        of which the earlier success had already overwritten -- and it runs only on that retry, so a
        candidate whose *first* attempt fails keeps whatever the design stage or a caller left there
        (the mixed-scale ranking rule reads those). ``isoform_coverage`` and ``conservation_score``
        stay either way: those are measurements the run made, not a score it claims.
        """
        candidate.composite_score = None
        candidate.weight_set_version = ""
        candidate.weight_vector = ""
        candidate.score_off_target = None
        candidate.score_target_accessibility = None
        candidate.score_asymmetry = None
        candidate.score_gc_content = None
        candidate.score_ago_start = None
        candidate.score_pos1_mismatch = None
        candidate.score_supp_13_16 = None

    def _generate_orf_report(self, orf_results: dict[str, Any], report_file: Path) -> DataFrame[ORFValidationSchema]:
        """Generate ORF validation report in tab-delimited format with schema validation.

        Returns:
            Validated DataFrame conforming to ORFValidationSchema
        """
        # Handle empty results case
        if not orf_results:
            logger.warning("No ORF results to report - creating empty report file")
            report_file.parent.mkdir(parents=True, exist_ok=True)
            # Create empty DataFrame with required columns for schema validation
            empty_df = pd.DataFrame(
                columns=[
                    "transcript_id",
                    "sequence_length",
                    "gc_content",
                    "orfs_found",
                    "has_valid_orf",
                    "longest_orf_start",
                    "longest_orf_end",
                    "longest_orf_length",
                    "longest_orf_frame",
                    "start_codon",
                    "stop_codon",
                    "orf_gc_content",
                    "utr5_length",
                    "utr3_length",
                    "predicted_sequence_type",
                ]
            )
            # Set correct dtypes to match schema - using Any types for nullable fields
            empty_df = empty_df.astype(
                {
                    "transcript_id": str,
                    "sequence_length": "Int64",
                    "gc_content": float,
                    "orfs_found": "Int64",
                    "has_valid_orf": bool,
                    "longest_orf_start": "object",
                    "longest_orf_end": "object",
                    "longest_orf_length": "object",
                    "longest_orf_frame": "object",
                    "start_codon": "object",
                    "stop_codon": "object",
                    "orf_gc_content": "object",
                    "utr5_length": "object",
                    "utr3_length": "object",
                    "predicted_sequence_type": "object",
                }
            )
            validated_df = ORFValidationSchema.validate(empty_df)
            validated_df.to_csv(report_file, sep="\t", index=False)
            return validated_df

        # Prepare data for DataFrame
        rows: list[dict[str, Any]] = []
        for transcript_id, analysis in orf_results.items():
            row_data: dict[str, Any] = {
                "transcript_id": transcript_id,
                "sequence_length": getattr(analysis, "sequence_length", None),
                "gc_content": getattr(analysis, "gc_content", None),
                "orfs_found": len(getattr(analysis, "orfs", []) or []),
                "has_valid_orf": getattr(analysis, "has_valid_orf", False),
                "utr5_length": getattr(analysis, "utr5_length", None),
                "utr3_length": getattr(analysis, "utr3_length", None),
                "predicted_sequence_type": getattr(
                    getattr(analysis, "sequence_type", None), "value", str(getattr(analysis, "sequence_type", ""))
                ),
            }

            if getattr(analysis, "longest_orf", None):
                orf = analysis.longest_orf
                row_data.update(
                    {
                        "longest_orf_start": orf.start_pos,
                        "longest_orf_end": orf.end_pos,
                        "longest_orf_length": orf.length,
                        "longest_orf_frame": orf.reading_frame,
                        "start_codon": orf.start_codon,
                        "stop_codon": orf.stop_codon,
                        "orf_gc_content": orf.gc_content,
                    }
                )
            else:
                row_data.update(
                    {
                        "longest_orf_start": None,
                        "longest_orf_end": None,
                        "longest_orf_length": None,
                        "longest_orf_frame": None,
                        "start_codon": None,
                        "stop_codon": None,
                        "orf_gc_content": None,
                    }
                )
            rows.append(row_data)

        # Create DataFrame and validate with pandera - let failures bubble up
        df = pd.DataFrame(rows)
        logger.debug(f"Validating ORF report DataFrame with {len(df)} rows")

        # Validate DataFrame with our validation middleware
        orf_validation = self.validation.validate_dataframe_output(df, "orf_validation")
        if not orf_validation.overall_result.is_valid:
            logger.warning(f"ORF DataFrame validation issues: {len(orf_validation.overall_result.errors)} errors")

        # Runtime validation with Pandera schema
        validated_df = ORFValidationSchema.validate(df)
        logger.info(f"ORF report schema validation passed for {len(validated_df)} transcripts")

        # Write validated DataFrame to file
        validated_df.to_csv(report_file, sep="\t", index=False)

        return validated_df

    def _summarize_transcripts(self, transcripts: list[TranscriptInfo]) -> dict[str, Any]:
        """Summarize transcript retrieval results."""
        return {
            "total_transcripts": len(transcripts),
            "transcript_types": list({t.transcript_type for t in transcripts}),
            "databases": list({t.database for t in transcripts}),
            "avg_length": (
                sum(t.length for t in transcripts if t.length is not None)
                / len([t for t in transcripts if t.length is not None])
                if any(t.length is not None for t in transcripts)
                else 0
            ),
        }

    def _summarize_orf_results(self, orf_results: dict[str, Any]) -> dict[str, Any]:
        """Summarize ORF validation results."""
        results = orf_results.get("results", {})
        valid_count = sum(1 for r in results.values() if r.has_valid_orf)

        return {
            "total_analyzed": len(results),
            "valid_orfs": valid_count,
            "validation_rate": valid_count / len(results) if results else 0,
        }

    def _summarize_design_results(self, design_results: DesignResult) -> dict[str, Any]:
        """Summarize siRNA design results, including the inputs design never covered.

        ``design_input_shortfalls`` is the only place a dropped transcript appears: ``base`` is
        derived from the DesignResults that were produced, and ``input_sequences`` still counts every
        transcript handed to design, so on its own it reports a target set the run never covered
        (#100). The count sits beside the map for the same reason ``repeat_flagged_count`` does --
        a reader comparing it against ``input_sequences`` sees the divergence without parsing reasons.
        """
        base = design_results.get_summary()
        total = design_results.total_candidates
        passed = design_results.filtered_candidates
        failed = max(0, total - passed)
        # Flagged, which is not the same as excluded once the repeat gate's action can be `warn`:
        # naming it *_excluded_count put two contradicting numbers for one quantity in the same JSON,
        # this one and selection_summary.repeat_excluded. This is the descriptive count; the selection
        # summary owns what was actually held out.
        repeat_flagged = sum(1 for c in design_results.candidates if c.repeat_flagged)
        base.update(
            {
                "pass_count": passed,
                "fail_count": failed,
                "repeat_flagged_count": repeat_flagged,
                "repeat_threshold_fraction": self.config.design_params.filters.max_repeat_transcript_fraction,
                "top_n_requested": self.config.top_n,
                "dirty_controls_added": getattr(self, "_dirty_controls_added", 0),
                "threads_used": self.config.num_threads,
                "design_input_shortfalls": dict(self._design_input_shortfalls),
                "design_input_dropped_count": len(self._design_input_shortfalls),
            }
        )
        return base


def _parse_env_flag(value: str | None) -> bool:
    """Parse a boolean-ish environment value."""
    if value is None:
        return False
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _load_zfn_sharding_overrides_from_env() -> dict[str, Any]:
    """Load optional JSON sharding overrides from ``SIRNAFORGE_ZFN_SHARDING_JSON``."""
    raw = os.getenv("SIRNAFORGE_ZFN_SHARDING_JSON")
    if not raw:
        return {}

    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as exc:
        logger.warning("Ignoring invalid SIRNAFORGE_ZFN_SHARDING_JSON: %s", exc)
        return {}

    if not isinstance(parsed, dict):
        logger.warning("Ignoring SIRNAFORGE_ZFN_SHARDING_JSON: expected object, got %s", type(parsed).__name__)
        return {}

    return cast(dict[str, Any], parsed)


def _normalize_sharding_override_value(overrides: dict[str, Any]) -> dict[str, Any]:
    """Normalize override aliases/units for ZFN sharding payloads."""
    normalized = dict(overrides)

    chunk_size_mb = normalized.pop("chunk_size_mb", None)
    if chunk_size_mb is not None and "chunk_size_bp" not in normalized:
        try:
            normalized["chunk_size_bp"] = int(float(chunk_size_mb) * 1_000_000)
        except (TypeError, ValueError):
            logger.warning("Ignoring invalid ZFN sharding chunk_size_mb: %r", chunk_size_mb)

    chromosomes = normalized.get("chromosomes")
    if isinstance(chromosomes, str):
        normalized["chromosomes"] = parse_csv(chromosomes)

    return normalized


def apply_zfn_runtime_overrides(
    zfn_design_params: ZFNDesignParameters,
    nextflow_config_overrides: dict[str, Any],
) -> ZFNDesignParameters:
    """Apply non-CLI ZFN sharding/runtime overrides and project Nextflow params.

        Starts from ``zfn_design_params.sharding`` (typed defaults are authoritative),
        merges optional JSON overrides from ``SIRNAFORGE_ZFN_SHARDING_JSON``, and
        mirrors the resolved sharding values into ``nextflow_config_overrides`` so
        the Nextflow route and direct Python route share the same effective config.

    The runtime search implementation remains generic and contig-aware, including
    chunk sharding on single-contig inputs when sharding is enabled.
    """
    merged_overrides: dict[str, Any] = {}

    if zfn_design_params.sharding.enabled:
        merged_overrides.update(zfn_design_params.sharding.model_dump(mode="python"))

    merged_overrides.update(_normalize_sharding_override_value(_load_zfn_sharding_overrides_from_env()))

    if merged_overrides:
        sharding = ZFNShardingConfig(**merged_overrides)
        zfn_design_params = zfn_design_params.model_copy(update={"sharding": sharding}, deep=True)

    # Optional Nextflow route can be enabled without CLI changes.
    if _parse_env_flag(os.getenv("SIRNAFORGE_ZFN_USE_NEXTFLOW")):
        nextflow_config_overrides["design_mode"] = "zfn"

    sharding_cfg = zfn_design_params.sharding
    nextflow_config_overrides.setdefault("zfn_sharding_enabled", sharding_cfg.enabled)
    nextflow_config_overrides.setdefault("zfn_shard_chunk_mb", max(1, sharding_cfg.chunk_size_bp // 1_000_000))
    nextflow_config_overrides.setdefault("zfn_shard_overlap_bp", sharding_cfg.overlap_bp)
    nextflow_config_overrides.setdefault("zfn_shard_chromosomes", ",".join(sharding_cfg.chromosomes))

    return zfn_design_params


# Convenience function for running complete workflow
async def run_sirna_workflow(
    gene_query: str,
    output_dir: str,
    input_fasta: str | None = None,
    database: str = "ensembl",
    design_mode: str = "sirna",
    top_n_candidates: int | None = None,
    screen_species: list[str] | None = None,
    query_species: str | None = None,
    transcriptome_indices: str | None = None,
    mirna_database: str = "mirgenedb",
    mirna_species: Sequence[str] | None = None,
    transcriptome_fasta: str | None = None,
    transcriptome_filter: str | None = None,
    transcriptome_selection: ReferenceSelection | None = None,
    ortholog_mapping_file: Path | str | None = None,
    resolved_policy: ResolvedRunPolicy | None = None,
    run_mode: str | None = None,
    policy_config: Path | str | None = None,
    filter_actions: Mapping[str, Any] | None = None,
    gc_min: float | None = None,
    gc_max: float | None = None,
    sirna_length: int | None = None,
    modification_pattern: str | None = None,
    overhang: str | None = None,
    zfn_design_params: ZFNDesignParameters | None = None,
    zfn_annotation: GenomicAnnotationConfig | None = None,
    check_off_targets: bool | None = None,
    # Variant targeting parameters
    variant_ids: list[str] | None = None,
    variant_vcf_file: Path | None = None,
    variant_mode: str = "avoid",
    variant_min_af: float = 0.01,
    variant_clinvar_filters: str = "Pathogenic,Likely pathogenic",
    variant_assembly: str = "GRCh38",
    log_file: str | None = None,
    write_json_summary: bool = True,
    num_threads: int | None = None,
    allow_transcriptome_with_input_fasta: bool = False,
    default_transcriptome_sources: Sequence[str] = DEFAULT_TRANSCRIPTOME_SOURCES,
    keep_nextflow_work: bool = False,
    nextflow_docker_image: str | None = None,
    max_hits: int | None = None,
    max_off_targets: int | None = None,
    min_asymmetry_score: float | None = None,
    max_paired_fraction: float | None = None,
    min_empirical_score: float | None = None,
    min_isoform_coverage: float | None = None,
    plfold_window: int | None = None,
    plfold_max_bp_span: int | None = None,
    accessibility_log_floor: float | None = None,
    **renamed: Any,
) -> dict[str, Any]:
    """Run complete siRNA design workflow.

    Args:
        gene_query: Gene name or ID to search for
        output_dir: Directory for output files
        input_fasta: Local path or remote URI to an input FASTA file
        database: Database to search (ensembl, refseq, gencode)
        design_mode: Design mode (sirna, mirna, or zfn)
        top_n_candidates: Cap on how many top-ranked candidates are reported (None = no cap, the
            default). Enumeration and screening always cover every candidate.
        screen_species: Species to screen against for off-target liabilities
        query_species: Organism the TARGET transcripts belong to. Defaults to the organism the
            gene-query database serves (human), which is also the species of the default
            transcriptome; set it when designing against an input FASTA from another organism.
        transcriptome_indices: Comma-separated species:/index_prefix transcriptome references the
            caller has already built. Resolved by the same resolver as the defaults, so an override
            and a default differ only in provenance (#99).
        mirna_database: miRNA reference database identifier
        mirna_species: miRNA reference species identifiers
        transcriptome_fasta: Path or URL to transcriptome FASTA for off-target analysis
        transcriptome_filter: Comma-separated filter names (protein_coding, canonical_only)
        transcriptome_selection: Pre-resolved transcriptome selection metadata
        ortholog_mapping_file: JSON mapping of query gene -> species -> orthologue gene IDs. Supply
            it to resolve cross-species orthology offline instead of calling Ensembl Compara.
        resolved_policy: A policy already resolved by ``config.run_policy.resolve_run_policy`` -- how
            the CLI passes its resolution, so the command line and this function cannot diverge.
            Supplying it together with any of the threshold arguments below is an error, because two
            resolutions of the same run could disagree.
        run_mode: design_only, exploratory or qualified. Defaults to qualified; ``check_off_targets``
            set False maps to design_only, as the legacy skip flag does.
        policy_config: JSON/TOML policy file; beats the built-in profile, loses to explicit values.
        filter_actions: ``filter_id -> off|warn|fail``. A filter set to off is not evaluated.
        gc_min: Minimum GC content percentage. None means unstated, so the profile applies.
        gc_max: Maximum GC content percentage. None means unstated, so the profile applies -- which
            for siRNA mode is 60.0.
        sirna_length: siRNA length in nucleotides (None = profile default)
        modification_pattern: Chemical modification pattern (None = profile default)
        overhang: Overhang sequence, dTdT for DNA or UU for RNA (None = profile default, and UU in
            miRNA design mode)
        zfn_design_params: Optional ZFN design parameters for ZFN mode workflow
        zfn_annotation: Optional genomic annotation config for ZFN off-target classification
        check_off_targets: Perform off-target analysis stage. None means unstated (screening on);
            False maps the run mode to design_only.
        variant_ids: List of variant identifiers (rsID, chr:pos:ref:alt, or HGVS) to target or avoid
        variant_vcf_file: Path to VCF file containing variants to target or avoid
        variant_mode: How to handle variants (avoid/target/both) - default is avoid
        variant_min_af: Minimum allele frequency threshold for variant filtering (default: 0.01)
        variant_clinvar_filters: Comma-separated ClinVar significance levels to include (default: Pathogenic,Likely pathogenic)
        variant_assembly: Reference genome assembly for variants (only GRCh38 supported)
        log_file: Path to centralized log file
        write_json_summary: Write logs/workflow_summary.json
        num_threads: Optional override for design parallelism
        allow_transcriptome_with_input_fasta: Opt in to resolving ``default_transcriptome_sources``
            when ``input_fasta`` is supplied (default: False). Left False, an input-FASTA run is
            design-only unless ``transcriptome_fasta`` names a reference explicitly: supplying your
            own sequences should never trigger a multi-gigabyte reference download you did not ask
            for. Set True to screen an input-FASTA run against the bundled defaults.
        default_transcriptome_sources: Ordered list of transcriptome identifiers evaluated by default
        keep_nextflow_work: Keep Nextflow work directory symlink in output
        nextflow_docker_image: Override Docker image used by the embedded Nextflow pipeline
        max_hits: Override the pipeline's per-candidate off-target hit cap (None keeps the pipeline's
            exhaustive default; set a lower value, e.g. 10000, to speed up large gene-family searches)
        max_off_targets: Override the genuine off-target ceiling that gates PASS vs
            EXCESS_OFF_TARGETS (None keeps OffTargetFilterCriteria's default of 15). Unlike
            max_hits this changes the verdict, not how many hits are recorded.
        min_asymmetry_score: Override the thermodynamic asymmetry floor gating LOW_ASYMMETRY
            (None keeps FilterCriteria's default of 0.65).
        max_paired_fraction: Override the guide self-structure ceiling gating EXCESS_PAIRING
            (None keeps FilterCriteria's default of 0.6).
        min_empirical_score: Override the empirical design-rule floor gating
            LOW_EMPIRICAL_SCORE (None keeps FilterCriteria's default).
        min_isoform_coverage: Opt into the protein-coding isoform coverage gate
            (LOW_ISOFORM_COVERAGE). None, the default, means no gate: coverage is reported on
            every candidate either way, and no floor has been calibrated against truth data.
        plfold_window: Override the RNAplfold averaging window W used by the target_accessibility
            term (None keeps the default of 150).
        plfold_max_bp_span: Override the RNAplfold maximum base-pair span L (None keeps 100).
        accessibility_log_floor: Override the log10-probability floor the accessibility term is
            normalised against (None keeps -5.0). Changing any of these three changes the numeric
            scale of composite_score, so results are not comparable with a default run.
        **renamed: Accepted only to refuse the ``genome_*`` names the #99 rename removed, quoting
            the replacement instead of "unexpected keyword argument".

    Returns:
        Dictionary with complete workflow results
    """
    refuse_renamed_arguments(renamed)
    # Resolve the policy exactly once, before anything creates a directory or fetches a reference.
    # An unstated argument is None, so the profile applies to it; a stated one wins. This is the
    # same call the CLI makes -- when the CLI has already made it, its result arrives as
    # resolved_policy and is used verbatim rather than being rebuilt from values.
    stated: dict[str, Any] = {
        key: value
        for key, value in (
            ("gc_min", gc_min),
            ("gc_max", gc_max),
            ("sirna_length", sirna_length),
            ("modification_pattern", modification_pattern),
            ("default_overhang", overhang),
            ("top_n", top_n_candidates),
            ("check_off_targets", check_off_targets),
            ("max_off_target_count", max_off_targets),
            ("min_asymmetry_score", min_asymmetry_score),
            ("max_paired_fraction", max_paired_fraction),
            ("min_empirical_score", min_empirical_score),
            ("min_isoform_coverage", min_isoform_coverage),
            ("plfold_window", plfold_window),
            ("plfold_max_bp_span", plfold_max_bp_span),
            ("accessibility_log_floor", accessibility_log_floor),
        )
        if value is not None
    }
    if modification_pattern is not None:
        stated["apply_modifications"] = modification_pattern.lower() != "none"
    if resolved_policy is not None:
        conflicting = sorted(set(stated) - {"apply_modifications"})
        if conflicting or run_mode or policy_config or filter_actions:
            raise RunPolicyError(
                f"resolved_policy was supplied together with {conflicting or 'run_mode/policy_config/filter_actions'}; "
                "a run must be resolved once, so pass either the policy or the individual settings"
            )
        policy = resolved_policy
    else:
        policy = resolve_run_policy(
            entry_point=EntryPoint.SCREENING_WORKFLOW,
            design_mode=design_mode,
            run_mode=run_mode,
            config_file=policy_config,
            stated=stated,
            filter_actions=filter_actions,
            query_species=query_species,
            screen_species=screen_species or (),
            # The same fact the CLI reports, derived from the same inputs: an input FASTA with no
            # transcriptome argument resolves no transcriptome reference, so this run cannot complete
            # that channel. Without it the API and the CLI resolved the same inputs differently, and
            # the API kept requiring evidence it could never obtain.
            transcriptome_reference_available=(
                False if (input_fasta and not transcriptome_fasta and not transcriptome_indices) else None
            ),
        )
    mode_enum = policy.design_mode
    design_params = policy.design_parameters
    check_off_targets = design_params.check_off_targets
    database_enum = DatabaseType(database.lower())

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    resolved_input: InputSource | None = None
    input_path: Path | None = None
    if input_fasta:
        inputs_dir = output_path / "inputs"
        resolved_input = resolve_input_source(input_fasta, inputs_dir)
        input_path = resolved_input.local_path

    if transcriptome_selection is None:
        input_spec = WorkflowInputSpec(
            input_fasta=input_fasta,
            transcriptome_argument=transcriptome_fasta,
            default_transcriptomes=default_transcriptome_sources,
            # check_off_targets=False means "do not screen", so there is nothing to resolve a
            # reference for. Hardcoding design_only=False here made the flag download and index
            # multi-gigabyte references before announcing that screening was skipped.
            design_only=not check_off_targets,
            allow_transcriptome_for_input_fasta=allow_transcriptome_with_input_fasta,
        )
        resolver = ReferencePolicyResolver(input_spec)
        transcriptome_selection = resolver.resolve_transcriptomes()

    # Configure variant targeting if specified
    variant_config_obj: VariantWorkflowConfig | None = None
    if variant_ids or variant_vcf_file:
        # Parse variant mode using helper that handles normalization
        variant_mode_enum = normalize_variant_mode(variant_mode)

        # Parse ClinVar filters
        clinvar_filters = parse_clinvar_filter_string(variant_clinvar_filters)

        variant_config_obj = VariantWorkflowConfig(
            variant_ids=variant_ids,
            vcf_file=Path(variant_vcf_file) if variant_vcf_file else None,
            variant_mode=variant_mode_enum,
            min_af=variant_min_af,
            clinvar_filter_levels=clinvar_filters,
            assembly=variant_assembly,
        )

    nextflow_config_overrides: dict[str, Any] = {}
    if nextflow_docker_image:
        nextflow_config_overrides["docker_image"] = nextflow_docker_image
    if max_hits is not None:
        nextflow_config_overrides["max_hits"] = max_hits

    # Build ZFN workflow config when in ZFN mode
    zfn_workflow_config: ZFNWorkflowConfig | None = None
    if mode_enum == DesignMode.ZFN and zfn_design_params is not None:
        zfn_design_params = apply_zfn_runtime_overrides(zfn_design_params, nextflow_config_overrides)
        zfn_workflow_config = ZFNWorkflowConfig(
            zfn_params=zfn_design_params,
            annotation=zfn_annotation,
        )

    config = WorkflowConfig(
        output_dir=output_path,
        gene_query=gene_query,
        input_fasta=input_path,
        database=database_enum,
        design_params=design_params,
        transcriptome_indices=transcriptome_indices,
        screen_species=screen_species or ["human", "rat", "rhesus"],
        query_species=query_species,
        mirna_database=mirna_database,
        mirna_species=mirna_species,
        transcriptome_fasta=transcriptome_fasta,
        transcriptome_filter=transcriptome_filter,
        transcriptome_selection=transcriptome_selection,
        ortholog_mapping_file=ortholog_mapping_file,
        log_file=log_file,
        write_json_summary=write_json_summary,
        num_threads=num_threads,
        input_source=resolved_input,
        keep_nextflow_work=keep_nextflow_work,
        variant_config=variant_config_obj,
        nextflow_config=nextflow_config_overrides,
        zfn_config=zfn_workflow_config,
        resolved_policy=policy,
    )

    # Run workflow
    workflow = SiRNAWorkflow(config)
    return await workflow.run_complete_workflow()


if __name__ == "__main__":
    # Example usage
    async def main() -> None:
        """Run example siRNA workflow."""
        with tempfile.TemporaryDirectory() as temp_dir:
            results = await run_sirna_workflow(gene_query="TP53", output_dir=temp_dir, top_n_candidates=20)
            print(f"Workflow completed: {results}")

    asyncio.run(main())


async def run_offtarget_only_workflow(
    input_candidates_fasta: str,
    output_dir: str,
    screen_species: list[str] | None = None,
    query_species: str | None = None,
    transcriptome_indices: str | None = None,
    mirna_database: str = "mirgenedb",
    mirna_species: Sequence[str] | None = None,
    transcriptome_fasta: str | None = None,
    transcriptome_filter: str | None = None,
    transcriptome_selection: ReferenceSelection | None = None,
    ortholog_mapping_file: Path | str | None = None,
    log_file: str | None = None,
    nextflow_docker_image: str | None = None,
    resolved_policy: ResolvedRunPolicy | None = None,
    run_mode: str | None = None,
    policy_config: Path | str | None = None,
    filter_actions: Mapping[str, Any] | None = None,
    **renamed: Any,
) -> dict[str, Any]:
    """Run off-target-only workflow for pre-designed siRNA candidates.

    This is a simplified workflow that only runs the off-target analysis stage
    without transcript retrieval, ORF validation, or siRNA design. It accepts
    pre-designed 21-nt siRNA guide sequences and runs comprehensive off-target
    analysis using the embedded Nextflow pipeline.

    Args:
        input_candidates_fasta: Path to FASTA file with 21-nt siRNA guide sequences
        output_dir: Directory for output files
        screen_species: Species to screen against for off-target liabilities
        query_species: Organism the input guides were designed against (defaults to human)
        transcriptome_indices: Comma-separated species:/index_prefix transcriptome references the
            caller has already built
        mirna_database: miRNA reference database identifier
        mirna_species: miRNA reference species identifiers
        transcriptome_fasta: Path or URL to transcriptome FASTA for off-target analysis
        transcriptome_filter: Comma-separated filter names (protein_coding, canonical_only)
        transcriptome_selection: Pre-resolved transcriptome selection metadata
        ortholog_mapping_file: JSON mapping of query gene -> species -> orthologue gene IDs, used
            instead of Ensembl Compara for cross-species orthology
        log_file: Path to centralized log file
        nextflow_docker_image: Override Docker image used by the embedded Nextflow pipeline
        resolved_policy: A policy already resolved by ``config.run_policy.resolve_run_policy``, which
            is where this path's off-target thresholds come from.
        run_mode: design_only, exploratory or qualified (default: qualified).
        policy_config: JSON/TOML policy file; beats the built-in profile, loses to explicit values.
        filter_actions: ``filter_id -> off|warn|fail``.
        **renamed: Accepted only to refuse the ``genome_*`` names the #99 rename removed, quoting
            the replacement instead of "unexpected keyword argument".

    Returns:
        Dictionary with off-target analysis results: the seven long-standing keys, plus (#100)
        ``selection_summary`` -- the resolved selection, counted by exclusion reason -- and
        ``written_files``. A run that could qualify nobody is a result, not a failure: it reports an
        empty shortlist and still returns normally.
    """
    refuse_renamed_arguments(renamed)
    # Resolved before the output tree is created, exactly as the other two entry points do it.
    policy = resolved_policy or resolve_run_policy(
        entry_point=EntryPoint.OFFTARGET_ONLY,
        run_mode=run_mode,
        config_file=policy_config,
        filter_actions=filter_actions,
        query_species=query_species,
        screen_species=screen_species or (),
    )

    console.print("\n🎯 [bold cyan]Starting Off-Target Analysis (Pre-Designed siRNAs)[/bold cyan]")
    console.print(f"Input Candidates: [yellow]{input_candidates_fasta}[/yellow]")
    console.print(f"Output Directory: [blue]{output_dir}[/blue]")

    start_time = time.perf_counter()

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Create output structure
    (output_path / "results").mkdir(exist_ok=True)
    (output_path / "logs").mkdir(exist_ok=True)

    # Parse input candidates
    input_fasta_path = Path(input_candidates_fasta)
    sequences = FastaUtils.read_fasta(input_fasta_path)

    if not sequences:
        raise ValueError("Input FASTA file is empty")

    console.print(f"📄 Loaded {len(sequences)} siRNA candidates from {input_fasta_path.name}")

    # Convert sequences to SiRNACandidate objects for off-target analysis
    # Calculate metrics for pre-designed candidates using the same methods as design workflow
    candidates: list[SiRNACandidate] = []
    for i, (header, seq) in enumerate(sequences):
        # Extract ID from header (first token)
        candidate_id = header.split()[0] if header else f"candidate_{i}"

        # The input guide sequence is the antisense strand (what targets the mRNA)
        # Generate the sense/passenger strand as the reverse complement
        guide_sequence = seq.upper()
        passenger_sequence = str(Seq(guide_sequence).reverse_complement())

        # Calculate GC content
        gc_count = guide_sequence.count("G") + guide_sequence.count("C")
        gc_content = (gc_count / len(guide_sequence)) * 100 if len(guide_sequence) > 0 else 0.0

        # Calculate thermodynamic properties
        asymmetry_score = 0.0
        duplex_stability = 0.0
        try:
            calc = ThermodynamicCalculator()

            # Create a temporary candidate for thermodynamic calculations
            temp_candidate = SiRNACandidate(
                id=candidate_id,
                transcript_id="pre_designed",
                position=1,
                guide_sequence=guide_sequence,
                passenger_sequence=passenger_sequence,
                length=len(guide_sequence),
                gc_content=gc_content,
                asymmetry_score=0.0,
                paired_fraction=0.0,
                duplex_stability=0.0,
                off_target_count=0,
                off_target_penalty=0.0,
                transcript_hit_count=0,
                transcript_hit_fraction=0.0,
                passes_filters=True,
            )

            # Calculate asymmetry score (5' vs 3' end stability)
            _, _, asymmetry_score = calc.calculate_asymmetry_score(temp_candidate)

            # Calculate duplex stability
            duplex_stability = calc.calculate_duplex_stability(guide_sequence, passenger_sequence)

        except Exception as e:
            # If thermodynamic calculations fail, use default values
            logger.warning(f"Failed to calculate thermodynamics for {candidate_id}: {e}")
            asymmetry_score = 0.0
            duplex_stability = 0.0

        # Create final candidate with computed metrics
        candidate = SiRNACandidate(
            id=candidate_id,
            transcript_id="pre_designed",  # Placeholder since these are pre-designed
            position=1,  # Must be >= 1 per validation
            guide_sequence=guide_sequence,
            passenger_sequence=passenger_sequence,  # Computed as reverse complement
            length=len(guide_sequence),
            gc_content=gc_content,  # Computed from guide sequence
            asymmetry_score=asymmetry_score,  # Computed thermodynamically
            paired_fraction=0.0,  # Not applicable for pre-designed guides
            duplex_stability=duplex_stability,  # Computed thermodynamically
            off_target_count=0,  # Will be populated by off-target analysis
            off_target_penalty=0.0,  # Will be populated by off-target analysis
            transcript_hit_count=0,  # Will be populated by off-target analysis
            transcript_hit_fraction=0.0,  # Will be populated by off-target analysis
            # design_score and composite_score stay None: a pre-designed guide has neither until
            # screening scores it, and 0.0 would read as a computed worst-possible candidate.
            passes_filters=True,  # Assume valid since user provided them
        )
        candidates.append(candidate)

    # Prepare candidates FASTA for off-target analysis
    candidates_fasta = output_path / "input_candidates.fasta"
    candidate_sequences = [(c.id, c.guide_sequence) for c in candidates]
    FastaUtils.save_sequences_fasta(candidate_sequences, candidates_fasta)

    console.print(f"📝 Prepared {len(candidates)} candidates for off-target analysis")

    # Set up Nextflow configuration
    nextflow_config: dict[str, Any] = {}
    if nextflow_docker_image:
        nextflow_config["docker_image"] = nextflow_docker_image

    # Resolve transcriptome policy
    if transcriptome_selection is None and transcriptome_fasta:
        input_spec = WorkflowInputSpec(
            input_fasta=None,
            transcriptome_argument=transcriptome_fasta,
            default_transcriptomes=DEFAULT_TRANSCRIPTOME_SOURCES,
            design_only=False,
        )
        resolver = ReferencePolicyResolver(input_spec)
        transcriptome_selection = resolver.resolve_transcriptomes()

    if transcriptome_selection is None:
        transcriptome_selection = ReferenceSelection.disabled("no transcriptome configured")

    # Create a minimal workflow config for off-target analysis
    workflow_config = WorkflowConfig(
        output_dir=output_path,
        gene_query="offtarget_only",  # Placeholder name
        input_fasta=None,
        database=DatabaseType.ENSEMBL,  # Not used, but required
        resolved_policy=policy,
        nextflow_config=nextflow_config,
        transcriptome_indices=transcriptome_indices,
        screen_species=screen_species or ["human", "rat", "rhesus"],
        query_species=query_species,
        mirna_database=mirna_database,
        mirna_species=mirna_species,
        transcriptome_fasta=transcriptome_fasta,
        transcriptome_filter=transcriptome_filter,
        transcriptome_selection=transcriptome_selection,
        ortholog_mapping_file=ortholog_mapping_file,
        log_file=log_file,
        write_json_summary=False,  # Skip JSON summary for off-target-only
    )

    # Create workflow instance
    workflow = SiRNAWorkflow(workflow_config)

    # Run off-target analysis
    with Progress(console=console) as progress:
        task = progress.add_task("[cyan]Running off-target analysis...", total=None)

        offtarget_results = await workflow.run_nextflow_offtarget_analysis(
            candidates=candidates,
            input_fasta=candidates_fasta,
        )

        progress.remove_task(task)

    # The same evidence -> gates -> selection -> exports contract the design entry point gets (#100).
    # A no-op when the screen integrated hits: each of these skips a candidate that already holds the
    # verdict it would write.
    workflow._gate_without_screening_evidence(candidates)
    workflow._apply_selection(candidates)

    written_files: dict[str, str] = {}
    try:
        written_files = workflow.write_offtarget_only_exports(candidates)
    except Exception as e:  # Do not fail the run for reporting extras, exactly as step6 does not
        logger.warning(f"Failed to write off-target-only candidate exports: {e}")

    total_time = max(0.0, time.perf_counter() - start_time)

    # Compile results. The seven original keys are unchanged; selection_summary and written_files are
    # additive, because a caller that only reads offtarget_summary must keep working.
    final_results: dict[str, Any] = {
        "workflow_type": "offtarget_only",
        "input_candidates": str(input_candidates_fasta),
        "candidate_count": len(candidates),
        "output_dir": str(output_path),
        "processing_time": total_time,
        "offtarget_summary": offtarget_results,
        # Same reference record as the full workflow publishes: which references resolved, over which
        # species, and what did not.
        "reference_summary": workflow._summarize_screening_references(),
        # Why the shortlist is the size it is -- including "nothing qualified", which for pre-designed
        # guides is a legitimate answer this command must be able to state without failing.
        "selection_summary": workflow._selection_summary,
        "written_files": written_files,
    }

    console.print(f"\n✅ [bold green]Off-target analysis completed in {total_time:.2f}s[/bold green]")
    console.print(f"📊 Results saved to: [blue]{output_path}[/blue]")
    if written_files:
        console.print(
            "   - Candidate tables: sirnaforge/ (candidates_all.csv, candidates_pass.csv, "
            "candidates_qualified.csv, candidates_provisional.csv, manifest.json)"
        )

    return final_results
