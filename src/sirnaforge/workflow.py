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
import hashlib
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
from typing import Any, cast

import pandas as pd
from Bio.Seq import Seq
from pandera.typing import DataFrame
from rich.console import Console
from rich.progress import Progress

from sirnaforge import __version__
from sirnaforge.config import (
    DEFAULT_TRANSCRIPTOME_SOURCES,
    ReferenceChoice,
    ReferencePolicyResolver,
    ReferenceSelection,
    WorkflowInputSpec,
)
from sirnaforge.core.design import MiRNADesigner, SiRNADesigner
from sirnaforge.core.hit_classification import (
    ClassificationContext,
    HitClass,
    HitClassCounts,
    classify_hit,
)
from sirnaforge.core.off_target import OffTargetAnalysisManager
from sirnaforge.core.repeat_detection import (
    DEFAULT_REPEAT_TRANSCRIPT_FRACTION,
    RepeatDetector,
    normalize_guide_sequence,
)
from sirnaforge.core.scoring import (
    COMPOSITE_TERMS,
    SCORING_WEIGHT_SET_VERSION,
    ScoringError,
    compute_composite,
    conservation_sub_score,
    isoform_coverage_sub_score,
    off_target_sub_score,
)
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.data.base import DatabaseType, FastaUtils, TranscriptInfo
from sirnaforge.data.gene_search import GeneSearcher
from sirnaforge.data.orf_analysis import ORFAnalyzer
from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcript_annotation import EnsemblTranscriptModelClient
from sirnaforge.data.transcript_index import TranscriptGeneIndex
from sirnaforge.data.transcriptome_manager import TranscriptomeManager
from sirnaforge.models.schemas import ORFValidationSchema, SiRNACandidateSchema
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    DesignResult,
    FilterCriteria,
    OffTargetFilterCriteria,
    SiRNACandidate,
    build_candidate_row,
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
from sirnaforge.utils.cache_utils import resolve_cache_subdir, stable_cache_key
from sirnaforge.utils.control_candidates import DIRTY_CONTROL_LABEL, inject_dirty_controls
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.modification_patterns import apply_modifications_to_candidate
from sirnaforge.utils.resource_resolver import InputSource, resolve_input_source
from sirnaforge.utils.species import is_human_species
from sirnaforge.validation import ValidationConfig, ValidationMiddleware
from sirnaforge.workflow_variant import (
    VariantWorkflowConfig,
    normalize_variant_mode,
    parse_clinvar_filter_string,
    resolve_workflow_variants,
)
from sirnaforge.zfn.design import ZFNDesigner

logger = get_logger(__name__)
console = Console(record=True, force_terminal=False, legacy_windows=True)


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
        genome_indices_override: str | None = None,
        genome_species: list[str] | None = None,
        mirna_database: str = "mirgenedb",
        mirna_species: Sequence[str] | None = None,
        transcriptome_fasta: str | None = None,
        transcriptome_filter: str | None = None,
        transcriptome_selection: ReferenceSelection | None = None,
        validation_config: ValidationConfig | None = None,
        log_file: str | None = None,
        write_json_summary: bool = True,
        num_threads: int | None = None,
        input_source: InputSource | None = None,
        keep_nextflow_work: bool = False,
        variant_config: VariantWorkflowConfig | None = None,
        zfn_config: ZFNWorkflowConfig | None = None,
    ):
        """Initialize workflow configuration."""
        self.output_dir = Path(output_dir)
        self.input_source = input_source
        self.zfn_config = zfn_config

        resolved_input = input_source.local_path if input_source else (Path(input_fasta) if input_fasta else None)
        self.input_fasta = resolved_input
        # Preserve the user-supplied gene_query as the logical label even when using an input FASTA
        self.gene_query = gene_query
        self.database = database
        self.design_params = design_params or DesignParameters()
        # single source of truth: number of candidates selected everywhere
        self.top_n = self.design_params.top_n
        self.nextflow_config: dict[str, Any] = dict(nextflow_config) if nextflow_config else {}

        override_species: list[str] | None = None
        if genome_indices_override:
            self.nextflow_config["genome_indices"] = genome_indices_override
            override_species = self._extract_species_from_indices(genome_indices_override)

        # miRNA genome species: used for miRNA database lookups, not genomic DNA alignment
        # Track whether species were explicitly requested (vs defaulting)
        species_explicitly_provided = genome_species is not None or override_species is not None
        default_mirna_genomes = genome_species or ["human", "rat", "rhesus"]
        if override_species:
            default_mirna_genomes = override_species

        # Normalize all species names to canonical form for consistent comparisons
        normalized_genomes = [normalize_species_name(s) for s in default_mirna_genomes]
        self.mirna_genome_species: list[str] = list(dict.fromkeys(normalized_genomes))
        self.species_explicitly_requested = species_explicitly_provided
        self.mirna_database = mirna_database
        # Preserve explicit miRNA species order (values already normalized by CLI helpers)
        if mirna_species:
            filtered_species = [value for value in mirna_species if value]
            self.mirna_species = list(dict.fromkeys(filtered_species))
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

    @staticmethod
    def _extract_species_from_indices(indices: str) -> list[str]:
        """Derive species list from comma-separated species:/index_prefix entries."""
        species: list[str] = []
        for token in indices.split(","):
            entry = token.strip()
            if not entry:
                continue
            head = entry.split(":", 1)[0].strip() if ":" in entry else entry
            if head and head not in species:
                species.append(head)
        return species


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
        if config.design_params.design_mode == DesignMode.ZFN:
            self.zfn_designer = ZFNDesigner()
        elif config.design_params.design_mode == DesignMode.MIRNA:
            self.sirnaforgeer = MiRNADesigner(config.design_params)
        else:
            self.sirnaforgeer = SiRNADesigner(config.design_params)

        self.results: dict[str, Any] = {}
        self._nextflow_cache_info: dict[str, Any] | None = None
        self._annotation_summary: dict[str, Any] = {}
        self._gene_transcript_ids: set[str] = set()
        self._query_gene_ids: set[str] = set()
        self._query_gene_symbols: set[str] = set()
        self._protein_coding_transcript_ids: set[str] = set()
        self._protein_coding_transcript_count: int = 0
        self._transcript_index = TranscriptGeneIndex()
        self._species_explicitly_requested: bool = False
        self._representative_to_candidates: dict[str, list[SiRNACandidate]] = {}
        self._candidate_id_to_representative: dict[str, str] = {}
        # Single authoritative query species, set once (not re-inferred per call site). The
        # rest of this file is already human-centric (is_human_species, treated_as_human), so
        # the first screened species is the query species; an explicit query-species CLI
        # option is out of scope for this change.
        self._query_species: str = self.config.mirna_genome_species[0] if self.config.mirna_genome_species else "human"
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
            self._gene_transcript_ids = {
                self._normalize_transcript_id(t.transcript_id) for t in transcripts if t.transcript_id
            }
            self._query_gene_ids = {self._normalize_transcript_id(t.gene_id) for t in transcripts if t.gene_id}
            self._query_gene_symbols = {t.gene_name.strip().upper() for t in transcripts if t.gene_name}
            if not self._query_gene_symbols and self.config.gene_query:
                self._query_gene_symbols = {self.config.gene_query.strip().upper()}

            # Record protein-coding transcript set for isoform coverage scoring
            protein_coding_transcripts = {
                self._normalize_transcript_id(t.transcript_id)
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
            "design_parameters": design_parameters,
            "repeat_summary": self._repeat_summary,
            "offtarget_summary": offtarget_results,
            "reference_summary": {
                "transcriptome": self.config.transcriptome_selection.to_metadata(),
            },
        }

        # Optionally save workflow summary JSON (store in logs/)
        if self.config.write_json_summary:
            summary_file = self.config.output_dir / "logs" / "workflow_summary.json"
            with summary_file.open("w") as f:
                json.dump(final_results, f, indent=2, default=str)

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
                    batch_results, batch_guide_mapping = fut.result()
                    results.extend(batch_results)
                    # Merge guide-to-transcript mappings
                    for guide_seq, transcript_set in batch_guide_mapping.items():
                        guide_to_transcripts.setdefault(guide_seq, set()).update(transcript_set)
                    # Advance progress by the number of transcripts in this batch
                    batch = futures[fut]
                    progress.advance(task, len(batch))
                except Exception as e:
                    batch = futures[fut]
                    batch_transcript_ids = [t.transcript_id for t in batch]
                    logger.exception(f"Design failed for transcript batch {batch_transcript_ids}: {e}")
                    progress.advance(task, len(batch))

        self._store_guide_to_transcripts(guide_to_transcripts)

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
        all_candidates.sort(key=lambda x: x.composite_score, reverse=True)
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
            normalize_guide_sequence(guide): frozenset(self._normalize_transcript_id(tid) for tid in tids)
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

    def _process_transcript_batch(self, batch: list[TranscriptInfo]) -> tuple[list[DesignResult], dict[str, set[str]]]:
        """Process a batch of transcripts and return results plus guide-to-transcript mapping.

        Args:
            batch: List of transcripts to process

        Returns:
            Tuple of (design_results, guide_to_transcripts_mapping)
        """
        results: list[DesignResult] = []
        guide_to_transcripts: dict[str, set[str]] = {}

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
                continue

        return results, guide_to_transcripts

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
        fetching a second, redundant copy just to find the query species' FASTA.
        """
        distinct_guides = {normalize_guide_sequence(c.guide_sequence) for c in candidates}

        query_cdna_fasta = self._species_cdna_fasta.get(self._query_species)
        if not query_cdna_fasta:
            console.print("⚠️  Query species cDNA reference not available; skipping repeat detection")
            return {
                "status": "skipped",
                "reason": "reference_unavailable",
                "query_species": self._query_species,
                "repeat_flagged_count": 0,
                "threshold_fraction": DEFAULT_REPEAT_TRANSCRIPT_FRACTION,
            }

        console.print(
            f"🔍 Scanning {len(distinct_guides)} distinct guides against {self._query_species} cDNA (~47s)..."
        )
        detector = RepeatDetector(threshold_fraction=DEFAULT_REPEAT_TRANSCRIPT_FRACTION)
        scan_result = detector.scan(distinct_guides, query_cdna_fasta)

        # Stamp repeat verdicts on candidates
        observations = scan_result.observations
        for candidate in candidates:
            SiRNADesigner.stamp_repeat_verdict(candidate, observations)

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
        """Step 6: Generate comprehensive reports."""
        # No user-facing top-candidates FASTA or text/json summaries are produced anymore.
        # We only keep canonical CSV outputs (ALL + PASS) for candidates. Off-target analysis
        # prepares its own internal FASTA input under off_target/.

        if self.config.design_params.apply_modifications:
            self._apply_modifications_to_results(design_results)

        base = self.config.output_dir / "sirnaforge"
        all_csv = base / "candidates_all.csv"
        pass_csv = base / "candidates_pass.csv"
        pass_fasta = base / "candidates_pass.fasta"
        report_file = self.config.output_dir / "orf_reports" / "orf_validation.txt"
        is_mirna_mode = self.config.design_params.design_mode == DesignMode.MIRNA

        try:
            # Single shared row-builder (models/sirna.py) so this CSV and DesignResult.save_csv
            # can never drift on which columns they emit (issue #80 F2).
            rows: list[dict[str, Any]] = [build_candidate_row(candidate) for candidate in design_results.candidates]

            if rows:
                all_df = pd.DataFrame(rows)
            else:
                template_cols = list(SiRNACandidateSchema.to_schema().columns.keys())
                all_df = pd.DataFrame(columns=template_cols)

            for col in ("seed_7mer_hits", "seed_8mer_hits"):
                if col in all_df.columns:
                    all_df[col] = all_df[col].astype("Int64")

            if "passes_filters" not in all_df.columns:
                all_df["passes_filters"] = pd.Series(dtype="object")

            validated_all = SiRNACandidateSchema.validate(all_df)

            if not is_mirna_mode:
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
            pass_df = validated_all[validated_all["passes_filters"] == "PASS"].copy()

            validated_all.to_csv(all_csv, index=False)
            pass_df.to_csv(pass_csv, index=False)

            if pass_df.empty:
                pass_fasta.unlink(missing_ok=True)
            else:
                try:
                    self._write_pass_candidates_fasta(pass_df, pass_fasta)
                except Exception as e:
                    logger.warning(f"Failed to write PASS candidates FASTA: {e}")

        except Exception as e:  # Do not fail workflow for reporting extras
            logger.warning(f"Failed to write all/pass CSVs: {e}")

        try:
            variant_links_path = self.config.output_dir / "logs" / "candidate_variants.json"
            self._write_candidate_variant_links(design_results.candidates, variant_links_path)
        except Exception as e:
            logger.warning(f"Failed to write candidate variant links: {e}")

        try:
            manifest = self._build_fair_manifest(
                all_csv=all_csv,
                pass_csv=pass_csv,
                pass_fasta=pass_fasta,
                orf_report=report_file,
            )
            manifest_path = base / "manifest.json"
            with manifest_path.open("w") as mf:
                json.dump(manifest, mf, indent=2)
        except Exception as e:
            logger.warning(f"Failed to write FAIR manifest: {e}")

        console.print("📋 Generated comprehensive reports and FAIR metadata")
        console.print("   - ORF validation report: orf_reports/")
        console.print("   - siRNA candidate CSVs: sirnaforge/ (candidates_all.csv, candidates_pass.csv)")
        console.print("   - siRNA candidate FASTA: sirnaforge/ (candidates_pass.fasta)")

    def _write_pass_candidates_fasta(self, pass_df: pd.DataFrame, output_path: Path) -> None:
        """Write passing candidates to FASTA format with simple headers.

        Args:
            pass_df: DataFrame containing passing candidates
            output_path: Path to write the FASTA file
        """
        try:
            sequences: list[tuple[str, str]] = []
            for _, row in pass_df.iterrows():
                # Create simple header with candidate ID and score
                header = f"{row['id']} score={row['composite_score']:.1f}"
                sequence = str(row["guide_sequence"])
                sequences.append((header, sequence))

            # Use FastaUtils to write the sequences
            FastaUtils.save_sequences_fasta(sequences, output_path)
            logger.info(f"Saved {len(sequences)} passing candidates to FASTA: {output_path}")

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

    def _file_hash_sha256(self, path: Path) -> str:
        """Return SHA-256 hash of a file for integrity (non-security) tracking."""
        h = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()

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
        orf_report: Path,
    ) -> dict[str, Any]:
        """Create a manifest JSON describing generated outputs (checksums, sizes, counts)."""
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
                "sha256": self._file_hash_sha256(p),
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
        add_file("orf_validation_report", orf_report, "tsv")

        # Scoring metadata with weight set version and active terms
        weights_dict = self.config.design_params.scoring.model_dump(mode="json")

        return {
            "tool": "sirnaforge",
            "tool_version": __version__,
            "gene_query": self.config.gene_query,
            "run_timestamp": now,
            # Dumped wholesale: hand-listing fields silently dropped thresholds
            # (min_asymmetry_score, max_poly_runs, ...) from the run record.
            "design_parameters": self.config.design_params.model_dump(mode="json"),
            "scoring": {
                "weight_set_version": SCORING_WEIGHT_SET_VERSION,
                "weights": weights_dict,
                "active_terms": list(COMPOSITE_TERMS),
            },
            "files": files,
        }

    async def step5_offtarget_analysis(self, design_results: DesignResult) -> dict[str, Any]:
        """Step 5: Detect repeat elements, then run off-target analysis via the Nextflow pipeline.

        Repeat detection runs here rather than as its own workflow step so it can reuse the
        transcriptome reference this step already materializes for screening, instead of
        fetching it a second time just to locate the query species' FASTA.
        """
        candidates_for_offtarget = self._select_candidates_for_offtarget(design_results)

        if not candidates_for_offtarget:
            console.print("⚠️  No candidates available for off-target analysis")
            self._repeat_summary = {
                "status": "skipped",
                "reason": "no_candidates",
                "repeat_flagged_count": 0,
                "threshold_fraction": DEFAULT_REPEAT_TRANSCRIPT_FRACTION,
            }
            return {"status": "skipped", "reason": "no_candidates"}

        # Prepare input files
        input_fasta = await self._prepare_offtarget_input(candidates_for_offtarget)

        # Materialize transcriptome references once; repeat detection and Nextflow both reuse
        # this instead of each fetching/indexing their own copy.
        additional_params: dict[str, Any] = dict(self.config.nextflow_config)
        has_transcriptome = await self._configure_transcriptome_inputs(additional_params)
        self._repeat_summary = self._run_repeat_detection(candidates_for_offtarget)
        # Exclude repeat-flagged candidates from top_candidates now, so the exclusion holds even
        # if screening below never runs (user-disabled, Nextflow unavailable/failed).
        self._apply_post_screen_ranking(design_results)

        # If user disabled off-target checking via design parameters, skip entirely
        if not getattr(self.config.design_params, "check_off_targets", True):
            console.print("⚠️  Off-target analysis skipped by user request")
            return {"status": "skipped", "reason": "user_disabled"}

        # Try Nextflow pipeline first. We do NOT run the simplistic sequence-based fallback
        # (it produces low-value results) when Nextflow is unavailable. Instead mark as skipped
        # so downstream steps/users can see the explicit reason.
        try:
            offtarget_result = await self._run_nextflow_offtarget_analysis(
                candidates_for_offtarget, input_fasta, additional_params, has_transcriptome
            )
            self._apply_post_screen_ranking(design_results)
            return offtarget_result
        except Exception as e:
            console.print(f"⚠️  Nextflow execution failed: {e}")
            logger.exception("Nextflow pipeline execution error")
            return {"status": "skipped", "reason": "nextflow_failed", "error": str(e)}

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
        """Re-rank candidates by their post-screen composite score.

        Screening activates the off_target, isoform_coverage and conservation terms, so the
        design-time order no longer reflects the final ranking. Re-sorts design_results.candidates
        in place (step6_generate_reports writes CSVs in this order) and rebuilds top_candidates
        from the viable subset (excluding repeat-flagged and failing candidates).

        Called once right after repeat detection (so the exclusion holds even if screening
        never runs) and again after successful screening (so post-screen scores take effect).
        Idempotent: re-sorting/re-filtering an already-ranked list is harmless.
        """
        design_results.candidates.sort(key=lambda c: c.composite_score, reverse=True)
        rankable = [c for c in design_results.candidates if not c.repeat_flagged and self._passes_filters(c)]
        design_results.top_candidates = rankable[: self.config.top_n]
        repeat_excluded = sum(1 for c in design_results.candidates if c.repeat_flagged)
        logger.info(f"Re-ranked {len(rankable)} candidates after screening (excluded {repeat_excluded} repeat-flagged)")

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

        FastaUtils.save_sequences_fasta(sequences, input_fasta)

        # Store the mappings for fan-out during integration. Duplicate ids are rejected
        # above, so this id->id lookup is O(1) and never needs to compare candidates by value.
        self._representative_to_candidates = representative_to_candidates
        self._candidate_id_to_representative = candidate_id_to_representative

        console.print(
            f"📝 Prepared off-target input: {len(candidates)} candidates → "
            f"{len(sequences)} distinct sequences (deduplication ratio: {len(candidates) / max(1, len(sequences)):.2f}x)"
        )
        return input_fasta

    async def _prepare_transcriptome_database(
        self, transcriptome_ref: str, filter_spec: list[str] | None = None
    ) -> dict[str, Any] | None:
        """Prepare transcriptome database from user-provided reference.

        Args:
            transcriptome_ref: Can be:
                - Pre-configured source name (e.g., 'ensembl_human_cdna')
                - Local file path
                - HTTP(S)/FTP URL
            filter_spec: Optional list of filter names (e.g., ['protein_coding', 'canonical_only'])

        Returns:
            Dictionary with 'fasta' and 'index' paths, or None if preparation failed
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

                species = manager.SOURCES[transcriptome_ref].species or "transcriptome"
                enriched_result: dict[str, Any] = {"species": species}
                enriched_result.update(raw_result)
                return enriched_result

            # Otherwise treat as custom path/URL
            logger.info(f"Processing custom transcriptome reference: {transcriptome_ref}")
            if filter_spec and len(filter_spec) > 0:
                logger.warning("Filtering is not supported for custom transcriptome paths; ignoring filters")
            raw_custom = manager.get_custom_transcriptome(transcriptome_ref, build_index=True)
            if raw_custom is None:
                return None

            enriched_custom: dict[str, Any] = {"species": "transcriptome"}
            enriched_custom.update(raw_custom)
            return enriched_custom

        except Exception as e:
            logger.exception(f"Failed to prepare transcriptome database from {transcriptome_ref}")
            console.print(f"⚠️  Transcriptome preparation error: {e}")
            return None

    async def _materialize_transcriptome_reference(self, choice: ReferenceChoice) -> tuple[str, str] | None:
        """Prepare a transcriptome reference for Nextflow usage."""
        if not choice.value:
            return None

        console.print(f"📚 Transcriptome reference: {choice.value} ({choice.state.value})")

        # Parse filter specification from config
        from sirnaforge.data.transcriptome_filter import get_filter_spec  # noqa: PLC0415

        filter_spec: list[str] | None = None
        if self.config.transcriptome_filter:
            try:
                filter_spec = get_filter_spec(self.config.transcriptome_filter)
                if filter_spec:
                    console.print(f"🔍 Applying transcriptome filters: {', '.join(filter_spec)}")
            except ValueError as exc:
                logger.error(f"Invalid transcriptome filter specification: {exc}")
                console.print(f"⚠️  Invalid filter specification: {exc}")
                # Continue without filters rather than failing
                filter_spec = None

        try:
            transcriptome_result = await self._prepare_transcriptome_database(choice.value, filter_spec)
        except Exception as exc:  # pragma: no cover - defensive logging path
            logger.exception("Failed to prepare transcriptome database")
            console.print(f"⚠️  Transcriptome preparation failed: {exc}")
            return None

        if not transcriptome_result or not transcriptome_result.get("fasta"):
            console.print("⚠️  Failed to prepare transcriptome database, continuing without it")
            return None

        species_value = transcriptome_result.get("species")
        raw_species = species_value if isinstance(species_value, str) else "transcriptome"
        # Normalize species name to canonical form (e.g., 'hsa' -> 'human', 'mmu' -> 'mouse')
        transcriptome_species = normalize_species_name(raw_species)

        # Build a transcript->gene index for EVERY screened species (not just human) so
        # orthologs can be recognized across species. This fixes the central defect where
        # ortholog recognition only ever worked for human.
        if not self._transcript_index.for_species(transcriptome_species):
            self._transcript_index.build(transcriptome_species, Path(transcriptome_result["fasta"]))

        # Stash the resolved FASTA so repeat detection can reuse it instead of fetching it again.
        self._species_cdna_fasta.setdefault(transcriptome_species, Path(transcriptome_result["fasta"]))

        # Use pre-built index if available (host has bwa-mem2), otherwise pass FASTA path
        # Nextflow will build the index in Docker if needed
        transcriptome_path = transcriptome_result.get("index") or transcriptome_result["fasta"]
        transcriptome_index = str(transcriptome_path)

        if transcriptome_result.get("index"):
            console.print(
                f"✨ Transcriptome database prepared: {transcriptome_result['fasta'].name} "
                f"(index: {transcriptome_result['index'].name})"
            )
        else:
            console.print(
                f"✨ Transcriptome database prepared: {transcriptome_result['fasta'].name} (Nextflow will build index)"
            )
        return transcriptome_species, transcriptome_index

    def _resolve_active_genome_species(self, params: Mapping[str, Any]) -> list[str]:
        """Filter genome species down to those with available indices."""
        requested = [species.strip() for species in self.config.mirna_genome_species if species.strip()]
        available: set[str] = set()
        for key in ("genome_indices", "genome_fastas"):
            raw_value = params.get(key) or self.config.nextflow_config.get(key)
            available.update(self._parse_species_entries(raw_value))

        if available:
            filtered = [species for species in requested if species in available]
            for species in sorted(available):
                if species not in filtered:
                    filtered.append(species)
            return filtered
        return requested

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
        genome_species: Sequence[str],
        additional_params: Mapping[str, Any],
        pipeline_revision: str,
    ) -> dict[str, Any]:
        """Configure cached work and home directories for Nextflow runs.

        Args:
            nf_config: Nextflow configuration
            genome_species: Species for miRNA genome lookups (used in cache key)
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
            "genome_species": sorted(genome_species),
            "additional_params": self._normalize_param_dict(additional_params),
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

    def _load_offtarget_aggregates(self, results_dir: Path) -> dict[str, Any]:
        """Load aggregated Nextflow summary JSON files when available."""
        aggregated: dict[str, Any] = {}
        results_path = Path(results_dir)
        search_roots: list[Path] = []
        agg_dir = results_path / "aggregated"
        if agg_dir.exists():
            search_roots.append(agg_dir)
        search_roots.append(results_path)

        def _load_json(filename: str) -> dict[str, Any] | None:
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

        transcriptome_summary = _load_json("combined_summary.json")
        if transcriptome_summary:
            aggregated["transcriptome"] = transcriptome_summary

        mirna_summary = _load_json("combined_mirna_summary.json")
        if mirna_summary:
            aggregated["mirna"] = mirna_summary

        return aggregated

    async def _configure_transcriptome_inputs(self, additional_params: dict[str, Any]) -> bool:
        """Prepare transcriptome inputs for Nextflow runs."""
        selection = self.config.transcriptome_selection
        if not selection.enabled:
            console.print(f"ℹ️  Transcriptome off-target disabled ({selection.disabled_reason})")
            return False
        if not self.config.transcriptome_references:
            return False

        prepared_entries: list[str] = []
        prepared_species: list[str] = []
        for choice in selection.choices:
            materialized = await self._materialize_transcriptome_reference(choice)
            if not materialized:
                continue
            transcriptome_species, transcriptome_index = materialized
            if transcriptome_species not in self.config.mirna_genome_species:
                self.config.mirna_genome_species.append(transcriptome_species)
            prepared_entries.append(f"{transcriptome_species}:{transcriptome_index}")
            prepared_species.append(transcriptome_species)

        if not prepared_entries:
            return False

        existing_indices = additional_params.get("transcriptome_indices")
        merged_entries = [token.strip() for token in existing_indices.split(",")] if existing_indices else []
        merged_entries = [entry for entry in merged_entries if entry]
        for entry in prepared_entries:
            if entry not in merged_entries:
                merged_entries.append(entry)
        additional_params["transcriptome_indices"] = ",".join(merged_entries)
        additional_params["transcriptome_species"] = ",".join(dict.fromkeys(prepared_species))
        return True

    def _log_nextflow_targets(
        self,
        active_species: Sequence[str],
        has_transcriptome: bool,
        additional_params: Mapping[str, Any],
    ) -> None:
        """Emit console updates about genome and transcriptome targets."""
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
            additional_params: Pre-configured Nextflow parameters (transcriptome inputs already
                materialized), or None to configure them here.
            has_transcriptome: Paired with additional_params; ignored when that is None.
        """
        if additional_params is None:
            additional_params = dict(self.config.nextflow_config)
            has_transcriptome = await self._configure_transcriptome_inputs(additional_params)
        has_transcriptome = bool(has_transcriptome)
        active_species = self._resolve_active_genome_species(additional_params)
        has_transcriptome = has_transcriptome or bool(additional_params.get("transcriptome_indices"))
        self._log_nextflow_targets(active_species, has_transcriptome, additional_params)

        if not active_species and not has_transcriptome:
            console.print("ℹ️  No transcriptome indices configured; skipping Nextflow run")
            return await self._basic_offtarget_analysis(candidates)

        runner, _ = self._setup_nextflow_runner(active_species, additional_params)

        if not self._validate_nextflow_environment(runner):
            return {"status": "skipped", "reason": "nextflow_unavailable"}

        # Execute pipeline
        console.print("🚀 Running embedded Nextflow off-target analysis...")
        nf_output_dir = self.config.output_dir / "off_target" / "results"

        results = await runner.run_offtarget_analysis(
            input_file=input_fasta,
            output_dir=nf_output_dir,
            genome_species=active_species,
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
        genome_species: Sequence[str],
        additional_params: Mapping[str, Any],
    ) -> tuple[NextflowRunner, dict[str, Any]]:
        """Configure Nextflow runner with user settings and cached workdirs.

        Args:
            genome_species: Species for miRNA genome lookups
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
            genome_species=genome_species,
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
        run_status = "completed"
        workflow_warnings: list[str] = []

        tx_summary = aggregated_views.get("transcriptome") if aggregated_views else None
        if tx_summary:
            missing_species = cast(list[str], tx_summary.get("missing_species") or [])
            if missing_species:
                run_status = "partial"
                warning_msg = (
                    "⚠️  No transcriptome alignment files were generated for: "
                    f"{', '.join(missing_species)}. This usually means the BWA-MEM2 indexing stage ran out of memory. "
                    "Increase Nextflow --max_memory (32GB+ recommended for human transcriptomes) or pre-build indices."
                )
                console.print(warning_msg)
                workflow_warnings.append(warning_msg)

        # Integrate off-target results into candidates with filtering
        filter_criteria = getattr(self.config.design_params, "offtarget_filters", None) or OffTargetFilterCriteria()
        updated_candidates, stats = self._integrate_offtarget_results(candidates, parsed, filter_criteria)
        self._log_offtarget_statistics(stats, aggregated_views, output_dir)

        # Map parsed results for return structure
        mapped = {}
        for c in updated_candidates:
            qid = c.id
            entry = parsed.get("results", {}).get(qid)
            if entry:
                mapped[qid] = {
                    "off_target_count": entry.get("off_target_count", 0),
                    "off_target_score": entry.get("off_target_score", 0.0),
                    "hits": entry.get("hits", []),
                }
            else:
                mapped[qid] = {"off_target_count": 0, "off_target_score": 0.0, "hits": []}

        return {
            "status": run_status,
            "method": "embedded_nextflow",
            "output_dir": str(output_dir),
            "results": mapped,
            "execution_metadata": results,
            "filtering_stats": stats,
            "aggregated": aggregated_views,
            "warnings": workflow_warnings,
        }

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
                ("failed_mirna_seed", "❌ {} failed: miRNA perfect seed matches"),
                ("failed_high_risk_mirna", "❌ {} failed: high-risk miRNA hits"),
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

        missing_on_target = stats.get("missing_on_target_hit", 0)
        self._log_missing_on_target(missing_on_target)

        if aggregated_views:
            tx_summary = aggregated_views.get("transcriptome")
            if tx_summary:
                species_counts = cast(dict[str, int], tx_summary.get("hits_per_species", {}) or {})
                human_hits = tx_summary.get("human_hits", 0)
                other_hits = tx_summary.get("other_species_hits", 0)
                console.print(f"   🧾 Aggregated transcriptome hits — human: {human_hits}, other: {other_hits}")
                if species_counts:
                    formatted = ", ".join(f"{k}: {v}" for k, v in sorted(species_counts.items()))
                    console.print(f"      per species: {formatted}")
                missing_species = cast(list[str], tx_summary.get("missing_species") or [])
                if missing_species:
                    console.print(
                        "      ⚠️ Transcriptome alignment files were missing for: "
                        f"{', '.join(missing_species)} (likely insufficient memory during BWA indexing)."
                    )

                species_analyzed = cast(list[str], tx_summary.get("species_analyzed", []) or [])
                zero_hit_species = [species for species in species_analyzed if species_counts.get(species, 0) == 0]
                if zero_hit_species and not missing_species:
                    console.print(f"      ℹ️ No transcriptome hits detected for: {', '.join(zero_hit_species)}")

            mirna_summary = aggregated_views.get("mirna")
            if mirna_summary:
                human_hits = mirna_summary.get("human_hits", 0)
                other_hits = mirna_summary.get("other_species_hits", 0)
                console.print(f"   🌱 Aggregated miRNA hits — human: {human_hits}, other: {other_hits}")

        trace_file = Path(output_dir) / "pipeline_info" / "execution_trace.txt"
        if trace_file.exists():
            console.print(f"   📘 Nextflow execution trace: {trace_file}")

    @staticmethod
    def _log_missing_on_target(missing_on_target: int) -> None:
        """Warn when no candidates confirmed a 0-mismatch hit against their own source transcript."""
        if missing_on_target:
            console.print(
                f"   ⚠️ {missing_on_target} candidate(s) had no confirmed 0-mismatch hit against their own "
                "source transcript (on-target self-match not found in the transcriptome index)"
            )

    async def _basic_offtarget_analysis(self, candidates: list[SiRNACandidate]) -> dict[str, Any]:
        """Fallback basic off-target analysis."""
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

        console.print(f"📊 Basic off-target analysis completed for {len(candidates)} candidates")
        return {"status": "completed", "method": "basic", "results": results, "aggregated": {}}

    async def _parse_nextflow_results(self, output_dir: Path) -> dict[str, Any]:  # noqa: PLR0912
        """Parse results from Nextflow off-target analysis.

        Parses BOTH genome/transcriptome AND miRNA results from their respective
        output directories and combines them into a single results structure for
        candidate filtering.
        """
        results: dict[str, dict[str, Any]] = {}

        if not output_dir.exists():
            return {"status": "missing", "method": "nextflow", "output_dir": str(output_dir), "results": results}

        # Check for combined genome/transcriptome results in aggregated subdirectory
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

        def _ingest_tsv(path: Path) -> bool:
            if not path.exists() or path.stat().st_size == 0:
                return False
            found = False
            with path.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    _ingest_row(row)
                    found = True
            return found

        def _ingest_json(path: Path) -> bool:
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
                found = True
            return found

        genome_hits_found = _ingest_tsv(_aggregate_path("combined_offtargets.tsv"))
        if not genome_hits_found:
            genome_hits_found = _ingest_json(_aggregate_path("combined_offtargets.json"))

        mirna_hits_found = _ingest_tsv(_aggregate_path("combined_mirna_hits.tsv"))
        if not mirna_hits_found:
            mirna_hits_found = _ingest_json(_aggregate_path("combined_mirna_hits.json"))

        if not genome_hits_found or not mirna_hits_found:
            genome_files: list[Path] = []
            mirna_files: list[Path] = []

            if not genome_hits_found:
                genome_dir = output_dir / "genome"
                if genome_dir.exists():
                    genome_files = list(genome_dir.glob("*_analysis.tsv"))

            if not mirna_hits_found:
                mirna_dir = output_dir / "mirna"
                if mirna_dir.exists():
                    mirna_files = list(mirna_dir.glob("*_analysis.tsv"))

            files: list[Path] = []
            if not genome_hits_found:
                files.extend(genome_files)
            if not mirna_hits_found:
                files.extend(mirna_files)

            if not files and not genome_hits_found and not mirna_hits_found:
                # Last resort: scan for any TSV files
                files = list(output_dir.glob("**/*_offtargets.tsv"))

            for fpath in files:
                _ingest_tsv(Path(fpath))

            if not mirna_hits_found:
                mirna_tsv = output_dir / "mirna" / "mirna_analysis.tsv"
                if _ingest_tsv(mirna_tsv):
                    logger.info(f"Parsing miRNA analysis results from {mirna_tsv}")
                    mirna_hits_found = True

        return {"status": "completed", "method": "nextflow", "output_dir": str(output_dir), "results": results}

    def _check_offtarget_filters(
        self,
        transcriptome_0mm: int,
        transcriptome_1mm: int,
        transcriptome_2mm: int,
        mirna_0mm_seed: int,
        mirna_high_risk: int,
        total_hits: int,
        genuine_off_target_count: int,
        filter_criteria: OffTargetFilterCriteria,
    ) -> tuple[bool, SiRNACandidate.FilterStatus | None]:
        """Check if candidate fails off-target filters.

        Returns:
            Tuple of (should_fail, fail_status enum or None)
        """
        # Define filter checks with their thresholds and enum members
        checks: list[tuple[int | None, int, SiRNACandidate.FilterStatus]] = [
            (
                filter_criteria.max_transcriptome_hits_0mm,
                transcriptome_0mm,
                SiRNACandidate.FilterStatus.TRANSCRIPTOME_PERFECT_MATCH,
            ),
            (
                filter_criteria.max_transcriptome_hits_1mm,
                transcriptome_1mm,
                SiRNACandidate.FilterStatus.TRANSCRIPTOME_1MM,
            ),
            (
                filter_criteria.max_transcriptome_hits_2mm,
                transcriptome_2mm,
                SiRNACandidate.FilterStatus.TRANSCRIPTOME_2MM,
            ),
            (
                filter_criteria.max_mirna_perfect_seed,
                mirna_0mm_seed,
                SiRNACandidate.FilterStatus.MIRNA_PERFECT_SEED,
            ),
            (
                filter_criteria.max_total_offtarget_hits,
                total_hits,
                SiRNACandidate.FilterStatus.TOTAL_OFFTARGETS,
            ),
            (
                filter_criteria.max_off_target_count,
                genuine_off_target_count,
                SiRNACandidate.FilterStatus.EXCESS_OFF_TARGETS,
            ),
        ]

        # Check all threshold-based filters
        for threshold, value, status in checks:
            if threshold is not None and value > threshold:
                return True, status

        # Check high-risk miRNA (boolean flag)
        if filter_criteria.fail_on_high_risk_mirna and mirna_high_risk > 0:
            return True, SiRNACandidate.FilterStatus.HIGH_RISK_MIRNA

        return False, None

    @staticmethod
    def _normalize_transcript_id(transcript_id: str) -> str:
        """Strip Ensembl-style version suffixes (e.g. ``.9``) for identity comparisons."""
        return re.sub(r"\.\d+$", "", transcript_id.strip())

    def _integrate_offtarget_results(  # noqa: PLR0912, C901
        self,
        candidates: list[SiRNACandidate],
        offtarget_data: dict[str, Any],
        filter_criteria: OffTargetFilterCriteria | None = None,
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

        Returns:
            Tuple of (updated candidates, statistics dict with hit class decomposition)
        """
        if not offtarget_data or offtarget_data.get("status") != "completed":
            logger.warning("No completed off-target data available; candidates keep design-time scores")
            return candidates, {}

        if filter_criteria is None:
            filter_criteria = OffTargetFilterCriteria()

        results = offtarget_data.get("results", {})

        # Single authoritative query species, set once in __init__ (see comment there).
        query_species = self._query_species

        # Conservation is keyed on the screened set, not on whether species were typed on the
        # CLI: the term goes inactive exactly when the screened set is query-species-only.
        requested_species = frozenset(normalize_species_name(s) for s in self.config.mirna_genome_species)
        requested_non_query = requested_species - {query_species}
        n_requested_non_query = len(requested_non_query)

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
        )

        # Fan out deduplicated results: each representative's results apply to all candidates sharing that sequence
        representative_results: dict[str, dict[str, Any]] = {}
        for repr_id, entry in results.items():
            representative_results[repr_id] = entry

        # Pre-seeded with every requested species (zero-filled) so a species that was screened
        # but produced no hits is distinguishable from one never requested (absent key). Species
        # seen on a hit but not requested (unexpected) still get a bucket via setdefault below.
        per_species: dict[str, dict[str, int]] = {
            species: {
                "on_target": 0,
                "ortholog": 0,
                "repeat": 0,
                "off_target": 0,
                "symbol_lookup_missing": 0,
                "species_index_missing": 0,
            }
            for species in requested_species
        }

        stats: dict[str, Any] = {
            "candidates_analyzed": len(candidates),
            "candidates_with_offtargets": 0,
            "hit_classes": {"on_target": 0, "ortholog": 0, "repeat": 0, "off_target": 0},
            "query_gene_transcripts_recognised": len(self._gene_transcript_ids),
            "ortholog_symbol_lookup_misses": 0,
            "species_index_misses": 0,
            "per_species": per_species,
            "failed_perfect_match": 0,
            "failed_transcriptome_1mm": 0,
            "failed_transcriptome_2mm": 0,
            "failed_mirna_seed": 0,
            "failed_high_risk_mirna": 0,
            "failed_excess_off_targets": 0,
            "human_transcriptome_hits": 0,
            "other_transcriptome_hits": 0,
            "human_mirna_hits": 0,
            "other_mirna_hits": 0,
            # User intent vs the default species list, recorded (not used to gate scoring).
            "species_explicitly_requested": self._species_explicitly_requested,
        }

        for candidate in candidates:
            candidate_id = candidate.id

            # O(1) lookup: candidate id -> its representative's id -> that representative's
            # results. Falls back to the candidate's own id when it was never deduplicated.
            repr_id = self._candidate_id_to_representative.get(candidate_id, candidate_id)
            offtarget_entry = representative_results.get(repr_id)

            # Mark as screened even if no hits
            candidate.off_target_screened = True

            if not offtarget_entry or not offtarget_entry.get("hits"):
                # No hits: compute post-screen score with zero off-targets
                self._score_candidate_post_screen(candidate, HitClassCounts(), n_requested_non_query)
                continue

            stats["candidates_with_offtargets"] += 1

            # Classify each transcriptome hit and aggregate miRNA hits (unchanged)
            hit_counts = HitClassCounts()
            # transcriptome_totals/_human are stratified SUBSETS of transcriptome_off_target_total
            # / transcriptome_human_total below (nm>=3 hits count toward the totals but land in
            # no bucket), so the two families never contradict each other.
            transcriptome_totals = {0: 0, 1: 0, 2: 0}
            transcriptome_human = {0: 0, 1: 0, 2: 0}
            transcriptome_off_target_total = 0
            transcriptome_human_total = 0
            transcriptome_seed_0mm = 0
            mirna_total = 0
            mirna_human_total = 0
            mirna_0mm_seed_total = 0
            mirna_human_0mm_seed = 0
            mirna_1mm_seed = 0
            mirna_high_risk_total = 0
            mirna_high_risk_human = 0

            for hit in offtarget_entry.get("hits", []):
                nm = int(hit.get("nm", 0))
                seed_mismatches = int(hit.get("seed_mismatches", 0))
                offtarget_score = float(hit.get("offtarget_score", 0.0))
                species_label = hit.get("species")
                species_is_human = is_human_species(species_label)

                is_mirna = "mirna_id" in hit or "database" in hit

                if is_mirna:
                    mirna_total += 1
                    if species_is_human:
                        mirna_human_total += 1
                    if seed_mismatches == 0:
                        mirna_0mm_seed_total += 1
                        if species_is_human:
                            mirna_human_0mm_seed += 1
                        if offtarget_score < 5.0:
                            mirna_high_risk_total += 1
                            if species_is_human:
                                mirna_high_risk_human += 1
                    elif seed_mismatches == 1:
                        mirna_1mm_seed += 1
                else:
                    # Classify transcriptome hit using the four-way classifier
                    classification = classify_hit(hit, candidate.guide_sequence, classification_context)
                    # A blank/missing species label belongs to the query species (see classifier).
                    hit_species = normalize_species_name(species_label) if species_label else query_species
                    species_bucket = per_species.setdefault(
                        hit_species,
                        {
                            "on_target": 0,
                            "ortholog": 0,
                            "repeat": 0,
                            "off_target": 0,
                            "symbol_lookup_missing": 0,
                            "species_index_missing": 0,
                        },
                    )

                    # Aggregate by class
                    if classification.hit_class == HitClass.ON_TARGET:
                        hit_counts.on_target += 1
                        species_bucket["on_target"] += 1
                    elif classification.hit_class == HitClass.ORTHOLOG:
                        hit_counts.ortholog += 1
                        species_bucket["ortholog"] += 1
                        if classification.matched_symbol:
                            # Track which species had ortholog hits
                            hit_counts.ortholog_species = frozenset(hit_counts.ortholog_species | {hit_species})
                    elif classification.hit_class == HitClass.REPEAT:
                        hit_counts.repeat += 1
                        species_bucket["repeat"] += 1
                    elif classification.hit_class == HitClass.OFF_TARGET:
                        hit_counts.off_target += 1
                        species_bucket["off_target"] += 1

                    # Track shortfall counters
                    if classification.symbol_lookup_missing:
                        hit_counts.symbol_lookup_missing += 1
                        species_bucket["symbol_lookup_missing"] += 1
                    if classification.species_index_missing:
                        hit_counts.no_species_index += 1
                        species_bucket["species_index_missing"] += 1

                    # Only genuine off-targets feed the mismatch-stratified counters. Letting
                    # on-target isoform hits through here would fail every guide on a
                    # multi-isoform gene against max_transcriptome_hits_0mm, which is the
                    # near-total kill switch fixed in 0.5.2.
                    if classification.hit_class is not HitClass.OFF_TARGET:
                        continue

                    # Every genuine off-target hit counts toward the totals regardless of nm, so
                    # transcriptome_hits_total agrees with off_target_count even when nm>=3 hits
                    # occur (the default exhaustive search no longer caps hits at nm<=2).
                    transcriptome_off_target_total += 1
                    treated_as_human = species_is_human or not species_label
                    if treated_as_human:
                        transcriptome_human_total += 1

                    if nm == 0:
                        transcriptome_totals[0] += 1
                        if treated_as_human:
                            transcriptome_human[0] += 1
                    elif nm == 1:
                        transcriptome_totals[1] += 1
                        if treated_as_human:
                            transcriptome_human[1] += 1
                    elif nm == 2:
                        transcriptome_totals[2] += 1
                        if treated_as_human:
                            transcriptome_human[2] += 1

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
            candidate.off_target_count = hit_counts.off_target  # Redefined: genuine off-targets only
            candidate.ortholog_species = ",".join(sorted(hit_counts.ortholog_species))

            # Legacy fields (still needed for reporting and miRNA filters)
            candidate.transcriptome_hits_total = transcriptome_total_hits
            candidate.transcriptome_hits_0mm = transcriptome_totals[0]
            candidate.transcriptome_hits_1mm = transcriptome_totals[1]
            candidate.transcriptome_hits_2mm = transcriptome_totals[2]
            candidate.transcriptome_hits_seed_0mm = transcriptome_seed_0mm
            candidate.on_target_confirmed = hit_counts.on_target > 0
            candidate.mirna_hits_total = mirna_total
            candidate.mirna_hits_0mm_seed = mirna_0mm_seed_total
            candidate.mirna_hits_1mm_seed = mirna_1mm_seed
            candidate.mirna_hits_high_risk = mirna_high_risk_total
            candidate.off_target_penalty = offtarget_entry.get("off_target_score", 0.0)

            # Update global stats
            stats["hit_classes"]["on_target"] += hit_counts.on_target
            stats["hit_classes"]["ortholog"] += hit_counts.ortholog
            stats["hit_classes"]["repeat"] += hit_counts.repeat
            stats["hit_classes"]["off_target"] += hit_counts.off_target
            stats["ortholog_symbol_lookup_misses"] += hit_counts.symbol_lookup_missing
            stats["species_index_misses"] += hit_counts.no_species_index
            stats["human_transcriptome_hits"] += human_transcriptome_hits
            stats["other_transcriptome_hits"] += transcriptome_total_hits - human_transcriptome_hits
            stats["human_mirna_hits"] += mirna_human_total
            stats["other_mirna_hits"] += mirna_total - mirna_human_total

            # Score candidate with post-screen terms
            self._score_candidate_post_screen(candidate, hit_counts, n_requested_non_query)

            # Apply filtering criteria
            human_total_hits_for_filters = human_transcriptome_hits + mirna_human_total
            should_fail, fail_status = self._check_offtarget_filters(
                transcriptome_human[0],
                transcriptome_human[1],
                transcriptome_human[2],
                mirna_human_0mm_seed,
                mirna_high_risk_human,
                human_total_hits_for_filters,
                hit_counts.off_target,
                filter_criteria,
            )

            if should_fail and fail_status:
                candidate.passes_filters = fail_status
                logger.info(f"Candidate {candidate_id} failed off-target filter: {fail_status.value}")

                # Update stat counters
                if fail_status == SiRNACandidate.FilterStatus.TRANSCRIPTOME_PERFECT_MATCH:
                    stats["failed_perfect_match"] += 1
                elif fail_status == SiRNACandidate.FilterStatus.TRANSCRIPTOME_1MM:
                    stats["failed_transcriptome_1mm"] += 1
                elif fail_status == SiRNACandidate.FilterStatus.TRANSCRIPTOME_2MM:
                    stats["failed_transcriptome_2mm"] += 1
                elif fail_status == SiRNACandidate.FilterStatus.MIRNA_PERFECT_SEED:
                    stats["failed_mirna_seed"] += 1
                elif fail_status == SiRNACandidate.FilterStatus.HIGH_RISK_MIRNA:
                    stats["failed_high_risk_mirna"] += 1
                elif fail_status == SiRNACandidate.FilterStatus.EXCESS_OFF_TARGETS:
                    stats["failed_excess_off_targets"] += 1

        # Re-ranking (excluding repeat-flagged candidates) happens in step5_offtarget_analysis,
        # where design_results is in scope to receive the reordered candidates/top_candidates.
        return candidates, stats

    def _score_candidate_post_screen(
        self, candidate: SiRNACandidate, hit_counts: HitClassCounts, n_requested_non_query: int
    ) -> None:
        """Compute post-screen composite score with the full term set.

        Args:
            candidate: Candidate to score
            hit_counts: Aggregated hit class counts
            n_requested_non_query: Number of non-query species the user requested
        """
        # Build features from design-time component scores (reuse existing sub-scores)
        features: dict[str, float] = {}
        cs = candidate.component_scores or {}

        # Design-time terms: reuse from component_scores, dropping any NaN
        for term in ("asymmetry", "gc_content", "accessibility", "empirical"):
            value = cs.get(term)
            if value is not None and not math.isnan(value):
                features[term] = value

        # The sub-score helpers raise on a numerator exceeding its denominator. Both numerators are
        # subsets of their denominators today, so neither is reachable -- but the guarantee rests on
        # invariants several call sites away, and the cost of one being broken later must not be an
        # aborted run after screening has already been paid for.
        try:
            features["off_target"] = off_target_sub_score(hit_counts.off_target)

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

            conservation = conservation_sub_score(len(hit_counts.ortholog_species), n_requested_non_query)
            if conservation is not None:
                features["conservation"] = conservation
                candidate.conservation_score = conservation

            result = compute_composite(features, self.config.design_params.scoring)
            candidate.composite_score = result.score
            candidate.weight_set_version = result.weight_set_version
            candidate.scored_after_screening = True
            # Write per-term contributions
            candidate.score_asymmetry = result.contributions.get("asymmetry")
            candidate.score_gc_content = result.contributions.get("gc_content")
            candidate.score_accessibility = result.contributions.get("accessibility")
            candidate.score_empirical = result.contributions.get("empirical")
            candidate.score_off_target = result.contributions.get("off_target")
            candidate.score_isoform_coverage = result.contributions.get("isoform_coverage")
            candidate.score_conservation = result.contributions.get("conservation")
        except (ScoringError, ValueError) as exc:
            logger.warning(f"Post-screen scoring failed for {candidate.id}: {exc}. Keeping design-time score.")
            candidate.scored_after_screening = False

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
        """Summarize siRNA design results."""
        base = design_results.get_summary()
        total = design_results.total_candidates
        passed = design_results.filtered_candidates
        failed = max(0, total - passed)
        repeat_excluded = sum(1 for c in design_results.candidates if c.repeat_flagged)
        base.update(
            {
                "pass_count": passed,
                "fail_count": failed,
                "repeat_excluded_count": repeat_excluded,
                "repeat_threshold_fraction": DEFAULT_REPEAT_TRANSCRIPT_FRACTION,
                "top_n_requested": self.config.top_n,
                "dirty_controls_added": getattr(self, "_dirty_controls_added", 0),
                "threads_used": self.config.num_threads,
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
        normalized["chromosomes"] = [token.strip() for token in chromosomes.split(",") if token.strip()]

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
    top_n_candidates: int = 20,
    genome_species: list[str] | None = None,
    genome_indices_override: str | None = None,
    mirna_database: str = "mirgenedb",
    mirna_species: Sequence[str] | None = None,
    transcriptome_fasta: str | None = None,
    transcriptome_filter: str | None = None,
    transcriptome_selection: ReferenceSelection | None = None,
    gc_min: float = 30.0,
    gc_max: float = 52.0,
    sirna_length: int = 21,
    modification_pattern: str = "standard_2ome",
    overhang: str = "dTdT",
    zfn_design_params: ZFNDesignParameters | None = None,
    zfn_annotation: GenomicAnnotationConfig | None = None,
    check_off_targets: bool = True,
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
    allow_transcriptome_with_input_fasta: bool = True,
    default_transcriptome_sources: Sequence[str] = DEFAULT_TRANSCRIPTOME_SOURCES,
    keep_nextflow_work: bool = False,
    nextflow_docker_image: str | None = None,
    max_hits: int | None = None,
) -> dict[str, Any]:
    """Run complete siRNA design workflow.

    Args:
        gene_query: Gene name or ID to search for
        output_dir: Directory for output files
        input_fasta: Local path or remote URI to an input FASTA file
        database: Database to search (ensembl, refseq, gencode)
        design_mode: Design mode (sirna, mirna, or zfn)
        top_n_candidates: Number of top candidates to generate
        genome_species: Species genomes for off-target analysis
        genome_indices_override: Comma-separated species:/index_prefix overrides for off-target analysis
        mirna_database: miRNA reference database identifier
        mirna_species: miRNA reference species identifiers
        transcriptome_fasta: Path or URL to transcriptome FASTA for off-target analysis
        transcriptome_filter: Comma-separated filter names (protein_coding, canonical_only)
        transcriptome_selection: Pre-resolved transcriptome selection metadata
        gc_min: Minimum GC content percentage
        gc_max: Maximum GC content percentage
        sirna_length: siRNA length in nucleotides
        modification_pattern: Chemical modification pattern
        overhang: Overhang sequence (dTdT for DNA, UU for RNA)
        zfn_design_params: Optional ZFN design parameters for ZFN mode workflow
        zfn_annotation: Optional genomic annotation config for ZFN off-target classification
        check_off_targets: Perform off-target analysis stage (default: True)
        variant_ids: List of variant identifiers (rsID, chr:pos:ref:alt, or HGVS) to target or avoid
        variant_vcf_file: Path to VCF file containing variants to target or avoid
        variant_mode: How to handle variants (avoid/target/both) - default is avoid
        variant_min_af: Minimum allele frequency threshold for variant filtering (default: 0.01)
        variant_clinvar_filters: Comma-separated ClinVar significance levels to include (default: Pathogenic,Likely pathogenic)
        variant_assembly: Reference genome assembly for variants (only GRCh38 supported)
        log_file: Path to centralized log file
        write_json_summary: Write logs/workflow_summary.json
        num_threads: Optional override for design parallelism
        allow_transcriptome_with_input_fasta: Allow default transcriptome analysis when using input FASTA
        default_transcriptome_sources: Ordered list of transcriptome identifiers evaluated by default
        keep_nextflow_work: Keep Nextflow work directory symlink in output
        nextflow_docker_image: Override Docker image used by the embedded Nextflow pipeline
        max_hits: Override the pipeline's per-candidate off-target hit cap (None keeps the pipeline's
            exhaustive default; set a lower value, e.g. 10000, to speed up large gene-family searches)

    Returns:
        Dictionary with complete workflow results
    """
    # Parse design mode
    try:
        mode_enum = DesignMode(design_mode.lower())
    except ValueError:
        mode_enum = DesignMode.SIRNA

    # Configure filter criteria
    filter_criteria = FilterCriteria(
        gc_min=gc_min,
        gc_max=gc_max,
    )

    # Configure workflow with modification parameters
    design_params = DesignParameters(
        design_mode=mode_enum,
        top_n=top_n_candidates,
        sirna_length=sirna_length,
        filters=filter_criteria,
        check_off_targets=check_off_targets,
        apply_modifications=modification_pattern.lower() != "none",
        modification_pattern=modification_pattern,
        default_overhang=overhang,
    )
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
            design_only=False,
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
        genome_indices_override=genome_indices_override,
        genome_species=genome_species or ["human", "rat", "rhesus"],
        mirna_database=mirna_database,
        mirna_species=mirna_species,
        transcriptome_fasta=transcriptome_fasta,
        transcriptome_filter=transcriptome_filter,
        transcriptome_selection=transcriptome_selection,
        log_file=log_file,
        write_json_summary=write_json_summary,
        num_threads=num_threads,
        input_source=resolved_input,
        keep_nextflow_work=keep_nextflow_work,
        variant_config=variant_config_obj,
        nextflow_config=nextflow_config_overrides,
        zfn_config=zfn_workflow_config,
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
    genome_species: list[str] | None = None,
    genome_indices_override: str | None = None,
    mirna_database: str = "mirgenedb",
    mirna_species: Sequence[str] | None = None,
    transcriptome_fasta: str | None = None,
    transcriptome_filter: str | None = None,
    transcriptome_selection: ReferenceSelection | None = None,
    log_file: str | None = None,
    nextflow_docker_image: str | None = None,
) -> dict[str, Any]:
    """Run off-target-only workflow for pre-designed siRNA candidates.

    This is a simplified workflow that only runs the off-target analysis stage
    without transcript retrieval, ORF validation, or siRNA design. It accepts
    pre-designed 21-nt siRNA guide sequences and runs comprehensive off-target
    analysis using the embedded Nextflow pipeline.

    Args:
        input_candidates_fasta: Path to FASTA file with 21-nt siRNA guide sequences
        output_dir: Directory for output files
        genome_species: Species genomes for off-target analysis
        genome_indices_override: Comma-separated species:/index_prefix overrides
        mirna_database: miRNA reference database identifier
        mirna_species: miRNA reference species identifiers
        transcriptome_fasta: Path or URL to transcriptome FASTA for off-target analysis
        transcriptome_filter: Comma-separated filter names (protein_coding, canonical_only)
        transcriptome_selection: Pre-resolved transcriptome selection metadata
        log_file: Path to centralized log file
        nextflow_docker_image: Override Docker image used by the embedded Nextflow pipeline

    Returns:
        Dictionary with off-target analysis results
    """
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
                composite_score=0.0,
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
            composite_score=0.0,  # Not computed for pre-designed guides
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
    if genome_indices_override:
        nextflow_config["genome_indices"] = genome_indices_override
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
        design_params=DesignParameters(),  # Minimal params
        nextflow_config=nextflow_config,
        genome_indices_override=genome_indices_override,
        genome_species=genome_species or ["human", "rat", "rhesus"],
        mirna_database=mirna_database,
        mirna_species=mirna_species,
        transcriptome_fasta=transcriptome_fasta,
        transcriptome_filter=transcriptome_filter,
        transcriptome_selection=transcriptome_selection,
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

    total_time = max(0.0, time.perf_counter() - start_time)

    # Compile results
    final_results: dict[str, Any] = {
        "workflow_type": "offtarget_only",
        "input_candidates": str(input_candidates_fasta),
        "candidate_count": len(candidates),
        "output_dir": str(output_path),
        "processing_time": total_time,
        "offtarget_summary": offtarget_results,
    }

    console.print(f"\n✅ [bold green]Off-target analysis completed in {total_time:.2f}s[/bold green]")
    console.print(f"📊 Results saved to: [blue]{output_path}[/blue]")

    return final_results
