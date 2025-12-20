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

from sirnaforge.config import (
    DEFAULT_TRANSCRIPTOME_SOURCES,
    ReferenceChoice,
    ReferencePolicyResolver,
    ReferenceSelection,
    WorkflowInputSpec,
)
from sirnaforge.core.design import MiRNADesigner, SiRNADesigner
from sirnaforge.core.off_target import OffTargetAnalysisManager
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.data.base import DatabaseType, FastaUtils, TranscriptInfo
from sirnaforge.data.gene_search import GeneSearcher
from sirnaforge.data.orf_analysis import ORFAnalyzer
from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcriptome_manager import TranscriptomeManager
from sirnaforge.models.schemas import ORFValidationSchema, SiRNACandidateSchema
from sirnaforge.models.sirna import (
    DesignMode,
    DesignParameters,
    DesignResult,
    FilterCriteria,
    OffTargetFilterCriteria,
    SiRNACandidate,
)
from sirnaforge.models.sirna import SiRNACandidate as _ModelSiRNACandidate
from sirnaforge.models.variant import VariantMode, VariantRecord
from sirnaforge.pipeline import NextflowConfig, NextflowRunner
from sirnaforge.utils.cache_utils import resolve_cache_subdir, stable_cache_key
from sirnaforge.utils.control_candidates import DIRTY_CONTROL_LABEL, inject_dirty_controls
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.modification_patterns import apply_modifications_to_candidate, get_modification_summary
from sirnaforge.utils.resource_resolver import InputSource, resolve_input_source
from sirnaforge.utils.species import is_human_species
from sirnaforge.validation import ValidationConfig, ValidationMiddleware
from sirnaforge.workflow_variant import VariantWorkflowConfig, parse_clinvar_filter_string, resolve_workflow_variants

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
    ):
        """Initialize workflow configuration."""
        self.output_dir = Path(output_dir)
        self.input_source = input_source

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
        default_mirna_genomes = genome_species or ["human", "rat", "rhesus"]
        if override_species:
            default_mirna_genomes = override_species

        # Normalize all species names to canonical form for consistent comparisons
        normalized_genomes = [normalize_species_name(s) for s in default_mirna_genomes]
        self.mirna_genome_species: list[str] = list(dict.fromkeys(normalized_genomes))
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
        (self.output_dir / "transcripts").mkdir(exist_ok=True)
        (self.output_dir / "orf_reports").mkdir(exist_ok=True)
        (self.output_dir / "sirnaforge").mkdir(exist_ok=True)
        (self.output_dir / "off_target").mkdir(exist_ok=True)
        (self.output_dir / "logs").mkdir(exist_ok=True)

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


class SiRNAWorkflow:
    """Main workflow orchestrator for siRNA design pipeline."""

    def __init__(self, config: WorkflowConfig):
        """Initialize the siRNA workflow orchestrator."""
        self.config = config
        self.gene_searcher = GeneSearcher()
        self.orf_analyzer = ORFAnalyzer()
        self.validation = ValidationMiddleware(config.validation_config)

        # Select designer based on design mode
        self.sirnaforgeer: SiRNADesigner
        if config.design_params.design_mode == DesignMode.MIRNA:
            self.sirnaforgeer = MiRNADesigner(config.design_params)
        else:
            self.sirnaforgeer = SiRNADesigner(config.design_params)

        self.results: dict[str, Any] = {}
        self._nextflow_cache_info: dict[str, Any] | None = None

    async def run_complete_workflow(self) -> dict[str, Any]:
        """Run the complete siRNA design workflow."""
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

            # Step 1.5: Variant Resolution (if configured)
            # Save resolved variants on the workflow instance for later use
            if self.config.variant_config and self.config.variant_config.has_variants:
                progress.update(main_task, description="[cyan]Resolving variants...")
                self.resolved_variants = await self.step1_5_resolve_variants(progress)
                progress.advance(main_task)
            else:
                # Skip variant resolution step
                self.resolved_variants = []
                progress.advance(main_task)

            # Step 2: ORF Validation
            progress.update(main_task, description="[cyan]Validating ORFs...")
            orf_results = await self.step2_validate_orfs(transcripts, progress)
            progress.advance(main_task)

            # Step 3: siRNAforge
            progress.update(main_task, description="[cyan]Designing siRNAs...")
            design_results = await self.step3_design_sirnas(transcripts, progress)
            progress.advance(main_task)

            # Step 4: Off-target Analysis (must run before reports to update candidate data)
            progress.update(main_task, description="[cyan]Running off-target analysis...")
            offtarget_results = await self.step5_offtarget_analysis(design_results)
            progress.advance(main_task)

            # Step 5: Generate Reports (after off-target analysis completes)
            progress.update(main_task, description="[cyan]Generating reports...")
            await self.step4_generate_reports(design_results)
            progress.advance(main_task)

        total_time = max(0.0, time.perf_counter() - start_time)

        # Compile final results
        # Serialize authoritative design parameters into the workflow summary.
        dp = self.config.design_params
        design_parameters: dict[str, Any] = {
            "top_n": dp.top_n,
            "sirna_length": dp.sirna_length,
            "filters": {
                "gc_min": dp.filters.gc_min,
                "gc_max": dp.filters.gc_max,
                "max_poly_runs": dp.filters.max_poly_runs,
                "max_paired_fraction": dp.filters.max_paired_fraction,
                "min_asymmetry_score": dp.filters.min_asymmetry_score,
            },
            "scoring": {
                "asymmetry": dp.scoring.asymmetry,
                "gc_content": dp.scoring.gc_content,
                "accessibility": dp.scoring.accessibility,
                "off_target": dp.scoring.off_target,
                "empirical": dp.scoring.empirical,
            },
            "avoid_snps": dp.avoid_snps,
            "check_off_targets": dp.check_off_targets,
            "predict_structure": dp.predict_structure,
            "snp_file": dp.snp_file,
            "genome_index": dp.genome_index,
        }

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
            "orf_summary": self._summarize_orf_results(orf_results),
            "design_summary": self._summarize_design_results(design_results),
            "design_parameters": design_parameters,
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

        return protein_transcripts

    async def step1_5_resolve_variants(self, progress: Progress) -> list[VariantRecord]:
        """Step 1.5: Resolve variants for targeting or avoidance.

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
        report_file = self.config.output_dir / "orf_reports" / f"{self.config.gene_query}_orf_validation.txt"
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
            design_result = self.sirnaforgeer.design_from_file(str(temp_fasta))
            added_controls = inject_dirty_controls(design_result)
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
        if added_controls:
            console.print(f"🧪 Added {len(added_controls)} {DIRTY_CONTROL_LABEL} candidates for signal verification")

        _ = self.validation.validate_design_results(combined)

        temp_fasta.unlink(missing_ok=True)
        console.print(f"🎯 Generated {len(combined.candidates)} siRNA candidates (threads={self.config.num_threads})")
        console.print(f"   Top {len(combined.top_candidates)} candidates selected for further analysis")
        return combined

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

    async def step4_generate_reports(self, design_results: DesignResult) -> None:  # noqa: C901, PLR0912
        """Step 4: Generate comprehensive reports."""
        # No user-facing top-candidates FASTA or text/json summaries are produced anymore.
        # We only keep canonical CSV outputs (ALL + PASS) for candidates. Off-target analysis
        # prepares its own internal FASTA input under off_target/.

        if self.config.design_params.apply_modifications:
            self._apply_modifications_to_results(design_results)

        base = self.config.output_dir / "sirnaforge"
        all_csv = base / f"{self.config.gene_query}_all.csv"
        pass_csv = base / f"{self.config.gene_query}_pass.csv"
        pass_fasta = base / f"{self.config.gene_query}_pass.fasta"
        report_file = self.config.output_dir / "orf_reports" / f"{self.config.gene_query}_orf_validation.txt"
        is_mirna_mode = self.config.design_params.design_mode == DesignMode.MIRNA

        try:
            rows: list[dict[str, Any]] = []
            for candidate in design_results.candidates:
                cs = getattr(candidate, "component_scores", {}) or {}
                mod_summary = get_modification_summary(candidate) if candidate.guide_metadata else {}

                pass_state = candidate.passes_filters
                if isinstance(pass_state, _ModelSiRNACandidate.FilterStatus):
                    normalized_pass_field: Any = pass_state.value
                else:
                    normalized_pass_field = pass_state

                def _maybe_attr(obj: Any, name: str, *, default: Any = None) -> Any:
                    return getattr(obj, name, default)

                rows.append(
                    {
                        "id": candidate.id,
                        "transcript_id": candidate.transcript_id,
                        "position": candidate.position,
                        "guide_sequence": candidate.guide_sequence,
                        "passenger_sequence": candidate.passenger_sequence,
                        "gc_content": candidate.gc_content,
                        "asymmetry_score": candidate.asymmetry_score,
                        "paired_fraction": candidate.paired_fraction,
                        "structure": getattr(candidate, "structure", None),
                        "mfe": getattr(candidate, "mfe", None),
                        "duplex_stability_dg": candidate.duplex_stability,
                        "duplex_stability_score": cs.get("duplex_stability_score"),
                        "dg_5p": cs.get("dg_5p"),
                        "dg_3p": cs.get("dg_3p"),
                        "delta_dg_end": cs.get("delta_dg_end"),
                        "melting_temp_c": cs.get("melting_temp_c"),
                        "off_target_count": candidate.off_target_count,
                        "off_target_penalty": candidate.off_target_penalty,
                        "transcriptome_hits_total": _maybe_attr(candidate, "transcriptome_hits_total", default=0),
                        "transcriptome_hits_0mm": _maybe_attr(candidate, "transcriptome_hits_0mm", default=0),
                        "transcriptome_hits_1mm": _maybe_attr(candidate, "transcriptome_hits_1mm", default=0),
                        "transcriptome_hits_2mm": _maybe_attr(candidate, "transcriptome_hits_2mm", default=0),
                        "transcriptome_hits_seed_0mm": _maybe_attr(candidate, "transcriptome_hits_seed_0mm", default=0),
                        "mirna_hits_total": _maybe_attr(candidate, "mirna_hits_total", default=0),
                        "mirna_hits_0mm_seed": _maybe_attr(candidate, "mirna_hits_0mm_seed", default=0),
                        "mirna_hits_1mm_seed": _maybe_attr(candidate, "mirna_hits_1mm_seed", default=0),
                        "mirna_hits_high_risk": _maybe_attr(candidate, "mirna_hits_high_risk", default=0),
                        "guide_pos1_base": _maybe_attr(candidate, "guide_pos1_base"),
                        "pos1_pairing_state": _maybe_attr(candidate, "pos1_pairing_state"),
                        "seed_class": _maybe_attr(candidate, "seed_class"),
                        "supp_13_16_score": _maybe_attr(candidate, "supp_13_16_score"),
                        "seed_7mer_hits": _maybe_attr(candidate, "seed_7mer_hits"),
                        "seed_8mer_hits": _maybe_attr(candidate, "seed_8mer_hits"),
                        "seed_hits_weighted": _maybe_attr(candidate, "seed_hits_weighted"),
                        "off_target_seed_risk_class": _maybe_attr(candidate, "off_target_seed_risk_class"),
                        "transcript_hit_count": candidate.transcript_hit_count,
                        "transcript_hit_fraction": candidate.transcript_hit_fraction,
                        "composite_score": candidate.composite_score,
                        "passes_filters": normalized_pass_field,
                        "guide_overhang": mod_summary.get("guide_overhang", ""),
                        "guide_modifications": mod_summary.get("guide_modifications", ""),
                        "passenger_overhang": mod_summary.get("passenger_overhang", ""),
                        "passenger_modifications": mod_summary.get("passenger_modifications", ""),
                    }
                )

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
        console.print("   - siRNA candidate CSVs: sirnaforge/ (<gene>_all.csv, <gene>_pass.csv)")
        console.print("   - siRNA candidate FASTA: sirnaforge/ (<gene>_pass.fasta)")

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

        return {
            "tool": "sirnaforge",
            "gene_query": self.config.gene_query,
            "run_timestamp": now,
            "design_parameters": {
                "top_n": self.config.top_n,
                "sirna_length": self.config.design_params.sirna_length,
                "gc_min": self.config.design_params.filters.gc_min,
                "gc_max": self.config.design_params.filters.gc_max,
            },
            "files": files,
        }

    async def step5_offtarget_analysis(self, design_results: DesignResult) -> dict[str, Any]:
        """Step 5: Run off-target analysis using embedded Nextflow pipeline."""
        candidates_for_offtarget = self._select_candidates_for_offtarget(design_results)

        if not candidates_for_offtarget:
            console.print("⚠️  No candidates available for off-target analysis")
            return {"status": "skipped", "reason": "no_candidates"}

        # Prepare input files
        input_fasta = await self._prepare_offtarget_input(candidates_for_offtarget)

        # If user disabled off-target checking via design parameters, skip entirely
        if not getattr(self.config.design_params, "check_off_targets", True):
            console.print("⚠️  Off-target analysis skipped by user request")
            return {"status": "skipped", "reason": "user_disabled"}

        # Try Nextflow pipeline first. We do NOT run the simplistic sequence-based fallback
        # (it produces low-value results) when Nextflow is unavailable. Instead mark as skipped
        # so downstream steps/users can see the explicit reason.
        try:
            return await self._run_nextflow_offtarget_analysis(candidates_for_offtarget, input_fasta)
        except Exception as e:
            console.print(f"⚠️  Nextflow execution failed: {e}")
            logger.exception("Nextflow pipeline execution error")
            return {"status": "skipped", "reason": "nextflow_failed", "error": str(e)}

    def _select_candidates_for_offtarget(self, design_results: DesignResult) -> list[SiRNACandidate]:
        """Return top-N candidates plus any dirty controls for off-target analysis.

        Off-target pipelines expect to see at least one "fails on purpose"
        control so we always tack the :func:`inject_dirty_controls` output onto
        the regular top-N selection. The controls remain true siRNA designs that
        merely failed QC, which makes them perfect sentinels for verifying that
        downstream aligners, Nextflow modules, and reports are actually running.
        """
        selected: list[SiRNACandidate] = list(design_results.top_candidates[: self.config.top_n])

        dirty_controls = [c for c in design_results.top_candidates if self._is_dirty_control_candidate(c)]
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

    async def _prepare_offtarget_input(self, candidates: list[SiRNACandidate]) -> Path:
        """Prepare FASTA input file for off-target analysis.

        ``candidates`` must already include any dirty controls (handled by
        :meth:`_select_candidates_for_offtarget`). We simply persist the ID and
        guide sequence for each entry so the generated
        ``off_target/input_candidates.fasta`` mirrors whatever the off-target
        runner receives.
        """
        input_fasta = self.config.output_dir / "off_target" / "input_candidates.fasta"
        sequences = [(f"{c.id}", c.guide_sequence) for c in candidates]
        FastaUtils.save_sequences_fasta(sequences, input_fasta)
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
    ) -> dict[str, Any]:
        """Run Nextflow-based off-target analysis."""
        additional_params: dict[str, Any] = dict(self.config.nextflow_config)
        has_transcriptome = await self._configure_transcriptome_inputs(additional_params)
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
        filter_criteria: OffTargetFilterCriteria,
    ) -> tuple[bool, str | None]:
        """Check if candidate fails off-target filters.

        Returns:
            Tuple of (should_fail, fail_reason)
        """
        # Define filter checks with their thresholds and messages
        checks = [
            (
                filter_criteria.max_transcriptome_hits_0mm,
                transcriptome_0mm,
                f"TRANSCRIPTOME_PERFECT_MATCH ({transcriptome_0mm} hits)",
            ),
            (
                filter_criteria.max_transcriptome_hits_1mm,
                transcriptome_1mm,
                f"TRANSCRIPTOME_1MM ({transcriptome_1mm} hits)",
            ),
            (
                filter_criteria.max_transcriptome_hits_2mm,
                transcriptome_2mm,
                f"TRANSCRIPTOME_2MM ({transcriptome_2mm} hits)",
            ),
            (filter_criteria.max_mirna_perfect_seed, mirna_0mm_seed, f"MIRNA_PERFECT_SEED ({mirna_0mm_seed} hits)"),
            (filter_criteria.max_total_offtarget_hits, total_hits, f"TOTAL_OFFTARGETS ({total_hits} hits)"),
        ]

        # Check all threshold-based filters
        for threshold, value, message in checks:
            if threshold is not None and value > threshold:
                return True, message

        # Check high-risk miRNA (boolean flag)
        if filter_criteria.fail_on_high_risk_mirna and mirna_high_risk > 0:
            return True, f"HIGH_RISK_MIRNA ({mirna_high_risk} hits)"

        return False, None

    def _integrate_offtarget_results(  # noqa: PLR0912
        self,
        candidates: list[SiRNACandidate],
        offtarget_data: dict[str, Any],
        filter_criteria: OffTargetFilterCriteria | None = None,
    ) -> tuple[list[SiRNACandidate], dict[str, int]]:
        """Integrate off-target analysis results into siRNA candidates.

        Args:
            candidates: List of siRNA candidates to update
            offtarget_data: Off-target results from Nextflow pipeline
            filter_criteria: Optional filtering criteria for off-targets

        Returns:
            Tuple of (updated candidates, statistics dict)
        """
        if not offtarget_data or offtarget_data.get("status") != "completed":
            logger.warning("No completed off-target data available for integration")
            return candidates, {}

        if filter_criteria is None:
            filter_criteria = OffTargetFilterCriteria()

        results = offtarget_data.get("results", {})
        stats = {
            "candidates_analyzed": len(candidates),
            "candidates_with_offtargets": 0,
            "failed_perfect_match": 0,
            "failed_transcriptome_1mm": 0,
            "failed_transcriptome_2mm": 0,
            "failed_mirna_seed": 0,
            "failed_high_risk_mirna": 0,
            "human_transcriptome_hits": 0,
            "other_transcriptome_hits": 0,
            "human_mirna_hits": 0,
            "other_mirna_hits": 0,
        }

        for candidate in candidates:
            candidate_id = candidate.id
            offtarget_entry = results.get(candidate_id, {})

            if not offtarget_entry or not offtarget_entry.get("hits"):
                continue

            stats["candidates_with_offtargets"] += 1

            # Parse hit details
            transcriptome_totals = {0: 0, 1: 0, 2: 0}
            transcriptome_human = {0: 0, 1: 0, 2: 0}
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

                # Check if this is a miRNA hit (has species/database/mirna_id fields)
                is_mirna = "mirna_id" in hit or "database" in hit

                if is_mirna:
                    mirna_total += 1
                    if species_is_human:
                        mirna_human_total += 1
                    if seed_mismatches == 0:
                        mirna_0mm_seed_total += 1
                        if species_is_human:
                            mirna_human_0mm_seed += 1
                        # High risk: perfect seed + low penalty score (likely strong binding)
                        if offtarget_score < 5.0:
                            mirna_high_risk_total += 1
                            if species_is_human:
                                mirna_high_risk_human += 1
                    elif seed_mismatches == 1:
                        mirna_1mm_seed += 1
                else:
                    # Transcriptome hit
                    treated_as_human = species_is_human or not species_label
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

            transcriptome_total_hits = sum(transcriptome_totals.values())
            human_transcriptome_hits = sum(transcriptome_human.values())

            # Update candidate fields
            candidate.transcriptome_hits_total = transcriptome_total_hits
            candidate.transcriptome_hits_0mm = transcriptome_totals[0]
            candidate.transcriptome_hits_1mm = transcriptome_totals[1]
            candidate.transcriptome_hits_2mm = transcriptome_totals[2]
            candidate.transcriptome_hits_seed_0mm = transcriptome_seed_0mm
            candidate.mirna_hits_total = mirna_total
            candidate.mirna_hits_0mm_seed = mirna_0mm_seed_total
            candidate.mirna_hits_1mm_seed = mirna_1mm_seed
            candidate.mirna_hits_high_risk = mirna_high_risk_total
            candidate.off_target_count = len(offtarget_entry.get("hits", []))
            candidate.off_target_penalty = offtarget_entry.get("off_target_score", 0.0)

            stats["human_transcriptome_hits"] += human_transcriptome_hits
            stats["other_transcriptome_hits"] += transcriptome_total_hits - human_transcriptome_hits
            stats["human_mirna_hits"] += mirna_human_total
            stats["other_mirna_hits"] += mirna_total - mirna_human_total

            # Apply filtering criteria
            human_total_hits_for_filters = human_transcriptome_hits + mirna_human_total
            should_fail, fail_reason = self._check_offtarget_filters(
                transcriptome_human[0],
                transcriptome_human[1],
                transcriptome_human[2],
                mirna_human_0mm_seed,
                mirna_high_risk_human,
                human_total_hits_for_filters,
                filter_criteria,
            )

            # Update filter status and stats if failed
            if should_fail and fail_reason:
                candidate.passes_filters = fail_reason  # type: ignore
                logger.info(f"Candidate {candidate_id} failed off-target filter: {fail_reason}")

                # Update appropriate stat counter
                if "PERFECT_MATCH" in fail_reason:
                    stats["failed_perfect_match"] += 1
                elif "1MM" in fail_reason:
                    stats["failed_transcriptome_1mm"] += 1
                elif "2MM" in fail_reason:
                    stats["failed_transcriptome_2mm"] += 1
                elif "MIRNA_PERFECT_SEED" in fail_reason:
                    stats["failed_mirna_seed"] += 1
                elif "HIGH_RISK_MIRNA" in fail_reason:
                    stats["failed_high_risk_mirna"] += 1

        return candidates, stats

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
        base.update(
            {
                "pass_count": passed,
                "fail_count": failed,
                "top_n_requested": self.config.top_n,
                "threads_used": self.config.num_threads,
            }
        )
        return base


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
) -> dict[str, Any]:
    """Run complete siRNA design workflow.

    Args:
        gene_query: Gene name or ID to search for
        output_dir: Directory for output files
        input_fasta: Local path or remote URI to an input FASTA file
        database: Database to search (ensembl, refseq, gencode)
        design_mode: Design mode (sirna or mirna)
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
        variant_ids: List of variant identifiers (rsID, chr:pos:ref:alt, or HGVS) to target or avoid
        variant_vcf_file: Path to VCF file containing variants to target or avoid
        variant_mode: How to handle variants (avoid/target/both) - default is avoid
        variant_min_af: Minimum allele frequency threshold for variant filtering (default: 0.01)
        variant_clinvar_filters: Comma-separated ClinVar significance levels to include (default: Pathogenic,Likely pathogenic)
        variant_assembly: Reference genome assembly for variants (only GRCh38 supported)
        log_file: Path to centralized log file
        write_json_summary: Write logs/workflow_summary.json
        num_threads: Optional override for design parallelism
        allow_transcriptome_with_input_fasta: Force transcriptome analysis even when using input FASTA
        default_transcriptome_sources: Ordered list of transcriptome identifiers evaluated by default
        keep_nextflow_work: Keep Nextflow work directory symlink in output

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
        # Parse variant mode
        try:
            variant_mode_enum = VariantMode(variant_mode.lower())
        except ValueError:
            variant_mode_enum = VariantMode.AVOID

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
            dg_5p, dg_3p, asymmetry_score = calc.calculate_asymmetry_score(temp_candidate)

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

        offtarget_results = await workflow._run_nextflow_offtarget_analysis(
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
