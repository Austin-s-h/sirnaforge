"""
siRNAforge Workflow Orchestrator

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
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from pandera.typing import DataFrame
from rich.console import Console
from rich.progress import Progress

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.core.off_target import OffTargetAnalysisManager
from sirnaforge.data.base import DatabaseType, FastaUtils, TranscriptInfo
from sirnaforge.data.gene_search import GeneSearcher
from sirnaforge.data.orf_analysis import ORFAnalyzer
from sirnaforge.models.schemas import ORFValidationSchema, SiRNACandidateSchema
from sirnaforge.models.sirna import DesignParameters, DesignResult, FilterCriteria, SiRNACandidate
from sirnaforge.models.sirna import SiRNACandidate as _ModelSiRNACandidate
from sirnaforge.pipeline import NextflowConfig, NextflowRunner
from sirnaforge.utils.logging_utils import get_logger
from sirnaforge.utils.modification_patterns import apply_modifications_to_candidate, get_modification_summary
from sirnaforge.utils.resource_resolver import InputSource, resolve_input_source
from sirnaforge.validation import ValidationConfig, ValidationMiddleware

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
        nextflow_config: dict | None = None,
        genome_species: list[str] | None = None,
        validation_config: ValidationConfig | None = None,
        log_file: str | None = None,
        write_json_summary: bool = True,
        num_threads: int | None = None,
        input_source: InputSource | None = None,
    ):
        self.output_dir = Path(output_dir)
        self.input_source = input_source

        resolved_input = input_source.local_path if input_source else (Path(input_fasta) if input_fasta else None)
        self.input_fasta = resolved_input
        # If input_fasta is provided, use its stem as the gene_query identifier
        self.gene_query = gene_query if resolved_input is None else resolved_input.stem
        self.database = database
        self.design_params = design_params or DesignParameters()
        # single source of truth: number of candidates selected everywhere
        self.top_n = self.design_params.top_n
        self.nextflow_config = nextflow_config or {}
        self.genome_species = genome_species or ["human", "rat", "rhesus"]
        self.validation_config = validation_config or ValidationConfig()
        self.log_file = log_file
        self.write_json_summary = write_json_summary
        # Parallelism for design stage (cap at 4 CPUs for better efficiency)
        requested_threads = num_threads if num_threads is not None else (os.cpu_count() or 4)
        self.num_threads = max(1, min(4, requested_threads))

        # Create output structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "transcripts").mkdir(exist_ok=True)
        (self.output_dir / "orf_reports").mkdir(exist_ok=True)
        (self.output_dir / "sirnaforge").mkdir(exist_ok=True)
        (self.output_dir / "off_target").mkdir(exist_ok=True)
        (self.output_dir / "logs").mkdir(exist_ok=True)


class SiRNAWorkflow:
    """Main workflow orchestrator for siRNA design pipeline."""

    def __init__(self, config: WorkflowConfig):
        self.config = config
        self.gene_searcher = GeneSearcher()
        self.orf_analyzer = ORFAnalyzer()
        self.validation = ValidationMiddleware(config.validation_config)
        self.sirnaforgeer = SiRNADesigner(config.design_params)
        self.results: dict = {}

    async def run_complete_workflow(self) -> dict:
        """Run the complete siRNA design workflow."""
        console.print("\n🧬 [bold cyan]Starting siRNAforge Workflow[/bold cyan]")
        console.print(f"Gene Query: [yellow]{self.config.gene_query}[/yellow]")
        console.print(f"Output Directory: [blue]{self.config.output_dir}[/blue]")

        start_time = time.time()

        # Validate input parameters (quiet: avoid verbose warnings in console)
        _ = self.validation.validate_input_parameters(self.config.design_params)

        with Progress(console=console) as progress:
            main_task = progress.add_task("[cyan]Overall Progress", total=5)

            # Step 1: Transcript Retrieval
            progress.update(main_task, description="[cyan]Retrieving transcripts...")
            transcripts = await self.step1_retrieve_transcripts(progress)
            progress.advance(main_task)

            # Step 2: ORF Validation
            progress.update(main_task, description="[cyan]Validating ORFs...")
            orf_results = await self.step2_validate_orfs(transcripts, progress)
            progress.advance(main_task)

            # Step 3: siRNAforge
            progress.update(main_task, description="[cyan]Designing siRNAs...")
            design_results = await self.step3_design_sirnas(transcripts, progress)
            progress.advance(main_task)

            # Step 4: Generate Reports
            progress.update(main_task, description="[cyan]Generating reports...")
            await self.step4_generate_reports(design_results)
            progress.advance(main_task)

            # Step 5: Off-target Analysis
            progress.update(main_task, description="[cyan]Running off-target analysis...")
            offtarget_results = await self.step5_offtarget_analysis(design_results)
            progress.advance(main_task)

        total_time = time.time() - start_time

        # Compile final results
        # Serialize authoritative design parameters into the workflow summary.
        dp = self.config.design_params
        design_parameters = {
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

        final_results = {
            "workflow_config": {
                "gene_query": self.config.gene_query,
                "database": self.config.database.value,
                "output_dir": str(self.config.output_dir),
                "processing_time": total_time,
            },
            "transcript_summary": self._summarize_transcripts(transcripts),
            "orf_summary": self._summarize_orf_results(orf_results),
            "design_summary": self._summarize_design_results(design_results),
            "design_parameters": design_parameters,
            "offtarget_summary": offtarget_results,
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

    async def step2_validate_orfs(self, transcripts: list[TranscriptInfo], progress: Progress) -> dict:
        """Step 2: Validate ORFs and generate validation report."""
        task = progress.add_task("[yellow]Analyzing ORFs...", total=len(transcripts) + 1)

        orf_results = {}
        valid_transcripts = []

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
            progress.advance(task)
            _ = self.validation.validate_design_results(design_result)
            temp_fasta.unlink(missing_ok=True)
            console.print(f"🎯 Generated {len(design_result.candidates)} siRNA candidates")
            console.print(f"   Top {len(design_result.top_candidates)} candidates selected for further analysis")
            return design_result

        # Parallel per-transcript path
        start = time.time()
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

        processing_time = time.time() - start
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
        )

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

        batches = []
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

        # Apply chemical modifications if enabled
        if self.config.design_params.apply_modifications:
            self._apply_modifications_to_results(design_results)

        # Emit candidate CSVs: always keep ALL and a filtered PASSING file
        try:
            # Import modification summary helper
            # Build DataFrame from all candidates with columns matching SiRNACandidateSchema
            rows = []
            for c in design_results.candidates:
                cs = getattr(c, "component_scores", {}) or {}

                # Get modification summary if modifications were applied
                mod_summary = get_modification_summary(c) if c.guide_metadata else {}

                rows.append(
                    {
                        "id": c.id,
                        "transcript_id": c.transcript_id,
                        "position": c.position,
                        "guide_sequence": c.guide_sequence,
                        "passenger_sequence": c.passenger_sequence,
                        "gc_content": c.gc_content,
                        "asymmetry_score": c.asymmetry_score,
                        "paired_fraction": c.paired_fraction,
                        # Thermodynamics
                        "structure": getattr(c, "structure", None),
                        "mfe": getattr(c, "mfe", None),
                        "duplex_stability_dg": c.duplex_stability,
                        "duplex_stability_score": cs.get("duplex_stability_score"),
                        "dg_5p": cs.get("dg_5p"),
                        "dg_3p": cs.get("dg_3p"),
                        "delta_dg_end": cs.get("delta_dg_end"),
                        "melting_temp_c": cs.get("melting_temp_c"),
                        "off_target_count": c.off_target_count,
                        "transcript_hit_count": c.transcript_hit_count,
                        "transcript_hit_fraction": c.transcript_hit_fraction,
                        "composite_score": c.composite_score,
                        "passes_filters": (
                            c.passes_filters.value if hasattr(c.passes_filters, "value") else c.passes_filters
                        ),
                        # Chemical modifications
                        "guide_overhang": mod_summary.get("guide_overhang", ""),
                        "guide_modifications": mod_summary.get("guide_modifications", ""),
                        "passenger_overhang": mod_summary.get("passenger_overhang", ""),
                        "passenger_modifications": mod_summary.get("passenger_modifications", ""),
                    }
                )

            all_df = pd.DataFrame(rows)
            # Validate with schema (will raise if invalid)
            validated_all = SiRNACandidateSchema.validate(all_df)

            # Note: design parameters are not appended as per-row columns anymore
            # to avoid cluttering the candidate CSVs. Full parameters are included
            # in the `workflow_summary.json` under the key `design_parameters`.

            # Normalize passes_filters values so legacy booleans are treated as 'PASS'/'FAIL'
            def _normalize_pass(v: Any) -> str:
                """Normalize various legacy pass values into 'PASS' or 'FAIL'.

                Accepts booleans, numeric 0/1, and common strings. For unknown
                values we fall back to truthiness (truthy -> 'PASS', falsy -> 'FAIL').
                Always returns a string to satisfy static typing checks.
                """
                result: str
                try:
                    # Handle booleans and numeric 0/1
                    if v is True or (isinstance(v, (int, float)) and v == 1):
                        result = "PASS"
                    elif v is False or (isinstance(v, (int, float)) and v == 0):
                        result = "FAIL"
                    elif isinstance(v, str):
                        vs = v.strip().upper()
                        if vs in {"PASS", "TRUE", "YES"}:
                            result = "PASS"
                        elif vs in {"FAIL", "FALSE", "NO"}:
                            result = "FAIL"
                        else:
                            # Unknown string: return uppercased form
                            result = vs
                    else:
                        # Fallback: use truthiness
                        result = "PASS" if bool(v) else "FAIL"
                except Exception:
                    # On unexpected errors, fail-safe to 'FAIL'
                    result = "FAIL"

                return result

            validated_all["passes_filters"] = validated_all["passes_filters"].apply(_normalize_pass)

            # Split into pass/fail and write ALL + PASS files (append later too for consistency)
            pass_df = validated_all[validated_all["passes_filters"] == "PASS"].copy()
            base = self.config.output_dir / "sirnaforge"
            out_all = base / f"{self.config.gene_query}_all.csv"
            out_pass = base / f"{self.config.gene_query}_pass.csv"
            validated_all.to_csv(out_all, index=False)
            pass_df.to_csv(out_pass, index=False)

            # Generate FASTA file for passing candidates
            try:
                out_pass_fasta = base / f"{self.config.gene_query}_pass.fasta"
                self._write_pass_candidates_fasta(pass_df, out_pass_fasta)
            except Exception as e:
                logger.warning(f"Failed to write PASS candidates FASTA: {e}")

        except Exception as e:  # Do not fail workflow for reporting extras
            logger.warning(f"Failed to write all/pass CSVs: {e}")

        # Generate machine-readable candidate summary JSON (concise)
        # No separate thermodynamic CSV or candidate summary JSON is written anymore.

        # Build a FAIR manifest with checksums and counts
        try:
            manifest = self._build_fair_manifest(
                all_csv=(self.config.output_dir / "sirnaforge" / f"{self.config.gene_query}_all.csv"),
                pass_csv=(self.config.output_dir / "sirnaforge" / f"{self.config.gene_query}_pass.csv"),
                pass_fasta=(self.config.output_dir / "sirnaforge" / f"{self.config.gene_query}_pass.fasta"),
                orf_report=(self.config.output_dir / "orf_reports" / f"{self.config.gene_query}_orf_validation.txt"),
            )
            manifest_path = self.config.output_dir / "sirnaforge" / "manifest.json"
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
            sequences = []
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
    ) -> dict:
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

    async def step5_offtarget_analysis(self, design_results: DesignResult) -> dict:
        """Step 5: Run off-target analysis using embedded Nextflow pipeline."""
        top_candidates = design_results.top_candidates[: self.config.top_n]

        if not top_candidates:
            console.print("⚠️  No candidates available for off-target analysis")
            return {"status": "skipped", "reason": "no_candidates"}

        # Prepare input files
        input_fasta = await self._prepare_offtarget_input(top_candidates)

        # If user disabled off-target checking via design parameters, skip entirely
        if not getattr(self.config.design_params, "check_off_targets", True):
            console.print("⚠️  Off-target analysis skipped by user request")
            return {"status": "skipped", "reason": "user_disabled"}

        # Try Nextflow pipeline first. We do NOT run the simplistic sequence-based fallback
        # (it produces low-value results) when Nextflow is unavailable. Instead mark as skipped
        # so downstream steps/users can see the explicit reason.
        try:
            return await self._run_nextflow_offtarget_analysis(top_candidates, input_fasta)
        except Exception as e:
            console.print(f"⚠️  Nextflow execution failed: {e}")
            logger.exception("Nextflow pipeline execution error")
            return {"status": "skipped", "reason": "nextflow_failed", "error": str(e)}

    async def _prepare_offtarget_input(self, candidates: list[SiRNACandidate]) -> Path:
        """Prepare FASTA input file for off-target analysis."""
        input_fasta = self.config.output_dir / "off_target" / "input_candidates.fasta"
        sequences = [(f"{c.id}", c.guide_sequence) for c in candidates]
        FastaUtils.save_sequences_fasta(sequences, input_fasta)
        return input_fasta

    async def _run_nextflow_offtarget_analysis(self, candidates: list[SiRNACandidate], input_fasta: Path) -> dict:
        """Run Nextflow-based off-target analysis."""
        # Configure Nextflow runner
        runner = self._setup_nextflow_runner()

        # Validate installation
        if not self._validate_nextflow_environment(runner):
            # Do not fall back to basic analysis here; report skipped instead
            return {"status": "skipped", "reason": "nextflow_unavailable"}

        # Execute pipeline
        console.print("🚀 Running embedded Nextflow off-target analysis...")
        nf_output_dir = self.config.output_dir / "off_target" / "results"

        results = await runner.run_offtarget_analysis(
            input_file=input_fasta,
            output_dir=nf_output_dir,
            genome_species=self.config.genome_species,
            additional_params=self.config.nextflow_config,
            show_progress=True,
        )

        if results["status"] == "completed":
            return await self._process_nextflow_results(candidates, nf_output_dir, results)

        console.print(f"❌ Nextflow pipeline failed: {results}")
        return await self._basic_offtarget_analysis(candidates)

    def _setup_nextflow_runner(self) -> NextflowRunner:
        """Configure Nextflow runner with user settings."""
        nf_config = NextflowConfig.for_production()
        if self.config.nextflow_config:
            for key, value in self.config.nextflow_config.items():
                setattr(nf_config, key, value)
        return NextflowRunner(nf_config)

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
        self, candidates: list[SiRNACandidate], output_dir: Path, results: dict
    ) -> dict:
        """Process and map Nextflow pipeline results to candidates."""
        console.print("✅ Nextflow pipeline completed successfully")
        parsed = await self._parse_nextflow_results(output_dir)

        # Map parsed results back to candidates
        mapped = {}
        for c in candidates:
            qid = c.id
            entry = parsed.get("results", {}).get(qid)
            if entry:
                mapped[qid] = {
                    "off_target_count": entry.get("off_target_count", 0),
                    "off_target_score": entry.get("off_target_score", 0.0),
                    "hits": entry.get("hits", []),
                }
                # Update candidate object fields
                self._update_candidate_scores(c, mapped[qid])
            else:
                mapped[qid] = {"off_target_count": 0, "off_target_score": 0.0, "hits": []}

        return {
            "status": "completed",
            "method": "embedded_nextflow",
            "output_dir": str(output_dir),
            "results": mapped,
            "execution_metadata": results,
        }

    def _update_candidate_scores(self, candidate: SiRNACandidate, result_data: dict) -> None:
        """Update candidate object with off-target analysis results."""
        try:
            candidate.off_target_count = result_data["off_target_count"]
            candidate.off_target_penalty = result_data["off_target_score"]
        except KeyError as e:
            logger.warning(f"Missing key in result data: {e}")
        except Exception as e:
            logger.error(f"Unexpected error updating candidate scores: {e}")

    async def _basic_offtarget_analysis(self, candidates: list[SiRNACandidate]) -> dict:
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
        return {"status": "completed", "method": "basic", "results": results}

    async def _parse_nextflow_results(self, output_dir: Path) -> dict:  # noqa: PLR0912
        """Parse results from Nextflow off-target analysis."""
        results: dict = {}

        if not output_dir.exists():
            return {"status": "missing", "method": "nextflow", "output_dir": str(output_dir), "results": results}

        # Prefer combined TSV if present
        combined_tsv = output_dir / "combined_offtargets.tsv"
        combined_json = output_dir / "combined_offtargets.json"

        def _ensure_row_key(d: dict, key: str, default: int = 0) -> None:
            if key not in d:
                d[key] = default

        if combined_tsv.exists():
            with combined_tsv.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    qname = row.get("qname") or row.get("query") or row.get("id")
                    if not qname:
                        continue
                    try:
                        score = float(row.get("offtarget_score") or 0)
                    except Exception:
                        score = 0.0

                    entry = results.setdefault(qname, {"off_target_count": 0, "off_target_score": 0.0, "hits": []})
                    entry["off_target_count"] += 1
                    # keep the maximum per-query off-target score as a conservative summary
                    entry["off_target_score"] = max(entry["off_target_score"], score)
                    entry["hits"].append(row)

        elif combined_json.exists():
            with combined_json.open() as fh:
                try:
                    data = json.load(fh)
                except Exception:
                    data = []
                for item in data:
                    qname = item.get("qname") or item.get("query") or item.get("id")
                    if not qname:
                        continue
                    try:
                        score = float(item.get("offtarget_score") or 0)
                    except Exception:
                        score = 0.0
                    entry = results.setdefault(qname, {"off_target_count": 0, "off_target_score": 0.0, "hits": []})
                    entry["off_target_count"] += 1
                    entry["off_target_score"] = max(entry["off_target_score"], score)
                    entry["hits"].append(item)

        else:
            # Fallback: scan for any per-species TSV files under output_dir
            files = list(output_dir.glob("**/*_offtargets.tsv"))
            for fpath in files:
                with Path(fpath).open() as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    for row in reader:
                        qname = row.get("qname") or row.get("query") or row.get("id")
                        if not qname:
                            continue
                        try:
                            score = float(row.get("offtarget_score") or 0)
                        except Exception:
                            score = 0.0
                        entry = results.setdefault(qname, {"off_target_count": 0, "off_target_score": 0.0, "hits": []})
                        entry["off_target_count"] += 1
                        entry["off_target_score"] = max(entry["off_target_score"], score)
                        entry["hits"].append(row)

        return {"status": "completed", "method": "nextflow", "output_dir": str(output_dir), "results": results}

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
        rows = []
        for transcript_id, analysis in orf_results.items():
            row_data = {
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

    def _generate_candidate_report(self, design_results: DesignResult, report_file: Path) -> None:  # noqa: ARG002
        """Deprecated: summary report disabled (kept for backward API compatibility)."""
        report_file.parent.mkdir(parents=True, exist_ok=True)
        with report_file.open("w") as f:
            f.write("# Candidate summary report generation disabled.\n")

    def _summarize_transcripts(self, transcripts: list[TranscriptInfo]) -> dict:
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

    def _summarize_orf_results(self, orf_results: dict) -> dict:
        """Summarize ORF validation results."""
        results = orf_results.get("results", {})
        valid_count = sum(1 for r in results.values() if r.has_valid_orf)

        return {
            "total_analyzed": len(results),
            "valid_orfs": valid_count,
            "validation_rate": valid_count / len(results) if results else 0,
        }

    def _summarize_design_results(self, design_results: DesignResult) -> dict:
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
    top_n_candidates: int = 20,
    genome_species: list[str] | None = None,
    gc_min: float = 30.0,
    gc_max: float = 52.0,
    sirna_length: int = 21,
    modification_pattern: str = "standard_2ome",
    overhang: str = "dTdT",
    log_file: str | None = None,
    write_json_summary: bool = True,
) -> dict:
    """
    Run complete siRNA design workflow.

    Args:
        gene_query: Gene name or ID to search for
        output_dir: Directory for output files
        database: Database to search (ensembl, refseq, gencode)
        top_n_candidates: Number of top candidates to generate
        genome_species: Species genomes for off-target analysis
        gc_min: Minimum GC content percentage
        gc_max: Maximum GC content percentage
        sirna_length: siRNA length in nucleotides

    Returns:
        Dictionary with complete workflow results
    """
    # Configure filter criteria
    filter_criteria = FilterCriteria(
        gc_min=gc_min,
        gc_max=gc_max,
    )

    # Configure workflow with modification parameters
    design_params = DesignParameters(
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

    config = WorkflowConfig(
        output_dir=output_path,
        gene_query=gene_query,
        input_fasta=input_path,
        database=database_enum,
        design_params=design_params,
        genome_species=genome_species or ["human", "rat", "rhesus"],
        log_file=log_file,
        write_json_summary=write_json_summary,
        input_source=resolved_input,
    )

    # Run workflow
    workflow = SiRNAWorkflow(config)
    return await workflow.run_complete_workflow()


if __name__ == "__main__":
    # Example usage
    async def main() -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            results = await run_sirna_workflow(gene_query="TP53", output_dir=temp_dir, top_n_candidates=20)
            print(f"Workflow completed: {results}")

    asyncio.run(main())
