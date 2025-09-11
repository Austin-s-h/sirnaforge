"""
siRNAforge Workflow Orchestrator

Coordinates the complete siRNA design pipeline:
1. Transcript retrieval and validation
2. ORF validation and reporting
3. siRNA candidate generation and scoring
4. Top-N candidate selection and reporting
5. Off-target analysis with Nextflow pipeline
"""

import asyncio
import csv
import json
import math
import tempfile
import time
from pathlib import Path
from typing import Optional

import pandas as pd
from rich.console import Console
from rich.progress import Progress

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.core.off_target import OffTargetAnalysisManager
from sirnaforge.data.base import DatabaseType, FastaUtils, TranscriptInfo
from sirnaforge.data.gene_search import GeneSearcher
from sirnaforge.data.orf_analysis import ORFAnalyzer
from sirnaforge.models.schemas import ORFValidationSchema
from sirnaforge.models.sirna import DesignParameters, DesignResult, FilterCriteria, SiRNACandidate
from sirnaforge.pipeline import NextflowConfig, NextflowRunner
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)
console = Console()


class WorkflowConfig:
    """Configuration for the complete siRNA design workflow."""

    def __init__(
        self,
        output_dir: Path,
        gene_query: str,
        input_fasta: Optional[Path] = None,
        database: DatabaseType = DatabaseType.ENSEMBL,
        design_params: Optional[DesignParameters] = None,
        top_n_for_offtarget: int = 10,
        nextflow_config: Optional[dict] = None,
        genome_species: Optional[list[str]] = None,
        log_file: Optional[str] = None,
    ):
        self.output_dir = Path(output_dir)
        # If input_fasta is provided, use its stem as the gene_query identifier
        self.input_fasta = Path(input_fasta) if input_fasta else None
        self.gene_query = gene_query if not self.input_fasta else self.input_fasta.stem
        self.database = database
        self.design_params = design_params or DesignParameters()
        self.top_n_for_offtarget = top_n_for_offtarget
        self.nextflow_config = nextflow_config or {}
        self.genome_species = genome_species or ["human", "rat", "rhesus"]
        self.log_file = log_file

        # Create output structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "transcripts").mkdir(exist_ok=True)
        (self.output_dir / "orf_reports").mkdir(exist_ok=True)
        (self.output_dir / "sirnaforge").mkdir(exist_ok=True)
        (self.output_dir / "off_target").mkdir(exist_ok=True)


class SiRNAWorkflow:
    """Main workflow orchestrator for siRNA design pipeline."""

    def __init__(self, config: WorkflowConfig):
        self.config = config
        self.gene_searcher = GeneSearcher()
        self.orf_analyzer = ORFAnalyzer()
        self.sirnaforgeer = SiRNADesigner(config.design_params)
        self.results: dict = {}

    async def run_complete_workflow(self) -> dict:
        """Run the complete siRNA design workflow."""
        console.print("\n🧬 [bold cyan]Starting siRNAforge Workflow[/bold cyan]")
        console.print(f"Gene Query: [yellow]{self.config.gene_query}[/yellow]")
        console.print(f"Output Directory: [blue]{self.config.output_dir}[/blue]")

        start_time = time.time()

        with Progress() as progress:
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
            "offtarget_summary": offtarget_results,
        }

        # Save workflow summary
        summary_file = self.config.output_dir / "workflow_summary.json"
        with summary_file.open("w") as f:
            json.dump(final_results, f, indent=2, default=str)

        # Save a small manifest that records invocation parameters and environment pointers
        manifest = {
            "gene_query": self.config.gene_query,
            "database": self.config.database.value,
            "output_dir": str(self.config.output_dir),
            "processing_time": total_time,
            "start_time": start_time,
            "end_time": time.time(),
            "tool_versions": final_results.get("design_summary", {}).get("tool_versions", {}),
            "log_file": self.config.log_file,  # Effective log file path if configured
        }

        manifest_file = self.config.output_dir / "workflow_manifest.json"
        with manifest_file.open("w") as mf:
            json.dump(manifest, mf, indent=2, default=str)

        console.print(f"\n✅ [bold green]Workflow completed in {total_time:.2f}s[/bold green]")
        console.print(f"📊 Results saved to: [blue]{self.config.output_dir}[/blue]")

        return final_results

    async def step1_retrieve_transcripts(self, progress: Progress) -> list[TranscriptInfo]:
        """Step 1: Retrieve and validate transcript sequences."""
        task = progress.add_task("[yellow]Fetching transcripts...", total=3)
        # If an input FASTA was provided, read sequences directly and create TranscriptInfo objects
        if self.config.input_fasta:
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
        """Step 3: Design siRNA candidates for valid transcripts."""
        task = progress.add_task("[yellow]Designing siRNAs...", total=2)

        # Create temporary FASTA file for siRNA design
        temp_fasta = self.config.output_dir / "transcripts" / "temp_for_design.fasta"
        sequences = [(f"{t.transcript_id}", t.sequence) for t in transcripts if t.sequence]

        FastaUtils.save_sequences_fasta(sequences, temp_fasta)
        progress.advance(task)

        # Run siRNA design
        design_result = self.sirnaforgeer.design_from_file(str(temp_fasta))
        progress.advance(task)

        # Save design results
        results_file = self.config.output_dir / "sirnaforge" / f"{self.config.gene_query}_sirna_results.csv"
        design_result.save_csv(str(results_file))

        # Clean up temp file
        temp_fasta.unlink()

        console.print(f"🎯 Generated {len(design_result.candidates)} siRNA candidates")
        console.print(f"   Top {len(design_result.top_candidates)} candidates selected for further analysis")

        return design_result

    async def step4_generate_reports(self, design_results: DesignResult) -> None:
        """Step 4: Generate comprehensive reports."""

        # Generate candidate summary report
        summary_file = self.config.output_dir / "sirnaforge" / f"{self.config.gene_query}_candidate_summary.txt"
        self._generate_candidate_report(design_results, summary_file)

        # Generate top candidates FASTA for off-target analysis
        top_candidates_fasta = self.config.output_dir / "sirnaforge" / f"{self.config.gene_query}_top_candidates.fasta"
        top_sequences = [
            (f"{c.id} score:{c.composite_score:.1f}", c.guide_sequence)
            for c in design_results.top_candidates[: self.config.top_n_for_offtarget]
        ]

        FastaUtils.save_sequences_fasta(top_sequences, top_candidates_fasta)

        console.print("📋 Generated comprehensive reports")
        console.print("   - ORF validation report: orf_reports/")
        console.print("   - siRNA candidate summary: sirnaforge/")
        console.print(f"   - Top {len(top_sequences)} candidates for off-target analysis")

    async def step5_offtarget_analysis(self, design_results: DesignResult) -> dict:
        """Step 5: Run off-target analysis using embedded Nextflow pipeline."""
        top_candidates = design_results.top_candidates[: self.config.top_n_for_offtarget]

        if not top_candidates:
            console.print("⚠️  No candidates available for off-target analysis")
            return {"status": "skipped", "reason": "no_candidates"}

        # Prepare input files
        input_fasta = await self._prepare_offtarget_input(top_candidates)

        # Try Nextflow pipeline first, fall back to basic analysis
        try:
            return await self._run_nextflow_offtarget_analysis(top_candidates, input_fasta)
        except Exception as e:
            console.print(f"⚠️  Nextflow execution failed: {e}")
            logger.exception("Nextflow pipeline execution error")
            return await self._basic_offtarget_analysis(top_candidates)

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
            return await self._basic_offtarget_analysis(candidates)

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
            console.print("⚠️  Nextflow not available, using basic analysis")
            return False
        if not validation["workflow_files"]:
            console.print("⚠️  Nextflow workflows not found, using basic analysis")
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

    def _generate_orf_report(self, orf_results: dict, report_file: Path) -> None:
        """Generate ORF validation report in tab-delimited format with schema validation."""
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
                ]
            )
            # Set correct dtypes to match schema
            empty_df = empty_df.astype(
                {
                    "transcript_id": str,
                    "sequence_length": "Int64",
                    "gc_content": float,
                    "orfs_found": "Int64",
                    "has_valid_orf": bool,
                    "longest_orf_start": "Int64",
                    "longest_orf_end": "Int64",
                    "longest_orf_length": "Int64",
                    "longest_orf_frame": "Int64",
                    "start_codon": str,
                    "stop_codon": str,
                    "orf_gc_content": float,
                }
            )
            validated_df = ORFValidationSchema.validate(empty_df)
            validated_df.to_csv(report_file, sep="\t", index=False)
            return

        # Prepare data for DataFrame
        rows = []
        for transcript_id, analysis in orf_results.items():
            row_data = {
                "transcript_id": transcript_id,
                "sequence_length": analysis.sequence_length,
                "gc_content": analysis.gc_content,
                "orfs_found": len(analysis.orfs),
                "has_valid_orf": analysis.has_valid_orf,
            }

            if analysis.longest_orf:
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
        validated_df = ORFValidationSchema.validate(df)
        logger.info(f"ORF report schema validation passed for {len(validated_df)} transcripts")

        # Write validated DataFrame to file
        validated_df.to_csv(report_file, sep="\t", index=False)

    def _generate_candidate_report(self, design_results: DesignResult, report_file: Path) -> None:
        """Generate siRNA candidate summary report."""
        with report_file.open("w") as f:
            f.write(f"siRNAforge Results for {self.config.gene_query}\n")
            f.write("=" * 60 + "\n\n")

            summary = design_results.get_summary()
            f.write("SUMMARY\n")
            f.write(f"  Input Sequences: {summary['input_sequences']}\n")
            f.write(f"  Total Candidates: {summary['total_candidates']}\n")
            f.write(f"  Filtered Candidates: {summary['filtered_candidates']}\n")
            f.write(f"  Top Candidates: {summary['top_candidates']}\n")
            f.write(f"  Best Score: {summary['best_score']:.1f}\n")
            f.write(f"  Processing Time: {summary['processing_time']}\n\n")

            f.write("TOP 10 CANDIDATES\n")
            f.write("-" * 40 + "\n")

            for i, candidate in enumerate(design_results.top_candidates[:10], 1):
                f.write(f"{i}. {candidate.id}\n")
                f.write(f"   Guide Sequence: {candidate.guide_sequence}\n")
                f.write(f"   Score: {candidate.composite_score:.1f}\n")
                f.write(f"   GC Content: {candidate.gc_content:.1f}%\n")
                f.write(f"   Position: {candidate.position}\n")
                f.write(f"   Passes Filters: {'✓' if candidate.passes_filters else '✗'}\n")
                if candidate.quality_issues:
                    f.write(f"   Issues: {'; '.join(candidate.quality_issues)}\n")
                f.write("\n")

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
        return design_results.get_summary()


# Convenience function for running complete workflow
async def run_sirna_workflow(
    gene_query: str,
    output_dir: str,
    input_fasta: Optional[str] = None,
    database: str = "ensembl",
    top_n_candidates: int = 20,
    top_n_offtarget: int = 10,
    genome_species: Optional[list[str]] = None,
    gc_min: float = 30.0,
    gc_max: float = 52.0,
    sirna_length: int = 21,
    log_file: Optional[str] = None,
) -> dict:
    """
    Run complete siRNA design workflow.

    Args:
        gene_query: Gene name or ID to search for
        output_dir: Directory for output files
        database: Database to search (ensembl, refseq, gencode)
        top_n_candidates: Number of top candidates to generate
        top_n_offtarget: Number of candidates for off-target analysis
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

    # Configure workflow
    design_params = DesignParameters(top_n=top_n_candidates, sirna_length=sirna_length, filters=filter_criteria)
    database_enum = DatabaseType(database.lower())

    config = WorkflowConfig(
        output_dir=Path(output_dir),
        gene_query=gene_query,
        input_fasta=Path(input_fasta) if input_fasta else None,
        database=database_enum,
        design_params=design_params,
        top_n_for_offtarget=top_n_offtarget,
        genome_species=genome_species or ["human", "rat", "rhesus"],
        log_file=log_file,
    )

    # Run workflow
    workflow = SiRNAWorkflow(config)
    return await workflow.run_complete_workflow()


if __name__ == "__main__":
    # Example usage
    async def main() -> None:
        # Use secure temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            results = await run_sirna_workflow(
                gene_query="TP53", output_dir=temp_dir, top_n_candidates=20, top_n_offtarget=10
            )
            print(f"Workflow completed: {results}")

    asyncio.run(main())
