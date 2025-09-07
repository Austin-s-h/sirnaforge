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
import subprocess
import time
from pathlib import Path
from typing import Optional

from rich.console import Console
from rich.progress import Progress

from sirnaforge.core.design import SiRNADesigner
from sirnaforge.core.off_target import OffTargetAnalyzer
from sirnaforge.data.base import DatabaseType, FastaUtils, TranscriptInfo
from sirnaforge.data.gene_search import GeneSearcher
from sirnaforge.data.orf_analysis import ORFAnalyzer
from sirnaforge.models.sirna import DesignParameters, DesignResult, FilterCriteria, SiRNACandidate
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)
console = Console()


class WorkflowConfig:
    """Configuration for the complete siRNA design workflow."""

    def __init__(
        self,
        output_dir: Path,
        gene_query: str,
        database: DatabaseType = DatabaseType.ENSEMBL,
        design_params: Optional[DesignParameters] = None,
        top_n_for_offtarget: int = 10,
        nextflow_config: Optional[dict] = None,
        genome_species: Optional[list[str]] = None,
    ):
        self.output_dir = Path(output_dir)
        self.gene_query = gene_query
        self.database = database
        self.design_params = design_params or DesignParameters()
        self.top_n_for_offtarget = top_n_for_offtarget
        self.nextflow_config = nextflow_config or {}
        self.genome_species = genome_species or ["human", "rat", "rhesus"]

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

        console.print(f"\n✅ [bold green]Workflow completed in {total_time:.2f}s[/bold green]")
        console.print(f"📊 Results saved to: [blue]{self.config.output_dir}[/blue]")

        return final_results

    async def step1_retrieve_transcripts(self, progress: Progress) -> list[TranscriptInfo]:
        """Step 1: Retrieve and validate transcript sequences."""
        task = progress.add_task("[yellow]Fetching transcripts...", total=3)

        # Search for gene
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
        """Step 5: Run off-target analysis using Nextflow pipeline."""

        top_candidates = design_results.top_candidates[: self.config.top_n_for_offtarget]

        if not top_candidates:
            console.print("⚠️  No candidates available for off-target analysis")
            return {"status": "skipped", "reason": "no_candidates"}

        # Create Nextflow input
        input_fasta = self.config.output_dir / "off_target" / "input_candidates.fasta"
        sequences = [(f"{c.id}", c.guide_sequence) for c in top_candidates]

        FastaUtils.save_sequences_fasta(sequences, input_fasta)

        # Prepare Nextflow command
        nf_script = self._find_nextflow_script()
        if not nf_script:
            console.print("⚠️  Nextflow off-target pipeline not found, using basic analysis")
            return await self._basic_offtarget_analysis(top_candidates)

        # Run Nextflow pipeline
        try:
            nf_work_dir = self.config.output_dir / "off_target" / "nextflow_work"
            nf_output_dir = self.config.output_dir / "off_target" / "results"

            cmd = [
                "nextflow",
                "run",
                str(nf_script),
                "--input",
                str(input_fasta),
                "--outdir",
                str(nf_output_dir),
                "--genome_species",
                ",".join(self.config.genome_species),
                "-w",
                str(nf_work_dir),
                "-resume",
            ]

            console.print("🚀 Running Nextflow off-target analysis...")
            console.print(f"   Command: {' '.join(cmd)}")

            result = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=1800)  # 30 min timeout

            if result.returncode == 0:
                console.print("✅ Nextflow pipeline completed successfully")
                parsed = await self._parse_nextflow_results(nf_output_dir)

                # Map parsed results back to the top candidates so classification can use them
                mapped = {}
                for c in top_candidates:
                    qid = c.id
                    entry = parsed.get("results", {}).get(qid)
                    if entry:
                        mapped[qid] = {
                            "off_target_count": entry.get("off_target_count", 0),
                            "off_target_score": entry.get("off_target_score", 0.0),
                            "hits": entry.get("hits", []),
                        }
                        # Update candidate object fields if available
                        try:
                            c.off_target_count = mapped[qid]["off_target_count"]
                            c.off_target_penalty = mapped[qid]["off_target_score"]
                        except Exception:
                            pass
                    else:
                        mapped[qid] = {"off_target_count": 0, "off_target_score": 0.0, "hits": []}

                return {
                    "status": "completed",
                    "method": "nextflow",
                    "output_dir": str(nf_output_dir),
                    "results": mapped,
                }
            console.print(f"❌ Nextflow pipeline failed: {result.stderr}")
            return await self._basic_offtarget_analysis(top_candidates)

        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            console.print(f"⚠️  Nextflow execution failed: {e}")
            return await self._basic_offtarget_analysis(top_candidates)

    def _find_nextflow_script(self) -> Optional[Path]:
        """Find the Nextflow off-target analysis script."""
        possible_locations = [
            Path(__file__).parent.parent.parent / "nextflow_pipeline" / "main.nf",
            Path.cwd() / "nextflow_pipeline" / "main.nf",
            Path(__file__).parent / "pipeline" / "offtarget_analysis.nf",
        ]

        for location in possible_locations:
            if location.exists():
                return location

        return None

    async def _basic_offtarget_analysis(self, candidates: list[SiRNACandidate]) -> dict:
        """Fallback basic off-target analysis."""
        analyzer = OffTargetAnalyzer()
        results = {}

        for candidate in candidates:
            off_target_count, penalty = analyzer.analyze_off_targets(candidate)
            score = analyzer.calculate_off_target_score(candidate)

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
        """Generate ORF validation report."""
        with report_file.open("w") as f:
            f.write(f"ORF Validation Report for {self.config.gene_query}\n")
            f.write("=" * 60 + "\n\n")

            valid_count = sum(1 for r in orf_results.values() if r.has_valid_orf)
            f.write(f"Summary: {valid_count}/{len(orf_results)} transcripts have valid ORFs\n\n")

            for transcript_id, analysis in orf_results.items():
                f.write(f"Transcript: {transcript_id}\n")
                f.write(f"  Sequence Length: {analysis.sequence_length} bp\n")
                f.write(f"  GC Content: {analysis.gc_content:.1f}%\n")
                f.write(f"  ORFs Found: {len(analysis.orfs)}\n")
                f.write(f"  Valid ORF: {'✓' if analysis.has_valid_orf else '✗'}\n")
                if analysis.longest_orf:
                    orf = analysis.longest_orf
                    f.write(f"  Longest ORF: {orf.start_pos}-{orf.end_pos} ({orf.length} bp)\n")
                f.write("\n")

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
            "avg_length": sum(t.length for t in transcripts if t.length is not None)
            / len([t for t in transcripts if t.length is not None])
            if any(t.length is not None for t in transcripts)
            else 0,
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
    database: str = "ensembl",
    top_n_candidates: int = 20,
    top_n_offtarget: int = 10,
    genome_species: Optional[list[str]] = None,
    gc_min: float = 30.0,
    gc_max: float = 52.0,
    sirna_length: int = 21,
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
        database=database_enum,
        design_params=design_params,
        top_n_for_offtarget=top_n_offtarget,
        genome_species=genome_species or ["human", "rat", "rhesus"],
    )

    # Run workflow
    workflow = SiRNAWorkflow(config)
    return await workflow.run_complete_workflow()


if __name__ == "__main__":
    # Example usage
    import asyncio

    async def main() -> None:
        results = await run_sirna_workflow(
            gene_query="TP53", output_dir="/tmp/sirna_workflow_test", top_n_candidates=20, top_n_offtarget=10
        )
        print(f"Workflow completed: {results}")

    asyncio.run(main())
