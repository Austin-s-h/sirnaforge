"""Modern CLI for siRNAforge using Typer and Rich."""

# Configure environment for ASCII compatibility before importing typer/rich
import os

os.environ.setdefault("FORCE_COLOR", "0")
os.environ.setdefault("NO_COLOR", "1")
os.environ.setdefault("TERM", "dumb")

import asyncio
import json
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Optional, TypeVar

import typer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

from sirnaforge.data.mirna_manager import MiRNADatabaseManager

# Monkey patch Rich console (best-effort). If this fails we continue with default behavior.
try:
    import rich.console

    original_init = rich.console.Console.__init__

    def patched_init(self: "rich.console.Console", *args: Any, **kwargs: Any) -> None:  # noqa: D401
        # Force simplified terminal capabilities for deterministic CI output.
        kwargs["legacy_windows"] = True
        kwargs["force_terminal"] = False
        original_init(self, *args, **kwargs)

    if not TYPE_CHECKING:  # Avoid confusing type checkers
        rich.console.Console.__init__ = patched_init  # type: ignore[attr-defined]
except (ImportError, AttributeError):  # Narrow exceptions
    # Silently ignore—output formatting will just be richer if available.
    # nosec B110 - acceptable silent fallback; not security relevant
    pass

from sirnaforge import __author__, __version__
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.data.base import DatabaseType
from sirnaforge.data.gene_search import (
    GeneSearcher,
    search_gene_sync,
    search_gene_with_fallback_sync,
    search_multiple_databases_sync,
)
from sirnaforge.models.sirna import DesignParameters, FilterCriteria
from sirnaforge.modifications import merge_metadata_into_fasta, parse_header
from sirnaforge.utils.logging_utils import configure_logging
from sirnaforge.workflow import run_sirna_workflow

app = typer.Typer(
    name="sirnaforge",
    help="siRNAforge - siRNA design toolkit for gene silencing",
    rich_markup_mode="rich",
)
# Configure console to use ASCII box characters for better compatibility
console = Console(force_terminal=False, legacy_windows=True)

# mypy-friendly alias for Typer command decorator
T = TypeVar("T", bound=Callable[..., object])
CommandDecorator = Callable[..., Callable[[T], T]]
app_command: CommandDecorator = app.command


def filter_transcripts(transcripts, include_types=None, exclude_types=None, canonical_only=False):  # type: ignore
    """Filter transcripts based on type and canonical status."""
    filtered = transcripts

    if canonical_only:
        filtered = [t for t in filtered if t.is_canonical]

    if include_types:
        filtered = [t for t in filtered if t.transcript_type in include_types]

    if exclude_types:
        filtered = [t for t in filtered if t.transcript_type not in exclude_types]

    return filtered


def extract_canonical_transcripts(transcripts, gene_name, output_dir=None):  # type: ignore
    """Extract canonical transcripts to a separate file."""
    canonical = [t for t in transcripts if t.is_canonical]

    if not canonical:
        return None, 0

    output_dir = Path.cwd() if output_dir is None else Path(output_dir)

    canonical_file = output_dir / f"{gene_name}_canonical.fasta"

    # Create a temporary searcher to use the save method
    # TODO: directly use fasta utils
    searcher = GeneSearcher()
    searcher.save_transcripts_fasta(canonical, canonical_file)

    return canonical_file, len(canonical)


@app_command()
def search(  # noqa: PLR0912
    query: str = typer.Argument(..., help="Gene ID, gene name, or transcript ID to search for"),
    output: Path = typer.Option(
        Path("transcripts.fasta"),
        "--output",
        "-o",
        help="Output FASTA file for transcript sequences",
    ),
    database: str = typer.Option(
        "ensembl",
        "--database",
        "-d",
        help="Database to search (ensembl, refseq, gencode)",
    ),
    all_databases: bool = typer.Option(
        False,
        "--all",
        "-a",
        help="Search all databases",
    ),
    fallback: bool = typer.Option(
        True,
        "--fallback/--no-fallback",
        help="Enable automatic fallback to other databases if access is blocked",
    ),
    no_sequence: bool = typer.Option(
        False,
        "--no-sequence",
        help="Skip sequence retrieval (metadata only)",
    ),
    canonical_only: bool = typer.Option(
        False,
        "--canonical-only",
        help="Extract only canonical isoforms",
    ),
    extract_canonical: bool = typer.Option(
        True,
        "--extract-canonical/--no-extract-canonical",
        help="Automatically extract canonical isoforms to separate file",
    ),
    transcript_types: str = typer.Option(
        "protein_coding,lncRNA",
        "--types",
        "-t",
        help="Comma-separated list of transcript types to include (e.g., protein_coding,lncRNA)",
    ),
    exclude_types: str = typer.Option(
        "nonsense_mediated_decay,retained_intron",
        "--exclude-types",
        help="Comma-separated list of transcript types to exclude",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
) -> None:
    """Search for gene transcripts and retrieve sequences."""

    try:
        # imports moved to top

        # Validate database choice
        try:
            db_type = DatabaseType(database.lower())
        except ValueError:
            console.print(f"❌ [red]Invalid database:[/red] {database}")
            console.print("Valid options: ensembl, refseq, gencode")
            raise typer.Exit(1)

        console.print(
            Panel.fit(
                f"🧬 [bold blue]Gene Search[/bold blue]\n"
                f"Query: [cyan]{query}[/cyan]\n"
                f"Database: [yellow]{database}[/yellow]\n"
                f"Fallback: [green]{'enabled' if fallback else 'disabled'}[/green]\n"
                f"Output: [cyan]{output}[/cyan]\n"
                f"Types: [green]{transcript_types}[/green]\n"
                f"Exclude: [red]{exclude_types}[/red]",
                title="Search Configuration",
            )
        )

        # Parse transcript types
        include_types = [t.strip() for t in transcript_types.split(",") if t.strip()] if transcript_types else []
        exclude_types_list = [t.strip() for t in exclude_types.split(",") if t.strip()] if exclude_types else []

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            if all_databases:
                progress.add_task("Searching all databases...", total=None)
                results = search_multiple_databases_sync(
                    query=query, databases=list(DatabaseType), include_sequence=not no_sequence
                )
            elif fallback:
                progress.add_task("Searching with fallback...", total=None)
                results = [search_gene_with_fallback_sync(query=query, include_sequence=not no_sequence)]
            else:
                progress.add_task(f"Searching {database}...", total=None)
                results = [search_gene_sync(query=query, database=db_type, include_sequence=not no_sequence)]

        # Display results
        successful_results = [r for r in results if r.success]

        if not successful_results:
            console.print(f"❌ [red]No results found for:[/red] {query}")
            for result in results:
                if result.error:
                    db_name = result.database.value if hasattr(result.database, "value") else str(result.database)
                    console.print(f"  {db_name}: {result.error}")
            raise typer.Exit(1)

        # Apply filtering to all transcripts
        all_transcripts = []
        for result in successful_results:
            filtered_transcripts = filter_transcripts(
                result.transcripts,
                include_types=include_types,
                exclude_types=exclude_types_list,
                canonical_only=canonical_only,
            )
            result.transcripts = filtered_transcripts  # Update result with filtered transcripts
            all_transcripts.extend(filtered_transcripts)

        # Check if filtering removed all transcripts
        if not all_transcripts:
            console.print("❌ [red]No transcripts found after applying filters[/red]")
            console.print(f"  Include types: {include_types}")
            console.print(f"  Exclude types: {exclude_types_list}")
            raise typer.Exit(1)

        # Display summary table
        summary_table = Table(title="📊 Search Results")
        summary_table.add_column("Database", style="cyan")
        summary_table.add_column("Gene ID", style="blue")
        summary_table.add_column("Gene Name", style="green")
        summary_table.add_column("Transcripts", style="yellow")
        summary_table.add_column("Status", style="magenta")

        for result in results:
            if result.success:
                gene_id = result.gene_info.gene_id if result.gene_info else "N/A"
                gene_name = result.gene_info.gene_name if result.gene_info else "N/A"
                transcript_count = len(result.transcripts)
                status = "✅ Success"
            else:
                gene_id = "N/A"
                gene_name = "N/A"
                transcript_count = 0
                status = f"❌ {result.error}"

            db_name = result.database.value if hasattr(result.database, "value") else str(result.database)
            summary_table.add_row(db_name, gene_id, gene_name, str(transcript_count), status)

        console.print(summary_table)

        # Show transcript details for successful results
        if successful_results and not no_sequence:
            transcript_table = Table(title="🧬 Transcript Details")
            transcript_table.add_column("Transcript ID", style="cyan")
            transcript_table.add_column("Database", style="blue")
            transcript_table.add_column("Type", style="green")
            transcript_table.add_column("Length", style="yellow")
            transcript_table.add_column("Canonical", style="magenta")

            for transcript in all_transcripts[:10]:  # Show first 10
                db_name = (
                    transcript.database.value if hasattr(transcript.database, "value") else str(transcript.database)
                )
                transcript_table.add_row(
                    transcript.transcript_id,
                    db_name,
                    transcript.transcript_type or "N/A",
                    str(transcript.length) if transcript.length else "N/A",
                    "✓" if transcript.is_canonical else "",
                )

            console.print(transcript_table)

            if len(all_transcripts) > 10:
                console.print(f"... and {len(all_transcripts) - 10} more transcripts")

        # Save sequences to FASTA if requested
        if not no_sequence and all_transcripts:
            searcher = GeneSearcher()
            transcripts_with_sequence = [t for t in all_transcripts if t.sequence]

            if transcripts_with_sequence:
                searcher.save_transcripts_fasta(transcripts_with_sequence, output)
                console.print(f"\n✅ [green]Saved {len(transcripts_with_sequence)} sequences to:[/green] {output}")

                # Extract canonical isoforms if requested
                if extract_canonical and not canonical_only:
                    gene_name = None
                    for result in successful_results:
                        if result.gene_info and result.gene_info.gene_name:
                            gene_name = result.gene_info.gene_name
                            break

                    if gene_name:
                        canonical_file, canonical_count = extract_canonical_transcripts(
                            transcripts_with_sequence, gene_name, output.parent
                        )
                        if canonical_file:
                            console.print(
                                f"📌 [blue]Extracted {canonical_count} canonical isoform(s) to:[/blue] {canonical_file}"
                            )
                        else:
                            console.print("ℹ️  [yellow]No canonical isoforms found[/yellow]")
            else:
                console.print("⚠️  [yellow]No sequences available to save[/yellow]")

        console.print(
            f"\n📊 [blue]Summary:[/blue] Found {len(all_transcripts)} transcripts from {len(successful_results)} database(s)"
        )

    except Exception as e:
        console.print(f"❌ [red]Search error:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


@app_command()
def workflow(
    gene_query: str = typer.Argument(..., help="Gene name or ID to analyze"),
    input_fasta: Optional[str] = typer.Option(
        None,
        "--input-fasta",
        help="Local path or remote URI to an input FASTA file (http/https/ftp)",
    ),
    output_dir: Path = typer.Option(
        Path("sirna_workflow_output"),
        "--output-dir",
        "-o",
        help="Output directory for all workflow results",
    ),
    database: str = typer.Option(
        "ensembl",
        "--database",
        "-d",
        help="Database to search (ensembl, refseq, gencode)",
    ),
    top_n_candidates: int = typer.Option(
        20,
        "--top-n",
        "-n",
        min=1,
        help="Number of top siRNA candidates to select (also used for off-target analysis)",
    ),
    genome_species: str = typer.Option(
        "human,rat,rhesus",
        "--genome-species",
        help="Comma-separated list of genome species for off-target analysis",
    ),
    gc_min: float = typer.Option(
        30.0,
        "--gc-min",
        min=0.0,
        max=100.0,
        help="Minimum GC content percentage",
    ),
    gc_max: float = typer.Option(
        52.0,
        "--gc-max",
        min=0.0,
        max=100.0,
        help="Maximum GC content percentage",
    ),
    sirna_length: int = typer.Option(
        21,
        "--length",
        "-l",
        min=19,
        max=23,
        help="siRNA length in nucleotides",
    ),
    modification_pattern: str = typer.Option(
        "standard_2ome",
        "--modifications",
        "-m",
        help="Chemical modification pattern (standard_2ome, minimal_terminal, maximal_stability, none)",
    ),
    overhang: str = typer.Option(
        "dTdT",
        "--overhang",
        help="Overhang sequence (dTdT for DNA, UU for RNA)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    log_file: Optional[Path] = typer.Option(
        None,
        "--log-file",
        help="Path to centralized log file (overrides SIRNAFORGE_LOG_FILE env)",
    ),
    json_summary: bool = typer.Option(
        True,
        "--json-summary/--no-json-summary",
        help="Write logs/workflow_summary.json (disable to skip JSON output)",
    ),
) -> None:
    """Run complete siRNA design workflow from gene query to off-target analysis."""

    if gc_min >= gc_max:
        console.print("❌ Error: gc-min must be less than gc-max", style="red")
        raise typer.Exit(1)

    # Parse genome species
    species_list = [s.strip() for s in genome_species.split(",") if s.strip()]

    input_descriptor = gene_query
    if input_fasta:
        input_descriptor = input_fasta if "://" in input_fasta else Path(input_fasta).name

    console.print(
        Panel.fit(
            f"🧬 [bold blue]Complete siRNA Workflow[/bold blue]\n"
            f"Gene Query: [cyan]{input_descriptor}[/cyan]\n"
            f"Database: [yellow]{database}[/yellow]\n"
            f"Output Directory: [cyan]{output_dir}[/cyan]\n"
            f"siRNA Length: [yellow]{sirna_length}[/yellow] nt\n"
            f"GC Range: [yellow]{gc_min:.1f}%-{gc_max:.1f}%[/yellow]\n"
            f"Top Candidates (used for off-target): [yellow]{top_n_candidates}[/yellow]\n"
            f"Genome Species: [green]{', '.join(species_list)}[/green]\n"
            f"Modifications: [magenta]{modification_pattern}[/magenta]\n"
            f"Overhang: [magenta]{overhang}[/magenta]",
            title="Workflow Configuration",
        )
    )

    try:
        # Run the complete workflow
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task("Running complete siRNA design workflow...", total=None)

            # Configure logging to file inside the output dir by default if not provided
            effective_log = str(log_file) if log_file else str(Path(output_dir) / "logs" / "sirnaforge.log")
            configure_logging(log_file=effective_log, level=os.getenv("SIRNAFORGE_LOG_LEVEL"))

            results = asyncio.run(
                run_sirna_workflow(
                    gene_query=gene_query,
                    input_fasta=input_fasta,
                    output_dir=str(output_dir),
                    database=database,
                    top_n_candidates=top_n_candidates,
                    genome_species=species_list,
                    gc_min=gc_min,
                    gc_max=gc_max,
                    sirna_length=sirna_length,
                    modification_pattern=modification_pattern,
                    overhang=overhang,
                    log_file=effective_log,
                    write_json_summary=json_summary,
                )
            )

            progress.remove_task(task)
        # TODO: simplify the printing and console logging summaries
        # Display results summary
        console.print("\n✅ [bold green]Workflow completed successfully![/bold green]")

        # Workflow summary
        summary_table = Table(title="📊 Workflow Results Summary")
        summary_table.add_column("Phase", style="cyan")
        summary_table.add_column("Status", style="green")
        summary_table.add_column("Details", style="white")

        results.get("workflow_config", {})
        transcript_summary = results.get("transcript_summary", {})
        design_summary = results.get("design_summary", {})
        offtarget_summary = results.get("offtarget_summary", {})

        summary_table.add_row(
            "📄 Transcript Retrieval",
            "✅ Complete",
            f"{transcript_summary.get('total_transcripts', 0)} transcripts from {database}",
        )

        summary_table.add_row(
            "🧬 siRNAforge", "✅ Complete", f"{design_summary.get('total_candidates', 0)} candidates generated"
        )

        summary_table.add_row(
            "🎯 Off-target Analysis",
            "✅ Complete" if offtarget_summary.get("status") == "completed" else "⚠️  Partial",
            f"Method: {offtarget_summary.get('method', 'basic')}",
        )

        console.print(summary_table)

        # Output locations
        console.print(f"\n📁 [bold]Results saved to:[/bold] [cyan]{output_dir}[/cyan]")
        console.print("📂 Key files:")
        console.print(f"   • Transcripts: [blue]transcripts/{gene_query}_transcripts.fasta[/blue]")
        console.print(f"   • siRNA candidates (ALL): [blue]sirnaforge/{gene_query}_all.csv[/blue]")
        console.print(f"   • siRNA candidates (PASS): [blue]sirnaforge/{gene_query}_pass.csv[/blue]")
        console.print("   • Off-target results: [blue]off_target/results/[/blue]")
        console.print("   • Console stream log: [blue]logs/workflow_stream.log[/blue]")
        if json_summary:
            console.print("   • Workflow summary: [blue]logs/workflow_summary.json[/blue]")

        if offtarget_summary.get("method") == "nextflow":
            console.print("   • Full off-target report: [blue]off_target/results/offtarget_report.html[/blue]")

    except Exception as e:
        console.print(f"❌ [red]Workflow error:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


@app_command()
def design(
    input_file: Path = typer.Argument(
        ...,
        help="Input FASTA file containing transcript sequences",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    output: Path = typer.Option(
        Path("sirna_results.tsv"),
        "--output",
        "-o",
        help="Output file for siRNA candidates",
    ),
    length: int = typer.Option(
        21,
        "--length",
        "-l",
        min=19,
        max=23,
        help="siRNA length in nucleotides",
    ),
    top_n: int = typer.Option(
        10,
        "--top-n",
        "-n",
        min=1,
        help="Number of top candidates to return",
    ),
    gc_min: float = typer.Option(
        30.0,
        "--gc-min",
        min=0.0,
        max=100.0,
        help="Minimum GC content percentage",
    ),
    gc_max: float = typer.Option(
        52.0,
        "--gc-max",
        min=0.0,
        max=100.0,
        help="Maximum GC content percentage",
    ),
    max_poly_runs: int = typer.Option(
        3,
        "--max-poly-runs",
        min=1,
        help="Maximum consecutive identical nucleotides",
    ),
    genome_index: Optional[Path] = typer.Option(
        None,
        "--genome-index",
        help="Genome index for off-target analysis",
    ),
    snp_file: Optional[Path] = typer.Option(
        None,
        "--snp-file",
        help="VCF file with SNPs to avoid",
    ),
    skip_structure: bool = typer.Option(
        False,
        "--skip-structure",
        help="Skip secondary structure prediction (faster)",
    ),
    skip_off_targets: bool = typer.Option(
        False,
        "--skip-off-targets",
        help="Skip off-target analysis (faster)",
    ),
    modification_pattern: str = typer.Option(
        "standard_2ome",
        "--modifications",
        "-m",
        help="Chemical modification pattern (standard_2ome, minimal_terminal, maximal_stability, none)",
    ),
    overhang: str = typer.Option(
        "dTdT",
        "--overhang",
        help="Overhang sequence (dTdT for DNA, UU for RNA)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
) -> None:
    """Design siRNA candidates from transcript sequences."""

    if gc_min >= gc_max:
        console.print("❌ Error: gc-min must be less than gc-max", style="red")
        raise typer.Exit(1)

    # Create parameters
    filters = FilterCriteria(
        gc_min=gc_min,
        gc_max=gc_max,
        max_poly_runs=max_poly_runs,
    )

    parameters = DesignParameters(
        sirna_length=length,
        top_n=top_n,
        filters=filters,
        predict_structure=not skip_structure,
        check_off_targets=not skip_off_targets,
        genome_index=str(genome_index) if genome_index else None,
        snp_file=str(snp_file) if snp_file else None,
        apply_modifications=modification_pattern.lower() != "none",
        modification_pattern=modification_pattern,
        default_overhang=overhang,
    )

    console.print(
        Panel.fit(
            f"🧬 [bold blue]siRNAforge Toolkit[/bold blue]\n"
            f"Input: [cyan]{input_file}[/cyan]\n"
            f"Output: [cyan]{output}[/cyan]\n"
            f"Length: [yellow]{length}[/yellow] nt\n"
            f"GC range: [yellow]{gc_min:.1f}%-{gc_max:.1f}%[/yellow]\n"
            f"Top candidates: [yellow]{top_n}[/yellow]\n"
            f"Modifications: [magenta]{modification_pattern}[/magenta]\n"
            f"Overhang: [magenta]{overhang}[/magenta]",
            title="Configuration",
        )
    )

    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            # Import here to avoid slow startup

            task1 = progress.add_task("Loading sequences...", total=None)
            designer = SiRNADesigner(parameters)

            progress.update(task1, description="Designing siRNAs...")
            result = designer.design_from_file(str(input_file))

            # Apply chemical modifications if enabled
            if parameters.apply_modifications:
                from sirnaforge.utils.modification_patterns import apply_modifications_to_candidate  # noqa: PLC0415

                progress.update(task1, description="Applying modifications...")
                for candidate in result.candidates:
                    apply_modifications_to_candidate(
                        candidate,
                        pattern_name=parameters.modification_pattern,
                        overhang=parameters.default_overhang,
                    )

            progress.update(task1, description="Saving results...")
            result.save_csv(str(output))

        # Display results summary
        summary = result.get_summary()

        table = Table(title="📊 Design Summary")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="yellow")

        for key, value in summary.items():
            table.add_row(key.replace("_", " ").title(), str(value))

        console.print(table)

        if result.top_candidates:
            console.print("\n🏆 [bold green]Top Candidates:[/bold green]")

            candidates_table = Table()
            candidates_table.add_column("ID", style="cyan")
            candidates_table.add_column("Transcript", style="blue")
            candidates_table.add_column("Position", style="yellow")
            candidates_table.add_column("Sequence", style="green")
            candidates_table.add_column("GC%", style="magenta")
            candidates_table.add_column("Hits", style="white")
            candidates_table.add_column("Hit %", style="white")
            candidates_table.add_column("Score", style="red")

            for candidate in result.top_candidates[:5]:  # Show top 5
                candidates_table.add_row(
                    candidate.id,
                    candidate.transcript_id,
                    str(candidate.position),
                    candidate.guide_sequence,
                    f"{candidate.gc_content:.1f}",
                    str(candidate.transcript_hit_count),
                    f"{candidate.transcript_hit_fraction * 100:.1f}%",
                    f"{candidate.composite_score:.1f}",
                )

            console.print(candidates_table)

        console.print(f"\n✅ [green]Results saved to:[/green] {output}")

    except Exception as e:
        console.print(f"❌ [red]Error during design:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


@app_command()
def validate(
    input_file: Path = typer.Argument(
        ...,
        help="FASTA file to validate",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
) -> None:
    """Validate input FASTA file format and content."""

    try:
        with console.status("Validating FASTA file..."):
            sequences = list(SeqIO.parse(input_file, "fasta"))

        if not sequences:
            console.print("❌ [red]No sequences found in FASTA file[/red]")
            raise typer.Exit(1)

        # Validation stats
        total_seqs = len(sequences)
        total_length = sum(len(seq) for seq in sequences)
        min_length = min(len(seq) for seq in sequences)
        max_length = max(len(seq) for seq in sequences)
        avg_length = total_length / total_seqs

        # Check for problematic sequences
        short_seqs = [seq for seq in sequences if len(seq) < 50]
        ambiguous_seqs = [seq for seq in sequences if "N" in str(seq.seq)]

        table = Table(title="📋 FASTA Validation Results")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="yellow")

        table.add_row("Total sequences", str(total_seqs))
        table.add_row("Total length", f"{total_length:,} nt")
        table.add_row("Average length", f"{avg_length:.1f} nt")
        table.add_row("Min length", f"{min_length} nt")
        table.add_row("Max length", f"{max_length} nt")
        table.add_row("Short sequences (<50 nt)", str(len(short_seqs)))
        table.add_row("Ambiguous sequences (with N)", str(len(ambiguous_seqs)))

        console.print(table)

        if short_seqs:
            console.print(f"⚠️  [yellow]{len(short_seqs)} sequences are shorter than 50 nt[/yellow]")

        if ambiguous_seqs:
            console.print(f"⚠️  [yellow]{len(ambiguous_seqs)} sequences contain ambiguous bases (N)[/yellow]")

        console.print("✅ [green]FASTA validation complete[/green]")

    except Exception as e:
        console.print(f"❌ [red]Validation error:[/red] {str(e)}")
        raise typer.Exit(1)


@app_command()
def version() -> None:
    """Show version information."""

    try:
        # Prefer Docker build-time APP_VERSION when the image is built with a VERSION arg
        app_version = os.environ.get("APP_VERSION") if "APP_VERSION" in os.environ else __version__
        console.print(
            Panel.fit(
                f"🧬 [bold blue]siRNAforge Toolkit[/bold blue]\n"
                f"Version: [yellow]{app_version}[/yellow]\n"
                f"Author: [cyan]{__author__}[/cyan]",
                title="Version Info",
            )
        )

    except ImportError:
        console.print("❌ [red]Could not determine version[/red]")
        raise typer.Exit(1)


@app_command()
def config() -> None:
    """Show default configuration parameters."""

    default_params = DesignParameters()

    console.print("[bold blue]Default Design Parameters:[/bold blue]\n")

    # Basic parameters
    console.print("[cyan]Basic Parameters:[/cyan]")
    console.print(f"  siRNA length: {default_params.sirna_length} nt")
    console.print(f"  Top candidates: {default_params.top_n}")

    # Filtering criteria
    console.print("\n[cyan]Filtering Criteria:[/cyan]")
    filters = default_params.filters
    console.print(f"  GC content: {filters.gc_min}% - {filters.gc_max}%")
    console.print(f"  Max poly runs: {filters.max_poly_runs}")
    console.print(f"  Max paired fraction: {filters.max_paired_fraction}")

    # Scoring weights
    console.print("\n[cyan]Scoring Weights:[/cyan]")
    scoring = default_params.scoring
    console.print(f"  Asymmetry: {scoring.asymmetry}")
    console.print(f"  GC content: {scoring.gc_content}")
    console.print(f"  Accessibility: {scoring.accessibility}")
    console.print(f"  Off-target: {scoring.off_target}")
    console.print(f"  Empirical: {scoring.empirical}")


@app_command()
def cache(
    clear: bool = typer.Option(False, "--clear", help="Clear all cached miRNA databases"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show what would be deleted without actually deleting"),
    info: bool = typer.Option(False, "--info", help="Show cache information"),
) -> None:
    """Manage miRNA database cache."""
    if not any([clear, dry_run, info]):
        console.print("❓ [yellow]No action specified. Use --clear, --dry-run, or --info[/yellow]")
        console.print("   Example: sirnaforge cache --info")
        return

    manager = MiRNADatabaseManager()

    if info or dry_run:
        cache_info = manager.cache_info()
        console.print("📊 [bold blue]Cache Information:[/bold blue]")
        console.print(f"  Directory: [cyan]{cache_info['cache_directory']}[/cyan]")
        console.print(f"  Files: [green]{cache_info['total_files']}[/green]")
        console.print(f"  Size: [yellow]{cache_info['total_size_mb']:.2f} MB[/yellow]")
        console.print(f"  TTL: [magenta]{cache_info['cache_ttl_days']} days[/magenta]")

    if dry_run:
        result = manager.clear_cache(confirm=False)
        console.print("\n🔍 [bold yellow]Clear Preview:[/bold yellow]")
        console.print(f"  Files to delete: [red]{result['files_deleted']}[/red]")
        console.print(f"  Size to free: [yellow]{result['size_freed_mb']:.2f} MB[/yellow]")
        console.print(f"  Status: [dim]{result['status']}[/dim]")

    elif clear:
        result = manager.clear_cache(confirm=True)
        console.print("🧹 [bold green]Cache Cleared:[/bold green]")
        console.print(f"  Files deleted: [red]{result['files_deleted']}[/red]")
        console.print(f"  Size freed: [yellow]{result['size_freed_mb']:.2f} MB[/yellow]")
        console.print(f"  Status: [green]{result['status']}[/green]")


# Create sequences subcommand group
sequences_app = typer.Typer(help="Manage siRNA sequences and metadata")
app.add_typer(sequences_app, name="sequences")
sequences_command: CommandDecorator = sequences_app.command


class SequencesShowError(RuntimeError):
    """Custom error for sequence display operations."""


def _load_fasta_records(input_file: Path) -> list[SeqRecord]:
    records = list(SeqIO.parse(input_file, "fasta"))
    if not records:
        raise SequencesShowError("No sequences found in file")
    return records


def _filter_records_by_id(records: list[SeqRecord], sequence_id: str) -> list[SeqRecord]:
    filtered = [record for record in records if record.id == sequence_id]
    if not filtered:
        raise SequencesShowError(f"Sequence ID '{sequence_id}' not found")
    return filtered


def _metadata_value_to_json(value: Any) -> Any:
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if hasattr(value, "value"):
        return value.value
    if isinstance(value, list):
        return [_metadata_value_to_json(item) for item in value]
    return value


def _records_to_json(records: list[SeqRecord]) -> str:
    payload = []
    for record in records:
        metadata = parse_header(record)
        payload.append({key: _metadata_value_to_json(val) for key, val in metadata.items()})
    return json.dumps(payload, indent=2)


def _summarize_modifications(metadata: dict[str, Any]) -> str:
    mods = metadata.get("chem_mods") or []
    summary = []
    for mod in mods:
        mod_type = getattr(mod, "type", str(mod))
        positions = getattr(mod, "positions", [])
        length = len(positions) if isinstance(positions, list) else positions
        summary.append(f"{mod_type}({length})")
    return ", ".join(summary)


def _print_records_fasta(records: list[SeqRecord]) -> None:
    for record in records:
        console.print(f">{record.description}")
        console.print(str(record.seq))


def _print_records_table(records: list[SeqRecord], input_file: Path) -> None:
    table = Table(title=f"📋 Sequences from {input_file.name}")
    table.add_column("ID", style="cyan")
    table.add_column("Sequence", style="green")
    table.add_column("Length", style="yellow")
    table.add_column("Target", style="blue")
    table.add_column("Role", style="magenta")
    table.add_column("Modifications", style="white")

    for record in records:
        metadata = parse_header(record)
        sequence = str(record.seq)
        role = metadata.get("strand_role")
        if isinstance(role, Enum):
            role_display = role.value
        elif isinstance(role, str):
            role_display = role
        else:
            role_display = ""
        mods_summary = _summarize_modifications(metadata)

        table.add_row(
            metadata.get("id", record.id),
            f"{sequence[:30]}..." if len(sequence) > 30 else sequence,
            str(len(sequence)),
            metadata.get("target_gene", ""),
            role_display,
            mods_summary,
        )

    console.print(table)
    console.print(f"\n📊 Total sequences: {len(records)}")


@sequences_command("show")
def sequences_show(
    input_file: Path = typer.Argument(
        ...,
        help="FASTA file to display",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    sequence_id: Optional[str] = typer.Option(
        None,
        "--id",
        help="Show only this sequence ID",
    ),
    format: str = typer.Option(
        "table",
        "--format",
        "-f",
        help="Output format (table, json, fasta)",
    ),
) -> None:
    """Show sequences with their metadata from FASTA file."""
    format_normalized = format.lower()
    try:
        records = _load_fasta_records(input_file)
        if sequence_id:
            records = _filter_records_by_id(records, sequence_id)

        format_handlers: dict[str, Callable[[list[SeqRecord]], None]] = {
            "json": lambda seqs: console.print(_records_to_json(seqs)),
            "fasta": _print_records_fasta,
            "table": lambda seqs: _print_records_table(seqs, input_file),
        }

        handler = format_handlers.get(format_normalized)
        if handler is None:
            raise SequencesShowError("Unsupported format. Choose from table, json, or fasta.")

        handler(records)

    except SequencesShowError as exc:
        console.print(f"❌ [red]{exc}[/red]")
        raise typer.Exit(1) from exc
    except Exception as exc:
        console.print(f"❌ [red]Error:[/red] {exc}")
        raise typer.Exit(1) from exc


@sequences_command("annotate")
def sequences_annotate(
    input_fasta: Path = typer.Argument(
        ...,
        help="Input FASTA file",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    metadata_json: Path = typer.Argument(
        ...,
        help="JSON file with metadata",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output FASTA file (default: <input>_annotated.fasta)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
) -> None:
    """Merge metadata from JSON into FASTA headers."""
    try:
        # Determine output path
        if output is None:
            output = input_fasta.parent / f"{input_fasta.stem}_annotated.fasta"

        output_path = output

        console.print(
            Panel.fit(
                f"🧬 [bold blue]Annotate Sequences[/bold blue]\n"
                f"Input FASTA: [cyan]{input_fasta}[/cyan]\n"
                f"Metadata JSON: [yellow]{metadata_json}[/yellow]\n"
                f"Output: [green]{output_path}[/green]",
                title="Configuration",
            )
        )

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            progress.add_task("Merging metadata into FASTA...", total=None)
            updated_count = merge_metadata_into_fasta(input_fasta, metadata_json, output_path)

        console.print("\n✅ [green]Success![/green]")
        console.print(f"   Updated {updated_count} sequences with metadata")
        console.print(f"   Output saved to: [cyan]{output_path}[/cyan]")

    except Exception as e:
        console.print(f"❌ [red]Error:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
