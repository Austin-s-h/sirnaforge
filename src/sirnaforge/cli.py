"""Modern CLI for siRNAforge using Typer and Rich."""

# Configure environment for ASCII compatibility before importing typer/rich
import os

os.environ.setdefault("FORCE_COLOR", "0")
os.environ.setdefault("NO_COLOR", "1")
os.environ.setdefault("TERM", "dumb")

import asyncio
import json
import logging
from collections.abc import Callable, Iterable, Mapping
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, Any, NoReturn, Protocol, TypeVar, cast

import typer
from Bio.SeqIO import parse as seqio_parse
from Bio.SeqRecord import SeqRecord
from click.core import ParameterSource
from pydantic import ValidationError
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

# Monkey patch Rich console (best-effort). If this fails we continue with default behavior.
try:
    import rich.console

    original_init = rich.console.Console.__init__

    def patched_init(self: "rich.console.Console", *args: Any, **kwargs: Any) -> None:
        """Force simplified terminal capabilities for deterministic CI output."""
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
from sirnaforge.config import (
    DEFAULT_MIRNA_CANONICAL_SPECIES,
    DEFAULT_TRANSCRIPTOME_SOURCES,
    ReferencePolicyResolver,
    WorkflowInputSpec,
    render_reference_selection_label,
)
from sirnaforge.config.run_policy import (
    EntryPoint,
    ResolvedRunPolicy,
    RunPolicyError,
    default_for,
    format_validation_error,
    mirna_preset_default_for,
    resolve_run_policy,
    switchable_filter_ids,
)
from sirnaforge.core.design import SiRNADesigner
from sirnaforge.data.base import DatabaseType, FastaUtils, TranscriptInfo
from sirnaforge.data.gene_search import (
    GeneSearcher,
    GeneSearchResult,
    search_gene_sync,
    search_gene_with_fallback_sync,
    search_multiple_databases_sync,
)
from sirnaforge.models.policy import RunMode
from sirnaforge.models.sirna import (
    DesignMode,
    ranking_score,
)
from sirnaforge.models.variant import VariantMode
from sirnaforge.models.zfn import (
    DimerMode,
    GenomicAnnotationConfig,
    ZFNAlgorithm,
    ZFNDefaultSubfingerMutationConstraint,
    ZFNDesignParameters,
    ZFNHalfSiteConstraints,
    ZFNMutationConstraints,
    ZFNMutationType,
    ZFNOverallMutationConstraint,
    ZFNSearchBackend,
    ZFNShardingConfig,
    ZFNSpacerConstraints,
    ZFNSubfingerMutationConstraint,
)
from sirnaforge.modifications import merge_metadata_into_fasta, parse_header
from sirnaforge.pipeline.nextflow.config import DEFAULT_SIRNAFORGE_DOCKER_IMAGE
from sirnaforge.reporting import ReportInputError, build_payload, write_report
from sirnaforge.utils.cli_inputs import extract_declared_species_from_indices, resolve_species_inputs
from sirnaforge.utils.logging_utils import configure_logging
from sirnaforge.utils.typed_decorators import command_decorator_typed
from sirnaforge.workflow import run_offtarget_only_workflow, run_sirna_workflow
from sirnaforge.zfn import emit_zfn_experimental_warning
from sirnaforge.zfn.nextflow_bridge import (
    aggregate_zfn_shard_results,
    make_zfn_shard_manifest,
    run_zfn_shard_search,
)
from sirnaforge.zfn.search import build_zfn_search_index

app = typer.Typer(
    name="sirnaforge",
    help="siRNAforge - siRNA design toolkit for gene silencing",
    rich_markup_mode="rich",
)
# Configure console to use ASCII box characters for better compatibility
console = Console(force_terminal=False, legacy_windows=True)

# mypy-friendly alias for Typer command decorator
app_command = command_decorator_typed(app.command)

DEFAULT_SPECIES_ARGUMENT = ",".join(DEFAULT_MIRNA_CANONICAL_SPECIES)
REMOTE_RESOURCE_SCHEMES = ("http://", "https://", "ftp://", "file://")

# Internal performance defaults for exhaustive ZFN search.
DEFAULT_ZFN_WINDOW_STRIDE = 1
DEFAULT_ZFN_TOP_N_SITES = 5000
DEFAULT_ZFN_REPORT_N_SITES = 200


def _offtarget_results_line(offtarget_summary: dict[str, Any], results_path: str) -> str:
    """Describe where off-target results landed, or why there are none.

    Reporting the path unconditionally meant a run whose screening never happened still
    pointed at a directory that does not exist, with the actual reason many lines earlier
    in the log.
    """
    if offtarget_summary.get("status") == "skipped":
        reason = offtarget_summary.get("reason") or "not run"
        return f"   • Off-target results: [yellow]not produced ({reason})[/yellow]"
    return f"   • Off-target results: [blue]{results_path}[/blue]"


def _autotune_zfn_sharding(
    cores_budget: int | None = None,
    search_backend: ZFNSearchBackend = ZFNSearchBackend.PYAHOCORASICK,
) -> ZFNShardingConfig:
    """Return internal sharding defaults tuned by backend and host CPU budget."""
    cpu_count = cores_budget if cores_budget is not None else (os.cpu_count() or 1)

    if search_backend == ZFNSearchBackend.FM_INDEX:
        # FM-index can over-fragment and over-allocate quickly on full-primary genomes.
        max_workers = min(4, max(1, cpu_count // 2 if cpu_count >= 4 else cpu_count))
        chunk_size_bp = 16_000_000 if cpu_count >= 8 else 20_000_000
    elif search_backend == ZFNSearchBackend.PYAHOCORASICK:
        # pyahocorasick is memory-sensitive on large fallback runs; keep a conservative
        # default worker count and treat workers primarily as a memory-control knob.
        # chunk_size_bp here acts as the large-contig fallback threshold, not the
        # primary scheduling unit.
        max_workers = min(2, max(1, cpu_count))
        chunk_size_bp = 8_000_000 if cpu_count >= 8 else 12_000_000
    else:
        # exhaustive_python remains baseline-oriented and follows the broader
        # CPU-parallel profile.
        max_workers = min(8, max(1, cpu_count))
        chunk_size_bp = 8_000_000 if cpu_count >= 8 else 12_000_000

    return ZFNShardingConfig(
        enabled=True,
        chunk_size_bp=chunk_size_bp,
        overlap_bp=50,
        chromosomes=[],
        max_workers=max_workers,
    )


class TranscriptLike(Protocol):
    """Minimal transcript-like interface used by CLI filters."""

    transcript_type: str | None
    is_canonical: bool


TTranscript = TypeVar("TTranscript", bound=TranscriptLike)


def _parse_fasta_seqrecords(input_file: Path) -> list[SeqRecord]:
    """Parse FASTA records with an explicit typed boundary for static checkers."""
    parse_fasta = cast(Callable[[str, str], Iterable[SeqRecord]], seqio_parse)
    return list(parse_fasta(str(input_file), "fasta"))


def filter_transcripts(
    transcripts: list[TTranscript],
    include_types: list[str] | None = None,
    exclude_types: list[str] | None = None,
    canonical_only: bool = False,
) -> list[TTranscript]:
    """Filter transcript records by type and canonical status.

    Args:
        transcripts: Iterable of transcript-like objects that expose
            ``transcript_type`` and ``is_canonical`` attributes.
        include_types: Optional iterable of transcript types to keep.
        exclude_types: Optional iterable of transcript types to drop.
        canonical_only: When True, keep only canonical isoforms.

    Returns:
        A list of transcripts that match the requested filters.
    """
    filtered: list[TTranscript] = transcripts

    if canonical_only:
        filtered = [t for t in filtered if t.is_canonical]

    if include_types:
        filtered = [t for t in filtered if t.transcript_type in include_types]

    if exclude_types:
        filtered = [t for t in filtered if t.transcript_type not in exclude_types]

    return filtered


def extract_canonical_transcripts(
    transcripts: list[TranscriptInfo],
    gene_name: str,
    output_dir: Path | str | None = None,
) -> tuple[Path | None, int]:
    """Write canonical isoforms to a separate FASTA file.

    Args:
        transcripts: Iterable of transcript-like objects (must expose
            ``is_canonical`` and sequence attributes used by the underlying
            save routine).
        gene_name: Name used to derive the output FASTA filename.
        output_dir: Directory to write the FASTA file into (defaults to CWD).

    Returns:
        A tuple of ``(canonical_fasta_path, count)`` where the path is None
        when no canonical isoforms are available.
    """
    canonical = [t for t in transcripts if t.is_canonical]

    if not canonical:
        return None, 0

    output_dir_path = Path.cwd() if output_dir is None else Path(output_dir)

    canonical_file = output_dir_path / f"{gene_name}_canonical.fasta"

    # Create a temporary searcher to use the save method
    # TODO: directly use fasta utils
    searcher = GeneSearcher()
    searcher.save_transcripts_fasta(canonical, canonical_file)

    return canonical_file, len(canonical)


def _option_was_stated(ctx: typer.Context, parameter: str) -> bool:
    """Whether the user actually supplied one option, read from Click's parameter source.

    This replaces comparing a value against the siRNA default. A value test cannot tell an omitted
    ``--gc-max`` from one typed as ``--gc-max 60``, so in miRNA mode an explicit 60 was silently
    rewritten to 52; the same sentinel governed ``--overhang`` and ``--modifications``. An option
    set through its environment variable counts as stated, because a user set it.
    """
    source = ctx.get_parameter_source(parameter)
    # Compared by name: Typer vendors its own ParameterSource, so an identity test against Click's
    # enum is a non-overlapping comparison that silently reads every option as unstated.
    return source is not None and source.name != ParameterSource.DEFAULT.name


def _stated_settings(ctx: typer.Context, candidates: Mapping[str, tuple[str, Any]]) -> dict[str, Any]:
    """Policy settings the user stated, keyed by setting name.

    ``candidates`` maps a policy setting key to ``(cli_parameter_name, value)``. The parameter name
    is carried explicitly because six options are not named after the setting they set (``--length``
    sets ``sirna_length``, ``--max-off-targets`` sets ``max_off_target_count``, ...) and looking the
    wrong name up in Click reports every one of them as unstated.
    """
    stated: dict[str, Any] = {}
    for key, (parameter, value) in candidates.items():
        if _option_was_stated(ctx, parameter) and value is not None:
            stated[key] = value
    return stated


def _parse_filter_actions(entries: Iterable[str]) -> dict[str, str]:
    """Parse repeatable ``--filter-action filter_id=off|warn|fail`` entries.

    Splitting only: the action itself is validated by
    :data:`~sirnaforge.config.run_policy.SELECTABLE_ACTIONS`, which accepts all three words and
    refuses an action the addressed gate could not apply -- ``off`` on a design-stage float
    threshold, or anything but ``off`` on a gate with no threshold, no screening evidence in this
    run, or no 0.7.1 reader at all.

    ``warn`` is selectable, and the resolver accepting it is not the same as every gate honouring
    it: four design-stage gates (``gc_content_min``, ``gc_content_max``, ``max_poly_runs``,
    ``max_repeat_transcript_fraction``) reject regardless, because their rejection is written
    outside ``record_filter_verdict``. The per-command ``--filter-action`` help names them.
    """
    actions: dict[str, str] = {}
    for raw in entries:
        token = raw.strip()
        if not token:
            continue
        if "=" not in token:
            raise RunPolicyError(f"--filter-action expects filter_id=action, got {raw!r}")
        filter_id, action = (part.strip() for part in token.split("=", 1))
        actions[filter_id] = action
    return actions


def _fail_with_config_error(message: str, *, logger: logging.Logger | None = None) -> NoReturn:
    """Report a configuration problem as a message and exit 1, never as a traceback.

    Cross-field Pydantic errors (``--plfold-max-bp-span 100 --plfold-window 40``) are raised when
    the model is constructed rather than when the option is parsed, so Typer's own validation never
    saw them and they reached the user as an unhandled ValidationError.
    """
    if logger is not None:
        logger.error("Invalid configuration: %s", message)
    console.print(f"❌ Error: {message}", style="red")
    raise typer.Exit(1)


def _resolve_policy_or_exit(
    *,
    entry_point: EntryPoint,
    logger: logging.Logger | None = None,
    **kwargs: Any,
) -> ResolvedRunPolicy:
    """Resolve the run policy, turning every configuration error into a message and exit 1."""
    try:
        return resolve_run_policy(entry_point=entry_point, **kwargs)
    except RunPolicyError as exc:
        _fail_with_config_error(str(exc), logger=logger)
    except ValidationError as exc:
        _fail_with_config_error(format_validation_error(exc), logger=logger)


def _parse_zfn_mutation_types(raw_types: str, raw_constraint: str) -> list[ZFNMutationType]:
    """Parse and normalize mutation type tokens."""
    aliases = {"mismatch": "substitution", "mismatches": "substitution"}
    mutation_types = [t.strip().lower() for t in raw_types.split(",") if t.strip()]
    if not mutation_types:
        raise ValueError(f"Invalid ZFN mutation types in '{raw_constraint}'. At least one type is required.")
    return [ZFNMutationType(aliases.get(value, value)) for value in mutation_types]


def _parse_zfn_mutation_constraints(
    raw_constraints: list[str],
) -> tuple[
    list[ZFNSubfingerMutationConstraint],
    ZFNDefaultSubfingerMutationConstraint | None,
    list[ZFNOverallMutationConstraint],
]:
    """Parse CLI ZFN mutation constraints in ``scope:max:type1,type2`` format.

    Scope can be:
      - ``<int>`` for explicit sub-finger index
      - ``*`` for default per-sub-finger budget
      - ``overall`` for global mutation budgets
    """
    per_subfinger: list[ZFNSubfingerMutationConstraint] = []
    default_subfinger: ZFNDefaultSubfingerMutationConstraint | None = None
    overall_constraints: list[ZFNOverallMutationConstraint] = []
    for raw in raw_constraints:
        parts = [part.strip() for part in raw.split(":", maxsplit=2)]
        if len(parts) != 3:
            raise ValueError(
                f"Invalid ZFN sub-finger mutation constraint '{raw}'. "
                "Expected format: scope:max_mutations:type1,type2 "
                "(scope: subfinger index, '*', or 'overall')"
            )

        scope_raw, max_mut_raw, types_raw = parts
        try:
            max_mutations = int(max_mut_raw)
        except ValueError as exc:
            raise ValueError(f"Invalid ZFN mutation counts in '{raw}'. max_mutations must be an integer.") from exc

        mutation_types = _parse_zfn_mutation_types(types_raw, raw)
        scope_normalized = scope_raw.lower()

        if scope_normalized == "*":
            default_subfinger = ZFNDefaultSubfingerMutationConstraint(
                max_mutations=max_mutations,
                mutation_types=mutation_types,
            )
            continue

        if scope_normalized == "overall":
            overall_constraints.append(
                ZFNOverallMutationConstraint(
                    max_mutations=max_mutations,
                    mutation_types=mutation_types,
                )
            )
            continue

        try:
            subfinger_index = int(scope_raw)
        except ValueError as exc:
            raise ValueError(
                f"Invalid ZFN mutation scope '{scope_raw}' in '{raw}'. "
                "Use subfinger index, '*' (default per-subfinger), or 'overall'."
            ) from exc

        per_subfinger.append(
            ZFNSubfingerMutationConstraint(
                subfinger_index=subfinger_index,
                max_mutations=max_mutations,
                mutation_types=mutation_types,
            )
        )

    return per_subfinger, default_subfinger, overall_constraints


def _build_zfn_design_configuration(  # noqa: PLR0912
    *,
    zfn_subfinger_mutation: list[str],
    zfn_max_mismatches_per_subfinger: int | None,
    zfn_max_substitutions_overall: int | None,
    zfn_left_half_site: str,
    zfn_right_half_site: str,
    zfn_search_space: str | None,
    zfn_search_space_index: str | None,
    zfn_search_backend: ZFNSearchBackend,
    zfn_algorithm: ZFNAlgorithm,
    zfn_dimer_mode: DimerMode,
    zfn_spacer_lengths: str,
    zfn_max_mismatches: int,
    zfn_window_stride: int | None = None,
    zfn_top_n_sites: int | None = None,
    zfn_report_n_sites: int | None = None,
    workflow_cores: int | None = None,
    zfn_annotation: str | None = None,
) -> tuple[
    ZFNDesignParameters,
    GenomicAnnotationConfig | None,
    list[ZFNSubfingerMutationConstraint],
    ZFNDefaultSubfingerMutationConstraint | None,
    list[ZFNOverallMutationConstraint],
]:
    """Parse and validate CLI options into a typed ZFN design configuration."""
    merged_zfn_constraints = list(zfn_subfinger_mutation)
    if zfn_max_mismatches_per_subfinger is not None:
        merged_zfn_constraints.append(f"*:{zfn_max_mismatches_per_subfinger}:mismatch")
    if zfn_max_substitutions_overall is not None:
        merged_zfn_constraints.append(f"overall:{zfn_max_substitutions_overall}:substitution")

    zfn_constraints, zfn_default_constraint, zfn_overall_constraints = _parse_zfn_mutation_constraints(
        merged_zfn_constraints
    )

    try:
        parsed_spacers = [int(s.strip()) for s in zfn_spacer_lengths.split(",")]
    except ValueError as exc:
        raise ValueError("--zfn-spacer-lengths must be comma-separated integers") from exc
    mutation_constraints: ZFNMutationConstraints | None = None
    if zfn_constraints or zfn_default_constraint or zfn_overall_constraints:
        mutation_constraints = ZFNMutationConstraints(
            subfinger_mutations=zfn_constraints,
            default_subfinger_mutation=zfn_default_constraint,
            overall_mutations=zfn_overall_constraints,
        )

    search_space_fasta: str | None = None
    search_space_reference: str | None = "ensembl_human_hg38_primary"
    if zfn_search_space:
        if Path(zfn_search_space).exists() or zfn_search_space.startswith(REMOTE_RESOURCE_SCHEMES):
            search_space_fasta = zfn_search_space
            search_space_reference = None
        else:
            search_space_reference = zfn_search_space

    annotation: GenomicAnnotationConfig | None = None
    if zfn_annotation:
        if zfn_annotation.startswith(REMOTE_RESOURCE_SCHEMES):
            annotation = GenomicAnnotationConfig(annotation_reference=zfn_annotation)
        else:
            annotation = GenomicAnnotationConfig(annotation_path=zfn_annotation)

    zfn_design_params = ZFNDesignParameters(
        left_half_site=zfn_left_half_site,
        right_half_site=zfn_right_half_site,
        search_space_reference=search_space_reference,
        search_space_fasta=search_space_fasta,
        search_space_index=zfn_search_space_index,
        search_backend=zfn_search_backend,
        algorithm=zfn_algorithm,
        dimer_mode=zfn_dimer_mode,
        spacer_constraints=ZFNSpacerConstraints(allowed_spacer_lengths=parsed_spacers),
        half_site_constraints=ZFNHalfSiteConstraints(
            max_mismatches=zfn_max_mismatches,
            window_stride=zfn_window_stride or DEFAULT_ZFN_WINDOW_STRIDE,
        ),
        top_n_sites=zfn_top_n_sites or DEFAULT_ZFN_TOP_N_SITES,
        report_n_sites=zfn_report_n_sites or DEFAULT_ZFN_REPORT_N_SITES,
        mutation_constraints=mutation_constraints,
        sharding=_autotune_zfn_sharding(workflow_cores, zfn_search_backend),
    )

    return zfn_design_params, annotation, zfn_constraints, zfn_default_constraint, zfn_overall_constraints


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
    """Search transcript references and optionally fetch sequences.

    This command queries Ensembl/RefSeq/Gencode (depending on flags) for a gene
    or transcript identifier. When sequences are fetched, it writes them to a
    FASTA file and can optionally also emit a canonical-only FASTA.
    """
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

        results: list[GeneSearchResult]
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
        successful_results: list[GeneSearchResult] = [r for r in results if r.success]

        if not successful_results:
            console.print(f"❌ [red]No results found for:[/red] {query}")
            for result in results:
                if result.error:
                    db_name = result.database.value
                    console.print(f"  {db_name}: {result.error}")
            raise typer.Exit(1)

        # Apply filtering to all transcripts
        all_transcripts: list[TranscriptInfo] = []
        for result in successful_results:
            filtered_transcripts: list[TranscriptInfo] = filter_transcripts(
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

            db_name = result.database.value
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
                db_name = transcript.database.value
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
            transcripts_with_sequence: list[TranscriptInfo] = [t for t in all_transcripts if t.sequence]

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
def workflow(  # noqa: PLR0912
    ctx: typer.Context,
    gene_query: str = typer.Argument(..., help="Gene name or ID to analyze"),
    run_mode: str | None = typer.Option(
        None,
        "--run-mode",
        help=(
            "How much screening evidence this run claims: design_only, exploratory or qualified. "
            "Defaults to qualified for this command; --skip-off-targets maps to design_only. "
            "Independent of --design-mode."
        ),
    ),
    policy_config: Path | None = typer.Option(
        None,
        "--policy-config",
        help=(
            "JSON or TOML file of policy settings. Beats the built-in profile and loses to options "
            "given on the command line. Unknown setting names are rejected, not ignored."
        ),
    ),
    filter_action: list[str] = typer.Option(
        [],
        "--filter-action",
        help=(
            "Set one filter's action: filter_id=off|warn|fail (repeatable). 'off' clears the gate's "
            "threshold, so it is not evaluated -- which is not the same as passing -- and works for "
            f"these gates: {', '.join(switchable_filter_ids())}. The seven design-stage gates read a "
            "threshold with no absent state and are refused rather than faked; widen the threshold "
            "instead. 'warn' keeps the gate measuring and stops it rejecting: the outcome is "
            "published in <filter_id>_verdict and <filter_id>_observed and passes_filters is left "
            "alone, so a 'fail' verdict on a retained candidate is a finding, not a contradiction. "
            "It is honoured by max_paired_fraction, min_asymmetry_score, min_empirical_score and "
            "every post-screen gate. KNOWN GAP: gc_content_min, gc_content_max and max_poly_runs "
            "decide during enumeration and drop the candidate from the design output before the "
            "action is read, and max_repeat_transcript_fraction stamps REPEAT_ELEMENT without "
            "reading it -- those four resolve to 'warn' in the manifest and still reject, so treat "
            "warn on them as unimplemented rather than as a label. An action can only turn a gate "
            "off, never on: 'warn' and 'fail' are refused for a gate that would be off anyway (no "
            "threshold, no screening evidence, or no 0.7.1 reader -- max_mirna_1mm_seed)."
        ),
    ),
    input_fasta: str | None = typer.Option(
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
    design_mode: str = typer.Option(
        "sirna",
        "--design-mode",
        help="Design mode: sirna (default), mirna (miRNA-biogenesis-aware), or zfn (EXPERIMENTAL)",
    ),
    zfn_subfinger_mutation: list[str] = typer.Option(
        [],
        "--zfn-subfinger-mutation",
        help=(
            "ZFN sub-finger mutation allowance. "
            "Repeatable format: scope:max_mutations:type1,type2. "
            "scope can be subfinger index (e.g. 2), '*' for default per-subfinger, "
            "or 'overall' for global budgets. "
            "Use 'mismatch' as a shorthand alias for 'substitution'."
        ),
    ),
    zfn_max_mismatches_per_subfinger: int | None = typer.Option(
        None,
        "--zfn-max-mismatches-per-subfinger",
        min=0,
        help="Convenience option equivalent to --zfn-subfinger-mutation '*:<N>:mismatch'.",
    ),
    zfn_max_substitutions_overall: int | None = typer.Option(
        None,
        "--zfn-max-substitutions-overall",
        min=0,
        help="Convenience option equivalent to --zfn-subfinger-mutation 'overall:<N>:substitution'.",
    ),
    # ── ZFN half-site and search-space inputs (required when --design-mode zfn) ──
    zfn_left_half_site: str | None = typer.Option(
        None,
        "--zfn-left-half-site",
        help="Left ZFN half-site sequence (9-18 bp, IUPAC allowed). Required for --design-mode zfn.",
    ),
    zfn_right_half_site: str | None = typer.Option(
        None,
        "--zfn-right-half-site",
        help="Right ZFN half-site sequence (9-18 bp, IUPAC allowed). Required for --design-mode zfn.",
    ),
    zfn_search_space: str | None = typer.Option(
        None,
        "--zfn-search-space",
        help=(
            "Genome reference key or local FASTA path for ZFN off-target search space. "
            "Built-in keys: ensembl_human_hg38_primary, ensembl_mouse_grcm39_primary, "
            "ensembl_rat_grcr8_toplevel, ensembl_macaque_mmul10_toplevel. "
            "Default: ensembl_human_hg38_primary when --design-mode zfn."
        ),
    ),
    zfn_search_space_index: str | None = typer.Option(
        None,
        "--zfn-search-space-index",
        help=(
            "Optional persisted search-space index bundle path for indexed ZFN backends "
            "(currently fm_index; fm_index is experimental on large references)."
        ),
    ),
    zfn_search_backend: ZFNSearchBackend = typer.Option(
        ZFNSearchBackend.PYAHOCORASICK,
        "--zfn-search-backend",
        help=(
            "Half-site search backend: pyahocorasick (default), "
            "exhaustive_python (baseline), or fm_index (experimental)."
        ),
    ),
    zfn_algorithm: ZFNAlgorithm = typer.Option(
        ZFNAlgorithm.ZFN_V2,
        "--zfn-algorithm",
        help="ZFN off-target scoring algorithm: homology, conserved_g, or zfn_v2 (default).",
    ),
    zfn_dimer_mode: DimerMode = typer.Option(
        DimerMode.HETERODIMER_ONLY,
        "--zfn-dimer-mode",
        help="Dimer mode: heterodimer_only (default) or include_homodimers.",
    ),
    zfn_spacer_lengths: str = typer.Option(
        "5,6,7",
        "--zfn-spacer-lengths",
        help="Comma-separated allowed spacer lengths between half-sites (default: 5,6,7).",
    ),
    zfn_max_mismatches: int = typer.Option(
        2,
        "--zfn-max-mismatches",
        min=0,
        max=6,
        help="Max mismatches per half-site in exhaustive genomic search (default: 2).",
    ),
    zfn_window_stride: int | None = typer.Option(
        None,
        "--zfn-window-stride",
        envvar="SIRNAFORGE_ZFN_WINDOW_STRIDE",
        min=1,
        max=50,
        hidden=True,
        help=("Internal tuning: sliding-window stride in bp for half-site scan (1 = fully exhaustive)."),
    ),
    zfn_top_n_sites: int | None = typer.Option(
        None,
        "--zfn-top-n-sites",
        envvar="SIRNAFORGE_ZFN_TOP_N_SITES",
        min=1,
        hidden=True,
        help="Internal tuning: maximum ranked off-target sites retained before candidate summarization.",
    ),
    zfn_report_n_sites: int | None = typer.Option(
        None,
        "--zfn-report-n-sites",
        envvar="SIRNAFORGE_ZFN_REPORT_N_SITES",
        min=1,
        hidden=True,
        help="Internal tuning: number of top ranked sites included in report outputs.",
    ),
    cores: int | None = typer.Option(
        None,
        "--cores",
        min=1,
        envvar="SIRNAFORGE_CORES",
        help=(
            "Total CPU core budget for workflow execution. ZFN sharding and workflow parallel stages derive from this."
        ),
    ),
    zfn_annotation: str | None = typer.Option(
        None,
        "--zfn-annotation",
        help="Optional GTF/GFF annotation file for ZFN off-target region classification.",
    ),
    top_n_candidates: int | None = typer.Option(
        None,
        "--top-n",
        "-n",
        min=1,
        help=(
            "Cap how many top-ranked candidates are reported (default: no cap, report all). "
            "Screening and enumeration always cover every candidate, so this only truncates the "
            "reported set -- leave it unset to keep the full design space."
        ),
    ),
    species: str = typer.Option(
        DEFAULT_SPECIES_ARGUMENT,
        "--species",
        help=(
            "Comma-separated canonical species identifiers. This single parameter drives "
            "all off-target analysis: miRNA database lookups (default: 7 species) and "
            "transcriptome fetching from Ensembl (default: 4 species). "
            "Override specific layers with --mirna-species or --transcriptome-fasta. "
            "Supported: human, mouse, macaque, rat, chicken, pig, rhesus"
        ),
    ),
    query_species: str | None = typer.Option(
        None,
        "--query-species",
        help=(
            "Organism the TARGET transcripts belong to, which decides which species' hits are "
            "on-target and whose alignment must succeed before candidates can be scored after "
            "screening. --species is an unordered set of references to screen AGAINST and never "
            "sets this. Defaults to the organism the gene-query database serves (human), which "
            "is also the species of the default transcriptome; set it when designing against an "
            "input FASTA from another organism."
        ),
    ),
    mirna_db: str = typer.Option(
        "mirgenedb",
        "--mirna-db",
        help="miRNA reference database to use for seed analysis",
    ),
    mirna_species: str | None = typer.Option(
        None,
        "--mirna-species",
        help=(
            "Override miRNA species identifiers (comma-separated). "
            "When omitted, automatically maps from --species. "
            "Use this for surgical control of miRNA database queries."
        ),
    ),
    transcriptome_fasta: str | None = typer.Option(
        None,
        "--transcriptome-fasta",
        help=(
            "Override or extend transcriptome references for off-target analysis. "
            "Accepts: local file, HTTP(S) URL, or pre-configured source (e.g., 'ensembl_human_cdna'). "
            "Prefix it with a species to state one -- 'mouse:/path/custom_cdna.fa' -- otherwise the "
            "species is read from the reference's own Ensembl cDNA headers. "
            "When omitted, automatically fetches Ensembl cDNA for species selected via --species. "
            "Custom FASTA files are cached and indexed automatically. "
            "Use this to add novel sequences (e.g., synthetic contigs) to the default set."
        ),
    ),
    transcriptome_filter: str | None = typer.Option(
        None,
        "--transcriptome-filter",
        help=(
            "Filter transcriptome to reduce size and memory requirements. "
            "Comma-separated filter names: 'protein_coding' (only protein-coding genes), "
            "'canonical_only' (only canonical isoforms). "
            "Example: --transcriptome-filter protein_coding,canonical_only. "
            "Filtered versions are cached separately with automatic indexing."
        ),
    ),
    transcriptome_indices: str | None = typer.Option(
        None,
        "--transcriptome-indices",
        "--offtarget-indices",
        help=(
            "Transcriptome BWA-MEM2 indices you have already built, as species:prefix entries: "
            "human:/abs/path/hs_cdna,mouse:/abs/path/mm_cdna. The species named here is that "
            "reference's species and is added to the screen; it does NOT suppress the references "
            "--species resolves, so pass --transcriptome-fasta (or --run-mode design_only) if you "
            "do not want the Ensembl defaults fetched as well. The cDNA FASTA the index was built "
            "from must be readable beside the prefix (the prefix itself, or <prefix>.fa) so hits "
            "can be resolved to genes."
        ),
    ),
    ortholog_mapping: Path | None = typer.Option(
        None,
        "--ortholog-mapping",
        help=(
            "JSON file of query gene -> species -> orthologue gene IDs, e.g. "
            '{"ENSG00000141510": {"mouse": ["ENSMUSG00000059552"]}}. '
            "Use it to classify cross-species hits offline instead of calling Ensembl Compara; "
            "a species the file omits falls back to the labelled gene-symbol heuristic."
        ),
    ),
    gc_min: float | None = typer.Option(
        None,
        "--gc-min",
        min=0.0,
        max=100.0,
        help=(
            f"Minimum GC content percentage (default: {default_for('gc_min')}; "
            f"--design-mode mirna defaults to {mirna_preset_default_for('gc_min')}). Any value in "
            "0-100 is supported."
        ),
    ),
    gc_max: float | None = typer.Option(
        None,
        "--gc-max",
        min=0.0,
        max=100.0,
        help=(
            f"Maximum GC content percentage (default: {default_for('gc_max')}; "
            f"--design-mode mirna defaults to {mirna_preset_default_for('gc_max')}). Any value in "
            "0-100 is supported, and stating it keeps it in every design mode."
        ),
    ),
    sirna_length: int | None = typer.Option(
        None,
        "--length",
        "-l",
        min=19,
        max=23,
        help=f"siRNA length in nucleotides (default: {default_for('sirna_length')})",
    ),
    modification_pattern: str | None = typer.Option(
        None,
        "--modifications",
        "-m",
        help=(
            "Chemical modification pattern: standard_2ome, minimal_terminal, maximal_stability, none "
            f"(default: {default_for('modification_pattern')})"
        ),
    ),
    overhang: str | None = typer.Option(
        None,
        "--overhang",
        help=(
            f"Overhang sequence, dTdT for DNA or UU for RNA (default: {default_for('default_overhang')}; "
            f"--design-mode mirna defaults to {mirna_preset_default_for('default_overhang')})"
        ),
    ),
    skip_off_targets: bool = typer.Option(
        False,
        "--skip-off-targets",
        help="Skip off-target analysis (faster). Maps to --run-mode design_only.",
    ),
    # Variant targeting parameters
    snp: list[str] = typer.Option(
        [],
        "--snp",
        help=(
            "Variant identifier(s) for SNP targeting/avoidance. "
            "Accepts rsID (rs12345), coordinate (chr17:7577121:G:A), or HGVS (NM_000546.6:c.215C>G). "
            "Can be specified multiple times. All variants must be on GRCh38 assembly."
        ),
    ),
    snp_file: Path | None = typer.Option(
        None,
        "--snp-file",
        help=(
            "VCF file containing variants for targeting/avoidance. "
            "Preferably bgzip-compressed with tabix index (.vcf.gz + .tbi) for performance. "
            "Variants are filtered by --min-af and --clinvar-filter-levels."
        ),
    ),
    variant_mode: VariantMode = typer.Option(
        VariantMode.AVOID,
        "--variant-mode",
        help=(
            "How to handle variants in siRNA design: "
            "'avoid' = exclude candidates overlapping variants (default), "
            "'target' = design siRNAs specifically targeting variant alleles, "
            "'both' = generate candidates for both reference and alternate alleles."
        ),
    ),
    min_af: float = typer.Option(
        0.01,
        "--min-af",
        min=0.0,
        max=1.0,
        help=(
            "Minimum allele frequency threshold for variant inclusion. "
            "Variants with AF below this value are excluded (default: 0.01 = 1%%)."
        ),
    ),
    clinvar_filter_levels: str = typer.Option(
        "Pathogenic,Likely pathogenic",
        "--clinvar-filter-levels",
        help=(
            "Comma-separated ClinVar clinical significance levels to include. "
            "Default: 'Pathogenic,Likely pathogenic'. "
            "Other options: 'Benign', 'Likely benign', 'Uncertain significance'."
        ),
    ),
    variant_assembly: str = typer.Option(
        "GRCh38",
        "--variant-assembly",
        help="Reference genome assembly for variants (only GRCh38 supported)",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    log_file: Path | None = typer.Option(
        None,
        "--log-file",
        help="Path to centralized log file (overrides SIRNAFORGE_LOG_FILE env)",
    ),
    nextflow_docker_image: str | None = typer.Option(
        None,
        "--nextflow-docker-image",
        envvar="SIRNAFORGE_NEXTFLOW_IMAGE",
        help=(f"Override the Docker image passed to Nextflow (default: {DEFAULT_SIRNAFORGE_DOCKER_IMAGE})"),
    ),
    max_hits: int | None = typer.Option(
        None,
        "--max-hits",
        min=1,
        help=(
            "Cap off-target hits retained per candidate per species (default: exhaustive, no cap). "
            "Set a lower value (e.g. 10000) to speed up analysis of large gene families at the cost of "
            "censoring per-species hit counts."
        ),
    ),
    max_off_targets: int | None = typer.Option(
        None,
        "--max-off-targets",
        min=0,
        help=(
            f"Reject a candidate above this many genuine off-target sites (default: "
            f"{default_for('max_off_target_count')}). Counts only hits left after on-target, "
            "ortholog and repeat classification. Unlike --max-hits this changes the "
            "PASS/EXCESS_OFF_TARGETS gate, not how many hits are recorded."
        ),
    ),
    min_asymmetry: float | None = typer.Option(
        None,
        "--min-asymmetry",
        min=0.3,
        max=1.0,
        help=(
            f"Thermodynamic asymmetry floor gating LOW_ASYMMETRY (default: "
            f"{default_for('min_asymmetry_score')}). The default has not been calibrated against "
            "measured potency; lower it to widen the candidate pool."
        ),
    ),
    max_paired_fraction: float | None = typer.Option(
        None,
        "--max-paired-fraction",
        min=0.0,
        max=1.0,
        help=(
            f"Guide self-structure ceiling gating EXCESS_PAIRING (default: "
            f"{default_for('max_paired_fraction')}). Quantised to 2k/L by the dot-bracket, so on a "
            "21-23mer only a few values are attainable and nearby thresholds behave identically."
        ),
    ),
    min_empirical: float | None = typer.Option(
        None,
        "--min-empirical",
        help=(
            "Empirical design-rule floor gating LOW_EMPIRICAL_SCORE. Applied to the empirical "
            "component score, whose attainable range is narrow."
        ),
    ),
    min_isoform_coverage: float | None = typer.Option(
        None,
        "--min-isoform-coverage",
        min=0.0,
        max=1.0,
        help=(
            "Opt into the protein-coding isoform coverage gate (LOW_ISOFORM_COVERAGE). Off by "
            "default: coverage is reported on every candidate either way, and no floor has been "
            "calibrated against truth data. Not a scoring term."
        ),
    ),
    plfold_window: int | None = typer.Option(
        None,
        "--plfold-window",
        min=20,
        max=1000,
        help=(
            f"RNAplfold averaging window W for target-site accessibility (default: "
            f"{default_for('plfold_window')}). Larger is mildly better on the benchmark; W=150-240 "
            "is the plateau. Changing it changes the scale of composite_score."
        ),
    ),
    plfold_max_bp_span: int | None = typer.Option(
        None,
        "--plfold-max-bp-span",
        min=10,
        max=1000,
        help=(
            f"RNAplfold maximum base-pair span L (default: {default_for('plfold_max_bp_span')}). "
            "Must not exceed --plfold-window."
        ),
    ),
    accessibility_log_floor: float | None = typer.Option(
        None,
        "--accessibility-log-floor",
        help=(
            f"log10 probability treated as zero target-site accessibility (default: "
            f"{default_for('accessibility_log_floor')}, ~99% of observed sites). Fixed, not "
            "per-transcript, so scores stay comparable between targets."
        ),
    ),
    json_summary: bool = typer.Option(
        True,
        "--json-summary/--no-json-summary",
        help="Write logs/workflow_summary.json (disable to skip JSON output)",
    ),
) -> None:
    """Run the end-to-end workflow: transcripts → siRNA design → off-target.

    This is the main orchestration command. It resolves transcriptome and miRNA
    reference policies, designs candidates, and then runs off-target analysis on
    the selected top candidates.
    """
    # Resolve and validate the biological policy BEFORE anything creates a directory or a log file:
    # an invalid bound or an incompatible mode must cost nothing.
    try:
        actions = _parse_filter_actions(filter_action)
    except RunPolicyError as exc:
        _fail_with_config_error(str(exc))
    stated = _stated_settings(
        ctx,
        {
            "gc_min": ("gc_min", gc_min),
            "gc_max": ("gc_max", gc_max),
            "sirna_length": ("sirna_length", sirna_length),
            "modification_pattern": ("modification_pattern", modification_pattern),
            "default_overhang": ("overhang", overhang),
            "top_n": ("top_n_candidates", top_n_candidates),
            "max_off_target_count": ("max_off_targets", max_off_targets),
            "min_asymmetry_score": ("min_asymmetry", min_asymmetry),
            "max_paired_fraction": ("max_paired_fraction", max_paired_fraction),
            "min_empirical_score": ("min_empirical", min_empirical),
            "min_isoform_coverage": ("min_isoform_coverage", min_isoform_coverage),
            "plfold_window": ("plfold_window", plfold_window),
            "plfold_max_bp_span": ("plfold_max_bp_span", plfold_max_bp_span),
            "accessibility_log_floor": ("accessibility_log_floor", accessibility_log_floor),
        },
    )
    if "modification_pattern" in stated:
        stated["apply_modifications"] = str(stated["modification_pattern"]).lower() != "none"
    policy = _resolve_policy_or_exit(
        entry_point=EntryPoint.SCREENING_WORKFLOW,
        # Passed only when the user typed it, so the manifest does not report a mode nobody chose.
        design_mode=design_mode if _option_was_stated(ctx, "design_mode") else None,
        run_mode=run_mode,
        config_file=policy_config,
        stated=stated,
        filter_actions=actions,
        legacy_skip_screening=skip_off_targets or None,
        query_species=query_species,
        screen_species=[value.strip() for value in species.split(",") if value.strip()],
    )
    mode_enum = policy.design_mode
    resolved_filters = policy.design_parameters.filters

    log_destination = Path(log_file) if log_file else output_dir / "logs" / "sirnaforge.log"
    log_destination.parent.mkdir(parents=True, exist_ok=True)
    configure_logging(level=os.getenv("SIRNAFORGE_LOG_LEVEL"), log_file=str(log_destination))
    effective_log = str(log_destination)
    logger = logging.getLogger(__name__)

    if mode_enum == DesignMode.ZFN:
        emit_zfn_experimental_warning(console)

    merged_zfn_constraints = list(zfn_subfinger_mutation)
    if zfn_max_mismatches_per_subfinger is not None:
        merged_zfn_constraints.append(f"*:{zfn_max_mismatches_per_subfinger}:mismatch")
    if zfn_max_substitutions_overall is not None:
        merged_zfn_constraints.append(f"overall:{zfn_max_substitutions_overall}:substitution")

    try:
        zfn_constraints, zfn_default_constraint, zfn_overall_constraints = _parse_zfn_mutation_constraints(
            merged_zfn_constraints
        )
    except ValueError as exc:
        logger.error("Invalid ZFN sub-finger mutation configuration: %s", exc)
        console.print(f"❌ Error: {exc}", style="red")
        raise typer.Exit(1)

    if mode_enum != DesignMode.ZFN and (zfn_constraints or zfn_default_constraint or zfn_overall_constraints):
        logger.error("ZFN constraints provided while design mode is %s", mode_enum.value)
        console.print("❌ Error: --zfn-subfinger-mutation requires --design-mode zfn", style="red")
        raise typer.Exit(1)

    # ── Assemble ZFNDesignParameters when mode is ZFN ──
    zfn_design_params: ZFNDesignParameters | None = None
    annotation: GenomicAnnotationConfig | None = None
    if mode_enum == DesignMode.ZFN:
        if not zfn_left_half_site or not zfn_right_half_site:
            console.print(
                "❌ Error: --zfn-left-half-site and --zfn-right-half-site are required for --design-mode zfn",
                style="red",
            )
            raise typer.Exit(1)

        try:
            zfn_design_params, annotation, zfn_constraints, zfn_default_constraint, zfn_overall_constraints = (
                _build_zfn_design_configuration(
                    zfn_subfinger_mutation=zfn_subfinger_mutation,
                    zfn_max_mismatches_per_subfinger=zfn_max_mismatches_per_subfinger,
                    zfn_max_substitutions_overall=zfn_max_substitutions_overall,
                    zfn_left_half_site=zfn_left_half_site,
                    zfn_right_half_site=zfn_right_half_site,
                    zfn_search_space=zfn_search_space,
                    zfn_search_space_index=zfn_search_space_index,
                    zfn_search_backend=zfn_search_backend,
                    zfn_algorithm=zfn_algorithm,
                    zfn_dimer_mode=zfn_dimer_mode,
                    zfn_spacer_lengths=zfn_spacer_lengths,
                    zfn_max_mismatches=zfn_max_mismatches,
                    zfn_window_stride=zfn_window_stride,
                    zfn_top_n_sites=zfn_top_n_sites,
                    zfn_report_n_sites=zfn_report_n_sites,
                    workflow_cores=cores,
                    zfn_annotation=zfn_annotation,
                )
            )
        except ValueError as exc:
            logger.error("Invalid ZFN configuration: %s", exc)
            console.print(f"❌ Error: {exc}", style="red")
            raise typer.Exit(1)
        except Exception as exc:
            logger.error("ZFN parameter validation failed: %s", exc)
            console.print(
                f"❌ Error: ZFN parameter validation failed: {exc}",
                style="red",
            )
            raise typer.Exit(1)

    try:
        resolved_species = resolve_species_inputs(species=species, mirna_db=mirna_db, mirna_species=mirna_species)
        override_species = extract_declared_species_from_indices(transcriptome_indices)
    except ValueError as exc:
        logger.error("Species resolution failed: %s", exc)
        console.print(f"❌ Error: {exc}", style="red")
        raise typer.Exit(1)

    source_normalized = resolved_species.source_normalized
    canonical_species = resolved_species.canonical_species
    species_list = resolved_species.screen_species
    mirna_species_list = resolved_species.mirna_species

    if not mirna_species_list:
        logger.error(
            "Failed to resolve miRNA species for species=%s mirna_db=%s mirna_overrides=%s",
            species,
            mirna_db,
            mirna_species,
        )
        console.print("❌ Error: failed to resolve miRNA species for selected inputs", style="red")
        raise typer.Exit(1)

    input_descriptor = gene_query
    if input_fasta:
        input_descriptor = input_fasta if "://" in input_fasta else Path(input_fasta).name

    # Resolve transcriptome policy once so downstream layers receive metadata.
    # allow_transcriptome_for_input_fasta stays False: --input-fasta means "design against MY
    # sequences", and auto-resolving DEFAULT_TRANSCRIPTOME_SOURCES there downloads and indexes
    # four multi-gigabyte cDNA references nobody asked for (it also contradicted the documented
    # design-only behaviour and timed out the toy container workflow). --transcriptome-fasta is
    # the explicit opt-in, and it accepts a bundled source name such as ensembl_human_cdna.
    # Key reference selection on the *resolved* run mode, not on the legacy flag: --run-mode
    # design_only and --skip-off-targets are the same run, and reading the raw flag here made them
    # publish two different reference_summary records for it.
    design_only_run = policy.run_mode is RunMode.DESIGN_ONLY
    transcriptome_spec = WorkflowInputSpec(
        input_fasta=input_fasta,
        transcriptome_argument=transcriptome_fasta,
        default_transcriptomes=DEFAULT_TRANSCRIPTOME_SOURCES,
        design_only=design_only_run,
        allow_transcriptome_for_input_fasta=False,
    )
    transcriptome_selection = ReferencePolicyResolver(transcriptome_spec).resolve_transcriptomes()
    transcriptome_label = render_reference_selection_label(transcriptome_selection)
    if input_fasta and not transcriptome_fasta and not design_only_run:
        console.print(
            "ℹ️  --input-fasta without --transcriptome-fasta: transcriptome off-target screening and "
            "repeat detection are disabled (design-only). Pass --transcriptome-fasta "
            "ensembl_human_cdna (or a path/URL) to screen against a reference.",
            style="yellow",
        )
    screen_species_for_workflow = override_species or species_list
    index_override_label = transcriptome_indices or "resolved from --species"
    nextflow_image_label = nextflow_docker_image or DEFAULT_SIRNAFORGE_DOCKER_IMAGE

    console.print(
        Panel.fit(
            f"🧬 [bold blue]Complete siRNA Workflow[/bold blue]\n"
            f"Design Mode: [cyan]{mode_enum.value}[/cyan]\n"
            f"Run Mode: [cyan]{policy.run_mode.value}[/cyan] "
            f"([magenta]{policy.profile.name} {policy.profile.version}, experimental[/magenta])\n"
            f"Gene Query: [cyan]{input_descriptor}[/cyan]\n"
            f"Database: [yellow]{database}[/yellow]\n"
            f"Output Directory: [cyan]{output_dir}[/cyan]\n"
            f"siRNA Length: [yellow]{policy.design_parameters.sirna_length}[/yellow] nt\n"
            f"GC Range: [yellow]{resolved_filters.gc_min:.1f}%-{resolved_filters.gc_max:.1f}%[/yellow]\n"
            f"Reported Candidates: [yellow]{policy.design_parameters.top_n if policy.design_parameters.top_n is not None else 'all (uncapped)'}[/yellow]\n"
            f"Species (canonical): [green]{', '.join(canonical_species)}[/green]\n"
            f"  ↳ miRNA Database ({source_normalized}): [green]{', '.join(mirna_species_list)}[/green]\n"
            f"  ↳ Transcriptome Reference: [green]{transcriptome_label}[/green]\n"
            f"  ↳ Prebuilt Index Override: [green]{index_override_label}[/green]\n"
            f"  ↳ Nextflow Docker Image: [green]{nextflow_image_label}[/green]\n"
            f"Modifications: [magenta]{policy.design_parameters.modification_pattern}[/magenta]\n"
            f"Overhang: [magenta]{policy.design_parameters.default_overhang}[/magenta]\n"
            f"ZFN Constraints: [magenta]{len(zfn_constraints)} explicit, "
            f"{1 if zfn_default_constraint else 0} default, {len(zfn_overall_constraints)} overall[/magenta]",
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

            results = asyncio.run(
                run_sirna_workflow(
                    gene_query=gene_query,
                    input_fasta=input_fasta,
                    output_dir=str(output_dir),
                    database=database,
                    # One resolved policy, so the CLI and a direct API call cannot diverge. Every
                    # threshold, the design mode and the run mode travel inside it.
                    resolved_policy=policy,
                    screen_species=screen_species_for_workflow,
                    query_species=query_species,
                    transcriptome_indices=transcriptome_indices,
                    mirna_database=source_normalized,
                    mirna_species=mirna_species_list,
                    transcriptome_fasta=transcriptome_fasta,
                    transcriptome_filter=transcriptome_filter,
                    transcriptome_selection=transcriptome_selection,
                    ortholog_mapping_file=ortholog_mapping,
                    zfn_design_params=zfn_design_params,
                    zfn_annotation=annotation if mode_enum == DesignMode.ZFN else None,
                    # Variant parameters
                    variant_ids=list(snp) if snp else None,
                    variant_vcf_file=snp_file,
                    variant_mode=variant_mode.value,
                    variant_min_af=min_af,
                    variant_clinvar_filters=clinvar_filter_levels,
                    variant_assembly=variant_assembly,
                    log_file=effective_log,
                    write_json_summary=json_summary,
                    num_threads=cores,
                    nextflow_docker_image=nextflow_docker_image,
                    max_hits=max_hits,
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

        if mode_enum == DesignMode.ZFN:
            summary_table.add_row(
                "ZFN Pair Search",
                "Complete",
                f"{results.get('off_target_sites', 0)} off-target sites",
            )
            summary_table.add_row(
                "ZFN Candidate Scoring",
                "Complete",
                f"{results.get('candidates', 0)} candidates",
            )
            summary_table.add_row(
                "Annotation",
                "Complete" if results.get("annotation_source") else "⚠️  Skipped",
                str(results.get("annotation_source") or "none"),
            )
        else:
            transcript_summary = results.get("transcript_summary", {})
            design_summary = results.get("design_summary", {})
            offtarget_summary = results.get("offtarget_summary", {})

            summary_table.add_row(
                "Transcript Retrieval",
                "Complete",
                f"{transcript_summary.get('total_transcripts', 0)} transcripts from {database}",
            )

            summary_table.add_row(
                "🧬 siRNAforge", "✅ Complete", f"{design_summary.get('total_candidates', 0)} candidates generated"
            )

            summary_table.add_row(
                "Off-target Analysis",
                "Complete" if offtarget_summary.get("status") == "completed" else "⚠️  Partial",
                f"Method: {offtarget_summary.get('method', 'basic')}",
            )

        console.print(summary_table)

        # Output locations
        console.print(f"\n📁 [bold]Results saved to:[/bold] [cyan]{output_dir}[/cyan]")
        console.print("📂 Key files:")
        if mode_enum == DesignMode.ZFN:
            console.print("   • Off-target sites: [blue]sirnaforge/offtarget_sites.csv[/blue]")
            console.print("   • Candidate summary: [blue]sirnaforge/candidate_summary.json[/blue]")
            console.print("   • Console stream log: [blue]logs/workflow_stream.log[/blue]")
            if json_summary:
                console.print("   • Workflow summary: [blue]logs/workflow_summary.json[/blue]")
        else:
            offtarget_summary = results.get("offtarget_summary", {})
            console.print(f"   • Transcripts: [blue]transcripts/{gene_query}_transcripts.fasta[/blue]")
            console.print("   • siRNA candidates (ALL): [blue]sirnaforge/candidates_all.csv[/blue]")
            console.print("   • siRNA candidates (PASS): [blue]sirnaforge/candidates_pass.csv[/blue]")
            if offtarget_summary:
                console.print(_offtarget_results_line(offtarget_summary, "off_target/results/"))
            console.print("   • Console stream log: [blue]logs/workflow_stream.log[/blue]")
            if json_summary:
                console.print("   • Workflow summary: [blue]logs/workflow_summary.json[/blue]")

            if offtarget_summary.get("method") == "nextflow":
                console.print("   • Full off-target report: [blue]off_target/results/offtarget_report.html[/blue]")

    except Exception as e:
        logger.exception("Workflow execution failed")
        console.print(f"❌ [red]Workflow error:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


@app_command()
def offtarget(  # noqa: PLR0912
    run_mode: str | None = typer.Option(
        None,
        "--run-mode",
        help=(
            "How much screening evidence this run claims: design_only, exploratory or qualified "
            "(default: qualified). design_only would screen nothing, so it is rejected here."
        ),
    ),
    policy_config: Path | None = typer.Option(
        None,
        "--policy-config",
        help="JSON or TOML file of policy settings; loses to options given on the command line.",
    ),
    filter_action: list[str] = typer.Option(
        [],
        "--filter-action",
        help=(
            "Set one filter's action: filter_id=off|warn|fail (repeatable). 'off' clears the gate's "
            f"threshold, so it is not evaluated, which is not the same as passing: {', '.join(switchable_filter_ids())}. "
            "The seven design-stage gates have no absent threshold and are refused rather than faked. "
            "'warn' keeps the gate measuring and stops it rejecting: the outcome is published in "
            "<filter_id>_verdict and <filter_id>_observed and passes_filters is left alone. Every gate "
            "this command applies -- the off-target gates and min_isoform_coverage -- honours it, "
            "because each one records through the shared verdict recorder. The design-stage gates do "
            "not run on pre-designed guides at all, so an action set on one is accepted and applies "
            "to nothing. min_isoform_coverage, max_transcriptome_seed_perfect and "
            "max_total_offtarget_hits ship with no threshold, so set one before asking for 'warn' or "
            "'fail'; max_mirna_1mm_seed is read by no 0.7.1 gate, so only 'off' is accepted for it."
        ),
    ),
    input_candidates_fasta: Path = typer.Option(
        ...,
        "--input-candidates-fasta",
        "-i",
        help="FASTA file containing pre-designed siRNA guide sequences (any length)",
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    output_dir: Path = typer.Option(
        Path("offtarget_output"),
        "--output-dir",
        "-o",
        help="Output directory for off-target analysis results",
    ),
    species: str = typer.Option(
        DEFAULT_SPECIES_ARGUMENT,
        "--species",
        help=(
            "Comma-separated canonical species identifiers for off-target analysis. "
            "Drives transcriptome fetching from Ensembl and miRNA database lookups. "
            "Supported: human, mouse, macaque, rat, chicken, pig, rhesus"
        ),
    ),
    query_species: str | None = typer.Option(
        None,
        "--query-species",
        help=(
            "Organism the supplied guides were designed against, which decides whose alignment "
            "must succeed before candidates can be scored after screening. --species is the set "
            "of references to screen AGAINST and never sets this. Defaults to human."
        ),
    ),
    mirna_db: str = typer.Option(
        "mirgenedb",
        "--mirna-db",
        help="miRNA reference database to use for seed analysis",
    ),
    mirna_species: str | None = typer.Option(
        None,
        "--mirna-species",
        help=("Override miRNA species identifiers (comma-separated). When omitted, automatically maps from --species."),
    ),
    transcriptome_fasta: str | None = typer.Option(
        None,
        "--transcriptome-fasta",
        help=(
            "Override or extend transcriptome references for off-target analysis. "
            "Accepts: local file, HTTP(S) URL, or pre-configured source (e.g., 'ensembl_human_cdna'). "
            "Prefix it with a species to state one -- 'mouse:/path/custom_cdna.fa' -- otherwise the "
            "species is read from the reference's own Ensembl cDNA headers."
        ),
    ),
    transcriptome_filter: str | None = typer.Option(
        None,
        "--transcriptome-filter",
        help=(
            "Filter transcriptome to reduce size and memory requirements. "
            "Comma-separated filter names: 'protein_coding', 'canonical_only'. "
            "Example: --transcriptome-filter protein_coding,canonical_only."
        ),
    ),
    transcriptome_indices: str | None = typer.Option(
        None,
        "--transcriptome-indices",
        "--offtarget-indices",
        help=(
            "Transcriptome BWA-MEM2 indices you have already built, as species:prefix entries: "
            "human:/abs/path/hs_cdna,mouse:/abs/path/mm_cdna. The species named here is added to "
            "the screen and does not suppress the references --species resolves. The cDNA FASTA the "
            "index was built from must be readable beside the prefix so hits can be resolved to genes."
        ),
    ),
    ortholog_mapping: Path | None = typer.Option(
        None,
        "--ortholog-mapping",
        help=(
            "JSON file of query gene -> species -> orthologue gene IDs, used to classify "
            "cross-species hits offline instead of calling Ensembl Compara."
        ),
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    log_file: Path | None = typer.Option(
        None,
        "--log-file",
        help="Path to centralized log file (overrides SIRNAFORGE_LOG_FILE env)",
    ),
    nextflow_docker_image: str | None = typer.Option(
        None,
        "--nextflow-docker-image",
        envvar="SIRNAFORGE_NEXTFLOW_IMAGE",
        help=(f"Override the Docker image used by Nextflow (default: {DEFAULT_SIRNAFORGE_DOCKER_IMAGE})"),
    ),
) -> None:
    """Run off-target analysis on pre-designed siRNA candidates.

    This command accepts a FASTA file containing pre-designed siRNA guide sequences
    of any length and runs comprehensive off-target analysis including:
    - Transcriptome alignment (BWA-MEM2)
    - miRNA seed match analysis
    - Off-target hit classification and scoring

    The embedded Nextflow pipeline is used for parallel processing across species.

    Notes:
        - ``--species`` drives transcriptome fetching and miRNA lookup.
        - ``--offtarget-indices`` can override the indices used for alignment
          using ``species:/abs/path/index_prefix`` entries.
    """
    # Policy first: an invalid threshold must not cost a FASTA read, a download or a directory.
    try:
        actions = _parse_filter_actions(filter_action)
    except RunPolicyError as exc:
        _fail_with_config_error(str(exc))
    policy = _resolve_policy_or_exit(
        entry_point=EntryPoint.OFFTARGET_ONLY,
        run_mode=run_mode,
        config_file=policy_config,
        filter_actions=actions,
        query_species=query_species,
        screen_species=[value.strip() for value in species.split(",") if value.strip()],
    )
    if policy.run_mode is RunMode.DESIGN_ONLY:
        _fail_with_config_error(
            "--run-mode design_only would screen nothing; this command exists to screen pre-designed guides"
        )

    # Validate input FASTA contains sequences (any length accepted)
    try:
        sequences = FastaUtils.read_fasta(input_candidates_fasta)

        if not sequences:
            console.print("❌ [red]Error:[/red] Input FASTA file is empty", style="red")
            raise typer.Exit(1)

        # Report sequence statistics without enforcing length constraints
        seq_lengths = [len(seq) for _, seq in sequences]
        min_len = min(seq_lengths)
        max_len = max(seq_lengths)

        if min_len == max_len:
            console.print(f"✅ Validated {len(sequences)} siRNA candidates (all {min_len} nt)")
        else:
            console.print(f"✅ Validated {len(sequences)} siRNA candidates ({min_len}-{max_len} nt)")

    except Exception as e:
        if isinstance(e, typer.Exit):
            raise
        console.print(f"❌ [red]Error validating input FASTA:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)

    try:
        resolved_species = resolve_species_inputs(species=species, mirna_db=mirna_db, mirna_species=mirna_species)
        override_species = extract_declared_species_from_indices(transcriptome_indices)
    except ValueError as exc:
        console.print(f"❌ Error: {exc}", style="red")
        raise typer.Exit(1)

    source_normalized = resolved_species.source_normalized
    canonical_species = resolved_species.canonical_species
    species_list = resolved_species.screen_species
    mirna_species_list = resolved_species.mirna_species

    if not mirna_species_list:
        console.print("❌ Error: failed to resolve miRNA species for selected inputs", style="red")
        raise typer.Exit(1)

    # Resolve transcriptome policy
    transcriptome_spec = WorkflowInputSpec(
        input_fasta=None,  # Not using input transcripts for off-target-only
        transcriptome_argument=transcriptome_fasta,
        default_transcriptomes=DEFAULT_TRANSCRIPTOME_SOURCES,
        design_only=False,
    )
    transcriptome_selection = ReferencePolicyResolver(transcriptome_spec).resolve_transcriptomes()
    transcriptome_label = render_reference_selection_label(transcriptome_selection)

    screen_species_for_workflow = override_species or species_list
    index_override_label = transcriptome_indices or "resolved from --species"
    nextflow_image_label = nextflow_docker_image or DEFAULT_SIRNAFORGE_DOCKER_IMAGE

    console.print(
        Panel.fit(
            f"🎯 [bold blue]Off-Target Analysis (Pre-Designed siRNAs)[/bold blue]\n"
            f"Input Candidates: [cyan]{input_candidates_fasta.name}[/cyan]\n"
            f"Candidate Count: [yellow]{len(sequences)}[/yellow]\n"
            f"Output Directory: [cyan]{output_dir}[/cyan]\n"
            f"Run Mode: [cyan]{policy.run_mode.value}[/cyan] "
            f"([magenta]{policy.profile.name} {policy.profile.version}, experimental[/magenta])\n"
            f"Species (canonical): [green]{', '.join(canonical_species)}[/green]\n"
            f"  ↳ miRNA Database ({source_normalized}): [green]{', '.join(mirna_species_list)}[/green]\n"
            f"  ↳ Transcriptome Reference: [green]{transcriptome_label}[/green]\n"
            f"  ↳ Prebuilt Index Override: [green]{index_override_label}[/green]\n"
            f"  ↳ Nextflow Docker Image: [green]{nextflow_image_label}[/green]",
            title="Off-Target Configuration",
        )
    )

    try:
        # Run off-target-only workflow
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task("Running off-target analysis...", total=None)

            # Configure logging
            effective_log = str(log_file) if log_file else str(Path(output_dir) / "logs" / "sirnaforge.log")
            configure_logging(log_file=effective_log, level=os.getenv("SIRNAFORGE_LOG_LEVEL"))

            # Run workflow
            results = asyncio.run(
                run_offtarget_only_workflow(
                    input_candidates_fasta=str(input_candidates_fasta),
                    output_dir=str(output_dir),
                    screen_species=screen_species_for_workflow,
                    query_species=query_species,
                    transcriptome_indices=transcriptome_indices,
                    mirna_database=source_normalized,
                    mirna_species=mirna_species_list,
                    transcriptome_fasta=transcriptome_fasta,
                    transcriptome_filter=transcriptome_filter,
                    transcriptome_selection=transcriptome_selection,
                    ortholog_mapping_file=ortholog_mapping,
                    log_file=effective_log,
                    nextflow_docker_image=nextflow_docker_image,
                    resolved_policy=policy,
                )
            )

            progress.remove_task(task)

        # Display results summary
        console.print("\n✅ [bold green]Off-target analysis completed successfully![/bold green]")

        offtarget_summary = results.get("offtarget_summary", {})

        summary_table = Table(title="📊 Off-Target Results Summary")
        summary_table.add_column("Metric", style="cyan")
        summary_table.add_column("Value", style="white")

        summary_table.add_row(
            "Status", "✅ Complete" if offtarget_summary.get("status") == "completed" else "⚠️ Partial"
        )
        summary_table.add_row("Method", offtarget_summary.get("method", "N/A"))
        summary_table.add_row("Candidates Analyzed", str(len(sequences)))

        console.print(summary_table)

        # Output locations
        console.print(f"\n📁 [bold]Results saved to:[/bold] [cyan]{output_dir}[/cyan]")
        console.print("📂 Key files:")
        console.print("   • Input candidates: [blue]input_candidates.fasta[/blue]")
        console.print(_offtarget_results_line(offtarget_summary, "results/"))
        console.print("   • Console log: [blue]logs/sirnaforge.log[/blue]")

        if offtarget_summary.get("method") == "embedded_nextflow":
            console.print("   • Full off-target report: [blue]results/offtarget_report.html[/blue]")

    except Exception as e:
        console.print(f"❌ [red]Off-target analysis error:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


@app_command()
def zfn(
    output_dir: Path = typer.Option(
        Path("sirna_zfn_output"),
        "--output-dir",
        "-o",
        help="Output directory for ZFN activity evaluation results",
    ),
    zfn_subfinger_mutation: list[str] = typer.Option(
        [],
        "--zfn-subfinger-mutation",
        help=(
            "ZFN sub-finger mutation allowance. "
            "Repeatable format: scope:max_mutations:type1,type2. "
            "scope can be subfinger index (e.g. 2), '*' for default per-subfinger, "
            "or 'overall' for global budgets. "
            "Use 'mismatch' as a shorthand alias for 'substitution'."
        ),
    ),
    zfn_max_mismatches_per_subfinger: int | None = typer.Option(
        None,
        "--zfn-max-mismatches-per-subfinger",
        min=0,
        help="Convenience option equivalent to --zfn-subfinger-mutation '*:<N>:mismatch'.",
    ),
    zfn_max_substitutions_overall: int | None = typer.Option(
        None,
        "--zfn-max-substitutions-overall",
        min=0,
        help="Convenience option equivalent to --zfn-subfinger-mutation 'overall:<N>:substitution'.",
    ),
    zfn_left_half_site: str = typer.Option(
        ...,
        "--zfn-left-half-site",
        help="Left ZFN half-site sequence (9-18 bp, IUPAC allowed).",
    ),
    zfn_right_half_site: str = typer.Option(
        ...,
        "--zfn-right-half-site",
        help="Right ZFN half-site sequence (9-18 bp, IUPAC allowed).",
    ),
    zfn_search_space: str | None = typer.Option(
        None,
        "--zfn-search-space",
        help=(
            "Genome reference key or local FASTA path for ZFN off-target search space. "
            "Built-in keys: ensembl_human_hg38_primary, ensembl_mouse_grcm39_primary, "
            "ensembl_rat_grcr8_toplevel, ensembl_macaque_mmul10_toplevel. "
            "Default: ensembl_human_hg38_primary."
        ),
    ),
    zfn_search_space_index: str | None = typer.Option(
        None,
        "--zfn-search-space-index",
        help=(
            "Optional persisted search-space index bundle path for indexed ZFN backends "
            "(currently fm_index; fm_index is experimental on large references)."
        ),
    ),
    zfn_search_backend: ZFNSearchBackend = typer.Option(
        ZFNSearchBackend.PYAHOCORASICK,
        "--zfn-search-backend",
        help=(
            "Half-site search backend: pyahocorasick (default), "
            "exhaustive_python (baseline), or fm_index (experimental)."
        ),
    ),
    zfn_algorithm: ZFNAlgorithm = typer.Option(
        ZFNAlgorithm.ZFN_V2,
        "--zfn-algorithm",
        help="ZFN off-target scoring algorithm: homology, conserved_g, or zfn_v2 (default).",
    ),
    zfn_dimer_mode: DimerMode = typer.Option(
        DimerMode.HETERODIMER_ONLY,
        "--zfn-dimer-mode",
        help="Dimer mode: heterodimer_only (default) or include_homodimers.",
    ),
    zfn_spacer_lengths: str = typer.Option(
        "5,6,7",
        "--zfn-spacer-lengths",
        help="Comma-separated allowed spacer lengths between half-sites (default: 5,6,7).",
    ),
    zfn_max_mismatches: int = typer.Option(
        2,
        "--zfn-max-mismatches",
        min=0,
        max=6,
        help="Max mismatches per half-site in exhaustive genomic search (default: 2).",
    ),
    zfn_window_stride: int | None = typer.Option(
        None,
        "--zfn-window-stride",
        envvar="SIRNAFORGE_ZFN_WINDOW_STRIDE",
        min=1,
        max=50,
        hidden=True,
        help=("Internal tuning: sliding-window stride in bp for half-site scan (1 = fully exhaustive)."),
    ),
    zfn_top_n_sites: int | None = typer.Option(
        None,
        "--zfn-top-n-sites",
        envvar="SIRNAFORGE_ZFN_TOP_N_SITES",
        min=1,
        hidden=True,
        help="Internal tuning: maximum ranked off-target sites retained before candidate summarization.",
    ),
    zfn_report_n_sites: int | None = typer.Option(
        None,
        "--zfn-report-n-sites",
        envvar="SIRNAFORGE_ZFN_REPORT_N_SITES",
        min=1,
        hidden=True,
        help="Internal tuning: number of top ranked sites included in report outputs.",
    ),
    cores: int | None = typer.Option(
        None,
        "--cores",
        min=1,
        envvar="SIRNAFORGE_CORES",
        help=(
            "Total CPU core budget for workflow execution. ZFN sharding and workflow parallel stages derive from this."
        ),
    ),
    zfn_annotation: str | None = typer.Option(
        None,
        "--zfn-annotation",
        help="Optional GTF/GFF annotation file for ZFN off-target region classification.",
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
    log_file: Path | None = typer.Option(
        None,
        "--log-file",
        help="Path to centralized log file (overrides SIRNAFORGE_LOG_FILE env)",
    ),
    nextflow_docker_image: str | None = typer.Option(
        None,
        "--nextflow-docker-image",
        envvar="SIRNAFORGE_NEXTFLOW_IMAGE",
        help=(f"Override the Docker image passed to Nextflow (default: {DEFAULT_SIRNAFORGE_DOCKER_IMAGE})"),
    ),
    json_summary: bool = typer.Option(
        True,
        "--json-summary/--no-json-summary",
        help="Write logs/workflow_summary.json (disable to skip JSON output)",
    ),
) -> None:
    """Evaluate a ZFN pair and run exhaustive genome-wide off-target search (EXPERIMENTAL)."""
    log_destination = Path(log_file) if log_file else output_dir / "logs" / "sirnaforge.log"
    log_destination.parent.mkdir(parents=True, exist_ok=True)
    configure_logging(level=os.getenv("SIRNAFORGE_LOG_LEVEL"), log_file=str(log_destination))

    emit_zfn_experimental_warning(console)

    try:
        zfn_design_params, annotation, zfn_constraints, zfn_default_constraint, zfn_overall_constraints = (
            _build_zfn_design_configuration(
                zfn_subfinger_mutation=zfn_subfinger_mutation,
                zfn_max_mismatches_per_subfinger=zfn_max_mismatches_per_subfinger,
                zfn_max_substitutions_overall=zfn_max_substitutions_overall,
                zfn_left_half_site=zfn_left_half_site,
                zfn_right_half_site=zfn_right_half_site,
                zfn_search_space=zfn_search_space,
                zfn_search_space_index=zfn_search_space_index,
                zfn_search_backend=zfn_search_backend,
                zfn_algorithm=zfn_algorithm,
                zfn_dimer_mode=zfn_dimer_mode,
                zfn_spacer_lengths=zfn_spacer_lengths,
                zfn_max_mismatches=zfn_max_mismatches,
                zfn_window_stride=zfn_window_stride,
                zfn_top_n_sites=zfn_top_n_sites,
                zfn_report_n_sites=zfn_report_n_sites,
                workflow_cores=cores,
                zfn_annotation=zfn_annotation,
            )
        )
    except ValueError as exc:
        console.print(f"❌ Error: {exc}", style="red")
        raise typer.Exit(1)
    except Exception as exc:
        console.print(f"❌ Error: ZFN parameter validation failed: {exc}", style="red")
        raise typer.Exit(1)

    console.print(
        Panel.fit(
            f"🧬 [bold blue]ZFN Activity Evaluation[/bold blue]\n"
            f"Left half-site: [yellow]{zfn_design_params.left_half_site}[/yellow]\n"
            f"Right half-site: [yellow]{zfn_design_params.right_half_site}[/yellow]\n"
            f"Algorithm: [cyan]{zfn_design_params.algorithm.value}[/cyan]\n"
            f"Search backend: [cyan]{zfn_design_params.search_backend.value}[/cyan]\n"
            f"Dimer mode: [cyan]{zfn_design_params.dimer_mode.value}[/cyan]\n"
            f"Spacer lengths: [cyan]{zfn_design_params.spacer_constraints.allowed_spacer_lengths}[/cyan]\n"
            f"Internal tuning: [cyan]stride={zfn_design_params.half_site_constraints.window_stride}, "
            f"top_n={zfn_design_params.top_n_sites}, report_n={zfn_design_params.report_n_sites}[/cyan]\n"
            f"ZFN Constraints: [magenta]{len(zfn_constraints)} explicit, "
            f"{1 if zfn_default_constraint else 0} default, {len(zfn_overall_constraints)} overall[/magenta]\n"
            f"Output Directory: [cyan]{output_dir}[/cyan]",
            title="ZFN Configuration",
        )
    )

    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task = progress.add_task("Running ZFN activity workflow...", total=None)
            asyncio.run(
                run_sirna_workflow(
                    gene_query="zfn",
                    output_dir=str(output_dir),
                    database="ensembl",
                    design_mode=DesignMode.ZFN.value,
                    zfn_design_params=zfn_design_params,
                    zfn_annotation=annotation,
                    log_file=str(log_destination),
                    num_threads=cores,
                    write_json_summary=json_summary,
                    nextflow_docker_image=nextflow_docker_image,
                )
            )
            progress.remove_task(task)

        console.print("\n✅ [bold green]ZFN workflow completed successfully![/bold green]")
        console.print(f"📁 [bold]Results saved to:[/bold] [cyan]{output_dir}[/cyan]")
        console.print("📂 Key files:")
        console.print("   • Off-target sites: [blue]sirnaforge/offtarget_sites.csv[/blue]")
        console.print("   • Candidate summary: [blue]sirnaforge/candidate_summary.json[/blue]")
        if json_summary:
            console.print("   • Workflow summary: [blue]logs/workflow_summary.json[/blue]")
    except Exception as e:
        console.print(f"❌ [red]ZFN workflow error:[/red] {str(e)}")
        if verbose:
            console.print_exception()
        raise typer.Exit(1)


@app_command()
def design(  # noqa: PLR0912
    ctx: typer.Context,
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
    design_mode: str = typer.Option(
        "sirna",
        "--design-mode",
        help="Design mode: sirna (default) or mirna (miRNA-biogenesis-aware). For ZFN use 'sirnaforge zfn'.",
    ),
    policy_config: Path | None = typer.Option(
        None,
        "--policy-config",
        help=(
            "JSON or TOML file of policy settings. Beats the built-in profile and loses to options "
            "given on the command line. Unknown setting names are rejected, not ignored."
        ),
    ),
    filter_action: list[str] = typer.Option(
        [],
        "--filter-action",
        help=(
            "Set one filter's action: filter_id=off|warn|fail (repeatable). 'off' clears the gate's "
            "threshold, so it is not evaluated, which is not the same as passing; it is refused for "
            "the seven design-stage gates, whose thresholds are plain floats with no absent state "
            "(widen the threshold instead), and it is a no-op on a post-screen gate, which this "
            "design_only command never evaluates anyway. 'warn' means record the verdict in "
            "<filter_id>_verdict and leave passes_filters alone. KNOWN GAP: this command builds the "
            "designer without handing it the resolved actions, so 'warn' here is recorded in the "
            "manifest and is inert -- every design-stage gate applies its declared default instead. "
            "Use 'sirnaforge workflow' if you need warn honoured, and note that even there "
            "gc_content_min, gc_content_max, max_poly_runs and max_repeat_transcript_fraction "
            "reject regardless of their action."
        ),
    ),
    length: int | None = typer.Option(
        None,
        "--length",
        "-l",
        min=19,
        max=23,
        help=f"siRNA length in nucleotides (default: {default_for('sirna_length')})",
    ),
    top_n: int | None = typer.Option(
        None,
        "--top-n",
        "-n",
        min=1,
        help=(
            "Cap how many top-ranked candidates are reported (default: no cap, report all). "
            "All candidates are generated and screened regardless."
        ),
    ),
    gc_min: float | None = typer.Option(
        None,
        "--gc-min",
        min=0.0,
        max=100.0,
        help=(
            f"Minimum GC content percentage (default: {default_for('gc_min')}; "
            f"--design-mode mirna defaults to {mirna_preset_default_for('gc_min')})"
        ),
    ),
    gc_max: float | None = typer.Option(
        None,
        "--gc-max",
        min=0.0,
        max=100.0,
        help=(
            f"Maximum GC content percentage (default: {default_for('gc_max')}; "
            f"--design-mode mirna defaults to {mirna_preset_default_for('gc_max')}). Stating it "
            "keeps it in every design mode."
        ),
    ),
    max_poly_runs: int | None = typer.Option(
        None,
        "--max-poly-runs",
        min=1,
        help=f"Maximum consecutive identical nucleotides (default: {default_for('max_poly_runs')})",
    ),
    min_asymmetry: float | None = typer.Option(
        None,
        "--min-asymmetry",
        min=0.3,
        max=1.0,
        help=(f"Thermodynamic asymmetry floor gating LOW_ASYMMETRY (default: {default_for('min_asymmetry_score')})"),
    ),
    max_paired_fraction: float | None = typer.Option(
        None,
        "--max-paired-fraction",
        min=0.0,
        max=1.0,
        help=(f"Guide self-structure ceiling gating EXCESS_PAIRING (default: {default_for('max_paired_fraction')})"),
    ),
    min_empirical: float | None = typer.Option(
        None,
        "--min-empirical",
        help=(
            "Empirical design-rule floor gating LOW_EMPIRICAL_SCORE "
            f"(default: {default_for('min_empirical_score')}, which rejects nothing)"
        ),
    ),
    plfold_window: int | None = typer.Option(
        None,
        "--plfold-window",
        min=20,
        max=1000,
        help=(f"RNAplfold averaging window W for target-site accessibility (default: {default_for('plfold_window')})"),
    ),
    plfold_max_bp_span: int | None = typer.Option(
        None,
        "--plfold-max-bp-span",
        min=10,
        max=1000,
        help=(
            f"RNAplfold maximum base-pair span L (default: {default_for('plfold_max_bp_span')}); "
            "must not exceed --plfold-window"
        ),
    ),
    accessibility_log_floor: float | None = typer.Option(
        None,
        "--accessibility-log-floor",
        help=(
            "log10 probability treated as zero target-site accessibility "
            f"(default: {default_for('accessibility_log_floor')})"
        ),
    ),
    genome_index: Path | None = typer.Option(
        None,
        "--genome-index",
        help="Genome index for off-target analysis",
    ),
    snp_file: Path | None = typer.Option(
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
    modification_pattern: str | None = typer.Option(
        None,
        "--modifications",
        "-m",
        help=(
            "Chemical modification pattern: standard_2ome, minimal_terminal, maximal_stability, none "
            f"(default: {default_for('modification_pattern')})"
        ),
    ),
    overhang: str | None = typer.Option(
        None,
        "--overhang",
        help=(
            f"Overhang sequence, dTdT for DNA or UU for RNA (default: {default_for('default_overhang')}; "
            f"--design-mode mirna defaults to {mirna_preset_default_for('default_overhang')})"
        ),
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        "-v",
        help="Enable verbose output",
    ),
) -> None:
    """Design siRNA candidates from a transcript FASTA file.

    Outputs a TSV/CSV-like table of candidates, optionally including secondary
    structure scoring, off-target checks, and chemical modification annotations.
    """
    try:
        actions = _parse_filter_actions(filter_action)
    except RunPolicyError as exc:
        _fail_with_config_error(str(exc))

    stated: dict[str, Any] = _stated_settings(
        ctx,
        {
            "gc_min": ("gc_min", gc_min),
            "gc_max": ("gc_max", gc_max),
            "max_poly_runs": ("max_poly_runs", max_poly_runs),
            "min_asymmetry_score": ("min_asymmetry", min_asymmetry),
            "max_paired_fraction": ("max_paired_fraction", max_paired_fraction),
            "min_empirical_score": ("min_empirical", min_empirical),
            "plfold_window": ("plfold_window", plfold_window),
            "plfold_max_bp_span": ("plfold_max_bp_span", plfold_max_bp_span),
            "accessibility_log_floor": ("accessibility_log_floor", accessibility_log_floor),
            "sirna_length": ("length", length),
            "top_n": ("top_n", top_n),
            "modification_pattern": ("modification_pattern", modification_pattern),
            "default_overhang": ("overhang", overhang),
        },
    )
    if skip_structure:
        stated["predict_structure"] = False
    if "modification_pattern" in stated:
        stated["apply_modifications"] = str(stated["modification_pattern"]).lower() != "none"
    policy = _resolve_policy_or_exit(
        entry_point=EntryPoint.DESIGN_COMMAND,
        design_mode=design_mode if _option_was_stated(ctx, "design_mode") else None,
        config_file=policy_config,
        stated=stated,
        filter_actions=actions,
        legacy_skip_screening=skip_off_targets or None,
        passthrough={
            "genome_index": str(genome_index) if genome_index else None,
            "snp_file": str(snp_file) if snp_file else None,
        },
    )
    mode_enum = policy.design_mode
    if mode_enum == DesignMode.ZFN:
        console.print(
            "❌ Error: --design-mode zfn is not supported in the 'design' command. Use 'sirnaforge zfn' instead.",
            style="red",
        )
        raise typer.Exit(1)

    parameters = policy.design_parameters
    filters = parameters.filters

    console.print(
        Panel.fit(
            f"🧬 [bold blue]siRNAforge Toolkit[/bold blue]\n"
            f"Design Mode: [cyan]{mode_enum.value}[/cyan]\n"
            f"Input: [cyan]{input_file}[/cyan]\n"
            f"Output: [cyan]{output}[/cyan]\n"
            f"Run Mode: [cyan]{policy.run_mode.value}[/cyan] "
            f"([magenta]{policy.profile.name} {policy.profile.version}, experimental[/magenta])\n"
            f"Length: [yellow]{parameters.sirna_length}[/yellow] nt\n"
            f"GC range: [yellow]{filters.gc_min:.1f}%-{filters.gc_max:.1f}%[/yellow]\n"
            f"Reported candidates: [yellow]{parameters.top_n if parameters.top_n is not None else 'all (uncapped)'}[/yellow]\n"
            f"Modifications: [magenta]{parameters.modification_pattern}[/magenta]\n"
            f"Overhang: [magenta]{parameters.default_overhang}[/magenta]",
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
            from sirnaforge.core.design import MiRNADesigner  # noqa: PLC0415

            task1 = progress.add_task("Loading sequences...", total=None)

            # Select designer based on design mode. The actions travel with the thresholds:
            # DesignParameters carries only the numbers, so without them a --filter-action reached the
            # manifest and was silently ignored by every design gate on this command.
            from sirnaforge.workflow import _policy_filter_actions  # noqa: PLC0415

            design_actions = _policy_filter_actions(policy)
            designer = (
                MiRNADesigner(parameters, filter_actions=design_actions)
                if mode_enum == DesignMode.MIRNA
                else SiRNADesigner(parameters, filter_actions=design_actions)
            )

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
                    f"{ranking_score(candidate):.1f}",
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
    """Validate a FASTA file and report basic statistics.

    This performs lightweight validation (parseable FASTA, presence of
    sequences, and common issues like short/ambiguous sequences).
    """
    try:
        with console.status("Validating FASTA file..."):
            sequences = _parse_fasta_seqrecords(input_file)

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
    """Show CLI version and author information."""
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
    """Print the parameter values a default run resolves to."""
    # Resolved rather than constructed: a run takes the profile's numbers, and gc_min is the one
    # place the profile knowingly differs from the model field default.
    policy = resolve_run_policy(entry_point=EntryPoint.SCREENING_WORKFLOW)
    default_params = policy.design_parameters

    console.print("[bold blue]Resolved Design Parameters:[/bold blue]")
    console.print(
        f"  [dim]profile {policy.profile.name} {policy.profile.version} "
        f"({'experimental' if policy.profile.experimental else 'calibrated'}), "
        f"run mode {policy.run_mode.value}[/dim]\n"
    )

    # Basic parameters
    console.print("[cyan]Basic Parameters:[/cyan]")
    console.print(f"  siRNA length: {default_params.sirna_length} nt")
    console.print(
        f"  Reported candidates: {default_params.top_n if default_params.top_n is not None else 'all (uncapped)'}"
    )

    # Filtering criteria
    console.print("\n[cyan]Filtering Criteria:[/cyan]")
    filters = default_params.filters
    console.print(f"  GC content: {filters.gc_min}% - {filters.gc_max}%")
    console.print(f"  Max poly runs: {filters.max_poly_runs}")
    console.print(f"  Max paired fraction: {filters.max_paired_fraction}")

    # Gates, with the three ways one can be in force or not
    console.print("\n[cyan]Gates:[/cyan]")
    for resolved in policy.filters:
        descriptor = resolved.descriptor
        threshold = "no threshold" if descriptor.threshold is None else descriptor.threshold
        console.print(
            f"  {resolved.filter_id}: {descriptor.comparator.value} {threshold} "
            f"([yellow]{descriptor.action.value}[/yellow], {descriptor.stage.value})"
        )

    # Scoring weights: one hand-authored vector per stage/mode, each summing to 1.0 and never
    # rescaled at runtime.
    console.print("\n[cyan]Scoring Weight Vectors:[/cyan]")
    for vector in default_params.scoring.all_vectors():
        console.print(f"  [bold]{vector.name}[/bold]")
        for term, weight in vector.as_mapping().items():
            console.print(f"    {term}: {weight}")
    console.print(
        "  [dim]Reported but not scored: empirical (min_empirical_score gate), isoform_coverage "
        "(optional gate), conservation, paired_fraction (EXCESS_PAIRING gate)[/dim]"
    )


@app_command()
def cache(
    clear: bool = typer.Option(False, "--clear", help="Clear all cached databases (miRNA + transcriptomes)"),
    clear_mirna: bool = typer.Option(False, "--clear-mirna", help="Clear only miRNA databases"),
    clear_transcriptome: bool = typer.Option(False, "--clear-transcriptome", help="Clear only transcriptomes"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show what would be deleted without actually deleting"),
    info: bool = typer.Option(False, "--info", help="Show cache information for all databases"),
) -> None:
    """Inspect and clear the unified reference cache.

    This command can display cache statistics and/or delete cached assets for
    miRNA databases and transcriptomes.
    """
    from sirnaforge.utils.unified_cache import UnifiedCacheManager  # noqa: PLC0415

    if not any([clear, clear_mirna, clear_transcriptome, dry_run, info]):
        console.print("❓ [yellow]No action specified. Use --info, --clear, or specific clear options[/yellow]")
        console.print("   Example: sirnaforge cache --info")
        console.print("   Example: sirnaforge cache --clear-transcriptome --dry-run")
        return

    manager = UnifiedCacheManager()

    if info or dry_run:
        # Display cache information using unified manager
        cache_info = manager.get_info()

        if "mirna" in cache_info:
            stats = cache_info["mirna"]
            console.print("\n📊 [bold blue]miRNA Database Cache:[/bold blue]")
            console.print(f"  Directory: [cyan]{stats['cache_directory']}[/cyan]")
            console.print(f"  Files: [green]{stats['total_files']}[/green]")
            console.print(f"  Size: [yellow]{stats['total_size_mb']:.2f} MB[/yellow]")
            console.print(f"  TTL: [magenta]{stats['cache_ttl_days']} days[/magenta]")

        if "transcriptome" in cache_info:
            stats = cache_info["transcriptome"]
            console.print("\n📚 [bold blue]Transcriptome Cache:[/bold blue]")
            console.print(f"  Directory: [cyan]{stats['cache_directory']}[/cyan]")
            console.print(f"  Files: [green]{stats['total_files']}[/green]")
            console.print(f"  Size: [yellow]{stats['total_size_mb']:.2f} MB[/yellow]")
            console.print(f"  TTL: [magenta]{stats['cache_ttl_days']} days[/magenta]")

        # Show total
        totals = manager.get_total_stats()
        console.print("\n📈 [bold cyan]Total Cache:[/bold cyan]")
        console.print(f"  Files: [green]{totals['total_files']}[/green]")
        console.print(f"  Size: [yellow]{totals['total_size_mb']:.2f} MB[/yellow]")

    if dry_run:
        console.print("\n🔍 [bold yellow]Clear Preview (dry run):[/bold yellow]")

        results = manager.clear(
            clear_mirna=clear or clear_mirna,
            clear_transcriptome=clear or clear_transcriptome,
            dry_run=True,
        )

        for component, result in results.items():
            console.print(f"\n  {component.title()}:")
            console.print(f"    Files to delete: [red]{result['files_deleted']}[/red]")
            console.print(f"    Size to free: [yellow]{result['size_freed_mb']:.2f} MB[/yellow]")

    elif clear or clear_mirna or clear_transcriptome:
        console.print("\n🧹 [bold green]Clearing Cache:[/bold green]")

        results = manager.clear(
            clear_mirna=clear or clear_mirna,
            clear_transcriptome=clear or clear_transcriptome,
            dry_run=False,
        )

        for component, result in results.items():
            console.print(f"\n  {component.title()}:")
            console.print(f"    Files deleted: [red]{result['files_deleted']}[/red]")
            console.print(f"    Size freed: [yellow]{result['size_freed_mb']:.2f} MB[/yellow]")
            console.print(f"    Status: [green]{result['status']}[/green]")


# Create sequences subcommand group
sequences_app = typer.Typer(help="Manage siRNA sequences and metadata")
app.add_typer(sequences_app, name="sequences")
sequences_command = command_decorator_typed(sequences_app.command)


class SequencesShowError(RuntimeError):
    """Raised when sequence display/formatting input is invalid."""


def _load_fasta_records(input_file: Path) -> list[SeqRecord]:
    """Load FASTA records from disk.

    Raises:
        SequencesShowError: If the file contains no records.
    """
    records = _parse_fasta_seqrecords(input_file)
    if not records:
        raise SequencesShowError("No sequences found in file")
    return records


def _filter_records_by_id(records: list[SeqRecord], sequence_id: str) -> list[SeqRecord]:
    """Filter FASTA records by record id.

    Raises:
        SequencesShowError: If no matching records are found.
    """
    filtered = [record for record in records if record.id == sequence_id]
    if not filtered:
        raise SequencesShowError(f"Sequence ID '{sequence_id}' not found")
    return filtered


def _metadata_value_to_json(value: Any) -> Any:
    """Convert parsed FASTA header metadata into JSON-serializable values."""
    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if hasattr(value, "value"):
        return value.value
    if isinstance(value, list):
        json_items: list[Any] = []
        items: list[object] = list(value)
        for item in items:
            json_items.append(_metadata_value_to_json(item))
        return json_items
    return value


def _records_to_json(records: list[SeqRecord]) -> str:
    """Render FASTA record header metadata as a JSON string."""
    payload: list[dict[str, Any]] = []
    for record in records:
        metadata = parse_header(record)
        payload.append({key: _metadata_value_to_json(val) for key, val in metadata.items()})
    return json.dumps(payload, indent=2)


def _summarize_modifications(metadata: dict[str, Any]) -> str:
    """Summarize chemical modifications from parsed header metadata."""
    raw_mods = metadata.get("chem_mods")
    mods: list[object] = list(raw_mods) if isinstance(raw_mods, list) else []
    summary: list[str] = []
    for mod in mods:
        mod_type = getattr(mod, "type", str(mod))
        positions_value = getattr(mod, "positions", [])
        if isinstance(positions_value, list):
            positions_list: list[object] = list(positions_value)
            length: int | Any = len(positions_list)
        else:
            length = positions_value
        summary.append(f"{mod_type}({length})")
    return ", ".join(summary)


def _print_records_fasta(records: list[SeqRecord]) -> None:
    """Print records as FASTA to stdout."""
    for record in records:
        console.print(f">{record.description}")
        console.print(str(record.seq))


def _print_records_table(records: list[SeqRecord], input_file: Path) -> None:
    """Print records as a Rich table with parsed header metadata."""
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
    sequence_id: str | None = typer.Option(
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
    """Show sequences from a FASTA file in table, JSON, or FASTA format.

    Use ``--id`` to select a single record. ``--format`` controls output:
    ``table`` (default), ``json`` (header metadata only), or ``fasta``.
    """
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
    output: Path | None = typer.Option(
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
    """Merge metadata from a JSON file into FASTA headers.

    The JSON is expected to conform to the project metadata schema used by the
    modification/annotation utilities.
    """
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


internal_app = typer.Typer(help="Internal workflow commands", hidden=True)
app.add_typer(internal_app, name="_internal", hidden=True)
internal_command = command_decorator_typed(internal_app.command)


@internal_command("zfn-make-shards")
def internal_zfn_make_shards(
    genome_fasta: Path = typer.Option(..., "--genome-fasta", exists=True, file_okay=True, dir_okay=False),
    left_half_site: str = typer.Option(..., "--left-half-site"),
    right_half_site: str = typer.Option(..., "--right-half-site"),
    spacer_lengths: str = typer.Option(..., "--spacer-lengths"),
    max_mismatches: int = typer.Option(..., "--max-mismatches"),
    sharding_enabled: str = typer.Option("true", "--sharding-enabled"),
    shard_chunk_mb: float = typer.Option(20.0, "--shard-chunk-mb"),
    shard_overlap_bp: int = typer.Option(50, "--shard-overlap-bp"),
    shard_chromosomes: str = typer.Option("", "--shard-chromosomes"),
    output: Path = typer.Option(Path("zfn_shards.tsv"), "--output"),
) -> None:
    """Build ZFN shard manifest for Nextflow execution."""
    make_zfn_shard_manifest(
        genome_fasta=genome_fasta,
        left_half_site=left_half_site,
        right_half_site=right_half_site,
        spacer_lengths=spacer_lengths,
        max_mismatches=max_mismatches,
        sharding_enabled=sharding_enabled,
        shard_chunk_mb=shard_chunk_mb,
        shard_overlap_bp=shard_overlap_bp,
        shard_chromosomes=shard_chromosomes,
        output_tsv=output,
    )


@internal_command("zfn-build-search-index")
def internal_zfn_build_search_index(
    genome_fasta: Path = typer.Option(..., "--genome-fasta", exists=True, file_okay=True, dir_okay=False),
    search_backend: ZFNSearchBackend = typer.Option(ZFNSearchBackend.FM_INDEX, "--search-backend"),
    output_dir: Path | None = typer.Option(None, "--output-dir"),
) -> None:
    """Build a persisted ZFN search-space index bundle for indexed backends."""
    summary = build_zfn_search_index(
        backend=search_backend,
        genome_fasta=genome_fasta,
        output_dir=output_dir,
    )
    console.print(json.dumps(summary, indent=2, sort_keys=True))


@internal_command("zfn-search-shard")
def internal_zfn_search_shard(
    shard_id: str = typer.Option(..., "--shard-id"),
    shard_chrom: str = typer.Option(..., "--shard-chrom"),
    scan_start_1: int = typer.Option(..., "--scan-start-1"),
    scan_end_1: int = typer.Option(..., "--scan-end-1"),
    core_start_1: int | None = typer.Option(
        None, "--core-start-1", help="Core window start (1-based). Defaults to scan-start-1."
    ),
    core_end_1: int | None = typer.Option(
        None, "--core-end-1", help="Core window end (1-based). Defaults to scan-end-1."
    ),
    shard_max_mismatches: int = typer.Option(..., "--shard-max-mismatches"),
    left_half_site: str = typer.Option(..., "--left-half-site"),
    right_half_site: str = typer.Option(..., "--right-half-site"),
    genome_fasta: Path = typer.Option(..., "--genome-fasta", exists=True, file_okay=True, dir_okay=False),
    search_backend: ZFNSearchBackend = typer.Option(ZFNSearchBackend.EXHAUSTIVE_PYTHON, "--search-backend"),
    search_space_index: Path | None = typer.Option(None, "--search-space-index"),
    algorithm: ZFNAlgorithm = typer.Option(..., "--algorithm"),
    dimer_mode: DimerMode = typer.Option(..., "--dimer-mode"),
    spacer_lengths: str = typer.Option(..., "--spacer-lengths"),
    annotation_file: Path | None = typer.Option(None, "--annotation-file"),
    output_sites_csv: Path = typer.Option(Path("offtarget_sites.csv"), "--output-sites-csv"),
    output_summary_json: Path = typer.Option(Path("candidate_summary.json"), "--output-summary-json"),
) -> None:
    """Run one shard-scoped ZFN search and emit shard artifacts."""
    run_zfn_shard_search(
        shard_id=shard_id,
        shard_chrom=shard_chrom,
        scan_start_1=scan_start_1,
        scan_end_1=scan_end_1,
        core_start_1=core_start_1,
        core_end_1=core_end_1,
        shard_max_mismatches=shard_max_mismatches,
        left_half_site=left_half_site,
        right_half_site=right_half_site,
        genome_fasta=genome_fasta,
        search_backend=search_backend,
        search_space_index=search_space_index,
        algorithm=algorithm,
        dimer_mode=dimer_mode,
        spacer_lengths=spacer_lengths,
        annotation_file=annotation_file,
        output_sites_csv=output_sites_csv,
        output_summary_json=output_summary_json,
    )


@internal_command("zfn-aggregate-shards")
def internal_zfn_aggregate_shards(
    shard_csv_glob: str = typer.Option("offtarget_sites_*.csv", "--shard-csv-glob"),
    output_sites_csv: Path = typer.Option(Path("offtarget_sites.csv"), "--output-sites-csv"),
    output_summary_json: Path = typer.Option(Path("candidate_summary.json"), "--output-summary-json"),
) -> None:
    """Aggregate shard-level ZFN outputs into final ranked outputs."""
    aggregate_zfn_shard_results(
        shard_csv_glob=shard_csv_glob,
        output_sites_csv=output_sites_csv,
        output_summary_json=output_summary_json,
    )


@app_command()
def report(
    run_dir: Path = typer.Argument(..., help="A completed run output directory"),
    output: Path = typer.Option(Path("report.html"), "--output", "-o", help="Where to write the report"),
) -> None:
    """Build a single self-contained HTML report over a finished run.

    Works on any finished or archived run: it reads candidates_all.csv, the aggregated hit tables and
    manifest.json, and needs no live pipeline state. The report classifies nothing and applies no
    threshold of its own -- gate descriptors come from the resolved policy and hit_class from the
    published table -- so it cannot disagree with the run it describes.
    """
    try:
        payload = build_payload(run_dir)
    except ReportInputError as exc:
        console.print(f"❌ [red]{exc}[/red]")
        raise typer.Exit(1) from exc

    written = write_report(payload, output)
    counts = payload.run["status_counts"]
    console.print(
        Panel.fit(
            f"📄 [bold blue]{written}[/bold blue]  ({written.stat().st_size / 1e6:.1f} MB, self-contained)\n"
            f"{payload.run['guides']:,} guides from {payload.run['candidate_rows']:,} candidate rows\n"
            f"[green]{counts['pass']} pass[/green] · "
            f"[yellow]{counts['warn']} pass with a warning[/yellow] · "
            f"[yellow]{counts['unknown']} not established[/yellow] · [red]{counts['fail']} fail[/red]",
            title=f"Report — {payload.run['gene_query']}",
        )
    )
    for caveat in payload.caveats:
        console.print(f"   ⚠️  {caveat}")


if __name__ == "__main__":
    app()
