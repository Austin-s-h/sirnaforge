"""Command-line entry points used by embedded Nextflow modules."""

import csv
import json
import shutil
from collections.abc import Sequence
from dataclasses import fields
from enum import Enum
from pathlib import Path
from typing import Any

from sirnaforge.config import DEFAULT_MIRNA_CANONICAL_SPECIES
from sirnaforge.core.off_target import (
    aggregate_mirna_results,
    aggregate_offtarget_results,
    build_bwa_index,
    parse_fasta_file,
    run_bwa_alignment_analysis,
    validate_index_files,
)
from sirnaforge.core.screening_evidence import (
    EVIDENCE_FILE_SUFFIX,
    EvidenceProducer,
    Reconciliation,
    build_plan,
    collect_evidence,
    guide_set_digest,
    reconcile,
    write_evidence,
    write_reconciliation,
)
from sirnaforge.core.transcript_seed import (
    ALL_SEED_CLASSES,
    SeedScanScope,
    SiteRegion,
    TranscriptSeedScanResult,
    TranscriptSeedSite,
    scan_transcript_seed_sites,
    transcript_seed_evidence_entry,
)
from sirnaforge.data.species_registry import normalize_species_name
from sirnaforge.data.transcript_index import TranscriptGeneIndex
from sirnaforge.models.evidence import EvidenceStatus, ObservedCounts, ScreeningEvidenceEntry, ScreeningPlan
from sirnaforge.models.policy import ScreeningChannel
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

DEFAULT_MIRNA_SPECIES_ARGUMENT = ",".join(DEFAULT_MIRNA_CANONICAL_SPECIES)

#: Default per-guide retention cap for the transcript-seed scan, mirroring ``workflow.py``'s
#: ``_TRANSCRIPT_SEED_SITE_CAP``. Not imported from there: ``workflow.py`` pulls in the whole design
#: stack, and this module is imported inside a Nextflow task that has no reason to. A test pins the
#: two equal so the two producers of ``<species>_transcript_seed_sites.tsv`` cannot drift apart.
#: A 7mer is expected roughly once per 16 kb, so an uncapped full-cDNA scan publishes thousands of
#: rows per guide; when the cap truncates, the scan reports CENSORED and every count it publishes is
#: declared a lower bound, so a ceiling read against it cannot claim to have been respected.
DEFAULT_TRANSCRIPT_SEED_SITE_CAP = 200

#: Guide-set digest stamped on a locally-derived fallback plan, used only when no ``evidence_plan``
#: file was threaded through a call. #100's real digests are 16-character sha256 prefixes, so this
#: sentinel never collides with one: a fallback plan entry can therefore only ever explain an
#: envelope that was NOT observed (reconciling to failed), never falsely match one that was.
_NO_PLAN_DIGEST_PLACEHOLDER = "no-evidence-plan-supplied"


def _reconcile_species(
    *,
    evidence_plan: str | None,
    transcriptome_species: Sequence[str] = (),
    mirna_species: Sequence[str] = (),
    transcript_seed_species: Sequence[str] = (),
    search_root: str | Path,
) -> Reconciliation:
    """Reconcile the species this call expected against whatever evidence was actually staged.

    Prefers the real plan threaded in via ``evidence_plan`` (#100's expected-plan contract, built by
    ``workflow.py`` before any reference can be dropped). With none supplied -- a bare Nextflow run,
    or a pipeline invocation that aborted before a plan file existed -- falls back to a plan built
    from this call's own species arguments, so a species with no evidence still produces a
    machine-readable failure record instead of vanishing with only a ``final_summary.txt`` behind.
    The fallback plan's digest is borrowed from whatever evidence was actually observed (a run
    screens one guide set at a time, by construction) so a real envelope still joins against it;
    with nothing observed at all there is no digest to borrow, and the placeholder is used purely
    to explain the absence -- it can never falsely match a real envelope written later.

    ``transcript_seed_species`` (#101) is treated exactly like the other two: a requested unit that
    published no envelope reconciles FAILED rather than vanishing. Without it in the fallback plan a
    seed envelope would arrive as ``unplanned`` -- reported, but with no way to notice the species
    that produced none -- so a species whose scan died would read as a channel nobody asked about.
    """
    observed = collect_evidence(search_root)
    if evidence_plan is not None:
        plan = ScreeningPlan.model_validate_json(Path(evidence_plan).read_text())
    else:
        digest = observed[0].entry.guide_set_digest if observed else _NO_PLAN_DIGEST_PLACEHOLDER
        plan = build_plan(
            guide_set_digest=digest,
            transcriptome=[(species, None) for species in transcriptome_species],
            mirna_species=mirna_species,
            transcript_seed=[(species, None) for species in transcript_seed_species],
        )
    return reconcile(plan, observed)


def _staged_channel_evidence(search_root: Path, destination: Path, *, channel: ScreeningChannel) -> list[Path]:
    """Every per-unit envelope for one channel under ``search_root``, excluding ``destination``'s own.

    Parameterised by channel (#101) because the envelope filename is
    ``<channel>_<species>_evidence.json``, so "the miRNA envelopes" and "the transcript-seed
    envelopes" differ only in that prefix -- and a glob that dropped the prefix would stage a third
    channel's evidence into the miRNA aggregate's directory, where its counts are not welcome.
    Resolved eagerly, and with anything already inside ``destination`` skipped, so a copy this call
    just made can never be re-found and copied onto itself.
    """
    pattern = f"{channel.value}_*{EVIDENCE_FILE_SUFFIX}"
    return [path for path in search_root.rglob(pattern) if path.parent != destination]


def _staged_mirna_evidence(search_root: Path, destination: Path) -> list[Path]:
    """Every miRNA per-unit evidence envelope under ``search_root``, excluding ``destination``'s own.

    The miRNA aggregate searches the results directory it is handed, so the envelopes have to be
    copied into it (#100). Named separately from :func:`_staged_channel_evidence` because this
    channel is the one with a copy destination at all.
    """
    return _staged_channel_evidence(search_root, destination, channel=ScreeningChannel.MIRNA_SEED)


def build_bwa_index_cli(fasta_file: str, species: str, output_dir: str = ".") -> dict[str, Any]:
    """Build a BWA-MEM2 index for one screening reference.

    Args:
        fasta_file: Path to input FASTA file
        species: Species identifier
        output_dir: Directory to write index files

    Returns:
        Dictionary with index prefix path
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Build index with species-specific prefix
    index_prefix = output_path / f"{species}_index"
    result_prefix = build_bwa_index(fasta_file=fasta_file, index_prefix=str(index_prefix))

    logger.info(f"Built BWA index for {species}: {result_prefix}")

    return {
        "species": species,
        "index_prefix": str(result_prefix),
        "index_files": list(output_path.glob(f"{species}_index*")),
    }


def offtarget_analysis_cli(
    species: str,
    index_prefix: str,
    candidates_file: str,
    output_dir: str = ".",
    max_hits: int | None = None,
    bwa_k: int = 12,
    bwa_T: int = 15,
    seed_start: int = 2,
    seed_end: int = 8,
) -> dict[str, Any]:
    """Align one species' reference for OFFTARGET_ANALYSIS, or record that it could not be aligned.

    A named index prefix that does not resolve to a usable BWA-MEM2 index used to reach the aligner,
    which logged a failure, returned no hits, and published a header-only table -- a completed screen
    with zero hits, indistinguishable from a clean one. The species is now published as an EMPTY
    analysis file plus a failed summary, which ``aggregate_offtarget_results`` reports as a per-species
    rejection rather than as a clean result. Either way, this function owns the failed-vs-completed
    decision, so it also writes the #100 ``transcriptome_<species>_evidence.json`` envelope on both
    branches -- directly here on the missing-index branch, and via ``run_bwa_alignment_analysis``'s
    ``evidence_dir`` on the aligned branch.

    Args:
        species: Species identifier this reference belongs to
        index_prefix: BWA-MEM2 index prefix as resolved for this species
        candidates_file: FASTA of candidate guides to screen
        output_dir: Directory to write ``<species>_analysis.tsv`` and ``<species>_summary.json``
        max_hits: Maximum hits per candidate (``None`` = exhaustive)
        bwa_k: BWA seed length
        bwa_T: BWA minimum score threshold
        seed_start: Seed region start (1-based)
        seed_end: Seed region end (1-based)

    Returns:
        Dictionary with the per-species status and the published file paths
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    analysis_file = output_path / f"{species}_analysis.tsv"
    summary_file = output_path / f"{species}_summary.json"

    if not validate_index_files(index_prefix, "bwa-mem2"):
        error = (
            f"No usable BWA-MEM2 index at prefix '{index_prefix}': the index files are missing, empty, "
            "or the path points at something that is not an index (a FASTA left behind by a failed "
            "index build looks exactly like this). No alignment was attempted."
        )
        logger.error(f"{species}: {error}")
        # Deliberately empty, not header-only: a header-only table is a screen that found nothing.
        analysis_file.write_bytes(b"")
        with summary_file.open("w") as handle:
            json.dump({"species": species, "status": "failed", "error": error}, handle, indent=2)
        # This branch owns the failed-vs-completed decision, so it is the one that must write the
        # #100 evidence envelope: a species whose index never resolved must never reconcile as
        # completed just because the aligner was never reached.
        digest = guide_set_digest(candidates_file)
        try:
            submitted_guides: int | None = len(parse_fasta_file(candidates_file))
        except (OSError, ValueError):
            submitted_guides = None
        write_evidence(
            output_dir,
            producer=EvidenceProducer.OFFTARGET_ANALYSIS,
            entry=ScreeningEvidenceEntry(
                channel=ScreeningChannel.TRANSCRIPTOME,
                species=species,
                guide_set_digest=digest,
                status=EvidenceStatus.FAILED,
                submitted_guide_digest=digest,
                submitted_guides=submitted_guides,
                detail=error,
            ),
        )
        return {"species": species, "status": "failed", "error": error, "analysis_file": str(analysis_file)}

    run_bwa_alignment_analysis(
        candidates_file=candidates_file,
        index_prefix=index_prefix,
        species=species,
        output_dir=output_dir,
        max_hits=max_hits,
        bwa_k=bwa_k,
        bwa_T=bwa_T,
        seed_start=seed_start,
        seed_end=seed_end,
        evidence_dir=output_dir,
    )

    # run_bwa_alignment_analysis names its outputs after the candidates file; the process declares
    # them per species.
    stem = Path(candidates_file).stem
    (output_path / f"{stem}_{species}_analysis.tsv").replace(analysis_file)
    (output_path / f"{stem}_{species}_summary.json").replace(summary_file)

    return {"species": species, "status": "completed", "analysis_file": str(analysis_file)}


#: Columns of ``<species>_transcript_seed_sites.tsv``, taken from the dataclass the artifact is a
#: serialisation of rather than restated. ``TranscriptSeedSiteSchema`` is ``strict=True``, so an extra
#: or renamed column is a contract break; deriving the header from
#: :class:`~sirnaforge.core.transcript_seed.TranscriptSeedSite` means the two producers of this file
#: (``workflow.py``'s in-process scan and this module) cannot disagree about its shape by accident.
TRANSCRIPT_SEED_SITE_COLUMNS: tuple[str, ...] = tuple(field.name for field in fields(TranscriptSeedSite))

#: Marks a filename as belonging to the transcript-seed channel. Used to keep those files out of the
#: transcriptome/miRNA staging below: ``<species>_transcript_seed_summary.json`` matches the
#: ``*_summary.json`` glob that feeds ``aggregate_offtarget_results``, so without this check a seed
#: site count could be folded into the miRNA hit counters (#101). A third channel's numbers are not
#: more of the first two's.
_TRANSCRIPT_SEED_ARTIFACT_MARKER = "transcript_seed"


def _is_transcript_seed_artifact(name: str) -> bool:
    """Whether a staged filename belongs to the transcript-seed channel rather than to the other two."""
    return _TRANSCRIPT_SEED_ARTIFACT_MARKER in name.lower()


def _transcript_seed_site_row(site: TranscriptSeedSite) -> dict[str, Any]:
    """One site as a published row, enums flattened to their declared string values."""
    row: dict[str, Any] = {}
    for field in fields(site):
        value = getattr(site, field.name)
        row[field.name] = value.value if isinstance(value, Enum) else value
    return row


def _unreadable_reference_result(
    *,
    species: str,
    cdna_fasta: str,
    submitted_guides: int,
    max_sites_per_guide: int | None,
    reference_id: str | None,
) -> TranscriptSeedScanResult:
    """A FAILED scan for a cDNA reference this task cannot read.

    The analogue of ``offtarget_analysis_cli``'s missing-index branch, and it exists for the same
    reason: the alternative is a scan that searched nothing, found nothing and published zero sites,
    which is byte-identical to a clean screen. FAILED keeps the (transcript_seed, species) pair out of
    ``completed_pairs``, so the three transcript-seed ceilings report UNKNOWN instead of passing.
    """
    return TranscriptSeedScanResult(
        sites=(),
        counts=ObservedCounts(),
        status=EvidenceStatus.FAILED,
        detail=(
            f"no readable cDNA reference at '{cdna_fasta}' for {species}: the file is missing or is not "
            "a file, so no transcript sequence was searched and no sites were screened"
        ),
        scope=SeedScanScope(
            species=species,
            region=SiteRegion.UNKNOWN,
            classes=ALL_SEED_CLASSES,
            max_sites_per_guide=max_sites_per_guide,
            reference_id=reference_id,
        ),
        submitted_guides=submitted_guides,
        processed_guides=0,
    )


def _transcript_seed_summary_document(result: TranscriptSeedScanResult) -> dict[str, Any]:
    """The published ``<species>_transcript_seed_summary.json`` payload.

    Deliberately the same shape ``workflow.py``'s in-process scan publishes under
    ``transcript_seed_summary``, so a reader does not have to know which producer ran. It carries the
    scope as well as the counts, because a site count without the region, the classes and the cap is
    unreadable: the same guide scores an order of magnitude higher against a full cDNA set than
    against a canonical-transcript one.
    """
    return {
        "species": result.scope.species,
        "status": result.status.value,
        "region": result.scope.region.value,
        "seed_classes": sorted(cls.value for cls in result.scope.classes),
        "max_sites_per_guide": result.scope.max_sites_per_guide,
        "reference_id": result.scope.reference_id,
        "sites": len(result.sites),
        "submitted_guides": result.submitted_guides,
        "processed_guides": result.processed_guides,
        "counts": result.counts.model_dump(mode="json"),
        "detail": result.detail,
    }


def transcript_seed_analysis_cli(
    species: str,
    cdna_fasta: str,
    candidates_file: str,
    output_dir: str = ".",
    region_scope: str = SiteRegion.FULL_CDNA.value,
    max_sites_per_guide: int | None = DEFAULT_TRANSCRIPT_SEED_SITE_CAP,
    reference_id: str | None = None,
) -> dict[str, Any]:
    """Scan one species' ALREADY-MATERIALISED cDNA for complementary seed sites (#101).

    The task side of the third liability channel. ``cdna_fasta`` is a required argument and is the
    reference the alignment channel already resolved and Nextflow already staged: nothing is
    downloaded here, nothing is indexed here, and there is no code path that could reach a reference
    resolver. That is the point of the channel's placement -- the expensive thing about a transcript
    scan is the reference, and this one has already been paid for.

    Separate from known-miRNA resemblance and sharing no code path with it. The miRNA scanner searches
    the guide's seed window FORWARD against a miRNA database, which is right for its own question;
    pointed at cDNA the same search finds guide-identical, passenger-orientation windows at the same
    expected rate as real sites, so the count-based null would look identical while every identity was
    wrong. :func:`~sirnaforge.core.transcript_seed.scan_transcript_seed_sites` searches the
    reverse complement instead.

    This function owns the failed-vs-completed decision, so it always publishes all three artifacts on
    every branch:

    * a scan that ran writes the sites table (header only when it found nothing -- a measurement),
      a summary and a COMPLETE (or CENSORED, when the per-guide cap truncated) envelope;
    * a scan that could not answer as asked -- an unrecognised region scope, a region no annotation in
      this repository can resolve, or an unreadable reference -- writes a DELIBERATELY EMPTY table
      (0 bytes, not a header: a header-only table is a scan that found nothing), a summary naming the
      reason and a FAILED envelope. ``workflow.py``'s parser skips a 0-byte table for exactly this
      reason, so the refusal cannot be read back as a measured zero.

    Args:
        species: Species key this reference belongs to. Canonicalised, so every output file is named
            after the canonical species and the output declarations in the module are globs.
        cdna_fasta: The staged cDNA FASTA to search. Also the source of the transcript→gene index, so
            the gene-level counts cost no second file.
        candidates_file: FASTA of deduplicated candidate guides; its record ids are the guide ids, so
            a site joins the alignment table on the same key.
        output_dir: Where the table, summary and envelope are written.
        region_scope: A :class:`~sirnaforge.core.transcript_seed.SiteRegion` value. Anything but
            ``full_cdna`` is refused with a reason naming the annotation that was absent, and an
            unrecognised string is refused as ``unknown`` rather than silently widened to the whole
            cDNA -- a UTR request answered over the full cDNA is a fabricated UTR screen.
        max_sites_per_guide: Retention cap per guide. 0 or below means uncapped, because a Nextflow
            ``val`` cannot carry Python's ``None`` and a cap of 0 would otherwise censor every site.
        reference_id: Reference identity for site provenance, when the caller knows it.

    Returns:
        The published summary document, with the paths of the three files it wrote.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    canonical_species = normalize_species_name(species)
    cap = max_sites_per_guide if (max_sites_per_guide or 0) > 0 else None

    guides = parse_fasta_file(candidates_file)
    digest = guide_set_digest(candidates_file)

    try:
        region = SiteRegion(region_scope)
    except ValueError:
        # SiteRegion.UNKNOWN, whose refusal reads "an unknown region scope was requested, which no
        # annotation can satisfy". A typo therefore produces a machine-readable failed unit rather
        # than either an exception with no envelope or a silent scan of something else.
        logger.error(f"{canonical_species}: unrecognised transcript-seed region scope '{region_scope}'")
        region = SiteRegion.UNKNOWN

    if not Path(cdna_fasta).is_file():
        result = _unreadable_reference_result(
            species=canonical_species,
            cdna_fasta=cdna_fasta,
            submitted_guides=len(guides),
            max_sites_per_guide=cap,
            reference_id=reference_id,
        )
    else:
        # Built from the SAME staged FASTA, so the gene-level counts need no annotation download. An
        # index built from headers that carry no gene leaves every site unresolved, which the scan
        # reports as unresolved_gene_sites with distinct_genes a lower bound -- never as zero genes.
        gene_index = TranscriptGeneIndex()
        gene_index.build(canonical_species, Path(cdna_fasta))
        result = scan_transcript_seed_sites(
            {guide_id: sequence for guide_id, sequence in guides.items() if sequence},
            cdna_fasta,
            species=canonical_species,
            region_scope=region,
            classes=ALL_SEED_CLASSES,
            max_sites_per_guide=cap,
            gene_index=gene_index,
            reference_id=reference_id,
        )

    sites_file = output_path / f"{canonical_species}_transcript_seed_sites.tsv"
    if result.status is EvidenceStatus.FAILED:
        logger.error(f"{canonical_species}: transcript-seed scan could not answer as scoped: {result.detail}")
        sites_file.write_bytes(b"")
    else:
        with sites_file.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=TRANSCRIPT_SEED_SITE_COLUMNS, delimiter="\t")
            writer.writeheader()
            for site in result.sites:
                writer.writerow(_transcript_seed_site_row(site))

    summary = _transcript_seed_summary_document(result)
    summary_file = output_path / f"{canonical_species}_transcript_seed_summary.json"
    with summary_file.open("w") as handle:
        json.dump(summary, handle, indent=2)

    evidence_file = write_evidence(
        output_dir,
        producer=EvidenceProducer.TRANSCRIPT_SEED_ANALYSIS,
        entry=transcript_seed_evidence_entry(result, guide_set_digest=digest),
    )

    return {
        **summary,
        "sites_file": str(sites_file),
        "summary_file": str(summary_file),
        "evidence_file": str(evidence_file),
    }


def _transcript_seed_report(search_root: Path, requested_species: Sequence[str]) -> dict[str, Any]:
    """What the transcript-seed channel published, read off its own envelopes and kept separate.

    Reported per species and NEVER summed across them (#101): each species is its own evidence unit
    with its own reference, and site counts scale with reference size, so one total across two
    references would be a number with no denominator. The three aggregation units -- sites, distinct
    transcripts, distinct genes -- are likewise carried through as they were observed rather than
    collapsed, because none is derivable from another.
    """
    by_species: dict[str, Any] = {}
    for envelope in collect_evidence(search_root):
        entry = envelope.entry
        if entry.channel is not ScreeningChannel.TRANSCRIPT_SEED:
            continue
        by_species[entry.species] = {
            "status": entry.status.value,
            "producer": envelope.producer.value,
            "reference_id": entry.reference_id,
            "counts": entry.counts.model_dump(mode="json"),
            "submitted_guides": entry.submitted_guides,
            "processed_guides": entry.processed_guides,
            "detail": entry.detail,
        }
    return {
        "requested_species": list(requested_species),
        "by_species": by_species,
        "envelopes_observed": len(by_species),
    }


def aggregate_results_cli(  # noqa: PLR0912
    transcriptome_species: str,
    output_dir: str = ".",
    mirna_db: str | None = None,
    mirna_species: str | None = None,
    analysis_files: list[str] | None = None,
    summary_files: list[str] | None = None,
    *,
    evidence_plan: str | None = None,
    transcript_seed_species: str | None = None,
) -> dict[str, Any]:
    """Aggregate off-target analysis results from multiple candidates and references.

    Args:
        transcriptome_species: Comma-separated list of species screened
        output_dir: Directory to write aggregated results
        mirna_db: The database that provided the reference
        mirna_species: The species code for the matching miRNA
        analysis_files: Optional explicit list of staged analysis files to aggregate
        summary_files: Optional explicit list of staged summary files to aggregate
        evidence_plan: Path to the serialized #100 ``ScreeningPlan`` (``workflow.py``'s
            ``screening_plan.json``), additive and keyword-only. Reconciliation against it is
            written unconditionally on every return path, including the zero-files branch, so a
            pipeline that aborted before producing a single analysis file still leaves a
            machine-readable ``evidence.json`` behind instead of only a ``final_summary.txt``.
        transcript_seed_species: Comma-separated species the transcript-seed channel was asked to
            scan (#101), additive and keyword-only. On a real pipeline run the plan arrives through
            ``evidence_plan`` and already names those units, because ``workflow.py`` builds it; this
            argument is what lets the *fallback* plan name them too, so a direct caller of this entry
            point is not left with a weaker contract for the third channel than for the first two.

    Returns:
        Dictionary with aggregation statistics. The ``transcript_seed`` section is always present and
        always separate: its counts are never folded into ``mirna``.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    species_list = [s.strip() for s in transcriptome_species.split(",") if s.strip()]
    seed_species_list = [s.strip() for s in (transcript_seed_species or "").split(",") if s.strip()]

    # Collect all analysis and summary files from current directory (Nextflow stages them)
    current_dir = Path()
    resolved_analysis_files = (
        [Path(path) for path in analysis_files]
        if analysis_files is not None
        else list(current_dir.rglob("*_analysis.tsv")) + list(current_dir.rglob("mirna_analysis.tsv"))
    )
    resolved_summary_files = (
        [Path(path) for path in summary_files]
        if summary_files is not None
        else list(current_dir.rglob("*_summary.json")) + list(current_dir.rglob("mirna_summary.json"))
    )

    # #101: the transcript-seed channel's tables are refused here by name, not merely left unstaged by
    # the subworkflow. ``<species>_transcript_seed_summary.json`` matches the ``*_summary.json`` glob
    # above, and ``species in f.name`` matches it again in the per-species staging below, so a seed
    # summary reaching this task by any route -- a hand-supplied file list, a publishDir collision, a
    # future channel wiring -- would have been fed to ``aggregate_offtarget_results`` and its site
    # count folded into the transcriptome and miRNA hit totals. A third channel's numbers are not more
    # of the first two's; the seed channel is reported through its own section below.
    resolved_analysis_files = [path for path in resolved_analysis_files if not _is_transcript_seed_artifact(path.name)]
    resolved_summary_files = [path for path in resolved_summary_files if not _is_transcript_seed_artifact(path.name)]
    transcript_seed = _transcript_seed_report(current_dir, seed_species_list)

    logger.info(f"Found {len(resolved_analysis_files)} analysis files and {len(resolved_summary_files)} summary files")

    if resolved_analysis_files or resolved_summary_files:
        # Create a temporary results directory structure for aggregation
        results_dir = Path("temp_results")
        results_dir.mkdir(exist_ok=True)

        # Organize files by species (extract from filename)
        for species in species_list:
            species_dir = results_dir / species
            species_dir.mkdir(exist_ok=True)

            # Copy relevant files to species directory
            for f in resolved_analysis_files:
                name_lower = f.name.lower()
                if name_lower.startswith("mirna"):
                    continue
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

            for f in resolved_summary_files:
                name_lower = f.name.lower()
                if name_lower.startswith("mirna"):
                    continue
                if species in f.name:
                    shutil.copy(f, species_dir / f.name)

        # Also persist miRNA batch results (if present) under dedicated directory
        mirna_results_dir = results_dir / "mirna"
        mirna_results_dir.mkdir(exist_ok=True)
        mirna_analysis_files: list[Path] = []
        for f in resolved_analysis_files:
            if "mirna" not in f.name.lower():
                continue
            dest_name = f.name
            if not dest_name.endswith("_mirna_analysis.tsv"):
                dest_name = "batch_mirna_analysis.tsv"
            dest_path = mirna_results_dir / dest_name
            shutil.copy(f, dest_path)
            mirna_analysis_files.append(dest_path)

        for f in resolved_summary_files:
            if "mirna" not in f.name.lower():
                continue
            dest_name = f.name
            if not dest_name.endswith("_mirna_summary.json"):
                dest_name = "batch_mirna_summary.json"
            dest_path = mirna_results_dir / dest_name
            shutil.copy(f, dest_path)

        # #100: stage the per-unit miRNA evidence envelopes beside those TSVs. Only the analysis and
        # summary files used to be copied here, while ``aggregate_mirna_results`` searched the very
        # directory it was handed for envelopes -- so it found none on any real run and
        # ``species_screened`` was decided by its no-envelope fallback rather than by the run. The
        # call below now names the real search root explicitly; this copy is what keeps the staged
        # directory a self-contained unit, so a caller that leaves ``evidence_root`` unset still
        # reconciles against real envelopes instead of against nothing. Only the miRNA channel's
        # envelopes are copied, which is all that function reads.
        for envelope in _staged_mirna_evidence(current_dir, mirna_results_dir):
            shutil.copy(envelope, mirna_results_dir / envelope.name)

        # Run aggregation using core function
        result_path = aggregate_offtarget_results(
            results_dir=str(results_dir), output_dir=output_dir, transcriptome_species=transcriptome_species
        )

        logger.info(f"Aggregation completed: {result_path}")

        mirna_summary: dict[str, Any] | None = None
        if mirna_analysis_files and (mirna_db or mirna_species):
            resolved_mirna_db = mirna_db or "mirgenedb"
            resolved_mirna_species = mirna_species or DEFAULT_MIRNA_SPECIES_ARGUMENT
            mirna_output = aggregate_mirna_results(
                results_dir=str(mirna_results_dir),
                output_dir=output_dir,
                mirna_db=resolved_mirna_db,
                mirna_species=resolved_mirna_species,
                # #100: the envelopes live where the tasks wrote them -- this task's own working
                # directory, staged by AGGREGATE_RESULTS -- not in the staging directory above, so
                # the search root is named explicitly rather than inherited from ``results_dir``.
                evidence_root=current_dir,
            )
            mirna_summary = {
                "mirna_db": resolved_mirna_db,
                "mirna_species": resolved_mirna_species,
                "output_dir": str(mirna_output),
                "analysis_files_processed": len(mirna_analysis_files),
            }

        reconciliation = _reconcile_species(
            evidence_plan=evidence_plan,
            transcriptome_species=species_list,
            transcript_seed_species=seed_species_list,
            search_root=current_dir,
        )
        reconciliation_file = write_reconciliation(output_dir, reconciliation)

        return {
            "status": "completed",
            "analysis_files_processed": len(resolved_analysis_files),
            "summary_files_processed": len(resolved_summary_files),
            "species": species_list,
            "output_dir": str(result_path),
            "mirna": mirna_summary,
            "transcript_seed": transcript_seed,
            "evidence": reconciliation.evidence.model_dump(mode="json"),
            "reconciliation_file": str(reconciliation_file),
        }

    # Create empty final summary
    final_summary = output_path / "final_summary.txt"
    with final_summary.open("w") as handle:
        handle.write("No analysis results found to aggregate\n")

    logger.warning("No files to aggregate")

    # Optional output presence is not completion evidence (#100): even though nothing was staged to
    # aggregate, every species this call expected to see still gets a reconciled -- here, failed --
    # entry, so a pipeline that aborted before a single analysis file existed leaves the same
    # machine-readable record behind as one that completed and legitimately found nothing.
    reconciliation = _reconcile_species(
        evidence_plan=evidence_plan,
        transcriptome_species=species_list,
        transcript_seed_species=seed_species_list,
        search_root=current_dir,
    )
    reconciliation_file = write_reconciliation(output_dir, reconciliation)

    return {
        "status": "empty",
        "analysis_files_processed": 0,
        "summary_files_processed": 0,
        "species": [],
        "output_dir": str(output_path),
        # Reported on this path too: a transcript-seed scan can succeed on a run where no alignment or
        # miRNA table was produced at all, and the section is what says so rather than the run looking
        # like nothing whatsoever happened (#101).
        "transcript_seed": transcript_seed,
        "evidence": reconciliation.evidence.model_dump(mode="json"),
        "reconciliation_file": str(reconciliation_file),
    }


def aggregate_mirna_results_cli(
    mirna_db: str,
    mirna_species: str,
    results_dir: str = ".",
    output_dir: str = ".",
    *,
    evidence_plan: str | None = None,
) -> dict[str, Any]:
    """Aggregate miRNA seed analysis results from multiple candidates.

    Args:
        mirna_db: miRNA database name used for analysis
        mirna_species: Comma-separated list of species analyzed
        results_dir: Directory containing individual miRNA results
        output_dir: Directory to write aggregated results
        evidence_plan: Path to the serialized #100 ``ScreeningPlan``, additive and keyword-only.
            ``aggregate_results.nf`` does not call this entry point, but it gets the same
            unconditional reconciliation treatment as ``aggregate_results_cli`` so a caller of it
            is not left with a weaker evidence contract than the transcriptome channel.

    Returns:
        Dictionary with aggregation statistics
    """
    logger.info(f"Aggregating miRNA results from {results_dir}")

    # Run aggregation using core function
    result_path = aggregate_mirna_results(
        results_dir=results_dir, output_dir=output_dir, mirna_db=mirna_db, mirna_species=mirna_species
    )

    # Load summary to get statistics
    summary_file = result_path / "combined_mirna_summary.json"
    stats = {}
    if summary_file.exists():
        with summary_file.open() as f:
            stats = json.load(f)

    logger.info(f"miRNA aggregation completed: {result_path}")

    reconciliation = _reconcile_species(
        evidence_plan=evidence_plan,
        mirna_species=[s.strip() for s in mirna_species.split(",") if s.strip()],
        search_root=results_dir,
    )
    reconciliation_file = write_reconciliation(output_dir, reconciliation)

    return {
        "status": "completed",
        "mirna_database": mirna_db,
        "species": mirna_species.split(","),
        "total_hits": stats.get("total_mirna_hits", 0),
        "candidates_analyzed": stats.get("total_candidates", 0),
        "output_dir": str(result_path),
        "evidence": reconciliation.evidence.model_dump(mode="json"),
        "reconciliation_file": str(reconciliation_file),
    }
