"""Workflow integration for variant targeting."""

import json
from pathlib import Path

from sirnaforge.data.variant_resolver import VariantResolver
from sirnaforge.models.variant import ClinVarSignificance, VariantMode, VariantRecord
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class VariantWorkflowConfig:
    """Configuration for variant targeting in the workflow."""

    def __init__(
        self,
        variant_ids: list[str] | None = None,
        vcf_file: Path | None = None,
        variant_mode: VariantMode = VariantMode.AVOID,
        min_af: float = 0.01,
        clinvar_filter_levels: list[ClinVarSignificance] | None = None,
        assembly: str = "GRCh38",
        cache_dir: Path | None = None,
    ):
        """Initialize variant workflow configuration.

        Args:
            variant_ids: List of variant identifiers (rsID, coordinate, HGVS)
            vcf_file: Path to VCF file with variants
            variant_mode: How to handle variants (avoid/target/both)
            min_af: Minimum allele frequency threshold
            clinvar_filter_levels: Allowed ClinVar significance levels
            assembly: Reference genome assembly (only GRCh38 supported)
            cache_dir: Cache directory for variant data
        """
        self.variant_ids = variant_ids or []
        self.vcf_file = vcf_file
        self.variant_mode = variant_mode
        self.min_af = min_af
        self.clinvar_filter_levels = clinvar_filter_levels or [
            ClinVarSignificance.PATHOGENIC,
            ClinVarSignificance.LIKELY_PATHOGENIC,
        ]
        self.assembly = assembly
        self.cache_dir = cache_dir

    @property
    def has_variants(self) -> bool:
        """Check if any variants are configured."""
        return len(self.variant_ids) > 0 or self.vcf_file is not None


async def resolve_workflow_variants(
    config: VariantWorkflowConfig,
    gene_name: str,
    output_dir: Path,
) -> list[VariantRecord]:
    """Resolve all variants for the workflow.

    This is Step 1.5 in the workflow: after gene selection, before siRNA design.

    Args:
        config: Variant workflow configuration
        gene_name: Gene name for logging/reporting
        output_dir: Output directory for reports

    Returns:
        List of resolved VariantRecords passing filters
    """
    if not config.has_variants:
        logger.info("No variants configured for workflow")
        return []

    logger.info(f"Resolving variants for {gene_name}")

    # Initialize resolver with variant mode for proper AF filtering
    resolver = VariantResolver(
        min_af=config.min_af,
        clinvar_filters=config.clinvar_filter_levels,
        assembly=config.assembly,
        cache_dir=config.cache_dir,
        variant_mode=config.variant_mode.value,  # Pass mode for population AF handling
    )

    resolved_variants: list[VariantRecord] = []

    # Resolve individual variant IDs
    if config.variant_ids:
        logger.info(f"Resolving {len(config.variant_ids)} variant identifiers")
        for variant_id in config.variant_ids:
            try:
                query = resolver.parse_identifier(variant_id)
                variant = await resolver.resolve_variant(query)
                if variant:
                    resolved_variants.append(variant)
                    logger.info(f"Resolved {variant_id} -> {variant.to_vcf_style()}")
                else:
                    logger.warning(f"Could not resolve or filtered out: {variant_id}")
            except Exception as e:
                logger.error(f"Error resolving {variant_id}: {e}")

    # Read variants from VCF file
    if config.vcf_file:
        logger.info(f"Reading variants from VCF: {config.vcf_file}")
        try:
            vcf_variants = resolver.read_vcf(config.vcf_file)
            resolved_variants.extend(vcf_variants)
            logger.info(f"Loaded {len(vcf_variants)} variants from VCF")
        except Exception as e:
            logger.error(f"Error reading VCF file: {e}")

    # Remove duplicates based on chr:pos:ref:alt
    unique_variants = _deduplicate_variants(resolved_variants)
    logger.info(f"Total unique variants after deduplication: {len(unique_variants)}")

    # Save variant resolution report
    if unique_variants:
        report_path = output_dir / "logs" / "resolved_variants.json"
        report_path.parent.mkdir(parents=True, exist_ok=True)
        _save_variant_report(unique_variants, report_path, gene_name, config)

    return unique_variants


def _deduplicate_variants(variants: list[VariantRecord]) -> list[VariantRecord]:
    """Remove duplicate variants based on chr:pos:ref:alt.

    When duplicates are found, keep the one with the highest priority source.

    Args:
        variants: List of variants to deduplicate

    Returns:
        Deduplicated list of variants
    """
    seen: dict[str, VariantRecord] = {}

    for variant in variants:
        key = variant.to_vcf_style()

        if key not in seen:
            seen[key] = variant
        else:
            # Keep variant with higher priority source
            existing = seen[key]
            existing_priority = existing.get_primary_source()
            new_priority = variant.get_primary_source()

            # Priority order: CLINVAR > ENSEMBL > DBSNP > LOCAL_VCF
            priority_order = ["clinvar", "ensembl", "dbsnp", "local-vcf"]
            if existing_priority and new_priority:
                existing_idx = (
                    priority_order.index(existing_priority.value)
                    if existing_priority.value in priority_order
                    else 999
                )
                new_idx = priority_order.index(new_priority.value) if new_priority.value in priority_order else 999

                if new_idx < existing_idx:
                    seen[key] = variant
                    logger.debug(f"Replaced {key} with higher priority source: {new_priority.value}")

    return list(seen.values())


def _save_variant_report(
    variants: list[VariantRecord],
    output_path: Path,
    gene_name: str,
    config: VariantWorkflowConfig,
) -> None:
    """Save variant resolution report to JSON.

    Args:
        variants: List of resolved variants
        output_path: Path to save report
        gene_name: Gene name
        config: Variant workflow configuration
    """
    report = {
        "gene": gene_name,
        "variant_mode": config.variant_mode.value,
        "filters": {
            "min_af": config.min_af,
            "clinvar_significance": [sig.value for sig in config.clinvar_filter_levels],
            "assembly": config.assembly,
        },
        "summary": {
            "total_variants": len(variants),
            "sources": _count_by_source(variants),
            "chromosomes": _count_by_chromosome(variants),
        },
        "variants": [
            {
                "id": v.id,
                "chr": v.chr,
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "sources": [s.value for s in v.sources],
                "primary_source": v.get_primary_source().value if v.get_primary_source() else None,
                "clinvar_significance": v.clinvar_significance.value if v.clinvar_significance else None,
                "af": v.af,
                "vcf_style": v.to_vcf_style(),
            }
            for v in variants
        ],
    }

    with output_path.open("w") as f:
        json.dump(report, f, indent=2)

    logger.info(f"Variant resolution report saved to {output_path}")


def _count_by_source(variants: list[VariantRecord]) -> dict[str, int]:
    """Count variants by primary source."""
    counts: dict[str, int] = {}
    for variant in variants:
        source = variant.get_primary_source()
        if source:
            counts[source.value] = counts.get(source.value, 0) + 1
    return counts


def _count_by_chromosome(variants: list[VariantRecord]) -> dict[str, int]:
    """Count variants by chromosome."""
    counts: dict[str, int] = {}
    for variant in variants:
        counts[variant.chr] = counts.get(variant.chr, 0) + 1
    return counts


def parse_clinvar_filter_string(filter_string: str) -> list[ClinVarSignificance]:
    """Parse comma-separated ClinVar filter string to enum list.

    Args:
        filter_string: Comma-separated string of significance levels

    Returns:
        List of ClinVarSignificance enums

    Raises:
        ValueError: If any significance level is invalid
    """
    levels = [s.strip() for s in filter_string.split(",")]
    result = []

    for level in levels:
        try:
            # Try exact match first
            sig = ClinVarSignificance(level)
            result.append(sig)
        except ValueError:
            # Try case-insensitive match
            for sig_enum in ClinVarSignificance:
                if sig_enum.value.lower() == level.lower():
                    result.append(sig_enum)
                    break
            else:
                raise ValueError(
                    f"Invalid ClinVar significance level: {level}. "
                    f"Valid options: {[s.value for s in ClinVarSignificance]}"
                )

    return result
