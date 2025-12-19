"""Variant resolution from multiple databases (ClinVar, Ensembl, dbSNP) with caching."""

import asyncio
import re
from pathlib import Path
from typing import Optional

import aiohttp
import pysam

from sirnaforge.data.variant_cache import VariantParquetCache
from sirnaforge.models.variant import (
    ClinVarSignificance,
    VariantQuery,
    VariantQueryType,
    VariantRecord,
    VariantSource,
)
from sirnaforge.utils.cache_utils import resolve_cache_subdir, stable_cache_key
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class VariantResolver:
    """Resolve variant identifiers to VariantRecord using multiple databases with priority ordering.

    Priority order: ClinVar > Ensembl > dbSNP
    Supports caching and local VCF files.
    Only supports GRCh38 assembly.
    """

    def __init__(
        self,
        min_af: float = 0.01,
        clinvar_filters: Optional[list[ClinVarSignificance]] = None,
        assembly: str = "GRCh38",
        source_priority: Optional[list[VariantSource]] = None,
        cache_dir: Optional[Path] = None,
        timeout: int = 30,
        variant_mode: Optional[str] = None,
    ):
        """Initialize variant resolver.

        Args:
            min_af: Minimum allele frequency threshold (default: 0.01)
            clinvar_filters: Allowed ClinVar significance levels (default: Pathogenic, Likely pathogenic)
            assembly: Reference genome assembly (only GRCh38 supported)
            source_priority: Source priority list (default: ClinVar > Ensembl > dbSNP)
            cache_dir: Cache directory for variant data
            timeout: HTTP request timeout in seconds
            variant_mode: Variant mode for AF filtering ('avoid', 'target', 'both').
                         In 'avoid' mode, uses max population AF if available to avoid SNPs
                         prevalent in any geographic group (e.g., >10% in one population).
        """
        if assembly != "GRCh38":
            raise ValueError(f"Only GRCh38 assembly is supported, got: {assembly}")

        self.min_af = min_af
        self.clinvar_filters = clinvar_filters or [
            ClinVarSignificance.PATHOGENIC,
            ClinVarSignificance.LIKELY_PATHOGENIC,
        ]
        self.assembly = assembly
        self.source_priority = source_priority or [
            VariantSource.CLINVAR,
            VariantSource.ENSEMBL,
            VariantSource.DBSNP,
        ]
        self.timeout = timeout
        
        # Import VariantMode for filtering
        from sirnaforge.models.variant import VariantMode as VM
        self.variant_mode = VM(variant_mode) if variant_mode else None

        # Initialize Parquet-based cache
        self.cache_dir = cache_dir or resolve_cache_subdir("variants")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Use Parquet cache
        self.cache = VariantParquetCache(self.cache_dir)
        logger.info(f"Variant cache at {self.cache_dir}")

    def parse_identifier(self, input_str: str) -> VariantQuery:
        """Parse variant identifier string into VariantQuery.

        Accepts:
        - rsID: rs12345
        - VCF-style coordinate: chr17:7577121:G:A or 17:7577121:G:A
        - HGVS: NM_000546.6:c.215C>G (basic support)

        Args:
            input_str: User-provided variant identifier

        Returns:
            VariantQuery with parsed components

        Raises:
            ValueError: If input format is not recognized
        """
        input_str = input_str.strip()

        # Check for rsID format
        rsid_match = re.match(r"^(rs\d+)$", input_str, re.IGNORECASE)
        if rsid_match:
            return VariantQuery(
                raw_input=input_str,
                query_type=VariantQueryType.RSID,
                rsid=rsid_match.group(1),
                assembly=self.assembly,
            )

        # Check for coordinate format: chr:pos:ref:alt or chrom:pos:ref:alt
        coord_match = re.match(r"^(chr)?(\d+|X|Y|MT?):(\d+):([ACGT]+):([ACGT]+)$", input_str, re.IGNORECASE)
        if coord_match:
            chr_prefix, chrom, pos, ref, alt = coord_match.groups()
            # Normalize chromosome to 'chrN' format
            normalized_chr = f"chr{chrom}" if not chr_prefix else f"chr{chrom}"
            return VariantQuery(
                raw_input=input_str,
                query_type=VariantQueryType.COORDINATE,
                chr=normalized_chr,
                pos=int(pos),
                ref=ref.upper(),
                alt=alt.upper(),
                assembly=self.assembly,
            )

        # Check for basic HGVS format (simplified detection)
        hgvs_match = re.match(r"^(NM_|NP_|ENST|ENSP|NC_)[\w.]+:[cgp]\.\S+$", input_str)
        if hgvs_match:
            return VariantQuery(
                raw_input=input_str,
                query_type=VariantQueryType.HGVS,
                hgvs=input_str,
                assembly=self.assembly,
            )

        raise ValueError(
            f"Unrecognized variant identifier format: {input_str}. "
            "Supported formats: rsID (rs12345), coordinate (chr17:7577121:G:A), HGVS (NM_000546.6:c.215C>G)"
        )

    async def resolve_variant(self, query: VariantQuery) -> Optional[VariantRecord]:
        """Resolve a variant query to a VariantRecord.

        Tries sources in priority order: ClinVar -> Ensembl -> dbSNP
        Applies AF and ClinVar filters.
        Uses cache when available.

        Args:
            query: Parsed variant query

        Returns:
            VariantRecord if found and passes filters, None otherwise
        """
        # Check cache first
        cache_key = self._get_cache_key(query)
        cached = self._get_from_cache(cache_key)
        if cached:
            logger.info(f"Cache hit for variant: {query.raw_input}")
            return cached

        # Try sources in priority order
        variant: Optional[VariantRecord] = None
        for source in self.source_priority:
            try:
                if source == VariantSource.CLINVAR:
                    variant = await self._query_clinvar(query)
                elif source == VariantSource.ENSEMBL:
                    variant = await self._query_ensembl(query)
                elif source == VariantSource.DBSNP:
                    variant = await self._query_dbsnp(query)

                if variant:
                    logger.info(f"Found variant in {source.value}: {query.raw_input}")
                    break
            except Exception as e:
                logger.warning(f"Error querying {source.value} for {query.raw_input}: {e}")
                continue

        if not variant:
            logger.warning(f"Variant not found in any source: {query.raw_input}")
            return None

        # Apply filters
        if not self._passes_filters(variant):
            logger.info(f"Variant filtered out: {query.raw_input}")
            return None

        # Cache the result
        self._put_to_cache(cache_key, variant)
        return variant

    async def _query_clinvar(self, query: VariantQuery) -> Optional[VariantRecord]:  # noqa: PLR0911
        """Query ClinVar for variant information.

        Uses ClinVar E-utilities API for GRCh38.

        Args:
            query: Parsed variant query

        Returns:
            VariantRecord if found, None otherwise
        """
        # ClinVar E-utilities base URL
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

        # Build query based on type
        search_term = None
        if query.query_type == VariantQueryType.RSID and query.rsid:
            search_term = query.rsid
        elif query.query_type == VariantQueryType.COORDINATE and query.chr and query.pos:
            # ClinVar coordinate search format
            chrom_num = query.chr.replace("chr", "")
            search_term = f"{chrom_num}[chr] AND {query.pos}[chrpos37]"  # Use chrpos38 for GRCh38
        else:
            return None

        if not search_term:
            return None

        try:
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=self.timeout)) as session:
                # Search ClinVar
                search_url = f"{base_url}/esearch.fcgi"
                params = {"db": "clinvar", "term": search_term, "retmode": "json", "retmax": 1}

                async with session.get(search_url, params=params) as response:
                    if response.status != 200:
                        logger.warning(f"ClinVar search failed: HTTP {response.status}")
                        return None

                    data = await response.json()
                    id_list = data.get("esearchresult", {}).get("idlist", [])

                    if not id_list:
                        return None

                    # Fetch variant details
                    clinvar_id = id_list[0]
                    fetch_url = f"{base_url}/esummary.fcgi"
                    params = {"db": "clinvar", "id": clinvar_id, "retmode": "json"}

                    async with session.get(fetch_url, params=params) as fetch_response:
                        if fetch_response.status != 200:
                            return None

                        summary_data = await fetch_response.json()
                        result = summary_data.get("result", {}).get(clinvar_id, {})

                        # Parse ClinVar response (simplified)
                        # Note: This is a basic implementation. Full ClinVar parsing is complex.
                        clinical_significance = result.get("clinical_significance", {}).get("description", "Unknown")

                        # Map to ClinVarSignificance enum
                        try:
                            sig_enum = ClinVarSignificance(clinical_significance)
                        except ValueError:
                            sig_enum = ClinVarSignificance.OTHER

                        # Extract coordinates if available
                        chr_val = query.chr or "unknown"
                        pos_val = query.pos or 0
                        ref_val = query.ref or ""
                        alt_val = query.alt or ""

                        return VariantRecord(
                            id=query.rsid or f"clinvar_{clinvar_id}",
                            chr=chr_val,
                            pos=pos_val,
                            ref=ref_val,
                            alt=alt_val,
                            assembly=self.assembly,
                            sources=[VariantSource.CLINVAR],
                            clinvar_significance=sig_enum,
                            provenance={"clinvar_id": clinvar_id, "source": "ClinVar E-utilities"},
                        )

        except Exception as e:
            logger.warning(f"ClinVar query failed: {e}")
            return None

    async def _query_ensembl(self, query: VariantQuery) -> Optional[VariantRecord]:  # noqa: PLR0911, PLR0912
        """Query Ensembl Variation API for variant information.

        Args:
            query: Parsed variant query

        Returns:
            VariantRecord if found, None otherwise
        """
        base_url = "https://rest.ensembl.org"

        try:
            async with aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=self.timeout)) as session:
                headers = {"Content-Type": "application/json"}

                # Build API endpoint based on query type
                if query.query_type == VariantQueryType.RSID and query.rsid:
                    url = f"{base_url}/variation/human/{query.rsid}"
                elif query.query_type == VariantQueryType.COORDINATE and query.chr and query.pos:
                    # Ensembl region endpoint
                    chrom = query.chr.replace("chr", "")
                    url = f"{base_url}/overlap/region/human/{chrom}:{query.pos}-{query.pos}?feature=variation"
                else:
                    return None

                async with session.get(url, headers=headers) as response:
                    if response.status != 200:
                        logger.debug(f"Ensembl API returned status {response.status} for {query.raw_input}")
                        return None

                    data = await response.json()

                    # Parse response based on endpoint
                    if query.query_type == VariantQueryType.RSID:
                        # Direct variant lookup response
                        if not data or "error" in data:
                            return None

                        # Extract allele frequency from gnomAD if available
                        populations = data.get("populations", [])
                        af = None
                        population_afs = {}
                        
                        for pop in populations:
                            pop_name = pop.get("population", "")
                            freq = pop.get("frequency")
                            
                            if freq is not None:
                                # Extract global gnomAD AF
                                if "gnomad" in pop_name.lower() and ":" not in pop_name:
                                    if af is None:
                                        af = freq
                                
                                # Extract population-specific AFs from gnomAD
                                # Format examples: "1000GENOMES:phase_3:AFR", "gnomAD:AFR"
                                if ":" in pop_name:
                                    parts = pop_name.split(":")
                                    # Get the last part which is usually the population code
                                    pop_code = parts[-1].upper()
                                    # Standard population codes: AFR, AMR, EAS, EUR, SAS, etc.
                                    if pop_code in ["AFR", "AMR", "EAS", "EUR", "SAS", "FIN", "ASJ", "OTH"]:
                                        # Keep the maximum AF if we see multiple sources for same population
                                        if pop_code not in population_afs or freq > population_afs[pop_code]:
                                            population_afs[pop_code] = freq

                        # Get first allele information
                        mappings = data.get("mappings", [])
                        if not mappings:
                            return None

                        mapping = mappings[0]
                        location = mapping.get("location", "")
                        match = re.match(r"(\d+|X|Y|MT):(\d+)-(\d+)", location)
                        if not match:
                            return None

                        chrom, start, end = match.groups()
                        allele_string = mapping.get("allele_string", "")
                        alleles = allele_string.split("/")
                        if len(alleles) < 2:
                            return None

                        return VariantRecord(
                            id=query.rsid,
                            chr=f"chr{chrom}",
                            pos=int(start),
                            ref=alleles[0],
                            alt=alleles[1],
                            assembly=self.assembly,
                            sources=[VariantSource.ENSEMBL],
                            af=af,
                            population_afs=population_afs,
                            annotations={"ensembl_data": data},
                            provenance={"source": "Ensembl Variation API"},
                        )

                    # Region overlap response - take first variant
                    if not data or not isinstance(data, list) or len(data) == 0:
                        return None

                    first_var = data[0]
                    var_id = first_var.get("id", "unknown")
                    chrom = first_var.get("seq_region_name", query.chr.replace("chr", ""))
                    pos = first_var.get("start", query.pos)
                    alleles = first_var.get("alleles", "").split("/")

                    if len(alleles) < 2:
                        return None

                    return VariantRecord(
                        id=var_id,
                        chr=f"chr{chrom}",
                        pos=pos,
                        ref=alleles[0],
                        alt=alleles[1],
                        assembly=self.assembly,
                        sources=[VariantSource.ENSEMBL],
                        annotations={"ensembl_data": first_var},
                        provenance={"source": "Ensembl Variation API"},
                    )

        except Exception as e:
            logger.warning(f"Ensembl query failed: {e}")
            return None

    async def _query_dbsnp(self, query: VariantQuery) -> Optional[VariantRecord]:
        """Query dbSNP for variant information.

        Note: This is a placeholder. Full dbSNP API integration would require
        additional implementation or use of local dbSNP VCF files.

        Args:
            query: Parsed variant query

        Returns:
            VariantRecord if found, None otherwise
        """
        # TODO: Implement dbSNP API query or VCF lookup
        logger.debug(f"dbSNP query not yet implemented for: {query.raw_input}")
        return None

    def _passes_filters(self, variant: VariantRecord) -> bool:
        """Check if variant passes AF and ClinVar filters.

        Args:
            variant: Variant record to check

        Returns:
            True if variant passes all filters
        """
        # Get effective AF based on variant mode
        if self.variant_mode:
            effective_af = variant.get_effective_af_for_mode(self.variant_mode)
        else:
            effective_af = variant.af

        # AF filter
        if effective_af is not None and effective_af < self.min_af:
            max_pop_af = variant.get_max_population_af()
            if max_pop_af and max_pop_af != effective_af:
                global_af_str = f"{variant.af:.4f}" if variant.af is not None else "N/A"
                logger.debug(
                    f"Variant filtered: effective AF {effective_af:.4f} < {self.min_af} "
                    f"(global AF: {global_af_str}, "
                    f"max population AF: {max_pop_af:.4f})"
                )
            else:
                logger.debug(f"Variant filtered: AF {effective_af:.4f} < {self.min_af}")
            return False

        # ClinVar significance filter
        if (
            VariantSource.CLINVAR in variant.sources
            and variant.clinvar_significance
            and variant.clinvar_significance not in self.clinvar_filters
        ):
            logger.debug(f"Variant filtered: ClinVar significance {variant.clinvar_significance} not in allowed list")
            return False

        return True

    def _get_cache_key(self, query: VariantQuery) -> str:
        """Generate cache key for variant query.

        Args:
            query: Parsed variant query

        Returns:
            Cache key string
        """
        cache_data = {
            "query": query.raw_input,
            "assembly": self.assembly,
            "min_af": self.min_af,
            "clinvar_filters": [f.value for f in self.clinvar_filters],
        }
        return stable_cache_key(cache_data)

    def _get_from_cache(self, cache_key: str) -> Optional[VariantRecord]:
        """Retrieve variant from cache.

        Args:
            cache_key: Cache key

        Returns:
            VariantRecord if found in cache, None otherwise
        """
        return self.cache.get(cache_key)

    def _put_to_cache(self, cache_key: str, variant: VariantRecord) -> None:
        """Store variant in cache.

        Args:
            cache_key: Cache key
            variant: Variant record to cache
        """
        self.cache.put(cache_key, variant)

    def read_vcf(self, vcf_path: Path) -> list[VariantRecord]:
        """Read variants from VCF file (supports bgzip+tabix).

        Args:
            vcf_path: Path to VCF file

        Returns:
            List of VariantRecord objects passing filters

        Raises:
            FileNotFoundError: If VCF file doesn't exist
        """
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")

        variants: list[VariantRecord] = []

        try:
            # Open VCF with pysam
            vcf = pysam.VariantFile(str(vcf_path))

            # Check if indexed
            if not vcf.index:
                logger.warning("VCF file is not indexed; performance may be slow for large files")

            for record in vcf:
                # Process each alternate allele
                for alt in record.alts or []:
                    # Extract AF if available
                    af = None
                    if "AF" in record.info:
                        af_values = record.info["AF"]
                        af = af_values[0] if isinstance(af_values, (list, tuple)) else af_values
                    elif "AC" in record.info and "AN" in record.info:
                        # Calculate AF from AC/AN
                        ac = record.info["AC"]
                        an = record.info["AN"]
                        ac_val = ac[0] if isinstance(ac, (list, tuple)) else ac
                        an_val = an[0] if isinstance(an, (list, tuple)) else an
                        if an_val > 0:
                            af = ac_val / an_val

                    # Normalize chromosome
                    chrom = record.chrom
                    if not chrom.startswith("chr"):
                        chrom = f"chr{chrom}"

                    variant = VariantRecord(
                        id=record.id or None,
                        chr=chrom,
                        pos=record.pos,
                        ref=record.ref,
                        alt=alt,
                        assembly=self.assembly,
                        sources=[VariantSource.LOCAL_VCF],
                        af=af,
                        annotations={"vcf_info": dict(record.info)},
                        provenance={"source": "local VCF", "file": str(vcf_path)},
                    )

                    # Apply filters
                    if self._passes_filters(variant):
                        variants.append(variant)

            vcf.close()

        except Exception as e:
            logger.error(f"Error reading VCF file {vcf_path}: {e}")
            raise

        logger.info(f"Loaded {len(variants)} variants from VCF file (after filtering)")
        return variants


def resolve_variant_sync(
    variant_id: str,
    min_af: float = 0.01,
    clinvar_filters: Optional[list[ClinVarSignificance]] = None,
    cache_dir: Optional[Path] = None,
    variant_mode: Optional[str] = None,
) -> Optional[VariantRecord]:
    """Synchronous wrapper for variant resolution.

    Args:
        variant_id: Variant identifier (rsID, coordinate, or HGVS)
        min_af: Minimum allele frequency
        clinvar_filters: Allowed ClinVar significance levels
        cache_dir: Cache directory
        variant_mode: Variant mode for AF filtering ('avoid', 'target', 'both')

    Returns:
        VariantRecord if found and passes filters, None otherwise
    """
    resolver = VariantResolver(
        min_af=min_af,
        clinvar_filters=clinvar_filters,
        cache_dir=cache_dir,
        variant_mode=variant_mode,
    )
    query = resolver.parse_identifier(variant_id)
    return asyncio.run(resolver.resolve_variant(query))
