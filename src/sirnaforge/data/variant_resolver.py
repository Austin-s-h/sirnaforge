"""Variant resolution from multiple databases (ClinVar, Ensembl, dbSNP) with caching."""

import asyncio
import re
from pathlib import Path
from typing import Any
from urllib.parse import quote, urlencode

import aiohttp
import pysam

from sirnaforge.data.variant_cache import VariantParquetCache
from sirnaforge.data.variant_resolution_helpers import (
    GRCH38_ACCESSIONS,
    VariantResolutionAmbiguity,
    exact_query_match,
    normalize_chromosome,
    unique_variant,
)
from sirnaforge.models.variant import (
    ClinVarSignificance,
    VariantQuery,
    VariantQueryType,
    VariantRecord,
    VariantSource,
)
from sirnaforge.models.variant import VariantMode as VM
from sirnaforge.utils.cache_utils import resolve_cache_subdir, stable_cache_key
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)

# Ensembl's gnomADe population codes are already the standard ones; REMAINING is the one that is not.
_ENSEMBL_POPULATION_CODES = {"REMAINING": "OTH"}


class VariantResolver:
    """Resolve variant identifiers to VariantRecord using multiple databases with priority ordering.

    Priority order: ClinVar > Ensembl > dbSNP
    Supports caching and local VCF files.
    Only supports GRCh38 assembly.
    """

    def __init__(
        self,
        min_af: float = 0.01,
        clinvar_filters: list[ClinVarSignificance] | None = None,
        assembly: str = "GRCh38",
        source_priority: list[VariantSource] | None = None,
        cache_dir: Path | None = None,
        timeout: int = 30,
        variant_mode: str | None = None,
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
        self._hgvs_normalizations: dict[str, VariantRecord] = {}

        # Import VariantMode for filtering

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
            normalized_chr = normalize_chromosome(chrom)
            if normalized_chr is None:
                raise ValueError(f"Unsupported human primary chromosome: {chrom}")
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

    async def resolve_variant(self, query: VariantQuery) -> VariantRecord | None:  # noqa: PLR0911, PLR0912
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
        if cached and exact_query_match(cached, query) and self._passes_filters(cached):
            logger.info(f"Cache hit for variant: {query.raw_input}")
            return cached

        try:
            normalized_query = await self._coordinate_query(query)
        except VariantResolutionAmbiguity as exc:
            logger.error("Variant query is ambiguous: %s", exc)
            return None
        if normalized_query is None:
            return None
        query = normalized_query

        # Try sources in priority order
        variant: VariantRecord | None = None
        for source in self.source_priority:
            try:
                if source == VariantSource.CLINVAR:
                    variant = await self._query_clinvar(query)
                elif source == VariantSource.ENSEMBL:
                    variant = await self._query_ensembl(query)
                elif source == VariantSource.DBSNP:
                    variant = await self._query_dbsnp(query)

                if variant and not exact_query_match(variant, query):
                    logger.warning("Ignoring %s response that does not match %s", source.value, query.raw_input)
                    variant = None
                if variant:
                    logger.info(f"Found variant in {source.value}: {query.raw_input}")
                    break
            except VariantResolutionAmbiguity as exc:
                logger.error("Variant query is ambiguous: %s", exc)
                return None
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

        if query.hgvs:
            variant.provenance["hgvs_normalization"] = {
                "input": query.hgvs,
                "source": "Ensembl VEP HGVS API",
                "resolved": variant.to_vcf_style(),
                "assembly": self.assembly,
            }

        # Cache the result
        self._put_to_cache(cache_key, variant)
        return variant

    async def _request_json(self, url: str, params: dict[str, str] | None = None) -> Any:
        """Read one bounded public API response; adapters validate biological identity."""
        async with (
            aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(total=self.timeout), trust_env=True) as session,
            session.get(url, params=params, headers={"Accept": "application/json"}) as response,
        ):
            response.raise_for_status()
            return await response.json()

    async def _coordinate_query(self, query: VariantQuery) -> VariantQuery | None:
        """Resolve HGVS once to an unambiguous GRCh38 allele before source selection."""
        if query.assembly != self.assembly:
            raise ValueError(f"Variant query assembly must be {self.assembly}")
        if query.query_type != VariantQueryType.HGVS:
            return query
        variant = await self._resolve_hgvs(query)
        if variant is None:
            return None
        return query.model_copy(
            update={
                "query_type": VariantQueryType.COORDINATE,
                "chr": variant.chr,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
            }
        )

    async def _resolve_hgvs(self, query: VariantQuery) -> VariantRecord | None:  # noqa: PLR0912
        """Use Ensembl's documented HGVS endpoint with ambiguity exposed, never hidden.

        https://rest.ensembl.org/documentation/info/vep_hgvs_get
        ``ambiguous_hgvs=1`` requests all protein-HGVS mappings so a protein
        description cannot silently select the server's most likely transcript.
        """
        if not query.hgvs:
            return None
        if query.hgvs in self._hgvs_normalizations:
            return self._hgvs_normalizations[query.hgvs].model_copy(deep=True)
        url = f"https://rest.ensembl.org/vep/human/hgvs/{quote(query.hgvs, safe='')}"
        try:
            data = await self._request_json(url, {"vcf_string": "1", "ambiguous_hgvs": "1"})
            records = []
            for item in data if isinstance(data, list) else []:
                if item.get("assembly_name") != self.assembly or item.get("input") != query.hgvs:
                    continue
                chrom = normalize_chromosome(str(item.get("seq_region_name", "")))
                if chrom is None:
                    continue
                vcf = item.get("vcf_string")
                if isinstance(vcf, str) and len(vcf.split("-")) == 4:
                    vcf_chr, vcf_pos, ref, alt = vcf.split("-")
                    if normalize_chromosome(vcf_chr) != chrom:
                        continue
                    pos = int(vcf_pos)
                else:
                    alleles = str(item.get("allele_string", "")).split("/")
                    if len(alleles) != 2 or item.get("strand") != 1:
                        continue
                    ref, alt = alleles
                    pos = int(item["start"])
                    if int(item["end"]) - pos + 1 != len(ref):
                        continue
                if not all(re.fullmatch(r"[ACGT]+", allele) for allele in (ref, alt)) or ref == alt:
                    continue
                records.append(
                    VariantRecord(
                        id=query.hgvs,
                        chr=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        sources=[VariantSource.ENSEMBL],
                        annotations={"vep_hgvs_data": item},
                        provenance={"source": "Ensembl VEP HGVS API", "query": query.hgvs, "endpoint": url},
                    )
                )
            result = unique_variant(records, query)
            if result is not None:
                self._hgvs_normalizations[query.hgvs] = result.model_copy(deep=True)
            else:
                logger.warning("HGVS %s has no verified GRCh38 nucleotide allele", query.hgvs)
            return result
        except VariantResolutionAmbiguity:
            raise
        except Exception as exc:
            logger.warning("HGVS normalization failed for %s: %s", query.hgvs, exc)
            return None

    async def _query_clinvar(self, query: VariantQuery) -> VariantRecord | None:  # noqa: PLR0912
        """Match every bounded ClinVar hit against explicit GRCh38 locus and alleles.

        Alleles come from a GRCh38 SequenceLocation or its matching canonical
        chromosome SPDI, never from the input's guessed reference/alternate.
        https://www.ncbi.nlm.nih.gov/clinvar/docs/programmatic_access/
        """
        normalized_query = await self._coordinate_query(query)
        if normalized_query is None:
            return None
        query = normalized_query
        if query.query_type == VariantQueryType.RSID and query.rsid:
            search_term = query.rsid.lower()
        elif query.query_type == VariantQueryType.COORDINATE and query.chr and query.pos:
            search_term = f"{query.chr.removeprefix('chr')}[chr] AND {query.pos}[chrpos38]"
        else:
            return None
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        try:
            search = await self._request_json(
                f"{base_url}/esearch.fcgi",
                {
                    "db": "clinvar",
                    "term": search_term,
                    "retmode": "json",
                    "retmax": "100",
                },
            )
            search_result = search.get("esearchresult", {})
            ids = search_result.get("idlist", [])
            if int(search_result.get("count", len(ids))) > len(ids):
                raise VariantResolutionAmbiguity("ClinVar search exceeded its complete 100-hit query budget")
            if not ids:
                return None
            summary = await self._request_json(
                f"{base_url}/esummary.fcgi",
                {
                    "db": "clinvar",
                    "id": ",".join(str(uid) for uid in ids),
                    "retmode": "json",
                },
            )
            if any(str(uid) not in summary.get("result", {}) for uid in ids):
                raise ValueError("ClinVar summary omitted one or more requested records")
            records = []
            for uid in ids:
                item = summary.get("result", {}).get(str(uid), {})
                significance = item.get("germline_classification", {}).get("description")
                try:
                    sig_enum = ClinVarSignificance(significance) if significance else None
                except ValueError:
                    sig_enum = ClinVarSignificance.OTHER
                for variation in item.get("variation_set") or []:
                    if query.rsid and not any(
                        str(xref.get("db_source", "")).lower() == "dbsnp"
                        and str(xref.get("db_id", "")).removeprefix("rs") == query.rsid.lower().removeprefix("rs")
                        for xref in variation.get("variation_xrefs", [])
                    ):
                        continue
                    for loc in variation.get("variation_loc") or []:
                        if loc.get("assembly_name") != self.assembly:
                            continue
                        chrom = normalize_chromosome(str(loc.get("chr", "")))
                        if chrom is None:
                            continue
                        pos = int(loc.get("start", 0))
                        ref, alt = loc.get("ref", ""), loc.get("alt", "")
                        if not ref or not alt:
                            spdi = str(variation.get("canonical_spdi", "")).split(":")
                            if len(spdi) != 4 or GRCH38_ACCESSIONS.get(spdi[0]) != chrom:
                                continue
                            accession, offset, ref, alt = spdi
                            if int(offset) + 1 != pos:
                                continue
                        if (
                            not all(
                                isinstance(allele, str) and re.fullmatch(r"[ACGT]+", allele) for allele in (ref, alt)
                            )
                            or ref == alt
                        ):
                            continue
                        if int(loc.get("stop") or pos + len(ref) - 1) != pos + len(ref) - 1:
                            continue
                        records.append(
                            VariantRecord(
                                id=query.rsid or f"clinvar_{uid}",
                                chr=chrom,
                                pos=pos,
                                ref=ref,
                                alt=alt,
                                sources=[VariantSource.CLINVAR],
                                clinvar_significance=sig_enum,
                                annotations={"clinvar_data": item},
                                provenance={
                                    "clinvar_id": str(uid),
                                    "source": "ClinVar E-utilities",
                                    "query": search_term,
                                },
                            )
                        )
            result = unique_variant(records, query)
            if result is not None:
                matching = [record for record in records if record.to_vcf_style() == result.to_vcf_style()]
                result.provenance["clinvar_ids"] = sorted({record.provenance["clinvar_id"] for record in matching})
                if len({record.clinvar_significance for record in matching}) > 1:
                    result.clinvar_significance = ClinVarSignificance.CONFLICTING
                    result.provenance["classification_conflict"] = [record.clinvar_significance for record in matching]
            return result
        except VariantResolutionAmbiguity:
            raise
        except Exception as exc:
            logger.warning("ClinVar query failed: %s", exc)
            return None

    @staticmethod
    def _ensembl_records(data: dict[str, Any], query: VariantQuery) -> list[VariantRecord]:
        """Parse every primary GRCh38 mapping and keep frequencies allele-specific."""
        records = []
        for mapping in data.get("mappings", []):
            if mapping.get("assembly_name") != "GRCh38":
                continue
            chrom = normalize_chromosome(str(mapping.get("seq_region_name", "")))
            alleles = str(mapping.get("allele_string", "")).split("/")
            if chrom is None or len(alleles) < 2 or mapping.get("strand") not in {-1, 1}:
                continue
            if not all(re.fullmatch(r"[ACGT]+", allele) for allele in alleles):
                continue
            start = int(mapping["start"])
            if int(mapping["end"]) - start + 1 != len(alleles[0]):
                continue
            for source_alt in alleles[1:]:
                ref, alt = alleles[0], source_alt
                if mapping["strand"] == -1:
                    ref, alt = (allele.translate(str.maketrans("ACGT", "TGCA"))[::-1] for allele in (ref, alt))
                if ref == alt:
                    continue
                frequencies = {}
                population_afs = {}
                for pop in data.get("populations") or []:
                    if not pop or pop.get("allele") != source_alt or pop.get("frequency") is None:
                        continue
                    freq = float(pop["frequency"])
                    if not 0 <= freq <= 1:
                        continue
                    name = str(pop.get("population", ""))
                    frequencies[name.lower()] = freq
                    if name.lower().startswith("gnomade:"):
                        code = name.rsplit(":", maxsplit=1)[-1].upper()
                        population_afs[_ENSEMBL_POPULATION_CODES.get(code, code)] = freq
                globals_by_priority = sorted(
                    (name for name in frequencies if name.endswith(":all") or name in {"global", "all"}),
                    key=lambda name: (0 if name == "gnomade:all" else 1, name),
                )
                global_source = next(iter(globals_by_priority), None)
                records.append(
                    VariantRecord(
                        id=data.get("name"),
                        chr=chrom,
                        pos=start,
                        ref=ref,
                        alt=alt,
                        sources=[VariantSource.ENSEMBL],
                        af=frequencies.get(global_source) if global_source else None,
                        population_afs=population_afs,
                        annotations={"ensembl_data": data},
                        provenance={"source": "Ensembl Variation API", "frequency_source": global_source},
                    )
                )
        return [record for record in records if exact_query_match(record, query)]

    async def _query_ensembl(self, query: VariantQuery) -> VariantRecord | None:  # noqa: PLR0911, PLR0912
        """Resolve all mappings and region hits without choosing the first alternate.

        https://rest.ensembl.org/documentation/info/variation_id
        """
        normalized_query = await self._coordinate_query(query)
        if normalized_query is None:
            return None
        query = normalized_query
        base_url = "https://rest.ensembl.org"
        try:
            if query.query_type == VariantQueryType.RSID and query.rsid:
                base_endpoint = f"{base_url}/variation/human/{query.rsid.lower()}"
                params = {"pops": "1"}
                data = await self._request_json(f"{base_endpoint}?{urlencode(params)}")
                aliases = {
                    str(data.get("name", "")).lower(),
                    *(str(alias).lower() for alias in data.get("synonyms") or []),
                }
                if query.rsid.lower() not in aliases:
                    return None
                return unique_variant(self._ensembl_records(data, query), query)
            if query.query_type != VariantQueryType.COORDINATE or not query.chr or not query.pos:
                return None
            chrom = query.chr.removeprefix("chr")
            data = await self._request_json(
                f"{base_url}/overlap/region/human/{chrom}:{query.pos}-{query.pos}", {"feature": "variation"}
            )
            records = []
            for item in data if isinstance(data, list) else []:
                alleles = item.get("alleles", [])
                if isinstance(alleles, str):
                    alleles = alleles.split("/")
                variation = {"name": item.get("id"), "mappings": [{**item, "allele_string": "/".join(alleles)}]}
                records.extend(self._ensembl_records(variation, query))
            result = unique_variant(records, query)
            if result is None:
                normalized = self._hgvs_normalizations.get(query.hgvs or "")
                return normalized.model_copy(deep=True) if normalized and exact_query_match(normalized, query) else None
            # Region-overlap records have no population frequencies; enrich only
            # the already-proven exact allele, never replace it with a neighbour.
            if result.id:
                try:
                    detail = await self._request_json(f"{base_url}/variation/human/{result.id}", {"pops": "1"})
                    enriched = unique_variant(self._ensembl_records(detail, query), query)
                    if enriched is not None:
                        return enriched
                except VariantResolutionAmbiguity:
                    raise
                except Exception as exc:
                    logger.warning("Ensembl exact-allele population enrichment unavailable: %s", exc)
            return result
        except VariantResolutionAmbiguity:
            raise
        except Exception as exc:
            logger.warning("Ensembl query failed: %s", exc)
            return None

    async def _query_dbsnp(self, query: VariantQuery) -> VariantRecord | None:  # noqa: PLR0912
        """Read NCBI RefSNP chromosome placements and match exact GRCh38 alleles.

        https://github.com/ncbi/dbsnp/tree/master/specs
        SPDI offsets are zero-based boundaries; VariantRecord positions are one-based.
        Empty SPDI alleles require reference anchoring and are not guessed here.
        """
        normalized_query = await self._coordinate_query(query)
        if normalized_query is None:
            return None
        query = normalized_query
        base_url = "https://api.ncbi.nlm.nih.gov/variation/v0"
        try:
            if query.query_type == VariantQueryType.RSID and query.rsid:
                ids = [query.rsid.lower().removeprefix("rs")]
            elif query.query_type == VariantQueryType.COORDINATE and query.chr and query.pos:
                accession = next((seq_id for seq_id, chrom in GRCH38_ACCESSIONS.items() if chrom == query.chr), None)
                if accession is None:
                    return None
                spdi = f"{accession}:{query.pos - 1}:{query.ref}:{query.alt}"
                lookup = await self._request_json(f"{base_url}/spdi/{quote(spdi, safe='')}/rsids")
                ids = lookup.get("data", {}).get("rsids", [])
            else:
                return None
            if len(ids) > 100:
                raise VariantResolutionAmbiguity("dbSNP exact-allele lookup exceeded its 100-record query budget")
            records = []
            for rsid in ids:
                data = await self._request_json(f"{base_url}/refsnp/{rsid}")
                if str(data.get("refsnp_id")) != str(rsid):
                    continue
                snapshot = data.get("primary_snapshot_data", {})
                annotations = snapshot.get("allele_annotations", [])
                for placement in snapshot.get("placements_with_allele", []):
                    accession = placement.get("seq_id")
                    chrom = GRCH38_ACCESSIONS.get(accession)
                    traits = placement.get("placement_annot", {}).get("seq_id_traits_by_assembly", [])
                    if chrom is None or not any(
                        (
                            trait.get("assembly_name") == self.assembly
                            or re.fullmatch(r"GRCh38\.p\d+", trait.get("assembly_name", ""))
                        )
                        and trait.get("is_top_level") is True
                        and trait.get("is_chromosome") is True
                        and not trait.get("is_alt")
                        and not trait.get("is_patch")
                        for trait in traits
                    ):
                        continue
                    for index, allele in enumerate(placement.get("alleles", [])):
                        spdi = allele.get("allele", {}).get("spdi", {})
                        ref, alt = spdi.get("deleted_sequence", ""), spdi.get("inserted_sequence", "")
                        if (
                            spdi.get("seq_id") != accession
                            or ref == alt
                            or not all(re.fullmatch(r"[ACGT]+", value) for value in (ref, alt))
                        ):
                            continue
                        frequencies = annotations[index].get("frequency", []) if index < len(annotations) else []
                        usable_frequencies = []
                        for frequency in frequencies:
                            total = frequency.get("total_count", 0)
                            count = frequency.get("allele_count", -1)
                            if isinstance(total, int) and isinstance(count, int) and total > 0 and 0 <= count <= total:
                                usable_frequencies.append({**frequency, "af": count / total})
                        usable_frequencies.sort(
                            key=lambda frequency: (
                                0 if "gnomad" in frequency.get("study_name", "").lower() else 1,
                                frequency.get("study_name", ""),
                            )
                        )
                        selected = next(iter(usable_frequencies), None)
                        records.append(
                            VariantRecord(
                                id=f"rs{rsid}",
                                chr=chrom,
                                pos=int(spdi["position"]) + 1,
                                ref=ref,
                                alt=alt,
                                sources=[VariantSource.DBSNP],
                                af=selected["af"] if selected else None,
                                annotations={"dbsnp_data": data, "dbsnp_frequencies": usable_frequencies},
                                provenance={
                                    "source": "NCBI RefSNP API",
                                    "sequence_accession": accession,
                                    "frequency_source": selected.get("study_name") if selected else None,
                                },
                            )
                        )
            return unique_variant(records, query)
        except VariantResolutionAmbiguity:
            raise
        except Exception as exc:
            logger.warning("dbSNP query failed: %s", exc)
            return None

    def _passes_filters(self, variant: VariantRecord) -> bool:
        """Check if variant passes AF and ClinVar filters.

        Args:
            variant: Variant record to check

        Returns:
            True if variant passes all filters
        """
        # Get effective AF based on variant mode
        effective_af = variant.get_effective_af_for_mode(self.variant_mode) if self.variant_mode else variant.af

        # AF filter - warn if no AF data but don't fail
        if effective_af is None:
            logger.warning(f"Variant {variant.to_vcf_style()} has no allele frequency data, skipping AF filter")
        elif effective_af < self.min_af:
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
        allowed_significances = set(self.clinvar_filters)
        if {ClinVarSignificance.PATHOGENIC, ClinVarSignificance.LIKELY_PATHOGENIC} <= allowed_significances:
            allowed_significances.add(ClinVarSignificance.PATHOGENIC_LIKELY_PATHOGENIC)
        if VariantSource.CLINVAR in variant.sources and variant.clinvar_significance not in allowed_significances:
            logger.debug(f"Variant filtered: ClinVar significance {variant.clinvar_significance} not in allowed list")
            return False

        return True

    def _get_cache_key(self, query: VariantQuery) -> str:
        """Generate cache key for variant query.

        The key must cover every setting that decides whether a variant is stored
        at all, because ``resolve_variant`` returns a cache hit without re-running
        ``_passes_filters``. ``variant_mode`` is one of those settings: avoid mode
        filters on the maximum population AF while target/both mode filters on the
        global AF, so at the same ``min_af`` the two modes admit different variants
        and must not share entries.

        Args:
            query: Parsed variant query

        Returns:
            Cache key string
        """
        cache_data: dict[str, Any] = {
            "resolver_schema": 2,
            "query": query.raw_input,
            "assembly": self.assembly,
            "min_af": self.min_af,
            "clinvar_filters": [f.value for f in self.clinvar_filters],
            "variant_mode": self.variant_mode.value if self.variant_mode else None,
            "source_priority": [source.value for source in self.source_priority],
        }
        return stable_cache_key(cache_data)

    def _get_from_cache(self, cache_key: str) -> VariantRecord | None:
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
                    af: float | None = None
                    if "AF" in record.info:
                        af_values = record.info["AF"]
                        af = af_values[0] if isinstance(af_values, list | tuple) else af_values
                    elif "AC" in record.info and "AN" in record.info:
                        # Calculate AF from AC/AN
                        ac = record.info["AC"]
                        an = record.info["AN"]
                        ac_val: float = ac[0] if isinstance(ac, list | tuple) else ac
                        an_val: float = an[0] if isinstance(an, list | tuple) else an
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
                        ref=str(record.ref),
                        alt=str(alt),
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
    clinvar_filters: list[ClinVarSignificance] | None = None,
    cache_dir: Path | None = None,
    variant_mode: str | None = None,
) -> VariantRecord | None:
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
