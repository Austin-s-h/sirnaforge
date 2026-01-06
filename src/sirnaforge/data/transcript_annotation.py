"""Transcript annotation providers using Ensembl REST and optional VEP enrichment.

This module provides clients for fetching genomic transcript annotations
(exon/CDS structure, coordinates, biotype) separate from sequence retrieval.

**Architecture Overview:**

- **EnsemblTranscriptModelClient**: Primary implementation using Ensembl REST API
- **VepConsequenceClient**: Optional enrichment client (placeholder for future development)

**Caching Strategy:**

Uses in-memory LRU cache with TTL rather than ReferenceManager's persistent file cache.
This design choice is intentional because:

1. **Data Size**: Annotation JSON responses are small (KB) vs. sequence files (GB)
2. **Volatility**: Annotations may update with new releases; TTL provides freshness
3. **Access Pattern**: High frequency, low latency requirements during workflow execution
4. **Scope**: Transient metadata enrichment vs. permanent reference datasets

The cache automatically evicts oldest entries when reaching max_cache_entries,
and entries expire after cache_ttl seconds.

**Relationship to GeneSearcher:**

- GeneSearcher: Discovers transcripts by gene name, fetches cDNA/protein sequences
- This module: Enriches known transcript IDs with genomic structural metadata
- Both can use Ensembl, but query different API endpoints for different purposes
- No redundancy: complementary data types that don't overlap
"""

import asyncio
from collections.abc import Mapping
from time import time
from typing import Any, cast

import aiohttp

from sirnaforge.config.reference_policy import ReferenceChoice
from sirnaforge.data.base import AbstractTranscriptAnnotationClient, DatabaseAccessError
from sirnaforge.models.transcript_annotation import Interval, TranscriptAnnotation, TranscriptAnnotationBundle
from sirnaforge.utils.logging_utils import get_logger

logger = get_logger(__name__)


class EnsemblTranscriptModelClient(AbstractTranscriptAnnotationClient):
    """Ensembl REST-based transcript annotation client.

    Retrieves transcript metadata including genomic coordinates, exon/CDS structure,
    and biotype information using Ensembl's public REST API.

    **API Endpoints Used:**

    1. **Lookup by ID** (`/lookup/id/:id?expand=1`):
       - Fetches detailed annotation for single transcript/gene ID
       - Returns exon coordinates, CDS intervals, biotype
       - Example: /lookup/id/ENST00000269305?expand=1

    2. **Overlap by Region** (`/overlap/region/:species/:region`):
       - Fetches all transcripts overlapping genomic region
       - Useful for region-based queries
       - Example: /overlap/region/human/17:7661779-7687550?feature=transcript

    **Caching Implementation:**

    - Cache key format: "id:{species}:{identifier}:{reference}" or "region:{species}:{region}:{reference}"
    - TTL: Configurable, default 1 hour (3600 seconds)
    - Eviction: LRU when max_cache_entries reached (default 1000)
    - Thread-safe: Single-process use only (workflow orchestration context)

    **Error Handling:**

    - 404: ID not found → added to unresolved list, no exception raised
    - 403/503: Server unavailable → DatabaseAccessError raised
    - Network errors: Wrapped in DatabaseAccessError with context
    - Timeout: Configurable via timeout parameter

    **Example Usage:**

        >>> client = EnsemblTranscriptModelClient()
        >>> reference = ReferenceChoice.explicit("GRCh38", reason="user-specified")
        >>> bundle = await client.fetch_by_ids(
        ...     ids=["ENST00000269305"],
        ...     species="human",
        ...     reference=reference
        ... )
        >>> print(f"Resolved: {bundle.resolved_count}, Unresolved: {bundle.unresolved_count}")
    """

    def __init__(
        self,
        timeout: int = 30,
        base_url: str = "https://rest.ensembl.org",
        cache_ttl: int = 3600,
        max_cache_entries: int = 1000,
    ):
        """Initialize Ensembl transcript annotation client.

        Args:
            timeout: Request timeout in seconds
            base_url: Ensembl REST API base URL
            cache_ttl: Cache time-to-live in seconds (default: 1 hour)
            max_cache_entries: Maximum number of cached entries (default: 1000)
        """
        super().__init__(timeout)
        self.base_url = base_url
        self.cache_ttl = cache_ttl
        self.max_cache_entries = max_cache_entries
        self._cache: dict[str, tuple[Any, float]] = {}

    def _get_cached(self, key: str) -> Any | None:
        """Retrieve item from cache if present and not expired.

        Args:
            key: Cache key (format: "type:species:identifier:reference")

        Returns:
            Cached value if present and not expired, None otherwise

        Side Effects:
            Removes expired entries from cache during lookup
        """
        if key not in self._cache:
            return None

        value, timestamp = self._cache[key]
        if time() - timestamp > self.cache_ttl:
            del self._cache[key]
            return None

        return value

    def _set_cache(self, key: str, value: Any) -> None:
        """Store item in cache with current timestamp.

        Implements simple LRU eviction: when cache reaches max_cache_entries,
        removes 10% of oldest entries by timestamp to make room.

        Args:
            key: Cache key (format: "type:species:identifier:reference")
            value: Value to cache (TranscriptAnnotation, dict, or None for unresolved)

        Side Effects:
            May evict up to 10% of oldest cache entries if at capacity
        """
        # Simple LRU eviction: remove oldest entries when cache is full
        if len(self._cache) >= self.max_cache_entries:
            # Remove 10% of oldest entries
            sorted_items = sorted(self._cache.items(), key=lambda x: x[1][1])
            remove_count = max(1, self.max_cache_entries // 10)
            for old_key, _ in sorted_items[:remove_count]:
                del self._cache[old_key]

        self._cache[key] = (value, time())

    async def fetch_by_ids(
        self, ids: list[str], *, species: str, reference: ReferenceChoice
    ) -> TranscriptAnnotationBundle:
        """Fetch transcript annotations by stable IDs using Ensembl lookup endpoint.

        Args:
            ids: List of transcript or gene IDs
            species: Species name (e.g., 'homo_sapiens', 'human')
            reference: Reference assembly/release choice

        Returns:
            TranscriptAnnotationBundle with resolved annotations
        """
        transcripts: dict[str, TranscriptAnnotation] = {}
        unresolved: list[str] = []

        # Normalize species name for Ensembl (convert 'human' to 'homo_sapiens')
        normalized_species = self._normalize_species(species)

        for identifier in ids:
            # Check cache first
            cache_key = f"id:{normalized_species}:{identifier}:{reference.value or 'default'}"
            cached = self._get_cached(cache_key)
            if cached is not None:
                if isinstance(cached, TranscriptAnnotation):
                    transcripts[identifier] = cached
                else:
                    unresolved.append(identifier)
                continue

            try:
                annotation = await self._fetch_annotation_by_id(identifier, normalized_species)
                if annotation:
                    # Update source metadata
                    annotation.provider = "ensembl_rest"
                    annotation.endpoint = f"{self.base_url}/lookup/id/{identifier}"
                    annotation.reference_choice = reference.value

                    transcripts[identifier] = annotation
                    self._set_cache(cache_key, annotation)
                else:
                    unresolved.append(identifier)
                    self._set_cache(cache_key, None)
            except DatabaseAccessError:
                logger.warning(f"Failed to fetch annotation for {identifier}")
                unresolved.append(identifier)
                self._set_cache(cache_key, None)

        return TranscriptAnnotationBundle(
            transcripts=transcripts,
            unresolved=unresolved,
            reference_choice=reference,
        )

    async def fetch_by_regions(
        self, regions: list[str], *, species: str, reference: ReferenceChoice
    ) -> TranscriptAnnotationBundle:
        """Fetch transcript annotations by genomic regions using Ensembl overlap endpoint.

        Args:
            regions: List of regions in format 'chr:start-end' (e.g., '17:7661779-7687550')
            species: Species name (e.g., 'homo_sapiens', 'human')
            reference: Reference assembly/release choice

        Returns:
            TranscriptAnnotationBundle with all transcripts overlapping regions
        """
        transcripts: dict[str, TranscriptAnnotation] = {}
        unresolved_regions: list[str] = []

        # Normalize species name
        normalized_species = self._normalize_species(species)

        for region in regions:
            # Check cache first
            cache_key = f"region:{normalized_species}:{region}:{reference.value or 'default'}"
            cached = self._get_cached(cache_key)
            if cached is not None:
                if isinstance(cached, dict):
                    transcripts.update(cached)
                continue

            try:
                region_transcripts = await self._fetch_annotations_by_region(region, normalized_species)

                # Update source metadata for all transcripts in region
                if region_transcripts:
                    for annotation in region_transcripts.values():
                        annotation.provider = "ensembl_rest"
                        annotation.endpoint = f"{self.base_url}/overlap/region/{normalized_species}/{region}"
                        annotation.reference_choice = reference.value

                    transcripts.update(region_transcripts)
                    self._set_cache(cache_key, region_transcripts)
                else:
                    self._set_cache(cache_key, {})
            except DatabaseAccessError:
                logger.warning(f"Failed to fetch annotations for region {region}")
                unresolved_regions.append(region)
                self._set_cache(cache_key, {})

        return TranscriptAnnotationBundle(
            transcripts=transcripts,
            unresolved=unresolved_regions,
            reference_choice=reference,
        )

    async def _fetch_annotation_by_id(self, identifier: str, species: str) -> TranscriptAnnotation | None:
        """Fetch a single transcript annotation by ID with expand=1.

        Makes HTTP GET request to Ensembl lookup/id endpoint with expand=1
        to retrieve detailed transcript information including exons and CDS.

        Args:
            identifier: Ensembl transcript or gene ID (e.g., "ENST00000269305", "TP53")
            species: Species identifier for Ensembl (e.g., "homo_sapiens")

        Returns:
            TranscriptAnnotation object if found, None if not found (404)

        Raises:
            DatabaseAccessError: For network errors, timeouts, or server unavailability (403/503)

        Note:
            The expand=1 parameter causes Ensembl to include Exon and Translation
            objects in the response, which are needed for structural annotation.
        """
        url = f"{self.base_url}/lookup/id/{identifier}?species={species}&expand=1"
        headers = {"Content-Type": "application/json"}

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, headers=headers) as response,
            ):
                if response.status == 404:
                    logger.debug(f"Transcript {identifier} not found in Ensembl")
                    return None

                if response.status == 200:
                    data = cast(dict[str, Any], await response.json())
                    annotation = self._parse_transcript_data(data)

                    # Avoid extra network calls by default. The lookup response usually includes
                    # `external_name` (gene symbol) already; only resolve via gene_id if missing.
                    if (not annotation.symbol) and annotation.gene_id:
                        gene_symbol, gene_interval = await self._fetch_gene_metadata(
                            annotation.gene_id,
                            species,
                            session,
                        )
                        if gene_symbol:
                            annotation.symbol = gene_symbol
                        if gene_interval:
                            annotation.gene_interval = gene_interval

                    return annotation

                if response.status in (403, 502, 503, 504):
                    raise DatabaseAccessError(f"HTTP {response.status}: Access denied or server unavailable", "Ensembl")

                logger.warning(f"Unexpected response status {response.status} for {identifier}")
                return None

        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "Ensembl") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "Ensembl") from e
        except (DatabaseAccessError, ValueError):
            raise
        except Exception as e:
            logger.exception(f"Unexpected error fetching annotation for {identifier}")
            raise DatabaseAccessError(f"Unexpected error: {e}", "Ensembl") from e

    async def _fetch_annotations_by_region(  # noqa: PLR0912
        self, region: str, species: str
    ) -> dict[str, TranscriptAnnotation]:
        """Fetch all transcript annotations overlapping a genomic region.

        Makes HTTP GET request to Ensembl overlap/region endpoint requesting
        transcript, exon, and CDS features. Groups features by transcript ID
        and constructs complete TranscriptAnnotation objects.

        Args:
            region: Genomic region in format "chr:start-end" (e.g., "17:7661779-7687550")
            species: Species identifier for Ensembl (e.g., "homo_sapiens")

        Returns:
            Dictionary mapping transcript IDs to TranscriptAnnotation objects.
            Empty dict if region not found or has no features.

        Raises:
            DatabaseAccessError: For network errors, timeouts, or server unavailability (403/503)

        Note:
            The response is a flat list of features (transcripts, exons, CDS) that
            must be grouped by transcript ID using the Parent field.
        """
        url = (
            f"{self.base_url}/overlap/region/{species}/{region}"
            "?feature=gene;feature=transcript;feature=exon;feature=cds"
        )
        headers = {"Content-Type": "application/json"}

        try:
            async with (
                aiohttp.ClientSession(timeout=aiohttp.ClientTimeout(self.timeout)) as session,
                session.get(url, headers=headers) as response,
            ):
                if response.status == 404:
                    logger.debug(f"Region {region} not found or no features")
                    return {}

                if response.status == 200:
                    data = cast(list[dict[str, Any]], await response.json())
                    annotations = self._parse_region_response(data)

                    # Only resolve gene metadata if the overlap response didn't provide it.
                    # This keeps region queries fast and avoids complicating mocked tests.
                    gene_ids = {a.gene_id for a in annotations.values() if a.gene_id and not a.symbol}
                    if gene_ids:
                        gene_symbol_map: dict[str, str] = {}
                        gene_interval_map: dict[str, Interval] = {}
                        for gene_id in gene_ids:
                            gene_symbol, gene_interval = await self._fetch_gene_metadata(
                                gene_id,
                                species,
                                session,
                            )
                            if gene_symbol:
                                gene_symbol_map[gene_id] = gene_symbol
                            if gene_interval:
                                gene_interval_map[gene_id] = gene_interval

                        for annotation in annotations.values():
                            if annotation.symbol is None:
                                resolved_symbol = gene_symbol_map.get(annotation.gene_id)
                                if resolved_symbol:
                                    annotation.symbol = resolved_symbol
                            resolved_interval = gene_interval_map.get(annotation.gene_id)
                            if resolved_interval:
                                annotation.gene_interval = resolved_interval

                    return annotations

                if response.status in (403, 502, 503, 504):
                    raise DatabaseAccessError(f"HTTP {response.status}: Access denied or server unavailable", "Ensembl")

                logger.warning(f"Unexpected response status {response.status} for region {region}")
                return {}

        except aiohttp.ClientConnectorError as e:
            raise DatabaseAccessError(f"Connection failed: {e}", "Ensembl") from e
        except asyncio.TimeoutError as e:
            raise DatabaseAccessError(f"Request timeout: {e}", "Ensembl") from e
        except (DatabaseAccessError, ValueError):
            raise
        except Exception as e:
            logger.exception(f"Unexpected error fetching region {region}")
            raise DatabaseAccessError(f"Unexpected error: {e}", "Ensembl") from e

    async def _fetch_gene_metadata(
        self,
        gene_id: str,
        species: str,
        session: aiohttp.ClientSession,
    ) -> tuple[str | None, Interval | None]:
        """Resolve gene symbol and coordinates for a gene stable ID via Ensembl lookup."""
        url = f"{self.base_url}/lookup/id/{gene_id}?species={species}"
        headers = {"Content-Type": "application/json"}

        async with session.get(url, headers=headers) as response:
            if response.status == 404:
                return None, None
            if response.status == 200:
                data = cast(dict[str, Any], await response.json())
                symbol = data.get("display_name") or data.get("external_name")
                interval: Interval | None = None
                seq_region_name = data.get("seq_region_name")
                start = data.get("start")
                end = data.get("end")
                strand = data.get("strand")
                if seq_region_name and start and end:
                    interval = Interval(
                        seq_region_name=str(seq_region_name),
                        start=int(start),
                        end=int(end),
                        strand=int(strand) if strand is not None else None,
                    )
                return (str(symbol) if symbol else None), interval
            if response.status in (403, 502, 503, 504):
                raise DatabaseAccessError(f"HTTP {response.status}: Access denied or server unavailable", "Ensembl")
            return None, None

    def _parse_transcript_data(self, data: Mapping[str, Any]) -> TranscriptAnnotation:
        """Parse Ensembl lookup/id response into TranscriptAnnotation.

        Extracts transcript metadata from Ensembl JSON response including:
        - Basic identifiers (transcript_id, gene_id, symbol, biotype)
        - Genomic coordinates (chr, start, end, strand)
        - Structural features (exons from Exon array, CDS from Translation object)

        Args:
            data: JSON response dict from Ensembl /lookup/id endpoint with expand=1

        Returns:
            TranscriptAnnotation object with all available fields populated

        Note:
            - Exons are parsed from the "Exon" array if present
            - CDS intervals are derived from "Translation" object if present
            - Provider/endpoint/reference metadata fields are set by caller
        """
        # Extract basic info
        transcript_id = str(data.get("id", ""))
        gene_id = str(data.get("Parent", "") or data.get("gene_id", ""))
        # Transcript-level `display_name` is typically gene+isoform (e.g., "TP53-201");
        # gene symbol is resolved separately from the parent gene record.
        symbol = data.get("external_name")
        biotype = data.get("biotype")

        # Genomic coordinates
        seq_region_name = str(data.get("seq_region_name", ""))
        start = int(data.get("start", 0))
        end = int(data.get("end", 0))
        strand = int(data.get("strand", 1))

        # Parse exons
        exons: list[Interval] = []
        if "Exon" in data and isinstance(data["Exon"], list):
            for exon_data in data["Exon"]:
                if isinstance(exon_data, dict):
                    exon = Interval(
                        seq_region_name=str(exon_data.get("seq_region_name", seq_region_name)),
                        start=int(exon_data.get("start", 0)),
                        end=int(exon_data.get("end", 0)),
                        strand=int(exon_data.get("strand", strand)),
                    )
                    exons.append(exon)

        # Parse CDS (coding sequence intervals)
        cds_intervals: list[Interval] = []
        if "Translation" in data and isinstance(data["Translation"], dict):
            translation = data["Translation"]
            if "start" in translation and "end" in translation:
                # For protein-coding transcripts, compute CDS from translation coordinates
                cds_start = int(translation.get("start", 0))
                cds_end = int(translation.get("end", 0))
                if cds_start > 0 and cds_end > 0:
                    cds_intervals.append(
                        Interval(
                            seq_region_name=seq_region_name,
                            start=cds_start,
                            end=cds_end,
                            strand=strand,
                        )
                    )

        return TranscriptAnnotation(
            transcript_id=transcript_id,
            gene_id=gene_id,
            symbol=symbol,
            biotype=biotype,
            seq_region_name=seq_region_name,
            start=start,
            end=end,
            strand=strand,
            exons=exons,
            cds=cds_intervals,
            provider="ensembl_rest",
            endpoint=None,  # Will be set by caller
            reference_choice=None,  # Will be set by caller
        )

    def _parse_region_response(self, features: list[dict[str, Any]]) -> dict[str, TranscriptAnnotation]:  # noqa: PLR0912
        """Parse Ensembl region overlap response into transcript annotations.

        The region endpoint returns a flat list of features (transcripts, exons, CDS).
        This method groups them by transcript ID using the Parent field and constructs
        complete TranscriptAnnotation objects.

        Args:
            features: List of feature dicts from /overlap/region endpoint, each with:
                - feature_type: "transcript", "exon", or "cds"
                - id: Feature stable ID
                - Parent: Parent transcript ID (for exons and CDS)
                - genomic coordinates and other metadata

        Returns:
            Dict mapping transcript IDs to complete TranscriptAnnotation objects

        Algorithm:
            1. First pass: Group features by transcript ID into transcript/exons/cds buckets
            2. Second pass: Build TranscriptAnnotation from each grouped set
            3. Provider/endpoint/reference metadata fields set by caller
        """
        # Group features by transcript ID. The API does not guarantee ordering,
        # so we allow exons/CDS to appear before transcript records.
        transcript_features: dict[str, dict[str, Any]] = {}
        gene_features: dict[str, dict[str, Any]] = {}

        def _ensure_bucket(transcript_id: str) -> dict[str, Any]:
            return transcript_features.setdefault(
                transcript_id,
                {
                    "transcript": None,
                    "exons": [],
                    "cds": [],
                },
            )

        for feature in features:
            feature_type = str(feature.get("feature_type", "") or "")

            if feature_type == "gene":
                gene_id = str(feature.get("id", "") or "")
                if gene_id:
                    gene_features[gene_id] = feature
                continue

            if feature_type == "transcript":
                transcript_id = str(feature.get("id", "") or "")
                if transcript_id:
                    _ensure_bucket(transcript_id)["transcript"] = feature
                continue

            if feature_type in {"exon", "cds"}:
                parent_id = str(feature.get("Parent", "") or "")
                if not parent_id:
                    continue
                bucket = _ensure_bucket(parent_id)
                bucket["exons" if feature_type == "exon" else "cds"].append(feature)
                continue

        # Build TranscriptAnnotation objects
        transcripts: dict[str, TranscriptAnnotation] = {}

        for transcript_id, grouped in transcript_features.items():
            transcript_data = grouped["transcript"]
            if not isinstance(transcript_data, dict):
                # Can't build a transcript annotation without a transcript record.
                continue

            # Basic info
            gene_id = str(transcript_data.get("Parent", "") or transcript_data.get("gene_id", ""))

            # Prefer the gene feature's symbol when available. In practice, transcript
            # features can carry a transcript name like "TP53-201" in `external_name`,
            # while the gene feature carries the canonical gene symbol "TP53".
            gene_feature = gene_features.get(gene_id) if gene_id else None
            gene_symbol = None
            gene_interval: Interval | None = None
            if isinstance(gene_feature, dict):
                gene_symbol = (
                    gene_feature.get("external_name")
                    or gene_feature.get("display_name")
                    or gene_feature.get("gene_name")
                )
                seq_region_name_g = gene_feature.get("seq_region_name")
                start_g = gene_feature.get("start")
                end_g = gene_feature.get("end")
                strand_g = gene_feature.get("strand")
                if seq_region_name_g and start_g and end_g:
                    gene_interval = Interval(
                        seq_region_name=str(seq_region_name_g),
                        start=int(start_g),
                        end=int(end_g),
                        strand=int(strand_g) if strand_g is not None else None,
                    )

            transcript_name = (
                transcript_data.get("external_name")
                or transcript_data.get("display_name")
                or transcript_data.get("gene_name")
            )
            symbol = gene_symbol or transcript_data.get("gene_name")
            if (not symbol) and isinstance(transcript_name, str) and "-" in transcript_name:
                # Common Ensembl convention: transcript name like "TP53-201".
                symbol = transcript_name.split("-", 1)[0]
            biotype = transcript_data.get("biotype")

            # Genomic coordinates
            seq_region_name = str(transcript_data.get("seq_region_name", ""))
            start = int(transcript_data.get("start", 0))
            end = int(transcript_data.get("end", 0))
            strand = int(transcript_data.get("strand", 1))

            # Build exon intervals
            exons: list[Interval] = []
            for exon_data in grouped["exons"]:
                exon = Interval(
                    seq_region_name=str(exon_data.get("seq_region_name", seq_region_name)),
                    start=int(exon_data.get("start", 0)),
                    end=int(exon_data.get("end", 0)),
                    strand=int(exon_data.get("strand", strand)),
                )
                exons.append(exon)

            # Build CDS intervals
            cds_intervals: list[Interval] = []
            for cds_data in grouped["cds"]:
                cds = Interval(
                    seq_region_name=str(cds_data.get("seq_region_name", seq_region_name)),
                    start=int(cds_data.get("start", 0)),
                    end=int(cds_data.get("end", 0)),
                    strand=int(cds_data.get("strand", strand)),
                )
                cds_intervals.append(cds)

            annotation = TranscriptAnnotation(
                transcript_id=transcript_id,
                gene_id=gene_id,
                symbol=symbol,
                biotype=biotype,
                seq_region_name=seq_region_name,
                start=start,
                end=end,
                strand=strand,
                gene_interval=gene_interval,
                exons=exons,
                cds=cds_intervals,
                provider="ensembl_rest",
                endpoint=None,  # Will be set by caller
                reference_choice=None,  # Will be set by caller
            )

            transcripts[transcript_id] = annotation

        return transcripts

    @staticmethod
    def _normalize_species(species: str) -> str:
        """Normalize species name for Ensembl API.

        Converts common species names and formats to Ensembl's expected format
        (lowercase with underscores, e.g., "homo_sapiens").

        Args:
            species: Species name in various formats:
                - Common name: "human", "mouse", "rat"
                - Scientific name: "Homo sapiens", "Mus musculus"
                - Ensembl format: "homo_sapiens" (returned as-is)

        Returns:
            Ensembl-compatible species identifier (e.g., "homo_sapiens")

        Examples:
            >>> _normalize_species("human")
            'homo_sapiens'
            >>> _normalize_species("Homo sapiens")
            'homo_sapiens'
            >>> _normalize_species("custom_species")
            'custom_species'
        """
        # Map common names to Ensembl species names
        species_map = {
            "human": "homo_sapiens",
            "mouse": "mus_musculus",
            "rat": "rattus_norvegicus",
            "rhesus": "macaca_mulatta",
            "macaque": "macaca_mulatta",
            "pig": "sus_scrofa",
            "chicken": "gallus_gallus",
        }

        normalized = species.lower().strip().replace(" ", "_")
        return species_map.get(normalized, normalized)


class VepConsequenceClient:
    """Optional VEP (Variant Effect Predictor) consequence enrichment client.

    Provides additional functional annotation for transcript variants.
    This is an optional enhancement and not required for base functionality.

    **Current Status: PLACEHOLDER**

    This client exists as a stub for future VEP integration. The `enrich_annotations`
    method currently returns the input bundle unchanged.

    **Future Implementation:**

    When activated (via config flag), this client will:
    1. Query Ensembl VEP REST API for consequence predictions
    2. Enrich TranscriptAnnotation objects with:
       - Variant consequence types (missense, nonsense, etc.)
       - Conservation scores
       - Regulatory feature overlaps
       - Population frequency data
    3. Maintain consistent caching strategy with EnsemblTranscriptModelClient

    **Design Rationale:**

    Separated from EnsemblTranscriptModelClient because:
    - VEP queries are expensive (rate-limited, slower)
    - Not all workflows need consequence predictions
    - Allows independent caching strategies
    - Can be enabled/disabled via configuration
    """

    def __init__(self, timeout: int = 30, base_url: str = "https://rest.ensembl.org"):
        """Initialize VEP client.

        Args:
            timeout: Request timeout in seconds
            base_url: Ensembl REST API base URL
        """
        self.timeout = timeout
        self.base_url = base_url
        logger.info("VEP enrichment client initialized (optional feature)")

    async def enrich_annotations(
        self,
        bundle: TranscriptAnnotationBundle,
        _species: str = "homo_sapiens",
    ) -> TranscriptAnnotationBundle:
        """Enrich transcript annotations with VEP consequence data.

        Args:
            bundle: Existing transcript annotation bundle
            species: Species name for VEP queries

        Returns:
            Enriched bundle (currently returns input unchanged - placeholder for future VEP integration)
        """
        # Placeholder for future VEP enrichment
        # This would query VEP REST API for consequence predictions
        logger.debug(f"VEP enrichment called for {bundle.resolved_count} transcripts (not yet implemented)")
        return bundle
