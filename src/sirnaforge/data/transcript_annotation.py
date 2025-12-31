"""Transcript annotation providers using Ensembl REST and optional VEP enrichment."""

import asyncio
from collections.abc import Mapping
from time import time
from typing import Any, Optional, cast

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

    def _get_cached(self, key: str) -> Optional[Any]:
        """Retrieve item from cache if present and not expired."""
        if key not in self._cache:
            return None

        value, timestamp = self._cache[key]
        if time() - timestamp > self.cache_ttl:
            del self._cache[key]
            return None

        return value

    def _set_cache(self, key: str, value: Any) -> None:
        """Store item in cache with current timestamp."""
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

    async def _fetch_annotation_by_id(self, identifier: str, species: str) -> Optional[TranscriptAnnotation]:
        """Fetch a single transcript annotation by ID with expand=1."""
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
                    return self._parse_transcript_data(data)

                if response.status in (403, 502, 503, 504):
                    raise DatabaseAccessError(
                        f"HTTP {response.status}: Access denied or server unavailable", "Ensembl"
                    )

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

    async def _fetch_annotations_by_region(
        self, region: str, species: str
    ) -> dict[str, TranscriptAnnotation]:
        """Fetch all transcript annotations overlapping a genomic region."""
        url = f"{self.base_url}/overlap/region/{species}/{region}?feature=transcript;feature=exon;feature=cds"
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
                    return self._parse_region_response(data)

                if response.status in (403, 502, 503, 504):
                    raise DatabaseAccessError(
                        f"HTTP {response.status}: Access denied or server unavailable", "Ensembl"
                    )

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

    def _parse_transcript_data(self, data: Mapping[str, Any]) -> TranscriptAnnotation:
        """Parse Ensembl lookup response into TranscriptAnnotation."""
        # Extract basic info
        transcript_id = str(data.get("id", ""))
        gene_id = str(data.get("Parent", "") or data.get("gene_id", ""))
        symbol = data.get("display_name") or data.get("external_name")
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

    def _parse_region_response(self, features: list[dict[str, Any]]) -> dict[str, TranscriptAnnotation]:
        """Parse Ensembl region overlap response into transcript annotations.

        The region endpoint returns a flat list of features (transcripts, exons, CDS).
        We need to group them by transcript ID and build complete annotations.
        """
        # Group features by transcript ID
        transcript_features: dict[str, dict[str, Any]] = {}

        for feature in features:
            feature_type = feature.get("feature_type", "")

            if feature_type == "transcript":
                transcript_id = str(feature.get("id", ""))
                if transcript_id:
                    transcript_features.setdefault(transcript_id, {
                        "transcript": feature,
                        "exons": [],
                        "cds": [],
                    })
            elif feature_type == "exon":
                parent_id = feature.get("Parent")
                if parent_id and parent_id in transcript_features:
                    transcript_features[parent_id]["exons"].append(feature)
            elif feature_type == "cds":
                parent_id = feature.get("Parent")
                if parent_id and parent_id in transcript_features:
                    transcript_features[parent_id]["cds"].append(feature)

        # Build TranscriptAnnotation objects
        transcripts: dict[str, TranscriptAnnotation] = {}

        for transcript_id, grouped in transcript_features.items():
            transcript_data = grouped["transcript"]

            # Basic info
            gene_id = str(transcript_data.get("Parent", "") or transcript_data.get("gene_id", ""))
            symbol = transcript_data.get("external_name") or transcript_data.get("gene_name")
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
        """Normalize species name for Ensembl API."""
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
