#!/usr/bin/env python3
"""miRNA Database Manager with multi-species support.

This module provides a clean interface for downloading, caching, and managing
miRNA databases from multiple sources (MirGeneDB, miRBase, TargetScan) with
automatic cache management and species-specific organization.
"""

import argparse
import hashlib
import json
import logging
import math
import os
import re
import tempfile
import time
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sirnaforge.utils.cache_utils import (
    ARTIFACT_MIRNA_COMBINED,
    fingerprint_inputs,
    is_artifact_stamp_current,
    write_artifact_stamp,
)

from .reference_manager import CacheMetadata, ReferenceManager, ReferenceSource
from .species_registry import (
    CANONICAL_SPECIES_ALIAS_MAP,
    CANONICAL_SPECIES_REGISTRY,
    MIRGENEDB_ALIAS_MAP,
    MIRGENEDB_SPECIES_TABLE,
)

logger = logging.getLogger(__name__)
# Plain "fasta" only. Biopython's "fasta-blast" reader assigns every record the *first*
# record's id/description, which silently collapses miRBase's ~270 species prefixes into
# one, and it does not exist at all below Biopython 1.85.
FASTA_PARSE_FORMAT = "fasta"
_MIRBASE_SPECIES_CODES = {
    "human": "hsa-",
    "mouse": "mmu-",
    "rat": "rno-",
    "zebrafish": "dre-",
    "fly": "dme-",
    "worm": "cel-",
    "chicken": "gga-",
    "dog": "cfa-",
    "pig": "ssc-",
    "rhesus": "mml-",
}


@dataclass
class MiRNASource(ReferenceSource):
    """miRNA-specific database source configuration.

    Inherits from ReferenceSource and can add miRNA-specific fields if needed.
    """

    pass


@dataclass(frozen=True)
class PinnedMiRNAReference:
    """Explicit local mature FASTA with independently supplied release and SHA-256.

    The checksum covers the original file, before species filtering. A pin never
    refreshes from an upstream endpoint and is not advertised as the latest release.
    """

    path: str | Path
    release: str
    sha256: str
    original_url: str | None = None

    def __post_init__(self) -> None:
        """Require explicit release identity and a full digest."""
        if not isinstance(self.release, str) or not self.release.strip():
            raise ValueError("Pinned miRNA reference requires a release")
        if not isinstance(self.sha256, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", self.sha256):
            raise ValueError("Pinned miRNA reference requires a SHA-256 checksum (64 hex characters)")
        if not isinstance(self.path, str | Path) or not str(self.path).strip():
            raise ValueError("Pinned miRNA reference requires a local FASTA path")
        if self.original_url is not None and not isinstance(self.original_url, str):
            raise ValueError("Pinned miRNA reference original_url must be a string")


def _atomic_write_text(destination: Path, content: str) -> None:
    """Install validated bytes without truncating an existing reference on failure."""
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=destination.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _build_mirbase_sources() -> dict[str, MiRNASource]:
    """Register mature miRBase sequences for every canonical screening species."""
    return {
        entry["genome"]: MiRNASource(
            name="mirbase_mature",
            url="https://www.mirbase.org/download/mature.fa",
            species=entry["genome"],
            format="fasta",
            compressed=False,
            description=f"All mature miRBase sequences, filtered for {species}",
        )
        for species, entry in CANONICAL_SPECIES_REGISTRY.items()
    }


def _build_combined_sources() -> dict[str, MiRNASource]:
    """Describe the explicit mature-only union of MirGeneDB and miRBase."""
    return {
        species: MiRNASource(
            name="mirgenedb+mirbase",
            url="combined:mirgenedb+mirbase",
            species=species,
            format="fasta",
            compressed=False,
            description="MirGeneDB mature plus full miRBase mature sequences (no star-only export)",
        )
        for species in _build_mirbase_sources()
    }


def _build_mirgenedb_sources() -> dict[str, MiRNASource]:
    """Construct MirGeneDB source map from the slug-based metadata table."""
    sources: dict[str, MiRNASource] = {}
    base_url = "https://www.mirgenedb.org/fasta/{slug}?mat=1"

    for slug, metadata in MIRGENEDB_SPECIES_TABLE.items():
        scientific_name = metadata.get("scientific_name", slug)
        taxonomy_id = metadata.get("taxonomy_id")
        common_name = metadata.get("common_name")
        description = f"MirGeneDB high-confidence miRNAs ({scientific_name}, NCBI:{taxonomy_id})"
        if common_name and common_name.lower() not in scientific_name.lower():
            description = f"{description} [{common_name}]"

        sources[slug] = MiRNASource(
            name="mirgenedb",
            url=base_url.format(slug=slug),
            species=slug,
            format="fasta",
            compressed=False,
            description=description,
        )

    return sources


class MiRNADatabaseManager(ReferenceManager[MiRNASource]):
    """Elegant miRNA database manager with caching and multi-species support."""

    # Database source configurations
    SOURCES = {
        "mirbase": _build_mirbase_sources(),
        "mirgenedb+mirbase": _build_combined_sources(),
        "mirbase_high_conf": {
            "human": MiRNASource(
                name="mirbase_mature_hc",
                url="https://www.mirbase.org/download/CURRENT/mature_high_conf.fa",
                species="human",
                format="fasta",
                compressed=False,
                description="miRBase high-confidence mature miRNA sequences (Homo sapiens - hsa)",
            ),
            "mouse": MiRNASource(
                name="mirbase_mature_hc",
                url="https://www.mirbase.org/download/CURRENT/mature_high_conf.fa",
                species="mouse",
                format="fasta",
                compressed=False,
                description="miRBase high-confidence mature miRNA sequences (Mus musculus - mmu)",
            ),
            "rat": MiRNASource(
                name="mirbase_mature_hc",
                url="https://www.mirbase.org/download/CURRENT/mature_high_conf.fa",
                species="rat",
                format="fasta",
                compressed=False,
                description="miRBase high-confidence mature miRNA sequences (Rattus norvegicus - rno)",
            ),
        },
        "mirbase_hairpin": {
            "human": MiRNASource(
                name="mirbase_hairpin",
                url="https://www.mirbase.org/download/hairpin.fa",
                species="human",
                format="fasta",
                compressed=False,
                description="miRBase hairpin precursor miRNA sequences (Homo sapiens - hsa)",
            ),
            "mouse": MiRNASource(
                name="mirbase_hairpin",
                url="https://www.mirbase.org/download/hairpin.fa",
                species="mouse",
                format="fasta",
                compressed=False,
                description="miRBase hairpin precursor miRNA sequences (Mus musculus - mmu)",
            ),
            "rat": MiRNASource(
                name="mirbase_hairpin",
                url="https://www.mirbase.org/download/hairpin.fa",
                species="rat",
                format="fasta",
                compressed=False,
                description="miRBase hairpin precursor miRNA sequences (Rattus norvegicus - rno)",
            ),
        },
        "mirgenedb": _build_mirgenedb_sources(),
        "targetscan": {
            "human": MiRNASource(
                name="targetscan",
                url="https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip",
                species="human",
                format="tsv",
                compressed=True,
                description="TargetScan miRNA family data",
            )
        },
    }

    @classmethod
    def get_available_sources(cls) -> list[str]:
        """Return sorted list of supported database sources."""
        return sorted(cls.SOURCES.keys())

    @classmethod
    def get_all_species(cls) -> list[str]:
        """Return sorted list of all species across sources."""
        species: set[str] = set()
        for species_map in cls.SOURCES.values():
            species.update(species_map.keys())
        return sorted(species)

    @classmethod
    def get_species_for_source(cls, source_name: str) -> list[str]:
        """Return sorted list of species supported by a given source."""
        return sorted(cls.SOURCES.get(source_name, {}).keys())

    @classmethod
    def get_species_aliases(cls, source_name: str) -> dict[str, list[str]]:
        """Return mapping of canonical species identifiers to their known aliases."""
        if source_name != "mirgenedb":
            return {species: [species] for species in cls.SOURCES.get(source_name, {})}

        aliases: dict[str, list[str]] = {}
        for slug, metadata in MIRGENEDB_SPECIES_TABLE.items():
            alias_values = metadata.get("aliases", [])
            aliases[slug] = sorted({slug, *[alias.lower() for alias in alias_values]})
        return aliases

    @classmethod
    def get_canonical_species(cls) -> list[str]:
        """Return sorted list of canonical species keys."""
        return sorted(CANONICAL_SPECIES_REGISTRY.keys())

    @classmethod
    def canonicalize_species_name(cls, species: str) -> str | None:
        """Normalize a raw species identifier to a canonical key."""
        if not species:
            return None
        return CANONICAL_SPECIES_ALIAS_MAP.get(species.lower())

    @classmethod
    def canonicalize_species_list(cls, species_list: Sequence[str]) -> list[str]:
        """Normalize a list of species identifiers to canonical keys, preserving order."""
        canonical: list[str] = []
        unknown: list[str] = []

        for raw_value in species_list:
            key = cls.canonicalize_species_name(raw_value)
            if key is None:
                unknown.append(raw_value)
                continue
            if key not in canonical:
                canonical.append(key)

        if unknown:
            raise ValueError(f"Unsupported species: {', '.join(unknown)}")

        return canonical

    @classmethod
    def get_genome_species_for_canonical(cls, canonical_species: Sequence[str]) -> list[str]:
        """Return genome species identifiers for canonical species keys.

        Note: These are used for miRNA database lookups, not genomic DNA alignment.
        The term 'genome' here refers to the organism's miRNA annotation set.
        """
        genome_species: list[str] = []
        for key in canonical_species:
            registry_entry = CANONICAL_SPECIES_REGISTRY.get(key)
            if not registry_entry:
                raise ValueError(f"Unknown canonical species '{key}'")
            genome_name = registry_entry.get("genome")
            if genome_name and genome_name not in genome_species:
                genome_species.append(genome_name)
        return genome_species

    @classmethod
    def get_mirna_slugs_for_canonical(cls, canonical_species: Sequence[str], source_name: str) -> list[str]:
        """Return normalized miRNA identifiers for canonical species."""
        slugs: list[str] = []
        for key in canonical_species:
            registry_entry = CANONICAL_SPECIES_REGISTRY.get(key)
            if not registry_entry:
                raise ValueError(f"Unknown canonical species '{key}'")

            if source_name == "mirgenedb":
                candidate = registry_entry.get("mirgenedb_slug")
                if not candidate:
                    raise ValueError(f"Canonical species '{key}' missing MirGeneDB slug mapping")
            else:
                candidate = registry_entry.get("genome")
                if not candidate:
                    raise ValueError(f"Canonical species '{key}' missing genome mapping")

            normalized = cls.normalize_species(source_name, candidate)
            if normalized is None:
                normalized = candidate

            if not cls.is_supported_species(source_name, normalized):
                raise ValueError(f"Species '{normalized}' not available for source '{source_name}'")
            if normalized not in slugs:
                slugs.append(normalized)

        return slugs

    @classmethod
    def get_supported_canonical_species_for_source(cls, source_name: str) -> list[str]:
        """Return canonical species supported by a given source."""
        supported: list[str] = []
        if source_name == "mirgenedb":
            for canonical_name, registry_entry in CANONICAL_SPECIES_REGISTRY.items():
                slug = registry_entry.get("mirgenedb_slug")
                if slug and slug in cls.SOURCES.get(source_name, {}):
                    supported.append(canonical_name)
        else:
            available_species = cls.SOURCES.get(source_name, {})
            for canonical_name, registry_entry in CANONICAL_SPECIES_REGISTRY.items():
                genome_name = registry_entry.get("genome")
                if genome_name and genome_name in available_species:
                    supported.append(canonical_name)

        return sorted(supported)

    @classmethod
    def resolve_species_selection(
        cls,
        requested_species: Sequence[str],
        source_name: str,
        mirna_overrides: Sequence[str] | None = None,
    ) -> dict[str, list[str]]:
        """Resolve canonical, genome, and miRNA identifiers for the requested species."""
        canonical_species = cls.canonicalize_species_list(requested_species)
        genome_species = cls.get_genome_species_for_canonical(canonical_species)

        if mirna_overrides:
            overridden: list[str] = []
            unknown: list[str] = []
            for raw_value in mirna_overrides:
                normalized = cls.normalize_species(source_name, raw_value)
                if normalized is None:
                    unknown.append(raw_value)
                    continue
                if not cls.is_supported_species(source_name, normalized):
                    unknown.append(raw_value)
                    continue
                if normalized not in overridden:
                    overridden.append(normalized)

            if unknown:
                raise ValueError(f"Unsupported miRNA species for source '{source_name}': {', '.join(unknown)}")
            mirna_species = overridden
        else:
            mirna_species = cls.get_mirna_slugs_for_canonical(canonical_species, source_name)

        return {
            "canonical": canonical_species,
            "genome": genome_species,  # miRNA genome annotations, not genomic DNA
            "mirna": mirna_species,
        }

    @classmethod
    def _normalize_mirgenedb_species(cls, species: str) -> str | None:
        if not species:
            return None
        return MIRGENEDB_ALIAS_MAP.get(species.lower())

    @classmethod
    def normalize_species(cls, source_name: str, species: str) -> str | None:
        """Normalize user-provided species identifiers to canonical keys."""
        if source_name == "mirgenedb":
            return cls._normalize_mirgenedb_species(species)
        if source_name in {"mirbase", "mirgenedb+mirbase"}:
            canonical = cls.canonicalize_species_name(species)
            if canonical:
                return str(CANONICAL_SPECIES_REGISTRY[canonical]["genome"])
        return species

    @classmethod
    def get_source_configuration(cls, source_name: str, species: str) -> MiRNASource | None:
        """Retrieve the MiRNASource configuration for a given source/species."""
        canonical_species = cls.normalize_species(source_name, species)
        if canonical_species is None:
            return None
        return cls.SOURCES.get(source_name, {}).get(canonical_species)

    @classmethod
    def is_supported_source(cls, source_name: str) -> bool:
        """Check if a source is supported."""
        return source_name in cls.SOURCES

    @classmethod
    def is_supported_species(cls, source_name: str, species: str) -> bool:
        """Check if a species is supported for the given source."""
        canonical_species = cls.normalize_species(source_name, species)
        if canonical_species is None:
            return False
        return canonical_species in cls.SOURCES.get(source_name, {})

    @classmethod
    def get_mirgenedb_species_metadata(cls) -> dict[str, dict[str, Any]]:
        """Expose the MirGeneDB species metadata table."""
        return {key: dict(value) for key, value in MIRGENEDB_SPECIES_TABLE.items()}

    def __init__(
        self,
        cache_dir: str | Path | None = None,
        cache_ttl_days: int = 30,
        *,
        pinned_references: dict[str, dict[str, PinnedMiRNAReference]] | None = None,
        download_timeout: int = 30,
        download_attempts: int = 2,
        retry_backoff_seconds: float = 0.25,
    ):
        """Initialize the miRNA database manager.

        Args:
            cache_dir: Directory for caching databases (default: ~/.cache/sirnaforge/mirna)
            cache_ttl_days: Cache time-to-live in days
            pinned_references: Explicit mature FASTA pins keyed by source and species.
            download_timeout: Socket timeout per live request, from 1 to 300 seconds.
            download_attempts: Maximum acquisition attempts, from 1 to 3.
            retry_backoff_seconds: Base exponential backoff, from 0 to 5 seconds.
        """
        if type(download_timeout) is not int or not 1 <= download_timeout <= 300:
            raise ValueError("miRNA download_timeout must be an integer from 1 to 300 seconds")
        if type(download_attempts) is not int or not 1 <= download_attempts <= 3:
            raise ValueError("miRNA download_attempts must be an integer from 1 to 3")
        if not math.isfinite(retry_backoff_seconds) or not 0 <= retry_backoff_seconds <= 5:
            raise ValueError("miRNA retry_backoff_seconds must be from 0 to 5 seconds")
        self.download_timeout = download_timeout
        self.download_attempts = download_attempts
        self.retry_backoff_seconds = retry_backoff_seconds
        super().__init__(cache_subdir="mirna", cache_dir=cache_dir, cache_ttl_days=cache_ttl_days)
        self.pinned_references: dict[str, dict[str, PinnedMiRNAReference]] = {}
        self._resolved_reference_metadata: dict[tuple[str, str], dict[str, Any]] = {}
        for source_name, species_pins in (pinned_references or {}).items():
            if source_name not in {"mirbase", "mirgenedb"}:
                raise ValueError(f"Pinned mature references require mirbase or mirgenedb, got {source_name!r}")
            for species, pin in species_pins.items():
                normalized = self.normalize_species(source_name, species)
                if normalized is None or not self.is_supported_species(source_name, normalized):
                    raise ValueError(f"Unsupported pinned miRNA source/species: {source_name}/{species}")
                if normalized in self.pinned_references.setdefault(source_name, {}):
                    raise ValueError(f"Duplicate pinned miRNA source/species: {source_name}/{normalized}")
                self.pinned_references[source_name][normalized] = pin

    def get_reference_metadata(self, source_name: str, species: str) -> dict[str, Any] | None:
        """Return provenance of a reference resolved by this manager in this run.

        Call after ``get_database``. None means acquisition failed or was never
        attempted. This never inspects a different process's or host's cache.
        """
        import copy  # noqa: PLC0415

        normalized = self.normalize_species(source_name, species)
        return copy.deepcopy(self._resolved_reference_metadata.get((source_name, normalized or species)))

    @staticmethod
    def _validate_mature_fasta(content: str) -> None:
        """Reject malformed payloads before they can replace a usable reference."""
        if (
            not content.strip().startswith(">")
            or "<" in content
            or ">" in "\n".join(line for line in content.splitlines() if not line.startswith(">"))
        ):
            raise ValueError("Expected plain mature FASTA, not HTML or malformed content")
        records = list(SeqIO.parse(StringIO(content), FASTA_PARSE_FORMAT))
        if not records or len({record.id for record in records}) != len(records):
            raise ValueError("Mature FASTA must contain records with unique identifiers")
        for record in records:
            if not record.seq or not re.fullmatch(r"[ACGTUNacgtun]+", str(record.seq)):
                raise ValueError(f"Invalid mature FASTA sequence: {record.id}")

    def _get_pinned_database(
        self, source_name: str, species: str, pin: PinnedMiRNAReference, force_refresh: bool
    ) -> Path | None:
        source = MiRNASource(
            name=source_name,
            url=f"pinned:sha256:{pin.sha256.lower()}:release:{pin.release}",
            species=species,
            description=f"Pinned mature {source_name} release {pin.release}",
        )
        key = source.cache_key()
        cache_file = self.cache_dir / f"{key}.fa"
        meta = self.metadata.get(key)
        # Explicit immutable pins do not expire by age. Integrity and producer
        # checks still apply; force_refresh explicitly reimports the local file.
        if (
            not force_refresh
            and meta is not None
            and self.is_producer_version_current(meta)
            and cache_file.is_file()
            and cache_file.stat().st_size > 0
            and self._compute_file_checksum(cache_file) == meta.checksum
            and meta.extra is not None
            and hashlib.sha256(cache_file.read_bytes()).hexdigest() == meta.extra.get("filtered_sha256")
        ):
            self._resolved_reference_metadata[source_name, species] = {**meta.extra, "acquisition_mode": "pinned_cache"}
            return cache_file
        try:
            raw = Path(pin.path).read_bytes()
            raw_sha256 = hashlib.sha256(raw).hexdigest()
            if raw_sha256 != pin.sha256.lower():
                raise ValueError("Pinned mature FASTA SHA-256 checksum mismatch")
            content = raw.decode("utf-8")
            self._validate_mature_fasta(content)
            if source_name == "mirbase":
                content = self._filter_species_sequences(content, species)
                if not all(
                    record.id.lower().startswith(_MIRBASE_SPECIES_CODES[species])
                    for record in SeqIO.parse(StringIO(content), FASTA_PARSE_FORMAT)
                ):
                    raise ValueError(f"Pinned miRBase FASTA contains mismatched species identifiers for {species}")
            # MirGeneDB exports are already species-specific. Verify record
            # IDs against the configured species instead of trusting the label.
            elif not all(
                record.id.lower().startswith(f"{species.lower()}-")
                for record in SeqIO.parse(StringIO(content), FASTA_PARSE_FORMAT)
            ):
                raise ValueError(f"Pinned MirGeneDB FASTA contains records outside species {species}")
            self._validate_mature_fasta(content)
            provenance = {
                "source": source_name,
                "release": pin.release,
                "original_url": pin.original_url,
                "acquisition_url": Path(pin.path).resolve().as_uri(),
                "raw_sha256": raw_sha256,
                "filtered_sha256": hashlib.sha256(content.encode("utf-8")).hexdigest(),
                "acquired_at": self._get_current_timestamp(),
                "species": species,
                "species_filtered": source_name == "mirbase",
                "acquisition_mode": "pinned_import",
            }
            _atomic_write_text(cache_file, content)
            self._record_cache_entry(key, source, cache_file, extra=provenance, persist=True)
            self._resolved_reference_metadata[source_name, species] = provenance
            return cache_file
        except (OSError, ValueError, UnicodeError) as exc:
            logger.error(
                "Pinned miRNA reference %s/%s unavailable: %s; supply the verified mature FASTA",
                source_name,
                species,
                exc,
            )
            return None

    def _metadata_from_dict(self, data: dict) -> CacheMetadata:
        """Create CacheMetadata from dictionary using MiRNASource."""
        return CacheMetadata.from_dict(data, source_class=MiRNASource)

    def _filter_species_sequences(self, fasta_content: str, species: str) -> str:
        """Filter FASTA content for specific species using miRBase three-letter codes.

        Returns the matching records as FASTA text, or an empty string when the species is
        absent from the payload (callers treat that as a download failure).

        Raises:
            ValueError: if no miRBase code is known for `species`.
        """
        # Species prefix mapping - simplified and clear
        species_codes = _MIRBASE_SPECIES_CODES

        code = species_codes.get(species)
        if not code:
            # Falling back to "return every species" would hand ~48k miRNAs from ~270 species
            # to seed screening under a single-species label, which looks like a successful run.
            # `get_database` only ever passes species it has already validated against SOURCES,
            # so reaching this is a programming error and must be loud rather than permissive.
            raise ValueError(
                f"No miRBase species code known for '{species}'; supported: {', '.join(sorted(species_codes))}"
            )

        # Some upstream endpoints occasionally return HTML-wrapped FASTA lines.
        normalized_content = re.sub(r"</?[^>\n]+>", "", fasta_content)
        # Block-style wrappers (a bare `<pre>` line around the whole payload) leave blank lines
        # ahead of the first record. Biopython 1.86's "fasta" reader treats anything before the
        # first ">" as a deprecated leading comment (BiopythonDeprecationWarning) and a future
        # release turns that into a ValueError, so drop the preamble instead of parsing it.
        lines = normalized_content.splitlines()
        first_record = next((index for index, line in enumerate(lines) if line.startswith(">")), None)
        normalized_content = "\n".join(lines[first_record:]) + "\n" if first_record is not None else ""

        input_handle = StringIO(normalized_content)
        output_handle = StringIO()
        filtered_count = 0

        for record in SeqIO.parse(input_handle, FASTA_PARSE_FORMAT):
            if code in record.description.lower():
                SeqIO.write(record, output_handle, "fasta")
                filtered_count += 1

        if filtered_count == 0:
            logger.error("No sequences found for species '%s' after filtering", species)

        logger.info(f"Filtered to {filtered_count} {species} sequences")
        return output_handle.getvalue().strip()

    def _download_available_source(self, source_name: str, source: MiRNASource) -> str | None:
        """Download supported references without substituting for a retired confidence subset."""
        if source_name == "mirbase_high_conf":
            logger.error(
                "miRBase no longer distributes the configured high-confidence FASTA. "
                "An existing valid cache can be reused; otherwise choose mirbase (all mature "
                "sequences), mirgenedb, or supply a reference FASTA explicitly."
            )
            return None
        for attempt in range(self.download_attempts):
            content = self._download_file(source, self.download_timeout)
            if content:
                return content
            if attempt + 1 < self.download_attempts:
                time.sleep(self.retry_backoff_seconds * (2**attempt))
        logger.error(
            "miRNA source %s/%s unavailable at %s after %s attempts (socket timeout %ss); "
            "supply a verified pinned mature FASTA through a reference manifest",
            source_name,
            source.species,
            source.url,
            self.download_attempts,
            self.download_timeout,
        )
        return None

    def get_database(self, source_name: str, species: str, force_refresh: bool = False) -> Path | None:
        """Resolve a single source or the explicit mature-only union.

        Args:
            source_name: Registered source name or mirgenedb+mirbase.
            species: Target organism identifier.
            force_refresh: Refresh both component references when combining.

        Returns:
            Validated reference path, or None if any required source failed.
        """
        normalized = self.normalize_species(source_name, species)
        self._resolved_reference_metadata.pop((source_name, normalized or species), None)
        if source_name != "mirgenedb+mirbase":
            return self._get_single_database(source_name, species, force_refresh)
        normalized = self.normalize_species(source_name, species)
        if normalized not in self.SOURCES[source_name]:
            logger.error("Unknown combined reference species '%s'", species)
            return None
        sources = ["mirgenedb", "mirbase"]
        if force_refresh:
            for component in sources:
                if self.get_database(component, normalized, force_refresh=True) is None:
                    return None
        return self.get_combined_database(sources, normalized)

    def _get_single_database(self, source_name: str, species: str, force_refresh: bool = False) -> Path | None:  # noqa: PLR0911, PLR0912
        """Get miRNA database, downloading and filtering if needed.

        Simplified caching: each species+source combination gets its own cache file.

        Args:
            source_name: Database source ("mirbase", "mirbase_high_conf", etc.)
            species: Species name ("human", "mouse", "rat")
            force_refresh: Force re-download even if cached

        Returns:
            Path to cached FASTA file, or None if failed
        """
        normalized_species = self.normalize_species(source_name, species)
        if normalized_species is None:
            logger.error("Unknown species '%s' for source '%s'", species, source_name)
            return None

        if source_name not in self.SOURCES or normalized_species not in self.SOURCES[source_name]:
            logger.error(f"Unknown source/species combination: {source_name}/{species}")
            return None

        source = self.SOURCES[source_name][normalized_species]
        pin = self.pinned_references.get(source_name, {}).get(normalized_species)
        if pin is not None:
            return self._get_pinned_database(source_name, normalized_species, pin, force_refresh)
        cache_key = source.cache_key()
        cache_file = self.cache_dir / f"{cache_key}.fa"

        # Check if we can use cached version
        if not force_refresh and self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached {source.name} ({source.species}): {cache_file}")
            meta = self.metadata[cache_key]
            self._resolved_reference_metadata[source_name, normalized_species] = {
                **(meta.extra or {}),
                "source": source_name,
                "release": (meta.extra or {}).get("release"),
                "original_url": source.url,
                "acquisition_url": source.url,
                "filtered_sha256": hashlib.sha256(cache_file.read_bytes()).hexdigest(),
                "acquired_at": meta.downloaded_at,
                "species": normalized_species,
                "acquisition_mode": "cached_reuse",
            }
            return cache_file

        def cleanup_cache() -> None:
            if cache_file.exists():
                cache_file.unlink(missing_ok=True)
            if cache_key in self.metadata:
                del self.metadata[cache_key]
                self._save_metadata()

        # Download and process
        logger.info(f"🔄 Downloading {source.name} ({source.species})...")

        content = self._download_available_source(source_name, source)
        failure = False

        if content is None:
            failure = True
        elif not content:
            logger.error("Downloaded content for %s/%s is empty", source_name, normalized_species)
            failure = True
        else:
            if source_name.startswith("mirbase"):
                logger.info(f"🔄 Filtering for {normalized_species}...")
                content = self._filter_species_sequences(content, normalized_species)
                if not content.strip():
                    logger.error(
                        "Filtered miRBase content for %s/%s produced no sequences; discarding cache update",
                        source_name,
                        normalized_species,
                    )
                    failure = True
            if not failure and not content.strip():
                logger.error(
                    "Downloaded content for %s/%s is empty after normalization",
                    source_name,
                    normalized_species,
                )
                failure = True

        if failure:
            # Failed acquisition is not permission to destroy a previous
            # reference. Return failure explicitly; no implicit stale fallback.
            return None

        assert content is not None
        try:
            self._validate_mature_fasta(content)
        except ValueError as exc:
            logger.error("Invalid miRNA reference %s/%s: %s", source_name, normalized_species, exc)
            return None

        # Save to cache
        try:
            _atomic_write_text(cache_file, content)
        except OSError as exc:
            logger.error("Cannot install miRNA reference %s/%s: %s", source_name, normalized_species, exc)
            return None

        if cache_file.stat().st_size == 0:
            logger.error(
                "Downloaded content for %s/%s is empty; removing cache file",
                source_name,
                normalized_species,
            )
            cleanup_cache()
            return None

        # Update metadata
        provenance = {
            "source": source_name,
            "release": None,
            "original_url": source.url,
            "acquisition_url": source.url,
            "filtered_sha256": hashlib.sha256(cache_file.read_bytes()).hexdigest(),
            "acquired_at": self._get_current_timestamp(),
            "species": normalized_species,
            "species_filtered": source_name.startswith("mirbase"),
            "acquisition_mode": "successful_refresh",
        }
        self._record_cache_entry(cache_key, source, cache_file, extra=provenance, persist=True)
        self._resolved_reference_metadata[source_name, normalized_species] = provenance

        logger.info(f"✅ Cached {source.name} ({source.species}): {cache_file} ({cache_file.stat().st_size:,} bytes)")
        return cache_file

    def _get_current_timestamp(self) -> str:
        """Get current timestamp in ISO format."""
        return datetime.now().isoformat()

    def _canonical_species_for_sources(self, sources: list[str], species: str) -> str:
        for source_name in sources:
            normalized = self.normalize_species(source_name, species)
            if normalized:
                return normalized
        return species

    def get_combined_database(self, sources: list[str], species: str, output_name: str | None = None) -> Path | None:
        """Combine multiple databases into a single file.

        Args:
            sources: List of source names to combine
            species: Target species
            output_name: Custom output filename (default: auto-generated)

        Returns:
            Path to combined FASTA file
        """
        canonical_species = self._canonical_species_for_sources(sources, species)
        combined_species = self.normalize_species("mirgenedb+mirbase", species) or species
        self._resolved_reference_metadata.pop(("mirgenedb+mirbase", combined_species), None)

        if output_name is None:
            output_name = f"combined_{canonical_species}_{'_'.join(sources)}.fa"

        combined_file = self.cache_dir / output_name
        source_files = []
        for source_name in sources:
            source_file = self.get_database(source_name, species)
            if source_file is None:
                logger.error(f"Failed to get {source_name} database")
                return None
            source_files.append(source_file)

        combined_provenance = {
            "source": "mirgenedb+mirbase",
            "species": combined_species,
            "acquisition_mode": "combined",
            "sources": [self.get_reference_metadata(source_name, species) for source_name in sources],
        }

        # The combined file is derived, so it inherits every defect of the FASTAs it was built
        # from -- and, being an ordinary file, it can also be truncated or corrupted *after* it was
        # written. Reuse is therefore gated on a stamp recorded at write time that pins three
        # things: the producer version, the exact input bytes it was combined from, and the
        # output's own size and digest, re-read from disk and subject to the cache TTL. The
        # mtime-only predecessor ("combined is newer than every source") had neither a TTL nor a
        # checksum of any kind, so a wrong, stale or half-written database was served forever.
        # Rebuilding unconditionally would also be correct but is the wrong default: it re-reads
        # and re-writes every source on every call for no information gained.
        input_digests = fingerprint_inputs(source_files)
        if combined_file.exists() and is_artifact_stamp_current(
            ARTIFACT_MIRNA_COMBINED,
            combined_file,
            inputs=input_digests,
            max_age_days=self.cache_ttl.days,
        ):
            logger.info(f"✅ Using existing combined database: {combined_file}")
            self._resolved_reference_metadata["mirgenedb+mirbase", combined_species] = {
                **combined_provenance,
                "filtered_sha256": hashlib.sha256(combined_file.read_bytes()).hexdigest(),
            }
            return combined_file

        logger.info(f"🔄 Combining {len(sources)} databases for {canonical_species}...")

        seen_sequences = set()
        total_sequences = 0

        with combined_file.open("w") as outfile:
            for source_name, source_file in zip(sources, source_files, strict=True):
                for record in SeqIO.parse(source_file, FASTA_PARSE_FORMAT):
                    sequence = str(record.seq).upper().replace("U", "T")
                    if sequence in seen_sequences:
                        continue
                    seen_sequences.add(sequence)

                    annotated_record = SeqRecord(
                        record.seq,
                        id=record.id,
                        name=record.name,
                        description=f"{record.description} [source:{source_name}]",
                    )
                    SeqIO.write(annotated_record, outfile, "fasta")
                    total_sequences += 1

        # Stamp the result *after* the bytes are on disk, so the recorded size/digest describe
        # what a later run will actually read back.
        write_artifact_stamp(
            ARTIFACT_MIRNA_COMBINED,
            combined_file,
            inputs=input_digests,
            extra={"sources": sources, "species": canonical_species, "sequences": total_sequences},
        )

        logger.info(f"✅ Combined database created: {combined_file} ({total_sequences} unique sequences)")
        self._resolved_reference_metadata["mirgenedb+mirbase", combined_species] = {
            **combined_provenance,
            "filtered_sha256": hashlib.sha256(combined_file.read_bytes()).hexdigest(),
        }
        return combined_file

    def list_available_databases(self) -> dict[str, dict[str, MiRNASource]]:
        """List all available database sources and species."""
        return self.SOURCES

    def clear_cache(self, confirm: bool = False) -> dict[str, Any]:
        """Clear miRNA cache using the shared reference-manager implementation."""
        return super().clear_cache(confirm=confirm)


def load_pinned_references(manifest_path: str | Path) -> dict[str, dict[str, PinnedMiRNAReference]]:
    """Read a version 1 local reference manifest, resolving paths relative to it.

    ``references`` is a nonempty list of source, species, release, sha256, path
    entries, optionally carrying original_url. This validates configuration only;
    acquisition verifies the actual FASTA bytes and never silently substitutes a source.
    """
    path = Path(manifest_path).resolve()
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ValueError(f"Cannot read miRNA reference manifest {path}: {exc}") from exc
    if (
        not isinstance(manifest, dict)
        or type(manifest.get("schema_version")) is not int
        or manifest["schema_version"] != 1
    ):
        raise ValueError("miRNA reference manifest requires schema_version: 1")
    entries = manifest.get("references")
    if not isinstance(entries, list) or not entries:
        raise ValueError("miRNA reference manifest requires a nonempty references list")
    result: dict[str, dict[str, PinnedMiRNAReference]] = {}
    required = {"source", "species", "release", "sha256", "path"}
    for entry in entries:
        if not isinstance(entry, dict) or not required.issubset(entry):
            raise ValueError("Each miRNA reference requires source, species, release, sha256 and path")
        if set(entry) - required - {"original_url"}:
            raise ValueError("Unknown miRNA reference manifest entry fields")
        if not all(isinstance(entry[key], str) and entry[key].strip() for key in required):
            raise ValueError("miRNA reference manifest fields must be nonempty strings")
        source = entry["source"]
        if source not in {"mirbase", "mirgenedb"}:
            raise ValueError(f"Unsupported pinned mature source: {source}")
        species = MiRNADatabaseManager.normalize_species(source, entry["species"])
        if species is None or not MiRNADatabaseManager.is_supported_species(source, species):
            raise ValueError(f"Unsupported pinned miRNA source/species: {source}/{entry['species']}")
        if species in result.setdefault(source, {}):
            raise ValueError(f"Duplicate pinned miRNA source/species: {source}/{species}")
        local_path = Path(entry["path"])
        if not local_path.is_absolute():
            local_path = path.parent / local_path
        result[source][species] = PinnedMiRNAReference(
            path=local_path.resolve(),
            release=entry["release"],
            sha256=entry["sha256"],
            original_url=entry.get("original_url"),
        )
    return result


def _create_parser() -> argparse.ArgumentParser:
    """Create and configure argument parser."""
    sources = MiRNADatabaseManager.get_available_sources()
    species_choices = MiRNADatabaseManager.get_all_species()

    parser = argparse.ArgumentParser(description="miRNA Database Manager")
    parser.add_argument(
        "--source",
        choices=sources,
        help="Database source",
    )
    parser.add_argument("--species", choices=species_choices, help="Target species")
    parser.add_argument(
        "--combine",
        nargs="+",
        choices=sources,
        help="Combine multiple sources",
    )
    parser.add_argument("--list", action="store_true", help="List available databases")
    parser.add_argument("--clean", action="store_true", help="Clean old cache files")
    parser.add_argument("--clear-cache", action="store_true", help="Clear all cache files")
    parser.add_argument(
        "--clear-cache-dry-run", action="store_true", help="Show what would be deleted without actually deleting"
    )
    parser.add_argument("--info", action="store_true", help="Show cache information")
    parser.add_argument("--force", action="store_true", help="Force refresh cached files")
    return parser


def main() -> None:  # noqa: PLR0912
    """CLI interface for the miRNA database manager."""
    args = _create_parser().parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    manager = MiRNADatabaseManager()

    if args.list:
        print("📋 Available miRNA databases:")
        for source_name, species_dict in manager.list_available_databases().items():
            print(f"\n🧬 {source_name}:")
            for species, source in species_dict.items():
                print(f"  • {species}: {source.description}")

    elif args.info:
        info = manager.cache_info()
        print("📊 Cache Information:")
        print(f"  Directory: {info['cache_directory']}")
        print(f"  Files: {info['total_files']}")
        print(f"  Size: {info['total_size_mb']:.1f} MB")
        print(f"  TTL: {info['cache_ttl_days']} days")

    elif args.clean:
        manager.clean_cache()

    elif args.clear_cache_dry_run:
        result = manager.clear_cache(confirm=False)
        print("🔍 Cache Clear Preview:")
        print(f"  Directory: {result['cache_directory']}")
        print(f"  Files to delete: {result['files_deleted']}")
        print(f"  Size to free: {result['size_freed_mb']:.2f} MB")
        print(f"  Status: {result['status']}")

    elif args.clear_cache:
        result = manager.clear_cache(confirm=True)
        print("🧹 Cache Cleared:")
        print(f"  Directory: {result['cache_directory']}")
        print(f"  Files deleted: {result['files_deleted']}")
        print(f"  Size freed: {result['size_freed_mb']:.2f} MB")
        print(f"  Status: {result['status']}")

    elif args.combine and args.species:
        output_file = manager.get_combined_database(args.combine, args.species)
        if output_file:
            print(f"✅ Combined database: {output_file}")
        else:
            print("❌ Failed to create combined database")

    elif args.source and args.species:
        output_file = manager.get_database(args.source, args.species, force_refresh=args.force)
        if output_file:
            print(f"✅ Database ready: {output_file}")
        else:
            print("❌ Failed to get database")

    else:
        _create_parser().print_help()


if __name__ == "__main__":
    main()
