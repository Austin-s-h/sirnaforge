#!/usr/bin/env python3
"""Elegant miRNA Database Manager with local caching and multi-species support.

This module provides a clean interface for downloading, caching, and managing
miRNA databases from multiple sources with automatic cache management and
species-specific organization.
"""

import argparse
import contextlib
import gzip
import hashlib
import html
import json
import logging
import os
import tempfile
import urllib.request
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Optional, Union

from .species_registry import (
    CANONICAL_SPECIES_ALIAS_MAP,
    CANONICAL_SPECIES_REGISTRY,
    MIRGENEDB_ALIAS_MAP,
    MIRGENEDB_SPECIES_TABLE,
)

logger = logging.getLogger(__name__)


@dataclass
class DatabaseSource:
    """Configuration for a miRNA database source."""

    name: str
    url: str
    species: str  # NCBI taxonomy ID or common name
    format: str  # "fasta", "json", "tsv"
    compressed: bool = False
    description: str = ""

    def cache_key(self) -> str:
        """Generate a unique cache key for this source."""
        content = f"{self.name}_{self.species}_{self.url}"
        return hashlib.md5(content.encode()).hexdigest()[:12]


@dataclass
class CacheMetadata:
    """Metadata for cached database files."""

    source: DatabaseSource
    downloaded_at: str
    file_size: int
    checksum: str
    version: str = "1.0"

    @classmethod
    def from_dict(cls, data: dict) -> "CacheMetadata":
        """Create CacheMetadata from dictionary."""
        source = DatabaseSource(**data["source"])
        return cls(
            source=source,
            downloaded_at=data["downloaded_at"],
            file_size=data["file_size"],
            checksum=data["checksum"],
            version=data.get("version", "1.0"),
        )

    def to_dict(self) -> dict:
        """Convert CacheMetadata to dictionary."""
        return asdict(self)


def _build_mirgenedb_sources() -> dict[str, DatabaseSource]:
    """Construct MirGeneDB source map from the slug-based metadata table."""
    sources: dict[str, DatabaseSource] = {}
    base_url = "https://www.mirgenedb.org/fasta/{slug}?mat=1"

    for slug, metadata in MIRGENEDB_SPECIES_TABLE.items():
        scientific_name = metadata.get("scientific_name", slug)
        taxonomy_id = metadata.get("taxonomy_id")
        common_name = metadata.get("common_name")
        description = f"MirGeneDB high-confidence miRNAs ({scientific_name}, NCBI:{taxonomy_id})"
        if common_name and common_name.lower() not in scientific_name.lower():
            description = f"{description} [{common_name}]"

        sources[slug] = DatabaseSource(
            name="mirgenedb",
            url=base_url.format(slug=slug),
            species=slug,
            format="fasta",
            compressed=False,
            description=description,
        )

    return sources


class MiRNADatabaseManager:
    """Elegant miRNA database manager with caching and multi-species support."""

    # Database source configurations
    SOURCES = {
        "mirbase": {
            "human": DatabaseSource(
                name="mirbase_mature",
                url="https://www.mirbase.org/download/CURRENT/mature.fa",
                species="human",
                format="fasta",
                compressed=False,
                description="miRBase mature miRNA sequences (all species, filtered for Homo sapiens - hsa)",
            ),
            "mouse": DatabaseSource(
                name="mirbase_mature",
                url="https://www.mirbase.org/download/CURRENT/mature.fa",
                species="mouse",
                format="fasta",
                compressed=False,
                description="miRBase mature miRNA sequences (all species, filtered for Mus musculus - mmu)",
            ),
            "rat": DatabaseSource(
                name="mirbase_mature",
                url="https://www.mirbase.org/download/CURRENT/mature.fa",
                species="rat",
                format="fasta",
                compressed=False,
                description="miRBase mature miRNA sequences (all species, filtered for Rattus norvegicus - rno)",
            ),
        },
        "mirbase_high_conf": {
            "human": DatabaseSource(
                name="mirbase_mature_hc",
                url="https://www.mirbase.org/download/CURRENT/mature_high_conf.fa",
                species="human",
                format="fasta",
                compressed=False,
                description="miRBase high-confidence mature miRNA sequences (Homo sapiens - hsa)",
            ),
            "mouse": DatabaseSource(
                name="mirbase_mature_hc",
                url="https://www.mirbase.org/download/CURRENT/mature_high_conf.fa",
                species="mouse",
                format="fasta",
                compressed=False,
                description="miRBase high-confidence mature miRNA sequences (Mus musculus - mmu)",
            ),
            "rat": DatabaseSource(
                name="mirbase_mature_hc",
                url="https://www.mirbase.org/download/CURRENT/mature_high_conf.fa",
                species="rat",
                format="fasta",
                compressed=False,
                description="miRBase high-confidence mature miRNA sequences (Rattus norvegicus - rno)",
            ),
        },
        "mirbase_hairpin": {
            "human": DatabaseSource(
                name="mirbase_hairpin",
                url="https://www.mirbase.org/download/CURRENT/hairpin.fa",
                species="human",
                format="fasta",
                compressed=False,
                description="miRBase hairpin precursor miRNA sequences (Homo sapiens - hsa)",
            ),
            "mouse": DatabaseSource(
                name="mirbase_hairpin",
                url="https://www.mirbase.org/download/CURRENT/hairpin.fa",
                species="mouse",
                format="fasta",
                compressed=False,
                description="miRBase hairpin precursor miRNA sequences (Mus musculus - mmu)",
            ),
            "rat": DatabaseSource(
                name="mirbase_hairpin",
                url="https://www.mirbase.org/download/CURRENT/hairpin.fa",
                species="rat",
                format="fasta",
                compressed=False,
                description="miRBase hairpin precursor miRNA sequences (Rattus norvegicus - rno)",
            ),
        },
        "mirgenedb": _build_mirgenedb_sources(),
        "targetscan": {
            "human": DatabaseSource(
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
    def canonicalize_species_name(cls, species: str) -> Optional[str]:
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
        """Return genome species identifiers for canonical species keys."""
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
        mirna_overrides: Optional[Sequence[str]] = None,
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
            "genome": genome_species,
            "mirna": mirna_species,
        }

    @classmethod
    def _normalize_mirgenedb_species(cls, species: str) -> Optional[str]:
        if not species:
            return None
        return MIRGENEDB_ALIAS_MAP.get(species.lower())

    @classmethod
    def normalize_species(cls, source_name: str, species: str) -> Optional[str]:
        """Normalize user-provided species identifiers to canonical keys."""
        if source_name == "mirgenedb":
            return cls._normalize_mirgenedb_species(species)
        return species

    @classmethod
    def get_source_configuration(cls, source_name: str, species: str) -> Optional[DatabaseSource]:
        """Retrieve the DatabaseSource configuration for a given source/species."""
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

    def __init__(self, cache_dir: Optional[Union[str, Path]] = None, cache_ttl_days: int = 30):
        """Initialize the miRNA database manager.

        Args:
            cache_dir: Directory for caching databases (default: ~/.cache/sirnaforge/mirna)
            cache_ttl_days: Cache time-to-live in days
        """
        fallback_locations: list[Path] = []
        env_override = os.getenv("SIRNAFORGE_CACHE_DIR")

        if cache_dir is not None:
            fallback_locations.append(Path(cache_dir))
        else:
            if env_override:
                fallback_locations.append(Path(env_override))

            xdg_cache = os.getenv("XDG_CACHE_HOME")
            if xdg_cache:
                fallback_locations.append(Path(xdg_cache) / "sirnaforge" / "mirna")

            fallback_locations.append(Path.home() / ".cache" / "sirnaforge" / "mirna")
            fallback_locations.append(Path.cwd() / ".sirnaforge_cache" / "mirna")
            fallback_locations.append(Path(tempfile.gettempdir()) / "sirnaforge" / "mirna")

        cache_path: Optional[Path] = None
        first_error: Optional[Exception] = None
        for candidate in fallback_locations:
            try:
                candidate.mkdir(parents=True, exist_ok=True)
                cache_path = candidate
                break
            except PermissionError as exc:
                logger.warning("Cannot create cache directory at %s; trying next option", candidate)
                if first_error is None:
                    first_error = exc

        if cache_path is None:
            if cache_dir is not None and first_error is not None:
                raise first_error
            raise RuntimeError("Cannot create a writable cache directory for miRNA databases")

        self.cache_dir = cache_path

        self.cache_ttl = timedelta(days=cache_ttl_days)
        self.metadata_file = self.cache_dir / "cache_metadata.json"

        self._load_metadata()

    def _load_metadata(self) -> None:
        """Load cache metadata from disk."""
        self.metadata: dict[str, CacheMetadata] = {}

        if self.metadata_file.exists():
            try:
                with self.metadata_file.open("r") as f:
                    data = json.load(f)
                    for key, meta_dict in data.items():
                        self.metadata[key] = CacheMetadata.from_dict(meta_dict)
            except (json.JSONDecodeError, KeyError) as e:
                logger.warning(f"Failed to load cache metadata: {e}")

    def _save_metadata(self) -> None:
        """Save cache metadata to disk."""
        try:
            data = {key: meta.to_dict() for key, meta in self.metadata.items()}
            with self.metadata_file.open("w") as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save cache metadata: {e}")

    def _compute_file_checksum(self, file_path: Path) -> str:
        """Compute MD5 checksum of a file."""
        hash_md5 = hashlib.md5()
        with file_path.open("rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def _is_cache_valid(self, cache_key: str) -> bool:
        """Check if cached data is still valid."""
        if cache_key not in self.metadata:
            return False

        meta = self.metadata[cache_key]
        cache_file = self.cache_dir / f"{cache_key}.fa"

        # Check if file exists
        if not cache_file.exists():
            return False

        # Reject zero-byte cache entries up front
        if cache_file.stat().st_size == 0:
            logger.warning("Cache file %s is empty; marking as invalid", cache_file)
            return False

        # Check TTL
        downloaded_at = datetime.fromisoformat(meta.downloaded_at)
        if datetime.now() - downloaded_at > self.cache_ttl:
            return False

        # Check file integrity
        if self._compute_file_checksum(cache_file) != meta.checksum:
            logger.warning(f"Cache file {cache_file} corrupted, will re-download")
            return False

        return True

    def _download_file(self, source: DatabaseSource) -> Optional[str]:
        """Download file from source URL and return as text."""
        try:
            logger.info(f"📥 Downloading {source.name} ({source.species}): {source.url}")

            request = urllib.request.Request(
                source.url,
                headers={
                    "User-Agent": "sirnaforge/1.0 (+https://github.com/austin-s-h/sirnaforge)",
                    "Accept": "text/plain,application/octet-stream",
                },
            )

            with urllib.request.urlopen(request, timeout=300) as response:
                data = response.read()

            if source.compressed and source.url.endswith(".gz"):
                data = gzip.decompress(data)
                # Note: Could add zip support here if needed

            # Try to decode as text
            try:
                content: str = data.decode("utf-8")
            except UnicodeDecodeError:
                try:
                    content = data.decode("latin-1")
                except UnicodeDecodeError:
                    logger.error(f"❌ Cannot decode {source.url} as text")
                    return None

            # Fix HTML entities (e.g., &gt; -> >, <br> -> newlines)
            content = html.unescape(content)
            content = content.replace("<br>", "\n").replace("<BR>", "\n")

            if not content.strip():
                logger.error("Received empty response from %s", source.url)
                return None

            logger.info(f"✅ Downloaded {len(content):,} characters")
            return content

        except Exception as e:
            logger.error(f"❌ Failed to download {source.url}: {e}")
            return None

    def _filter_species_sequences(self, fasta_content: str, species: str) -> str:
        """Filter FASTA content for specific species using miRBase three-letter codes."""
        # Species prefix mapping - simplified and clear
        species_codes = {
            "human": "hsa-",
            "mouse": "mmu-",
            "rat": "rno-",
            "zebrafish": "dre-",
            "fly": "dme-",
            "worm": "cel-",
            "chicken": "gga-",
            "dog": "cfa-",
            "pig": "ssc-",
        }

        code = species_codes.get(species)
        if not code:
            logger.warning(f"Unknown species '{species}', returning all sequences")
            return fasta_content

        # Simple filtering: include header+sequence pairs where header contains species code
        filtered_lines = []
        current_header = None

        for raw_line in fasta_content.split("\n"):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # New header - check if it matches our species
                if code in line:
                    current_header = line
                    filtered_lines.append(line)
                else:
                    current_header = None
            elif current_header and line:
                # Sequence line for a matching header
                filtered_lines.append(line)

        filtered_count = len([line for line in filtered_lines if line.startswith(">")])
        if filtered_count == 0:
            logger.error("No sequences found for species '%s' after filtering", species)

        logger.info(f"Filtered to {filtered_count} {species} sequences")
        return "\n".join(filtered_lines)

    def get_database(self, source_name: str, species: str, force_refresh: bool = False) -> Optional[Path]:
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
        cache_key = source.cache_key()
        cache_file = self.cache_dir / f"{cache_key}.fa"

        # Check if we can use cached version
        if not force_refresh and self._is_cache_valid(cache_key):
            logger.info(f"✅ Using cached {source.name} ({source.species}): {cache_file}")
            return cache_file

        def cleanup_cache() -> None:
            if cache_file.exists():
                cache_file.unlink(missing_ok=True)
            if cache_key in self.metadata:
                del self.metadata[cache_key]
                self._save_metadata()

        # Download and process
        logger.info(f"🔄 Downloading {source.name} ({source.species})...")

        content: Optional[str] = self._download_file(source)
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
            cleanup_cache()
            return None

        assert content is not None

        # Save to cache
        with cache_file.open("w", encoding="utf-8") as f:
            f.write(content)

        if cache_file.stat().st_size == 0:
            logger.error(
                "Downloaded content for %s/%s is empty; removing cache file",
                source_name,
                normalized_species,
            )
            cleanup_cache()
            return None

        # Update metadata
        checksum = self._compute_file_checksum(cache_file)
        self.metadata[cache_key] = CacheMetadata(
            source=source,
            downloaded_at=datetime.now().isoformat(),
            file_size=cache_file.stat().st_size,
            checksum=checksum,
        )
        self._save_metadata()

        logger.info(f"✅ Cached {source.name} ({source.species}): {cache_file} ({cache_file.stat().st_size:,} bytes)")
        return cache_file

    def _canonical_species_for_sources(self, sources: list[str], species: str) -> str:
        for source_name in sources:
            normalized = self.normalize_species(source_name, species)
            if normalized:
                return normalized
        return species

    def get_combined_database(
        self, sources: list[str], species: str, output_name: Optional[str] = None
    ) -> Optional[Path]:
        """Combine multiple databases into a single file.

        Args:
            sources: List of source names to combine
            species: Target species
            output_name: Custom output filename (default: auto-generated)

        Returns:
            Path to combined FASTA file
        """
        canonical_species = self._canonical_species_for_sources(sources, species)

        if output_name is None:
            output_name = f"combined_{canonical_species}_{'_'.join(sources)}.fa"

        combined_file = self.cache_dir / output_name

        # Check if we need to regenerate
        source_files = []
        for source_name in sources:
            source_file = self.get_database(source_name, species)
            if source_file is None:
                logger.error(f"Failed to get {source_name} database")
                return None
            source_files.append(source_file)

        # Check if combined file is newer than all sources
        if combined_file.exists():
            combined_mtime = combined_file.stat().st_mtime
            if all(source_file.stat().st_mtime <= combined_mtime for source_file in source_files):
                logger.info(f"✅ Using existing combined database: {combined_file}")
                return combined_file

            # Combine databases
            logger.info(f"🔄 Combining {len(sources)} databases for {canonical_species}...")

        seen_sequences = set()
        total_sequences = 0

        with combined_file.open("w") as outfile:
            for source_file in source_files:
                source_name = "unknown"
                for src_name in sources:
                    if src_name in source_file.name:
                        source_name = src_name
                        break

                with source_file.open("r") as infile:
                    header = None
                    for line in infile:
                        line_content = line.strip()
                        if line_content.startswith(">"):
                            header = f"{line_content} [source:{source_name or 'unknown'}]"
                        elif line_content and header:
                            seq_upper = line_content.upper()
                            if seq_upper not in seen_sequences:
                                outfile.write(f"{header}\n{line_content}\n")
                                seen_sequences.add(seq_upper)
                                total_sequences += 1
                            header = None

        logger.info(f"✅ Combined database created: {combined_file} ({total_sequences} unique sequences)")
        return combined_file

    def list_available_databases(self) -> dict[str, dict[str, DatabaseSource]]:
        """List all available database sources and species."""
        return self.SOURCES

    def clean_cache(self, older_than_days: Optional[int] = None) -> None:
        """Clean old cache files (both raw and filtered).

        Args:
            older_than_days: Remove files older than this (default: use TTL)
        """
        if older_than_days is None:
            older_than_days = self.cache_ttl.days

        cutoff = datetime.now() - timedelta(days=older_than_days)
        removed_count = 0

        for cache_key in list(self.metadata.keys()):
            meta = self.metadata[cache_key]
            downloaded_at = datetime.fromisoformat(meta.downloaded_at)

            if downloaded_at < cutoff:
                cache_file = self.cache_dir / f"{cache_key}.fa"
                if cache_file.exists():
                    cache_file.unlink()
                    removed_count += 1

                del self.metadata[cache_key]

        if removed_count > 0:
            self._save_metadata()
            logger.info(f"🧹 Cleaned {removed_count} old cache files")
        else:
            logger.info("🧹 No old cache files to clean")

    def cache_info(self) -> dict[str, Any]:
        """Get information about the current cache state."""
        total_files = len(list(self.cache_dir.glob("*.fa")))
        total_size = sum(f.stat().st_size for f in self.cache_dir.glob("*.fa"))

        return {
            "cache_directory": str(self.cache_dir),
            "total_files": total_files,
            "total_size_mb": total_size / (1024 * 1024),
            "cache_ttl_days": self.cache_ttl.days,
            "cached_databases": list(self.metadata.keys()),
        }

    def clear_cache(self, confirm: bool = False) -> dict[str, Any]:
        """Clear the miRNA cache directory.

        Args:
            confirm: If True, actually delete files. If False, just return what would be deleted.

        Returns:
            Dictionary with information about files deleted or that would be deleted.
        """
        if not self.cache_dir.exists():
            return {
                "cache_directory": str(self.cache_dir),
                "files_deleted": 0,
                "size_freed_mb": 0.0,
                "status": "Cache directory does not exist",
            }

        cache_files = list(self.cache_dir.glob("*.fa"))
        json_files = list(self.cache_dir.glob("*.json"))
        all_files = cache_files + json_files

        total_size = sum(f.stat().st_size for f in all_files if f.exists())

        result = {
            "cache_directory": str(self.cache_dir),
            "files_deleted": len(all_files),
            "size_freed_mb": total_size / (1024 * 1024),
            "status": "Would delete" if not confirm else "Deleted",
        }

        if confirm:
            # Actually delete the files
            for file_path in all_files:
                with contextlib.suppress(FileNotFoundError):
                    file_path.unlink()
            result["status"] = "Cache cleared successfully"
        else:
            result["status"] = f"Would delete {len(all_files)} files ({total_size / (1024 * 1024):.2f} MB)"

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
