"""Shared CLI input parsing and validation helpers.

These helpers are intentionally small and side-effect free so that Typer/Rich
commands can stay focused on orchestration and UX.
"""

from __future__ import annotations

from dataclasses import dataclass

from sirnaforge.data.mirna_manager import MiRNADatabaseManager


def parse_csv(value: str) -> list[str]:
    """Split a comma-separated string into normalized non-empty tokens."""
    return [token.strip() for token in value.split(",") if token.strip()]


def parse_required_csv(value: str, *, error_message: str) -> list[str]:
    """Parse a required CSV argument, raising ValueError when empty."""
    tokens = parse_csv(value)
    if not tokens:
        raise ValueError(error_message)
    return tokens


def parse_optional_csv(value: str | None, *, error_message: str) -> list[str] | None:
    """Parse an optional CSV argument.

    Returns None when value is None, otherwise returns a list of tokens.
    Raises ValueError when the provided string contains no usable tokens.
    """
    if value is None:
        return None
    tokens = parse_csv(value)
    if not tokens:
        raise ValueError(error_message)
    return tokens


def extract_override_species_from_offtarget_indices(offtarget_indices: str | None) -> list[str] | None:
    """Extract unique species tokens from an offtarget indices override string.

    Validates the expected ``species:/index_prefix`` format for each entry.
    Returns None when no override is provided.
    """
    if not offtarget_indices:
        return None

    entries = parse_csv(offtarget_indices)
    bad_entries = [entry for entry in entries if ":" not in entry]
    if bad_entries:
        raise ValueError("--offtarget-indices entries must be in species:/index_prefix form")

    override_species: list[str] = []
    for entry in entries:
        species_token = entry.split(":", 1)[0].strip() or entry
        if species_token and species_token not in override_species:
            override_species.append(species_token)
    return override_species


@dataclass(frozen=True)
class SpeciesResolution:
    """Resolved species identifiers for genome, canonical, and miRNA scopes."""

    source_normalized: str
    canonical_species: list[str]
    genome_species: list[str]
    mirna_species: list[str]


def resolve_species_inputs(*, species: str, mirna_db: str, mirna_species: str | None) -> SpeciesResolution:
    """Validate species + miRNA arguments and return resolved identifiers.

    Raises ValueError with a user-facing message when inputs are invalid.
    """
    requested_species = parse_required_csv(species, error_message="at least one species must be provided")
    source_normalized = mirna_db.lower()

    if not MiRNADatabaseManager.is_supported_source(source_normalized):
        valid_sources = ", ".join(MiRNADatabaseManager.get_available_sources())
        raise ValueError(f"unknown miRNA database '{mirna_db}'. Supported sources: {valid_sources}")

    mirna_overrides = parse_optional_csv(
        mirna_species,
        error_message="--mirna-species override must contain at least one value",
    )

    try:
        species_resolution = MiRNADatabaseManager.resolve_species_selection(
            requested_species,
            source_normalized,
            mirna_overrides=mirna_overrides,
        )
    except ValueError as exc:
        supported = MiRNADatabaseManager.get_supported_canonical_species_for_source(source_normalized)
        raise ValueError(f"{exc}. Supported canonical species: {', '.join(supported)}") from exc

    return SpeciesResolution(
        source_normalized=source_normalized,
        canonical_species=species_resolution["canonical"],
        genome_species=species_resolution["genome"],
        mirna_species=species_resolution["mirna"],
    )
