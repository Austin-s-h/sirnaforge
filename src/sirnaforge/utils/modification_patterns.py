"""Utility functions for applying chemical modification patterns to siRNA candidates.

This module provides functions to apply standard modification patterns to siRNA
candidates during the design workflow, enabling automated annotation of chemical
modifications for downstream synthesis and analysis.
"""

from typing import TYPE_CHECKING, Optional

from sirnaforge.models.modifications import (
    ChemicalModification,
    ConfirmationStatus,
    Provenance,
    SourceType,
    StrandMetadata,
)

if TYPE_CHECKING:
    from sirnaforge.models.sirna import SiRNACandidate


def apply_standard_2ome_pattern(sequence: str) -> list[ChemicalModification]:
    """Apply standard alternating 2'-O-methyl pattern.

    This is the industry-standard pattern providing balanced nuclease
    resistance and RISC loading efficiency.

    Args:
        sequence: RNA sequence to modify

    Returns:
        List containing one ChemicalModification with alternating positions
    """
    # Alternating positions (1, 3, 5, 7, ...)
    positions = [i for i in range(1, len(sequence) + 1) if i % 2 == 1]
    return [ChemicalModification(type="2OMe", positions=positions)]


def apply_minimal_terminal_pattern(sequence: str) -> list[ChemicalModification]:
    """Apply minimal terminal modifications for cost-effective protection.

    Modifies only the 3' terminal positions to provide basic nuclease
    resistance while minimizing synthesis cost.

    Args:
        sequence: RNA sequence to modify

    Returns:
        List containing one ChemicalModification with terminal positions
    """
    seq_len = len(sequence)
    # Last 3 positions for 3' terminal protection
    return [ChemicalModification(type="2OMe", positions=[seq_len - 2, seq_len - 1, seq_len])]


def apply_maximal_stability_pattern(sequence: str) -> list[ChemicalModification]:
    """Apply maximal stability pattern for in vivo applications.

    Fully modified pattern similar to FDA-approved therapeutics, providing
    maximum nuclease resistance and extended serum half-life.

    Args:
        sequence: RNA sequence to modify

    Returns:
        List containing ChemicalModifications (2OMe on all positions, PS at terminals)
    """
    seq_len = len(sequence)
    all_positions = list(range(1, seq_len + 1))

    # Full 2'-O-methyl modification
    modifications = [ChemicalModification(type="2OMe", positions=all_positions)]

    # Add phosphorothioate linkages at terminal dinucleotides
    # PS linkages are between nucleotides, so we mark positions involved
    terminal_ps = [1, 2, seq_len - 1, seq_len]
    modifications.append(ChemicalModification(type="PS", positions=terminal_ps))

    return modifications


def get_modification_pattern(pattern_name: str, sequence: str) -> list[ChemicalModification]:
    """Get modification pattern by name.

    Args:
        pattern_name: Name of the pattern (standard_2ome, minimal_terminal, maximal_stability, none)
        sequence: RNA sequence to apply pattern to

    Returns:
        List of ChemicalModification objects

    Raises:
        ValueError: If pattern_name is not recognized
    """
    pattern_name = pattern_name.lower()

    if pattern_name == "standard_2ome":
        return apply_standard_2ome_pattern(sequence)
    if pattern_name == "minimal_terminal":
        return apply_minimal_terminal_pattern(sequence)
    if pattern_name == "maximal_stability":
        return apply_maximal_stability_pattern(sequence)
    if pattern_name == "none":
        return []
    raise ValueError(
        f"Unknown modification pattern: {pattern_name}. "
        "Valid options: standard_2ome, minimal_terminal, maximal_stability, none"
    )


def apply_modifications_to_candidate(
    candidate: "SiRNACandidate",
    pattern_name: str = "standard_2ome",
    overhang: str = "dTdT",
    target_gene: Optional[str] = None,
) -> "SiRNACandidate":
    """Apply chemical modifications to a siRNA candidate.

    This function annotates both guide and passenger strands with the specified
    modification pattern and overhang, updating the candidate's metadata fields.

    Args:
        candidate: SiRNACandidate to annotate
        pattern_name: Modification pattern to apply (default: standard_2ome)
        overhang: Overhang sequence (default: dTdT)
        target_gene: Optional target gene name for metadata

    Returns:
        Updated SiRNACandidate with modification metadata
    """
    # Get modification pattern
    guide_mods = get_modification_pattern(pattern_name, candidate.guide_sequence)
    passenger_mods = get_modification_pattern(pattern_name, candidate.passenger_sequence)

    # Create provenance for designed sequences
    provenance = Provenance(
        source_type=SourceType.DESIGNED,
        identifier=f"sirnaforge_{candidate.id}",
        url="https://github.com/Austin-s-h/sirnaforge",
    )

    # Create guide strand metadata
    guide_metadata = StrandMetadata(
        id=f"{candidate.id}_guide",
        sequence=candidate.guide_sequence,
        overhang=overhang,
        chem_mods=guide_mods,
        provenance=provenance,
        confirmation_status=ConfirmationStatus.PENDING,
        notes=f"Guide strand with {pattern_name} modifications" + (f" targeting {target_gene}" if target_gene else ""),
    )

    # Create passenger strand metadata
    passenger_metadata = StrandMetadata(
        id=f"{candidate.id}_passenger",
        sequence=candidate.passenger_sequence,
        overhang=overhang,
        chem_mods=passenger_mods,
        provenance=provenance,
        confirmation_status=ConfirmationStatus.PENDING,
        notes=f"Passenger strand with {pattern_name} modifications",
    )

    # Update candidate with metadata
    candidate.guide_metadata = guide_metadata
    candidate.passenger_metadata = passenger_metadata

    return candidate


def get_modification_summary(candidate: "SiRNACandidate") -> dict[str, str]:
    """Get a summary of modifications for a candidate.

    Args:
        candidate: SiRNACandidate with modification metadata

    Returns:
        Dictionary with modification summary info
    """
    summary = {
        "guide_overhang": "",
        "guide_modifications": "",
        "passenger_overhang": "",
        "passenger_modifications": "",
    }

    if candidate.guide_metadata:
        summary["guide_overhang"] = candidate.guide_metadata.overhang or ""

        mods = []
        for mod in candidate.guide_metadata.chem_mods:
            mods.append(f"{mod.type}({len(mod.positions)})")
        summary["guide_modifications"] = "+".join(mods) if mods else "none"

    if candidate.passenger_metadata:
        summary["passenger_overhang"] = candidate.passenger_metadata.overhang or ""

        mods = []
        for mod in candidate.passenger_metadata.chem_mods:
            mods.append(f"{mod.type}({len(mod.positions)})")
        summary["passenger_modifications"] = "+".join(mods) if mods else "none"

    return summary
