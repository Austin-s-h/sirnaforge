"""Pydantic models for siRNA design data structures."""

from .modifications import (
    ChemicalModification,
    ConfirmationStatus,
    Provenance,
    SequenceRecord,
    SourceType,
    StrandMetadata,
    StrandRole,
)
from .off_target import (
    AggregatedMiRNASummary,
    AggregatedOffTargetSummary,
    AlignmentStrand,
    AnalysisMode,
    AnalysisSummary,
    MiRNAHit,
    MiRNASummary,
    OffTargetHit,
)
from .sirna import (
    DesignMode,
    DesignParameters,
    DesignResult,
    FilterCriteria,
    MiRNADesignConfig,
    ScoringWeights,
    SequenceType,
    SiRNACandidate,
)

__all__ = [
    # siRNA design models
    "DesignMode",
    "DesignParameters",
    "DesignResult",
    "FilterCriteria",
    "MiRNADesignConfig",
    "ScoringWeights",
    "SequenceType",
    "SiRNACandidate",
    # Chemical modification models
    "ChemicalModification",
    "ConfirmationStatus",
    "Provenance",
    "SequenceRecord",
    "SourceType",
    "StrandMetadata",
    "StrandRole",
    # Off-target analysis models
    "OffTargetHit",
    "MiRNAHit",
    "AnalysisSummary",
    "MiRNASummary",
    "AggregatedOffTargetSummary",
    "AggregatedMiRNASummary",
    "AlignmentStrand",
    "AnalysisMode",
]
