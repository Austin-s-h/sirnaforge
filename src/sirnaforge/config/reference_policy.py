"""Reference/default resolution utilities for workflow inputs."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum

DEFAULT_TRANSCRIPTOME_SOURCES: tuple[str, ...] = (
    "ensembl_human_cdna",
    "ensembl_mouse_cdna",
    "ensembl_rat_cdna",
    "ensembl_macaque_cdna",
)
DEFAULT_TRANSCRIPTOME_SOURCE = DEFAULT_TRANSCRIPTOME_SOURCES[0]

DEFAULT_MIRNA_SOURCE = "mirgenedb"
DEFAULT_MIRNA_CANONICAL_SPECIES: tuple[str, ...] = (
    "chicken",
    "pig",
    "rat",
    "mouse",
    "human",
    "rhesus",
    "macaque",
)


class ReferenceState(str, Enum):
    """Describe how a reference input was selected."""

    EXPLICIT = "explicit"
    DEFAULT = "default"
    DISABLED = "disabled"


@dataclass(frozen=True)
class ReferenceChoice:
    """Normalized representation of a resolved reference input."""

    value: str | None
    state: ReferenceState
    reason: str

    @property
    def enabled(self) -> bool:
        """Return True when a usable reference has been selected."""
        return self.value is not None and self.state is not ReferenceState.DISABLED

    @staticmethod
    def explicit(value: str, reason: str = "user-provided") -> ReferenceChoice:
        """Create an explicit user-selected reference choice."""
        return ReferenceChoice(value=value, state=ReferenceState.EXPLICIT, reason=reason)

    @staticmethod
    def default(value: str, reason: str = "auto-selected") -> ReferenceChoice:
        """Create a default-sourced reference choice."""
        return ReferenceChoice(value=value, state=ReferenceState.DEFAULT, reason=reason)

    @staticmethod
    def disabled(reason: str) -> ReferenceChoice:
        """Create a disabled reference choice with context."""
        return ReferenceChoice(value=None, state=ReferenceState.DISABLED, reason=reason)

    def to_metadata(self) -> dict[str, str | None | bool]:
        """Return a serializable snapshot for logs/JSON summaries."""
        return {
            "value": self.value,
            "state": self.state.value,
            "reason": self.reason,
            "enabled": self.enabled,
        }


@dataclass(frozen=True)
class ReferenceSelection:
    """Container describing zero or more resolved references."""

    choices: tuple[ReferenceChoice, ...] = field(default_factory=tuple)
    disabled_reason: str | None = None

    @property
    def enabled(self) -> bool:
        """Return True when at least one reference is configured."""
        return bool(self.choices)

    @staticmethod
    def disabled(reason: str) -> ReferenceSelection:
        """Create a disabled selection with a descriptive reason."""
        return ReferenceSelection(choices=(), disabled_reason=reason)

    def to_metadata(self) -> dict[str, object]:
        """Render selection metadata for logging."""
        return {
            "enabled": self.enabled,
            "disabled_reason": self.disabled_reason,
            "choices": [choice.to_metadata() for choice in self.choices],
        }


@dataclass(frozen=True)
class WorkflowInputSpec:
    """Raw workflow inputs prior to policy resolution."""

    input_fasta: str | None = None
    transcriptome_argument: str | None = None
    default_transcriptomes: Sequence[str] = field(default_factory=lambda: DEFAULT_TRANSCRIPTOME_SOURCES)
    design_only: bool = False
    allow_transcriptome_for_input_fasta: bool = False


class ReferencePolicyResolver:
    """Resolve workflow defaults while preserving intent metadata."""

    def __init__(self, spec: WorkflowInputSpec):
        """Create a resolver for a specific workflow input specification."""
        self.spec = spec

    def resolve_transcriptomes(self) -> ReferenceSelection:
        """Return one or more transcriptome references."""
        if self.spec.design_only:
            return ReferenceSelection.disabled("design-only mode requested")

        transcription_arg = (self.spec.transcriptome_argument or "").strip()
        if transcription_arg:
            choice = ReferenceChoice.explicit(transcription_arg, reason="explicit transcriptome override")
            return ReferenceSelection(choices=(choice,))

        if self.spec.input_fasta and not self.spec.allow_transcriptome_for_input_fasta:
            return ReferenceSelection.disabled("input FASTA runs default to design-only mode")

        defaults = tuple(ref.strip() for ref in self.spec.default_transcriptomes if ref and ref.strip())
        if defaults:
            choices = tuple(ReferenceChoice.default(ref, reason="auto transcriptome default") for ref in defaults)
            return ReferenceSelection(choices=choices)

        return ReferenceSelection.disabled("no transcriptome defaults configured")


def render_reference_selection_label(selection: ReferenceSelection) -> str:
    """Render a stable, human-readable label for CLI/config summaries."""
    if selection.enabled:
        rendered_choices = [f"{choice.value} ({choice.state.value})" for choice in selection.choices if choice.value]
        return ", ".join(rendered_choices)

    reason = selection.disabled_reason or "not available"
    return f"disabled ({reason})"


__all__ = [
    "DEFAULT_TRANSCRIPTOME_SOURCE",
    "DEFAULT_TRANSCRIPTOME_SOURCES",
    "DEFAULT_MIRNA_SOURCE",
    "DEFAULT_MIRNA_CANONICAL_SPECIES",
    "ReferenceChoice",
    "ReferenceSelection",
    "ReferencePolicyResolver",
    "ReferenceState",
    "WorkflowInputSpec",
    "render_reference_selection_label",
]
