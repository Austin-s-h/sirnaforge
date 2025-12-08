"""Unit tests for reference/default policy resolution."""

from sirnaforge.config import (
    DEFAULT_TRANSCRIPTOME_SOURCES,
    ReferencePolicyResolver,
    ReferenceState,
    WorkflowInputSpec,
)


def test_resolver_uses_default_when_no_overrides() -> None:
    """Resolver should return the bundled defaults when no inputs override them."""
    spec = WorkflowInputSpec()
    resolver = ReferencePolicyResolver(spec)
    selection = resolver.resolve_transcriptomes()

    assert selection.enabled is True
    assert len(selection.choices) == len(DEFAULT_TRANSCRIPTOME_SOURCES)
    values = [choice.value for choice in selection.choices]
    assert values[0] == DEFAULT_TRANSCRIPTOME_SOURCES[0]
    assert all(choice.state is ReferenceState.DEFAULT for choice in selection.choices)


def test_resolver_respects_explicit_argument() -> None:
    """Explicit transcriptome arguments take precedence over defaults."""
    spec = WorkflowInputSpec(transcriptome_argument="/data/ref.fa")
    resolver = ReferencePolicyResolver(spec)
    selection = resolver.resolve_transcriptomes()

    assert selection.enabled is True
    assert len(selection.choices) == 1
    choice = selection.choices[0]
    assert choice.value == "/data/ref.fa"
    assert choice.state is ReferenceState.EXPLICIT


def test_resolver_disables_when_input_fasta_without_override() -> None:
    """Input FASTA without transcriptome override disables transcriptome analysis."""
    spec = WorkflowInputSpec(input_fasta="custom.fasta")
    resolver = ReferencePolicyResolver(spec)
    selection = resolver.resolve_transcriptomes()

    assert selection.enabled is False
    assert selection.disabled_reason is not None
    assert "design-only" in selection.disabled_reason


def test_resolver_can_allow_transcriptome_with_input_fasta() -> None:
    """Policy can force transcriptome usage even when input FASTA is supplied."""
    spec = WorkflowInputSpec(
        input_fasta="custom.fasta",
        default_transcriptomes=("ensembl_mouse_cdna",),
        allow_transcriptome_for_input_fasta=True,
    )
    resolver = ReferencePolicyResolver(spec)
    selection = resolver.resolve_transcriptomes()

    assert selection.enabled is True
    assert len(selection.choices) == 1
    choice = selection.choices[0]
    assert choice.value == "ensembl_mouse_cdna"
    assert choice.state is ReferenceState.DEFAULT
