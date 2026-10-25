"""Typed inputs and outputs for transcript-seed site discovery."""

from enum import Enum

from pydantic import BaseModel, ConfigDict, Field


class TranscriptSeedClass(str, Enum):
    """Canonical seed-site geometry, named for the complementary target site."""

    SIXMER = "6mer"
    SEVENMER_A1 = "7mer-A1"
    SEVENMER_M8 = "7mer-m8"
    EIGHTMER = "8mer"


class TranscriptSequenceScope(str, Enum):
    """Sequence region searched by the transcript-seed scanner."""

    FULL_CDNA = "full_cdna"
    THREE_PRIME_UTR = "three_prime_utr"


class TranscriptSeedTranscript(BaseModel):
    """Transcript sequence plus optional explicitly annotated 3′UTR sequence.

    ``utr3_sequence=None`` means its annotation is unavailable. An empty string means the
    annotation is available and explicitly contains no 3′UTR bases.
    """

    transcript_id: str = Field(min_length=1)
    gene_id: str | None = Field(default=None, min_length=1)
    cdna_sequence: str = Field(min_length=1)
    utr3_sequence: str | None = None
    utr3_start: int | None = Field(default=None, ge=1, description="1-based start in transcript-oriented cDNA")

    model_config = ConfigDict(frozen=True)


class TranscriptSeedSite(BaseModel):
    """One canonical target-site match in transcript-oriented coordinates.

    Sites are unique by guide, transcript, scope, and the target interval pairing to guide
    positions 2–7. When overlapping motif classes share that core anchor, only the strongest
    class is retained. ``start``/``end`` describe its full annotated motif and the core fields
    make the site identity explicit.
    """

    guide_id: str
    transcript_id: str
    gene_id: str | None
    scope: TranscriptSequenceScope
    site_class: TranscriptSeedClass
    start: int = Field(ge=1, description="1-based inclusive transcript coordinate")
    end: int = Field(ge=1, description="1-based inclusive transcript coordinate")
    seed_core_start: int = Field(ge=1, description="1-based start of the shared guide positions 2–7 anchor")
    seed_core_end: int = Field(ge=1, description="1-based end of the shared guide positions 2–7 anchor")
    target_sequence: str

    model_config = ConfigDict(frozen=True)


class TranscriptSeedScanResult(BaseModel):
    """Sites and coverage counters from one scope-specific scan."""

    scope: TranscriptSequenceScope
    sites: tuple[TranscriptSeedSite, ...]
    site_count: int = Field(ge=0)
    transcript_count: int = Field(ge=0, description="Distinct transcripts with at least one site")
    gene_count: int = Field(ge=0, description="Distinct resolved genes with at least one site")
    unresolved_gene_count: int = Field(ge=0, description="Transcripts with a site but no gene identifier")
    transcripts_scanned: int = Field(ge=0)
    transcripts_unavailable: tuple[str, ...] = ()

    model_config = ConfigDict(frozen=True)
