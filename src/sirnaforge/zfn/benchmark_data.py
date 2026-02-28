"""Typed ingestion helpers for CCR5 ZFN benchmark tables (PROGNOS supplementary extracts)."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, cast

HalfOrientation = Literal["L", "R"]


@dataclass(frozen=True)
class MatchType:
    """Structured representation of PROGNOS match type, e.g. ``L-5-R``."""

    left: HalfOrientation
    spacer_len: int
    right: HalfOrientation


@dataclass(frozen=True)
class CCR5S10VisibleRow:
    """One visible row from the extracted S10 CCR5 off-target validation table."""

    closest_gene: str
    match_type: MatchType
    chrom: str
    pos_hg19: int
    plus_half_site: str
    minus_half_site: str
    empty_indels: int | None
    empty_total: int | None
    active_indels: int | None
    active_total: int | None
    active_mutation_freq_percent: float | None
    p_value: float | None
    notes: str
    sequencing_failure: bool


@dataclass(frozen=True)
class CCR5S11HomologyRow:
    """One visible row from the S11 homology ranking extract."""

    homology_rank: int
    t_mismatches: int
    plus_mismatches: int
    minus_mismatches: int
    interrogated_by: str
    closest_gene: str
    match_type: MatchType
    chrom: str
    pos_hg19: int
    plus_half_site: str
    minus_half_site: str


def parse_match_type(raw: str) -> MatchType:
    """Parse a match-type token such as ``L-5-R`` into structured fields."""
    match = re.fullmatch(r"\s*([LR])-(\d+)-([LR])\s*", raw)
    if match is None:
        raise ValueError(f"Invalid match type: {raw}")
    left = cast(HalfOrientation, match.group(1))
    right = cast(HalfOrientation, match.group(3))
    return MatchType(left=left, spacer_len=int(match.group(2)), right=right)


def parse_hg19_coordinate(raw: str) -> tuple[str, int]:
    """Parse hg19 coordinates in ``chrN:POS`` or compact ``chrNPOS`` format."""
    match = re.fullmatch(r"\s*(chr[0-9A-Za-z]+):(\d+)\s*", raw)
    if match is not None:
        return match.group(1), int(match.group(2))

    compact = raw.strip()
    if not compact.startswith("chr"):
        raise ValueError(f"Invalid hg19 coordinate: {raw}")

    suffix = compact[3:]
    prefixes = [
        "MT",
        "22",
        "21",
        "20",
        "19",
        "18",
        "17",
        "16",
        "15",
        "14",
        "13",
        "12",
        "11",
        "10",
        "9",
        "8",
        "7",
        "6",
        "5",
        "4",
        "3",
        "2",
        "1",
        "X",
        "Y",
        "M",
    ]
    for prefix in prefixes:
        if suffix.startswith(prefix):
            pos = suffix[len(prefix) :]
            if pos and pos.isdigit():
                return f"chr{prefix}", int(pos)

    raise ValueError(f"Invalid hg19 coordinate: {raw}")


def _parse_optional_int(raw: str) -> int | None:
    value = raw.strip()
    if value in {"", "?", "N/A", "Sequencing Failure"} or "?" in value:
        return None
    normalized = value.replace(" ", "")
    if not normalized.isdigit():
        return None
    return int(normalized)


def _parse_optional_percent(raw: str) -> float | None:
    value = raw.strip()
    if value in {"", "?", "N/A", "Sequencing Failure"}:
        return None
    if value.endswith("%"):
        try:
            return float(value[:-1])
        except ValueError:
            return None
    try:
        return float(value)
    except ValueError:
        return None


def _parse_optional_float(raw: str) -> float | None:
    value = raw.strip()
    if value in {"", "?", "N/A", "Sequencing Failure"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def load_ccr5_s10_visible_rows(path: str | Path) -> list[CCR5S10VisibleRow]:
    """Load visible S10 rows from CSV into typed records."""
    rows: list[CCR5S10VisibleRow] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["Closest gene"] is None:
                continue
            chrom, pos = parse_hg19_coordinate(row["hg19 coordinate"])
            sequencing_failure = row["Empty indels"].strip() == "Sequencing Failure"
            rows.append(
                CCR5S10VisibleRow(
                    closest_gene=row["Closest gene"].strip(),
                    match_type=parse_match_type(row["Match type"]),
                    chrom=chrom,
                    pos_hg19=pos,
                    plus_half_site=row["(+) half-site"].strip(),
                    minus_half_site=row["(−) half-site"].strip(),
                    empty_indels=_parse_optional_int(row["Empty indels"]),
                    empty_total=_parse_optional_int(row["Empty total"]),
                    active_indels=_parse_optional_int(row["Active indels"]),
                    active_total=_parse_optional_int(row["Active total"]),
                    active_mutation_freq_percent=_parse_optional_percent(row["Active mutation freq"]),
                    p_value=_parse_optional_float(row["p-value"]),
                    notes=row["Notes"].strip(),
                    sequencing_failure=sequencing_failure,
                )
            )
    return rows


def load_ccr5_s11_homology_rows(path: str | Path) -> list[CCR5S11HomologyRow]:
    """Load visible S11 homology rows from CSV into typed records."""
    rows: list[CCR5S11HomologyRow] = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["Homology rank"] is None:
                continue
            chrom, pos = parse_hg19_coordinate(row["hg19 coordinate"])
            rows.append(
                CCR5S11HomologyRow(
                    homology_rank=int(row["Homology rank"]),
                    t_mismatches=int(row["T mism"]),
                    plus_mismatches=int(row["+ mism"]),
                    minus_mismatches=int(row["− mism"]),
                    interrogated_by=row["Interrogated by"].strip(),
                    closest_gene=row["Closest gene"].strip(),
                    match_type=parse_match_type(row["Match type"]),
                    chrom=chrom,
                    pos_hg19=pos,
                    plus_half_site=row["(+) half-site"].strip(),
                    minus_half_site=row["(−) half-site"].strip(),
                )
            )
    return rows


__all__ = [
    "CCR5S10VisibleRow",
    "CCR5S11HomologyRow",
    "MatchType",
    "load_ccr5_s10_visible_rows",
    "load_ccr5_s11_homology_rows",
    "parse_hg19_coordinate",
    "parse_match_type",
]
