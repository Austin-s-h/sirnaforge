"""Exhaustive ZFN off-target search implementation for provided half-sites."""

from __future__ import annotations

import fnmatch
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import TypedDict

from sirnaforge.data.annotation_manager import AnnotationManager
from sirnaforge.data.genome_manager import GenomeManager
from sirnaforge.models.zfn import (
    DimerMode,
    GenomicAnnotationConfig,
    IUPACMode,
    MatchOrientation,
    Strand,
    ZFNAlgorithm,
    ZFNDesignParameters,
    ZFNOffTargetSite,
)
from sirnaforge.utils.fasta import load_fasta_sequences
from sirnaforge.utils.logging_utils import get_logger

from .annotation import GTFZFNAnnotationProvider
from .interfaces import ZFNAnnotationProvider
from .rank import rank_sites

logger = get_logger(__name__)

IUPAC_MAP: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "N": {"A", "C", "G", "T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "W": {"A", "T"},
    "S": {"C", "G"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
}

_RC_TRANS = str.maketrans("ACGTNRYWSKMBDHV", "TGCANYRWSMKVHDB")


@dataclass(slots=True)
class _HalfSiteHit:
    """Internal representation of a half-site near-match hit."""

    kind: str
    chrom: str
    start0: int
    end0: int
    strand: Strand
    query: str
    observed: str
    mismatches: int
    seed_mismatches: int
    mismatch_positions: list[int]
    aligned: str


@dataclass(slots=True)
class _ShardSpec:
    """Internal chromosome/chunk shard specification."""

    shard_id: str
    chrom: str
    core_start0: int
    core_end0: int
    scan_start0: int
    scan_end0: int


class _WindowMatch(TypedDict):
    """Typed match payload returned by one half-site window evaluation."""

    mismatches: int
    seed_mismatches: int
    mismatch_positions: list[int]
    aligned: str


class ExhaustiveZFNOffTargetSearcher:
    """Exhaustive sliding-window off-target search for a provided ZFN pair."""

    def __init__(self, annotation_provider: ZFNAnnotationProvider | None = None) -> None:
        """Initialize searcher with optional annotation provider."""
        self.annotation_provider = annotation_provider or GTFZFNAnnotationProvider()

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Return reverse complement for one uppercase DNA/IUPAC sequence."""
        return seq.translate(_RC_TRANS)[::-1]

    def search(
        self,
        params: ZFNDesignParameters,
        annotation: GenomicAnnotationConfig | None = None,
    ) -> list[ZFNOffTargetSite]:
        """Search all predicted cut sites with explicit mismatch + spacer constraints."""
        phase_timings: dict[str, float] = {}

        phase_start = time.perf_counter()
        fasta_path = self._resolve_search_space_fasta(params)
        if annotation is not None:
            self._resolve_annotation_path(annotation)
        phase_timings["resolve_inputs_s"] = time.perf_counter() - phase_start

        phase_start = time.perf_counter()
        chrom_sequences = self._load_fasta(fasta_path)
        phase_timings["load_fasta_s"] = time.perf_counter() - phase_start

        phase_start = time.perf_counter()
        shard_specs = self._build_shard_specs(chrom_sequences, params)
        phase_timings["build_shards_s"] = time.perf_counter() - phase_start
        all_sites: list[ZFNOffTargetSite] = []

        workers = min(params.sharding.max_workers, len(shard_specs)) if shard_specs else 1
        logger.info("Starting ZFN shard search: %s shard(s), workers=%s", len(shard_specs), workers)
        if workers == 1 and len(shard_specs) > 8:
            logger.info(
                "ZFN search is currently using one shard worker; runtime tuning remains internal and is auto-selected"
            )

        started_at = time.perf_counter()
        progress_interval = max(1, len(shard_specs) // 20) if shard_specs else 1
        if workers > 1 and len(shard_specs) > 1:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(
                        self._search_shard,
                        shard=shard,
                        chrom_sequences=chrom_sequences,
                        params=params,
                    ): shard.shard_id
                    for shard in shard_specs
                }
                for completed_count, future in enumerate(as_completed(futures), start=1):
                    all_sites.extend(future.result())
                    if completed_count == len(shard_specs) or completed_count % progress_interval == 0:
                        elapsed = time.perf_counter() - started_at
                        logger.info(
                            "ZFN shard progress: %s/%s complete (%.1f%%) after %.1fs",
                            completed_count,
                            len(shard_specs),
                            (completed_count * 100.0) / max(1, len(shard_specs)),
                            elapsed,
                        )
        else:
            for completed_count, shard in enumerate(shard_specs, start=1):
                all_sites.extend(self._search_shard(shard=shard, chrom_sequences=chrom_sequences, params=params))
                if completed_count == len(shard_specs) or completed_count % progress_interval == 0:
                    elapsed = time.perf_counter() - started_at
                    logger.info(
                        "ZFN shard progress: %s/%s complete (%.1f%%) after %.1fs",
                        completed_count,
                        len(shard_specs),
                        (completed_count * 100.0) / max(1, len(shard_specs)),
                        elapsed,
                    )
        phase_timings["search_shards_s"] = time.perf_counter() - started_at

        phase_start = time.perf_counter()
        deduped = self._dedupe_sites(all_sites)
        phase_timings["dedupe_s"] = time.perf_counter() - phase_start

        if annotation and self.annotation_provider:
            phase_start = time.perf_counter()
            deduped = self.annotation_provider.annotate(deduped, annotation)
            phase_timings["annotate_s"] = time.perf_counter() - phase_start

        phase_start = time.perf_counter()
        ranked = rank_sites(deduped, params)
        phase_timings["rank_s"] = time.perf_counter() - phase_start

        phase_start = time.perf_counter()
        if len(ranked) > params.top_n_sites:
            ranked = ranked[: params.top_n_sites]
        phase_timings["truncate_s"] = time.perf_counter() - phase_start

        logger.info(
            "ZFN search phase timings: %s",
            json.dumps(phase_timings, sort_keys=True),
        )
        logger.info(
            "ZFN search counts: shards=%s raw_sites=%s deduped_sites=%s ranked_sites=%s top_n=%s",
            len(shard_specs),
            len(all_sites),
            len(deduped),
            len(ranked),
            params.top_n_sites,
        )

        return ranked

    def _build_shard_specs(self, chrom_sequences: dict[str, str], params: ZFNDesignParameters) -> list[_ShardSpec]:
        """Build chromosome/chunk shard specs with overlap safety guarantees."""
        sharding = params.sharding
        selected_contigs = self._resolve_target_contigs(chrom_sequences, sharding.chromosomes)
        if not selected_contigs:
            raise ValueError("No matching chromosomes available for configured ZFN sharding filter")

        sharding_active = sharding.enabled

        if not sharding_active:
            return [
                _ShardSpec(
                    shard_id=f"{chrom}:1-{len(chrom_sequences[chrom])}",
                    chrom=chrom,
                    core_start0=0,
                    core_end0=len(chrom_sequences[chrom]),
                    scan_start0=0,
                    scan_end0=len(chrom_sequences[chrom]),
                )
                for chrom in selected_contigs
            ]

        max_spacer = max(params.spacer_constraints.allowed_spacer_lengths)
        required_overlap = max(
            50,
            max(len(params.left_half_site), params.half_site_constraints.max_len)
            + max_spacer
            + max(len(params.right_half_site), params.half_site_constraints.max_len),
        )
        overlap = max(sharding.overlap_bp, required_overlap)

        shards: list[_ShardSpec] = []
        for chrom in selected_contigs:
            chrom_len = len(chrom_sequences[chrom])
            chunk_size = max(1, sharding.chunk_size_bp)

            chunk_start0 = 0
            while chunk_start0 < chrom_len:
                chunk_end0 = min(chrom_len, chunk_start0 + chunk_size)
                scan_start0 = max(0, chunk_start0 - overlap)
                scan_end0 = min(chrom_len, chunk_end0 + overlap)

                shards.append(
                    _ShardSpec(
                        shard_id=f"{chrom}:{chunk_start0 + 1}-{chunk_end0}",
                        chrom=chrom,
                        core_start0=chunk_start0,
                        core_end0=chunk_end0,
                        scan_start0=scan_start0,
                        scan_end0=scan_end0,
                    )
                )
                chunk_start0 = chunk_end0

        logger.info(
            "ZFN sharding enabled: %s shards across %s contigs (chunk=%sbp, overlap=%sbp)",
            len(shards),
            len(selected_contigs),
            sharding.chunk_size_bp,
            overlap,
        )
        return shards

    def _resolve_target_contigs(self, chrom_sequences: dict[str, str], requested: list[str]) -> list[str]:
        """Resolve requested chromosome filters against loaded FASTA contigs."""
        if not requested:
            return list(chrom_sequences.keys())

        tokens = [token.strip() for token in requested if token.strip()]
        if not tokens:
            return list(chrom_sequences.keys())

        resolved: list[str] = []
        for chrom in chrom_sequences:
            normalized = self._normalize_chrom(chrom)
            if any(self._chrom_matches_token(chrom, normalized, token) for token in tokens):
                resolved.append(chrom)
        return resolved

    def _chrom_matches_token(self, chrom: str, normalized: str, token: str) -> bool:
        """Return whether a contig matches one filter token."""
        raw_token = token.strip().lower()
        normalized_token = self._normalize_chrom(token)
        result = False

        if raw_token in {"*", "all"}:
            result = True
        elif raw_token in {"autosomes", "autosome"}:
            result = normalized.isdigit() and 1 <= int(normalized) <= 22
        elif raw_token in {"sex", "gonosomes"}:
            result = normalized in {"x", "y"}
        elif raw_token in {"mito", "mitochondrial", "mtdna"}:
            result = normalized in {"m", "mt"}
        elif "-" in normalized_token:
            start, end = normalized_token.split("-", 1)
            if start.isdigit() and end.isdigit() and normalized.isdigit():
                lower = int(start)
                upper = int(end)
                if lower > upper:
                    lower, upper = upper, lower
                value = int(normalized)
                result = lower <= value <= upper
        elif any(ch in token for ch in "*?[]"):
            chrom_lc = chrom.lower()
            result = fnmatch.fnmatch(chrom_lc, raw_token) or fnmatch.fnmatch(normalized, normalized_token)
        else:
            result = normalized == normalized_token or chrom.lower() == raw_token

        return result

    def _normalize_chrom(self, chrom: str) -> str:
        """Normalize chromosome labels for robust alias matching (chr3 == 3)."""
        token = chrom.strip().lower()
        if token.startswith("chr"):
            token = token[3:]
        return token

    def _search_shard(
        self,
        shard: _ShardSpec,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
    ) -> list[ZFNOffTargetSite]:
        """Run one shard search and return local paired sites."""
        left_hits = self._scan_half_site(
            kind="L",
            query=params.left_half_site,
            chrom_sequences=chrom_sequences,
            params=params,
            target_chrom=shard.chrom,
            region_start0=shard.scan_start0,
            region_end0=shard.scan_end0,
        )
        right_hits = self._scan_half_site(
            kind="R",
            query=params.right_half_site,
            chrom_sequences=chrom_sequences,
            params=params,
            target_chrom=shard.chrom,
            region_start0=shard.scan_start0,
            region_end0=shard.scan_end0,
        )

        return self._pair_hits(left_hits, right_hits, chrom_sequences, params)

    def _dedupe_sites(self, sites: list[ZFNOffTargetSite]) -> list[ZFNOffTargetSite]:
        """Dedupe by coordinates+orientation and keep the highest-scoring site."""
        deduped: dict[tuple[str, int, int, MatchOrientation], ZFNOffTargetSite] = {}
        for site in sites:
            key = (site.chrom, site.start_1based, site.end_1based, site.orientation)
            prev = deduped.get(key)
            if prev is None or site.score > prev.score:
                deduped[key] = site
        return list(deduped.values())

    def _resolve_search_space_fasta(self, params: ZFNDesignParameters) -> Path:
        """Resolve search-space FASTA via explicit input or cache-managed sources."""
        manager = GenomeManager()

        if params.search_space_fasta:
            custom = manager.get_custom_genome(params.search_space_fasta, build_index=False)
            if custom is None or "fasta" not in custom:
                raise ValueError(f"Unable to resolve search_space_fasta: {params.search_space_fasta}")
            resolved = Path(custom["fasta"])
            logger.info(f"Resolved ZFN search space from explicit FASTA: {resolved}")
            return resolved

        source = params.search_space_reference
        if not source:
            raise ValueError("Provide either search_space_fasta or search_space_reference for ZFN off-target search")

        cached = manager.get_genome(source_name=source, build_index=False)
        if cached is None or "fasta" not in cached:
            raise ValueError(f"Unable to resolve search_space_reference via cache: {source}")

        resolved = Path(cached["fasta"])
        logger.info(f"Resolved ZFN search space from cache source '{source}': {resolved}")
        return resolved

    def _resolve_annotation_path(self, annotation: GenomicAnnotationConfig) -> Path | None:
        """Resolve annotation path for future annotation backends."""
        manager = AnnotationManager(cache_dir=annotation.cache_dir)

        input_location = annotation.annotation_path or annotation.annotation_reference
        if not input_location:
            return None

        resolved = manager.get_custom_annotation(input_location)
        if resolved is None or not resolved.exists():
            logger.warning("Unable to resolve annotation input: %s", input_location)
            return None

        annotation.annotation_path = str(resolved)
        return resolved

    def _load_fasta(self, fasta_path: Path) -> dict[str, str]:
        """Load contig/chromosome sequences from FASTA."""
        return load_fasta_sequences(fasta_path)

    def _scan_half_site(
        self,
        kind: str,
        query: str,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        target_chrom: str | None = None,
        region_start0: int | None = None,
        region_end0: int | None = None,
    ) -> list[_HalfSiteHit]:
        """Exhaustively scan all chromosomes and both strands for one half-site."""
        query = query.upper()
        qlen = len(query)
        seed_positions = self._seed_positions(kind, qlen, params.half_site_constraints.seed_len_from_fokI)

        hits: list[_HalfSiteHit] = []
        for chrom, chrom_seq in chrom_sequences.items():
            if target_chrom is not None and chrom != target_chrom:
                continue
            if len(chrom_seq) < qlen:
                continue

            scan_start0 = 0 if region_start0 is None else max(0, region_start0)
            scan_end0 = len(chrom_seq) if region_end0 is None else min(len(chrom_seq), region_end0)
            if scan_end0 - scan_start0 < qlen:
                continue

            stride = params.half_site_constraints.window_stride
            for i in range(scan_start0, scan_end0 - qlen + 1, stride):
                window_plus = chrom_seq[i : i + qlen]

                plus_match = self._evaluate_window(
                    kind=kind,
                    query=query,
                    observed=window_plus,
                    mode=params.half_site_constraints.iupac_mode,
                    seed_len=params.half_site_constraints.seed_len_from_fokI,
                    max_mm=params.half_site_constraints.max_mismatches,
                    max_seed_mm=params.half_site_constraints.seed_max_mismatches,
                    seed_positions=seed_positions,
                )
                if plus_match is not None:
                    hits.append(
                        _HalfSiteHit(
                            kind=kind,
                            chrom=chrom,
                            start0=i,
                            end0=i + qlen,
                            strand=Strand.PLUS,
                            query=query,
                            observed=window_plus,
                            mismatches=plus_match["mismatches"],
                            seed_mismatches=plus_match["seed_mismatches"],
                            mismatch_positions=plus_match["mismatch_positions"],
                            aligned=plus_match["aligned"],
                        )
                    )

                window_minus_oriented = self._reverse_complement(window_plus)
                minus_match = self._evaluate_window(
                    kind=kind,
                    query=query,
                    observed=window_minus_oriented,
                    mode=params.half_site_constraints.iupac_mode,
                    seed_len=params.half_site_constraints.seed_len_from_fokI,
                    max_mm=params.half_site_constraints.max_mismatches,
                    max_seed_mm=params.half_site_constraints.seed_max_mismatches,
                    seed_positions=seed_positions,
                )
                if minus_match is not None:
                    hits.append(
                        _HalfSiteHit(
                            kind=kind,
                            chrom=chrom,
                            start0=i,
                            end0=i + qlen,
                            strand=Strand.MINUS,
                            query=query,
                            observed=window_minus_oriented,
                            mismatches=minus_match["mismatches"],
                            seed_mismatches=minus_match["seed_mismatches"],
                            mismatch_positions=minus_match["mismatch_positions"],
                            aligned=minus_match["aligned"],
                        )
                    )

        return hits

    def _evaluate_window(
        self,
        kind: str,
        query: str,
        observed: str,
        mode: IUPACMode,
        seed_len: int | None,
        max_mm: int,
        max_seed_mm: int | None,
        expanded_query: str | None = None,
        seed_positions: set[int] | None = None,
    ) -> _WindowMatch | None:
        """Evaluate one window against one query under IUPAC + seed constraints."""
        query_for_match = expanded_query if (expanded_query and mode == IUPACMode.EXPAND_IUPAC) else query

        mismatches = 0
        seed_mismatches = 0
        mismatch_positions: list[int] = []
        aligned_chars: list[str] = []

        active_seed_positions = seed_positions
        if active_seed_positions is None:
            active_seed_positions = self._seed_positions(kind, len(query_for_match), seed_len)

        for idx, (q, o) in enumerate(zip(query_for_match, observed, strict=False)):
            is_match = self._base_match(q, o, mode)
            if is_match:
                aligned_chars.append(o)
            else:
                mismatches += 1
                mismatch_positions.append(idx)
                aligned_chars.append(o.lower())
                if idx in active_seed_positions:
                    seed_mismatches += 1
                if mismatches > max_mm:
                    return None

        if max_seed_mm is not None and seed_len is not None and seed_mismatches > max_seed_mm:
            return None

        return {
            "mismatches": mismatches,
            "seed_mismatches": seed_mismatches,
            "mismatch_positions": mismatch_positions,
            "aligned": "".join(aligned_chars),
        }

    def _seed_positions(self, kind: str, seq_len: int, seed_len: int | None) -> set[int]:
        """Return seed positions nearest FokI according to half-site side."""
        if seed_len is None or seed_len <= 0:
            return set()
        effective = min(seed_len, seq_len)

        if kind == "L":
            # Left half-site FokI-proximal side is toward 3' end in canonical L...R representation.
            return set(range(seq_len - effective, seq_len))

        # Right half-site FokI-proximal side is toward 5' end in canonical L...R representation.
        return set(range(0, effective))

    def _base_match(self, query_base: str, observed_base: str, mode: IUPACMode) -> bool:
        """Return whether one query base matches one observed base under configured IUPAC mode."""
        query_base = query_base.upper()
        observed_base = observed_base.upper()

        if mode == IUPACMode.NONE:
            return query_base == observed_base

        allowed = IUPAC_MAP.get(query_base, {query_base})
        return observed_base in allowed

    @staticmethod
    def _index_hits_by_pos(hits: list[_HalfSiteHit]) -> dict[str, dict[int, list[_HalfSiteHit]]]:
        """Build a two-level index: chrom → start0 → [hits].

        Allows O(1) lookup of hits at an expected genomic position when
        probing for valid spacer-distance partners.
        """
        index: dict[str, dict[int, list[_HalfSiteHit]]] = {}
        for hit in hits:
            index.setdefault(hit.chrom, {}).setdefault(hit.start0, []).append(hit)
        return index

    def _pair_hits(  # noqa: PLR0912
        self,
        left_hits: list[_HalfSiteHit],
        right_hits: list[_HalfSiteHit],
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
    ) -> list[ZFNOffTargetSite]:
        """Join half-site hits into full predicted cut sites by spacer and orientation rules."""
        sites: list[ZFNOffTargetSite] = []
        allowed_spacers = set(params.spacer_constraints.allowed_spacer_lengths)

        # Build coordinate indexes for O(hits × |spacers|) pairing instead of O(n²).
        # Index hits by (chrom, start0) so we can probe expected partner positions directly.
        right_by_pos = self._index_hits_by_pos(right_hits)
        left_by_pos = self._index_hits_by_pos(left_hits)

        # Case A: left_hit is leftward of right_hit (canonical L...R geometry)
        for left_hit in left_hits:
            chrom_rights = right_by_pos.get(left_hit.chrom)
            if not chrom_rights:
                continue
            for spacer in allowed_spacers:
                for right_hit in chrom_rights.get(left_hit.end0 + spacer, []):
                    site = self._build_site(
                        first=left_hit,
                        second=right_hit,
                        chrom_seq=chrom_sequences[left_hit.chrom],
                        allowed_spacers=allowed_spacers,
                        params=params,
                    )
                    if site is not None:
                        sites.append(site)

        # Case B: right_hit is leftward of left_hit (inverted R...L geometry)
        for right_hit in right_hits:
            chrom_lefts = left_by_pos.get(right_hit.chrom)
            if not chrom_lefts:
                continue
            for spacer in allowed_spacers:
                for left_hit in chrom_lefts.get(right_hit.end0 + spacer, []):
                    site = self._build_site(
                        first=left_hit,
                        second=right_hit,
                        chrom_seq=chrom_sequences[right_hit.chrom],
                        allowed_spacers=allowed_spacers,
                        params=params,
                    )
                    if site is not None:
                        sites.append(site)

        if params.dimer_mode == DimerMode.INCLUDE_HOMODIMERS:
            # Homodimer left×left: reuse the already-indexed left_by_pos
            for left_hit in left_hits:
                chrom_lefts = left_by_pos.get(left_hit.chrom)
                if not chrom_lefts:
                    continue
                for spacer in allowed_spacers:
                    for second in chrom_lefts.get(left_hit.end0 + spacer, []):
                        site = self._build_site(
                            first=left_hit,
                            second=second,
                            chrom_seq=chrom_sequences[left_hit.chrom],
                            allowed_spacers=allowed_spacers,
                            params=params,
                        )
                        if site is not None:
                            sites.append(site)

            # Homodimer right×right: reuse the already-indexed right_by_pos
            for right_hit in right_hits:
                chrom_rights = right_by_pos.get(right_hit.chrom)
                if not chrom_rights:
                    continue
                for spacer in allowed_spacers:
                    for second in chrom_rights.get(right_hit.end0 + spacer, []):
                        site = self._build_site(
                            first=right_hit,
                            second=second,
                            chrom_seq=chrom_sequences[right_hit.chrom],
                            allowed_spacers=allowed_spacers,
                            params=params,
                        )
                        if site is not None:
                            sites.append(site)

        deduped: dict[tuple[str, int, int, MatchOrientation], ZFNOffTargetSite] = {}
        for site in sites:
            key = (site.chrom, site.start_1based, site.end_1based, site.orientation)
            prev = deduped.get(key)
            if prev is None or site.score > prev.score:
                deduped[key] = site
        return list(deduped.values())

    def _build_site(
        self,
        first: _HalfSiteHit,
        second: _HalfSiteHit,
        chrom_seq: str,
        allowed_spacers: set[int],
        params: ZFNDesignParameters,
    ) -> ZFNOffTargetSite | None:
        """Build one site model from two half-site hits if all pairing constraints pass."""
        leftmost, rightmost = (first, second) if first.start0 <= second.start0 else (second, first)

        if params.spacer_constraints.require_opposite_strands and first.strand == second.strand:
            return None

        spacer_len = rightmost.start0 - leftmost.end0
        if spacer_len < 0 or spacer_len not in allowed_spacers:
            return None

        orientation = MatchOrientation(f"{leftmost.kind}...{rightmost.kind}")
        if params.dimer_mode == DimerMode.HETERODIMER_ONLY and orientation not in {
            MatchOrientation.LR,
            MatchOrientation.RL,
        }:
            return None

        site_start_0 = leftmost.start0
        site_end_0 = rightmost.end0
        sequence = chrom_seq[site_start_0:site_end_0]

        if first.kind == "L":
            left_hit, right_hit = first, second
        elif second.kind == "L":
            left_hit, right_hit = second, first
        else:
            left_hit, right_hit = leftmost, rightmost

        score, components = self._score_site_with_components(left_hit, right_hit, params.algorithm)
        total_mismatches = left_hit.mismatches + right_hit.mismatches
        site_id = (
            f"{left_hit.chrom}:{site_start_0 + 1}-{site_end_0}:"
            f"{orientation.value}:mm{left_hit.mismatches}+{right_hit.mismatches}"
        )

        return ZFNOffTargetSite(
            site_id=site_id,
            chrom=left_hit.chrom,
            start_1based=site_start_0 + 1,
            end_1based=site_end_0,
            strand=leftmost.strand,
            orientation=orientation,
            spacer_len=spacer_len,
            sequence=sequence,
            left_mismatches=left_hit.mismatches,
            right_mismatches=right_hit.mismatches,
            left_seed_mismatches=left_hit.seed_mismatches,
            right_seed_mismatches=right_hit.seed_mismatches,
            left_mismatch_positions=left_hit.mismatch_positions,
            right_mismatch_positions=right_hit.mismatch_positions,
            total_mismatches=total_mismatches,
            score=score,
            score_components=components,
            dimer_compatible=True,
            region="unknown",
            nearest_gene=None,
            left_aligned=left_hit.aligned,
            right_aligned=right_hit.aligned,
        )

    def _score_site(self, left_hit: _HalfSiteHit, right_hit: _HalfSiteHit, algorithm: ZFNAlgorithm) -> float:
        """Score one paired site using the selected algorithm family."""
        return self._score_site_with_components(left_hit, right_hit, algorithm)[0]

    def _score_site_with_components(
        self,
        left_hit: _HalfSiteHit,
        right_hit: _HalfSiteHit,
        algorithm: ZFNAlgorithm,
    ) -> tuple[float, dict[str, float]]:
        """Score one paired site and return explainable score components."""
        if algorithm == ZFNAlgorithm.HOMOLOGY:
            return self._score_homology_components(left_hit, right_hit)
        if algorithm == ZFNAlgorithm.CONSERVED_G:
            return self._score_conserved_g_components(left_hit, right_hit)
        return self._score_zfn_v2_components(left_hit, right_hit)

    def _score_homology(self, left_hit: _HalfSiteHit, right_hit: _HalfSiteHit) -> float:
        """Mismatch-count score in [0,100]."""
        return self._score_homology_components(left_hit, right_hit)[0]

    def _score_homology_components(
        self, left_hit: _HalfSiteHit, right_hit: _HalfSiteHit
    ) -> tuple[float, dict[str, float]]:
        """Mismatch-count score in [0,100]."""
        mismatch_penalty = (left_hit.mismatches + right_hit.mismatches) * 12.0
        seed_penalty = (left_hit.seed_mismatches + right_hit.seed_mismatches) * 8.0
        score = max(0.0, min(100.0, 100.0 - mismatch_penalty - seed_penalty))
        return score, {
            "algorithm_weighted_penalty": mismatch_penalty + seed_penalty,
            "mismatch_penalty": mismatch_penalty,
            "seed_penalty": seed_penalty,
        }

    def _score_conserved_g(self, left_hit: _HalfSiteHit, right_hit: _HalfSiteHit) -> float:
        """Homology score plus a heuristic bonus for conserved intended G contacts."""
        return self._score_conserved_g_components(left_hit, right_hit)[0]

    def _score_conserved_g_components(
        self, left_hit: _HalfSiteHit, right_hit: _HalfSiteHit
    ) -> tuple[float, dict[str, float]]:
        """Homology score plus a heuristic bonus for conserved intended G contacts."""
        base, base_components = self._score_homology_components(left_hit, right_hit)
        total_g_positions = left_hit.query.count("G") + right_hit.query.count("G")
        if total_g_positions == 0:
            base_components["conserved_g_bonus"] = 0.0
            return base, base_components

        conserved_left = sum(
            1 for q, o in zip(left_hit.query, left_hit.observed, strict=False) if q == "G" and o == "G"
        )
        conserved_right = sum(
            1 for q, o in zip(right_hit.query, right_hit.observed, strict=False) if q == "G" and o == "G"
        )
        conserved_ratio = (conserved_left + conserved_right) / total_g_positions
        bonus = conserved_ratio * 10.0
        score = max(0.0, min(100.0, base + bonus))
        base_components["conserved_g_bonus"] = bonus
        base_components["algorithm_weighted_penalty"] = max(
            0.0,
            base_components["mismatch_penalty"] + base_components["seed_penalty"] - bonus,
        )
        return score, base_components

    def _score_zfn_v2(self, left_hit: _HalfSiteHit, right_hit: _HalfSiteHit) -> float:
        """Finger-aware, polarity-weighted, compensatory score in [0,100]."""
        return self._score_zfn_v2_components(left_hit, right_hit)[0]

    def _score_zfn_v2_components(
        self, left_hit: _HalfSiteHit, right_hit: _HalfSiteHit
    ) -> tuple[float, dict[str, float]]:
        """Finger-aware, polarity-weighted, compensatory score in [0,100]."""
        left_penalty = self._half_site_v2_penalty(left_hit, is_left=True)
        right_penalty = self._half_site_v2_penalty(right_hit, is_left=False)

        # Compensation heuristic: unbalanced penalties partially compensate each other.
        compensation_bonus = min(6.0, abs(left_penalty - right_penalty) * 0.25)
        total_penalty = max(0.0, left_penalty + right_penalty - compensation_bonus)

        score = max(0.0, min(100.0, 100.0 - total_penalty))
        return score, {
            "algorithm_weighted_penalty": total_penalty,
            "left_penalty": left_penalty,
            "right_penalty": right_penalty,
            "compensation_bonus": compensation_bonus,
        }

    def _half_site_v2_penalty(self, hit: _HalfSiteHit, is_left: bool) -> float:
        """Compute one half-site penalty with finger and polarity weighting."""
        if not hit.mismatch_positions:
            return 0.0

        seq_len = len(hit.query)
        finger_len = 3
        penalties: list[float] = []

        for pos in hit.mismatch_positions:
            dist_from_foki = (seq_len - 1) - pos if is_left else pos

            polarity_weight = 1.0 + (0.08 * max(0, 5 - dist_from_foki))
            finger_idx = pos // finger_len
            finger_weight = 1.0 + (0.05 * finger_idx)
            penalties.append(9.0 * polarity_weight * finger_weight)

        return sum(penalties)
