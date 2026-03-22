"""Exhaustive ZFN off-target search implementation for provided half-sites."""

from __future__ import annotations

import fnmatch
import importlib
import json
import pickle
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import Any, Protocol, TypedDict, cast

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
    ZFNSearchBackend,
)
from sirnaforge.utils.cache_utils import stable_cache_key
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


class _HalfSiteScanEngine(Protocol):
    """Internal backend contract for half-site scan engines."""

    backend: ZFNSearchBackend

    def scan_half_site(
        self,
        kind: str,
        query: str,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        target_chrom: str | None = None,
        region_start0: int | None = None,
        region_end0: int | None = None,
    ) -> list[_HalfSiteHit]:
        """Return verified half-site hits for one query across the requested region."""
        ...


class _SearchBackendUnavailableError(RuntimeError):
    """Raised when a configured search backend dependency is unavailable."""


@dataclass(frozen=True, slots=True)
class _MultiFMIndexState:
    """Loaded multi-document FM index plus contig-to-doc mapping."""

    index: Any
    doc_ids_by_chrom: dict[str, int]


class _SearchIndexManifest(TypedDict):
    """Persisted search-space index manifest payload."""

    format_version: int
    backend: str
    backend_class: str
    genome_fasta: str
    artifact: str
    contig_names: list[str]
    contig_lengths: dict[str, int]


_SEARCH_INDEX_MANIFEST = "manifest.json"
_FM_INDEX_ARTIFACT = "multifm_index.pkl"
_SEARCH_INDEX_FORMAT_VERSION = 1


def _load_fm_index_module() -> Any:
    """Import the fm_index package with consistent backend error handling."""
    try:
        return importlib.import_module("fm_index")
    except ImportError as exc:
        raise _SearchBackendUnavailableError(
            "ZFN search backend 'fm_index' is unavailable because the installed environment is missing "
            "the required 'fm-index' package"
        ) from exc


def _resolve_search_index_bundle(index_path: str | Path) -> Path:
    """Resolve a search-space index bundle path to an absolute directory."""
    bundle_dir = Path(index_path).expanduser().resolve()
    if not bundle_dir.exists():
        raise ValueError(f"Configured search_space_index does not exist: {bundle_dir}")
    if not bundle_dir.is_dir():
        raise ValueError(f"Configured search_space_index must be a directory bundle: {bundle_dir}")
    return bundle_dir


def _load_search_index_manifest(bundle_dir: Path) -> _SearchIndexManifest:
    """Load and parse one persisted search-index manifest."""
    manifest_path = bundle_dir / _SEARCH_INDEX_MANIFEST
    if not manifest_path.exists():
        raise ValueError(f"search_space_index bundle is missing {_SEARCH_INDEX_MANIFEST}: {bundle_dir}")

    with manifest_path.open(encoding="utf-8") as handle:
        return cast(_SearchIndexManifest, json.load(handle))


def _chrom_sequence_lengths(chrom_sequences: dict[str, str]) -> dict[str, int]:
    """Return contig lengths for one loaded FASTA dictionary."""
    return {chrom: len(seq) for chrom, seq in chrom_sequences.items()}


def _resolve_cached_genome_fasta(genome_fasta: Path) -> Path:
    """Resolve one genome FASTA through the shared genome cache pipeline."""
    manager = GenomeManager()
    cached = manager.get_custom_genome(genome_fasta, build_index=False)
    if cached is None or "fasta" not in cached:
        raise ValueError(f"Unable to resolve genome FASTA for search-index build: {genome_fasta}")
    return Path(cached["fasta"])


def _default_search_index_bundle_dir(genome_fasta: Path, backend: ZFNSearchBackend) -> Path:
    """Return the canonical cache-backed bundle directory for one backend and cached FASTA."""
    bundle_key = stable_cache_key(
        {
            "backend": backend.value,
            "fasta": str(genome_fasta),
            "format_version": _SEARCH_INDEX_FORMAT_VERSION,
        }
    )
    return genome_fasta.parent / f"{bundle_key}_{backend.value}_index"


def build_zfn_search_index(
    *,
    backend: ZFNSearchBackend,
    genome_fasta: Path,
    output_dir: Path | None = None,
) -> dict[str, Any]:
    """Build a persisted search-space index bundle for indexed ZFN backends."""
    resolved_fasta = _resolve_cached_genome_fasta(genome_fasta)
    output_dir = (
        _default_search_index_bundle_dir(resolved_fasta, backend)
        if output_dir is None
        else output_dir.expanduser().resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    chrom_sequences = load_fasta_sequences(resolved_fasta)
    if not chrom_sequences:
        raise ValueError(f"No contigs found in genome FASTA: {resolved_fasta}")

    if backend != ZFNSearchBackend.FM_INDEX:
        raise ValueError(
            "Persisted search-space indexes are currently implemented only for the 'fm_index' backend; "
            "pyahocorasick remains query-derived and exhaustive_python is scan-only"
        )

    fm_index = _load_fm_index_module()
    contig_names = list(chrom_sequences.keys())
    multi_index = fm_index.MultiFMIndex([chrom_sequences[name] for name in contig_names])
    artifact_path = output_dir / _FM_INDEX_ARTIFACT
    with artifact_path.open("wb") as handle:
        pickle.dump(multi_index, handle)

    manifest: _SearchIndexManifest = {
        "format_version": _SEARCH_INDEX_FORMAT_VERSION,
        "backend": backend.value,
        "backend_class": "MultiFMIndex",
        "genome_fasta": str(resolved_fasta),
        "artifact": _FM_INDEX_ARTIFACT,
        "contig_names": contig_names,
        "contig_lengths": _chrom_sequence_lengths(chrom_sequences),
    }
    with (output_dir / _SEARCH_INDEX_MANIFEST).open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    return {
        "backend": backend.value,
        "bundle_dir": str(output_dir),
        "artifact": str(artifact_path),
        "genome_fasta": str(resolved_fasta),
        "contigs": len(contig_names),
        "bases": sum(len(seq) for seq in chrom_sequences.values()),
    }


def _reverse_complement_seq(seq: str) -> str:
    """Return reverse complement for one uppercase DNA/IUPAC sequence."""
    return seq.translate(_RC_TRANS)[::-1]


def _seed_positions_for(kind: str, seq_len: int, seed_len: int | None) -> set[int]:
    """Return seed positions nearest FokI according to half-site side."""
    if seed_len is None or seed_len <= 0:
        return set()
    effective = min(seed_len, seq_len)

    if kind == "L":
        return set(range(seq_len - effective, seq_len))

    return set(range(0, effective))


def _base_match_for(query_base: str, observed_base: str, mode: IUPACMode) -> bool:
    """Return whether one query base matches one observed base under configured IUPAC mode."""
    query_base = query_base.upper()
    observed_base = observed_base.upper()

    if mode == IUPACMode.NONE:
        return query_base == observed_base

    allowed = IUPAC_MAP.get(query_base, {query_base})
    return observed_base in allowed


def _evaluate_window_match(
    *,
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
    """Evaluate one window against one query under IUPAC and seed constraints."""
    query_for_match = expanded_query if (expanded_query and mode == IUPACMode.EXPAND_IUPAC) else query

    mismatches = 0
    seed_mismatches = 0
    mismatch_positions: list[int] = []
    aligned_chars: list[str] = []

    active_seed_positions = seed_positions
    if active_seed_positions is None:
        active_seed_positions = _seed_positions_for(kind, len(query_for_match), seed_len)

    for idx, (q, o) in enumerate(zip(query_for_match, observed, strict=False)):
        is_match = _base_match_for(q, o, mode)
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


class _BaseHalfSiteScanEngine:
    """Shared helpers for concrete half-site scan backends."""

    _OBSERVED_ALPHABET = tuple(IUPAC_MAP.keys())

    def __init__(self, owner: ExhaustiveZFNOffTargetSearcher) -> None:
        self._owner = owner
        self._pattern_cache: dict[tuple[Any, ...], tuple[str, ...]] = {}
        self._pattern_lock = Lock()

    def _pattern_cache_key(self, kind: str, query: str, params: ZFNDesignParameters) -> tuple[Any, ...]:
        constraints = params.half_site_constraints
        return (
            kind,
            query,
            constraints.iupac_mode,
            constraints.max_mismatches,
            constraints.seed_len_from_fokI,
            constraints.seed_max_mismatches,
        )

    def _candidate_patterns(self, kind: str, query: str, params: ZFNDesignParameters) -> tuple[str, ...]:
        key = self._pattern_cache_key(kind, query, params)
        cached = self._pattern_cache.get(key)
        if cached is not None:
            return cached

        with self._pattern_lock:
            cached = self._pattern_cache.get(key)
            if cached is not None:
                return cached

            constraints = params.half_site_constraints
            seed_positions = _seed_positions_for(kind, len(query), constraints.seed_len_from_fokI)
            max_seed_mm = constraints.seed_max_mismatches
            mode = constraints.iupac_mode
            results: set[str] = set()
            builder: list[str] = []

            def walk(idx: int, mismatches: int, seed_mismatches: int) -> None:
                if idx == len(query):
                    results.add("".join(builder))
                    return

                allowed = {query[idx]} if mode == IUPACMode.NONE else IUPAC_MAP.get(query[idx], {query[idx]})
                in_seed = idx in seed_positions
                for base in self._OBSERVED_ALPHABET:
                    delta = 0 if base in allowed else 1
                    next_mismatches = mismatches + delta
                    if next_mismatches > constraints.max_mismatches:
                        continue

                    next_seed_mismatches = seed_mismatches + (delta if in_seed else 0)
                    if (
                        max_seed_mm is not None
                        and constraints.seed_len_from_fokI is not None
                        and next_seed_mismatches > max_seed_mm
                    ):
                        continue

                    builder.append(base)
                    walk(idx + 1, next_mismatches, next_seed_mismatches)
                    builder.pop()

            walk(0, 0, 0)
            materialized = tuple(sorted(results))
            self._pattern_cache[key] = materialized
            return materialized

    def _build_verified_hit(
        self,
        *,
        kind: str,
        chrom: str,
        start0: int,
        strand: Strand,
        query: str,
        observed_oriented: str,
        params: ZFNDesignParameters,
    ) -> _HalfSiteHit | None:
        match = _evaluate_window_match(
            kind=kind,
            query=query,
            observed=observed_oriented,
            mode=params.half_site_constraints.iupac_mode,
            seed_len=params.half_site_constraints.seed_len_from_fokI,
            max_mm=params.half_site_constraints.max_mismatches,
            max_seed_mm=params.half_site_constraints.seed_max_mismatches,
        )
        if match is None:
            return None

        return _HalfSiteHit(
            kind=kind,
            chrom=chrom,
            start0=start0,
            end0=start0 + len(query),
            strand=strand,
            query=query,
            observed=observed_oriented,
            mismatches=match["mismatches"],
            seed_mismatches=match["seed_mismatches"],
            mismatch_positions=match["mismatch_positions"],
            aligned=match["aligned"],
        )


class _ExhaustivePythonHalfSiteScanEngine(_BaseHalfSiteScanEngine):
    """Baseline exhaustive Python sliding-window scan engine."""

    backend = ZFNSearchBackend.EXHAUSTIVE_PYTHON

    def scan_half_site(
        self,
        kind: str,
        query: str,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        target_chrom: str | None = None,
        region_start0: int | None = None,
        region_end0: int | None = None,
    ) -> list[_HalfSiteHit]:
        query = query.upper()
        qlen = len(query)
        seed_positions = _seed_positions_for(kind, qlen, params.half_site_constraints.seed_len_from_fokI)

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

                plus_match = _evaluate_window_match(
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

                window_minus_oriented = _reverse_complement_seq(window_plus)
                minus_match = _evaluate_window_match(
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


class _PyAhoCorasickHalfSiteScanEngine(_BaseHalfSiteScanEngine):
    """Automaton-backed exact candidate search with shared hit verification."""

    backend = ZFNSearchBackend.PYAHOCORASICK

    def __init__(self, owner: ExhaustiveZFNOffTargetSearcher) -> None:
        super().__init__(owner)
        self._automaton_cache: dict[tuple[Any, ...], Any] = {}
        self._automaton_lock = Lock()

    @staticmethod
    def _load_module() -> Any:
        try:
            return importlib.import_module("ahocorasick")
        except ImportError as exc:
            raise _SearchBackendUnavailableError(
                "ZFN search backend 'pyahocorasick' is unavailable because the installed environment is missing "
                "the required 'pyahocorasick' package"
            ) from exc

    def _automaton(self, kind: str, query: str, params: ZFNDesignParameters) -> Any:
        key = self._pattern_cache_key(kind, query, params)
        cached = self._automaton_cache.get(key)
        if cached is not None:
            return cached

        with self._automaton_lock:
            cached = self._automaton_cache.get(key)
            if cached is not None:
                return cached

            ahocorasick = self._load_module()
            payload_map: dict[str, list[tuple[Strand, str]]] = {}
            for oriented_pattern in self._candidate_patterns(kind, query, params):
                payload_map.setdefault(oriented_pattern, []).append((Strand.PLUS, oriented_pattern))
                minus_genomic_pattern = _reverse_complement_seq(oriented_pattern)
                payload_map.setdefault(minus_genomic_pattern, []).append((Strand.MINUS, oriented_pattern))

            automaton = ahocorasick.Automaton()
            for pattern, payloads in payload_map.items():
                automaton.add_word(pattern, tuple(payloads))
            automaton.make_automaton()
            self._automaton_cache[key] = automaton
            return automaton

    def _scan_with_prebuilt_index(
        self,
        *,
        kind: str,
        query: str,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        target_chrom: str | None,
        region_start0: int | None,
        region_end0: int | None,
    ) -> list[_HalfSiteHit] | None:
        """Aho-Corasick does not currently persist search-space bundles because its automaton is query-derived."""
        del kind, query, chrom_sequences, target_chrom, region_start0, region_end0
        if params.search_space_index:
            logger.info(
                "ZFN backend '%s' received search_space_index=%s, but pyahocorasick does not reuse persisted "
                "search-space bundles because its automaton is derived from the query patterns; falling back to "
                "an in-memory automaton build",
                self.backend.value,
                params.search_space_index,
            )
        return None

    def scan_half_site(
        self,
        kind: str,
        query: str,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        target_chrom: str | None = None,
        region_start0: int | None = None,
        region_end0: int | None = None,
    ) -> list[_HalfSiteHit]:
        query = query.upper()
        qlen = len(query)
        prebuilt_hits = self._scan_with_prebuilt_index(
            kind=kind,
            query=query,
            chrom_sequences=chrom_sequences,
            params=params,
            target_chrom=target_chrom,
            region_start0=region_start0,
            region_end0=region_end0,
        )
        if prebuilt_hits is not None:
            return prebuilt_hits
        automaton = self._automaton(kind, query, params)
        hits: list[_HalfSiteHit] = []

        for chrom, chrom_seq in chrom_sequences.items():
            if target_chrom is not None and chrom != target_chrom:
                continue

            scan_start0 = 0 if region_start0 is None else max(0, region_start0)
            scan_end0 = len(chrom_seq) if region_end0 is None else min(len(chrom_seq), region_end0)
            if scan_end0 - scan_start0 < qlen:
                continue

            for end_index, payloads in automaton.iter(chrom_seq, scan_start0, scan_end0):
                start0 = end_index - qlen + 1
                if start0 < scan_start0 or start0 + qlen > scan_end0:
                    continue
                for strand, observed_oriented in payloads:
                    hit = self._build_verified_hit(
                        kind=kind,
                        chrom=chrom,
                        start0=start0,
                        strand=strand,
                        query=query,
                        observed_oriented=observed_oriented,
                        params=params,
                    )
                    if hit is not None:
                        hits.append(hit)

        return hits


class _FMIndexHalfSiteScanEngine(_BaseHalfSiteScanEngine):
    """FM-index-backed exact candidate search with shared hit verification."""

    backend = ZFNSearchBackend.FM_INDEX

    def __init__(self, owner: ExhaustiveZFNOffTargetSearcher) -> None:
        super().__init__(owner)
        self._multi_index_cache: dict[tuple[Any, ...], _MultiFMIndexState] = {}
        self._multi_index_lock = Lock()

    @staticmethod
    def _load_module() -> Any:
        return _load_fm_index_module()

    @staticmethod
    def _multi_index_from_sequences(chrom_sequences: dict[str, str]) -> _MultiFMIndexState:
        """Build one live MultiFMIndex for the loaded contigs."""
        fm_index = _load_fm_index_module()
        contig_names = list(chrom_sequences.keys())
        index = fm_index.MultiFMIndex([chrom_sequences[name] for name in contig_names])
        return _MultiFMIndexState(
            index=index,
            doc_ids_by_chrom={chrom: doc_id for doc_id, chrom in enumerate(contig_names)},
        )

    def _live_multi_index(self, chrom_sequences: dict[str, str]) -> _MultiFMIndexState:
        """Return a cached live MultiFMIndex for one loaded FASTA snapshot."""
        key = ("live", tuple((chrom, id(seq)) for chrom, seq in chrom_sequences.items()))
        cached = self._multi_index_cache.get(key)
        if cached is not None:
            return cached

        with self._multi_index_lock:
            cached = self._multi_index_cache.get(key)
            if cached is not None:
                return cached

            state = self._multi_index_from_sequences(chrom_sequences)
            self._multi_index_cache[key] = state
            return state

    def _persisted_multi_index(self, bundle_dir: Path, chrom_sequences: dict[str, str]) -> _MultiFMIndexState:
        """Load and validate one persisted MultiFMIndex bundle."""
        key = ("persisted", str(bundle_dir))
        cached = self._multi_index_cache.get(key)
        if cached is not None:
            return cached

        with self._multi_index_lock:
            cached = self._multi_index_cache.get(key)
            if cached is not None:
                return cached

            manifest = _load_search_index_manifest(bundle_dir)
            if manifest.get("backend") != self.backend.value:
                raise ValueError(
                    f"search_space_index backend mismatch: expected '{self.backend.value}', got '{manifest.get('backend')}'"
                )
            if manifest.get("format_version") != _SEARCH_INDEX_FORMAT_VERSION:
                raise ValueError(
                    "search_space_index format version mismatch: "
                    f"expected {_SEARCH_INDEX_FORMAT_VERSION}, got {manifest.get('format_version')}"
                )

            contig_names = manifest["contig_names"]
            contig_lengths = manifest["contig_lengths"]

            live_lengths = _chrom_sequence_lengths(chrom_sequences)
            if list(chrom_sequences.keys()) != contig_names:
                raise ValueError("search_space_index contig ordering does not match the resolved FASTA")
            if live_lengths != contig_lengths:
                raise ValueError("search_space_index contig lengths do not match the resolved FASTA")

            artifact = manifest["artifact"]
            artifact_path = bundle_dir / artifact
            if not artifact_path.exists():
                raise ValueError(f"search_space_index artifact does not exist: {artifact_path}")

            with artifact_path.open("rb") as handle:
                index = pickle.load(handle)

            state = _MultiFMIndexState(
                index=index,
                doc_ids_by_chrom={chrom: doc_id for doc_id, chrom in enumerate(contig_names)},
            )
            self._multi_index_cache[key] = state
            return state

    def _scan_with_multi_index_state(
        self,
        *,
        state: _MultiFMIndexState,
        kind: str,
        query: str,
        params: ZFNDesignParameters,
        target_chrom: str | None,
        region_start0: int | None,
        region_end0: int | None,
    ) -> list[_HalfSiteHit]:
        """Scan one query against a shared MultiFMIndex state."""
        qlen = len(query)
        patterns = self._candidate_patterns(kind, query, params)
        scan_hits: list[_HalfSiteHit] = []

        selected_chroms = [target_chrom] if target_chrom is not None else list(state.doc_ids_by_chrom.keys())

        for chrom in selected_chroms:
            doc_id = state.doc_ids_by_chrom.get(chrom)
            if doc_id is None:
                continue

            scan_start = 0 if region_start0 is None else max(0, region_start0)
            scan_end = None if region_end0 is None else region_end0

            for oriented_pattern in patterns:
                for start0 in state.index.iter_locate(oriented_pattern, doc_id=doc_id):
                    if start0 < scan_start or (scan_end is not None and start0 + qlen > scan_end):
                        continue
                    hit = self._build_verified_hit(
                        kind=kind,
                        chrom=chrom,
                        start0=start0,
                        strand=Strand.PLUS,
                        query=query,
                        observed_oriented=oriented_pattern,
                        params=params,
                    )
                    if hit is not None:
                        scan_hits.append(hit)

                minus_genomic_pattern = _reverse_complement_seq(oriented_pattern)
                for start0 in state.index.iter_locate(minus_genomic_pattern, doc_id=doc_id):
                    if start0 < scan_start or (scan_end is not None and start0 + qlen > scan_end):
                        continue
                    hit = self._build_verified_hit(
                        kind=kind,
                        chrom=chrom,
                        start0=start0,
                        strand=Strand.MINUS,
                        query=query,
                        observed_oriented=oriented_pattern,
                        params=params,
                    )
                    if hit is not None:
                        scan_hits.append(hit)

        return scan_hits

    def _scan_with_prebuilt_index(
        self,
        *,
        kind: str,
        query: str,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        target_chrom: str | None,
        region_start0: int | None,
        region_end0: int | None,
    ) -> list[_HalfSiteHit] | None:
        """Load and query one persisted MultiFMIndex bundle when configured."""
        if not params.search_space_index:
            return None

        bundle_dir = _resolve_search_index_bundle(params.search_space_index)
        state = self._persisted_multi_index(bundle_dir, chrom_sequences)
        return self._scan_with_multi_index_state(
            state=state,
            kind=kind,
            query=query,
            params=params,
            target_chrom=target_chrom,
            region_start0=region_start0,
            region_end0=region_end0,
        )

    def scan_half_site(
        self,
        kind: str,
        query: str,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        target_chrom: str | None = None,
        region_start0: int | None = None,
        region_end0: int | None = None,
    ) -> list[_HalfSiteHit]:
        query = query.upper()
        prebuilt_hits = self._scan_with_prebuilt_index(
            kind=kind,
            query=query,
            chrom_sequences=chrom_sequences,
            params=params,
            target_chrom=target_chrom,
            region_start0=region_start0,
            region_end0=region_end0,
        )
        if prebuilt_hits is not None:
            return prebuilt_hits

        state = self._live_multi_index(chrom_sequences)
        return self._scan_with_multi_index_state(
            state=state,
            kind=kind,
            query=query,
            params=params,
            target_chrom=target_chrom,
            region_start0=region_start0,
            region_end0=region_end0,
        )


def _normalize_chrom_token(chrom: str) -> str:
    """Normalize chromosome labels for robust alias matching (chr3 == 3)."""
    token = chrom.strip().lower()
    if token.startswith("chr"):
        token = token[3:]
    return token


def _chrom_matches_filter(chrom: str, token: str) -> bool:
    """Return whether a contig matches one filter token."""
    raw_token = token.strip().lower()
    normalized = _normalize_chrom_token(chrom)
    normalized_token = _normalize_chrom_token(token)
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


def resolve_target_contigs(contig_names: list[str], requested: list[str]) -> list[str]:
    """Resolve requested chromosome filters against loaded contig names."""
    if not requested:
        return list(contig_names)

    tokens = [token.strip() for token in requested if token.strip()]
    if not tokens:
        return list(contig_names)

    return [chrom for chrom in contig_names if any(_chrom_matches_filter(chrom, token) for token in tokens)]


def build_zfn_shard_specs(contig_lengths: dict[str, int], params: ZFNDesignParameters) -> list[_ShardSpec]:
    """Build chromosome/chunk shard specs from contig lengths.

    This is the authoritative shard planning logic shared by the direct Python
    searcher and external orchestration layers such as Nextflow.
    """
    sharding = params.sharding
    selected_contigs = resolve_target_contigs(list(contig_lengths.keys()), sharding.chromosomes)
    if not selected_contigs:
        raise ValueError("No matching chromosomes available for configured ZFN sharding filter")

    if not sharding.enabled:
        return [
            _ShardSpec(
                shard_id=f"{chrom}:1-{contig_lengths[chrom]}",
                chrom=chrom,
                core_start0=0,
                core_end0=contig_lengths[chrom],
                scan_start0=0,
                scan_end0=contig_lengths[chrom],
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
        chrom_len = contig_lengths[chrom]
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


class ExhaustiveZFNOffTargetSearcher:
    """Exhaustive sliding-window off-target search for a provided ZFN pair."""

    _DEFAULT_MEMORY_RESERVE_GB = 2
    _DEFAULT_PER_WORKER_GB = 0.5

    def __init__(self, annotation_provider: ZFNAnnotationProvider | None = None) -> None:
        """Initialize searcher with optional annotation provider."""
        self.annotation_provider = annotation_provider or GTFZFNAnnotationProvider()
        self._scan_engines: dict[ZFNSearchBackend, _HalfSiteScanEngine] = {}

    def _scan_engine_for(self, params: ZFNDesignParameters) -> _HalfSiteScanEngine:
        """Return the configured half-site scan engine for this search request."""
        backend = params.search_backend
        cached = self._scan_engines.get(backend)
        if cached is not None:
            return cached

        if backend == ZFNSearchBackend.EXHAUSTIVE_PYTHON:
            engine: _HalfSiteScanEngine = _ExhaustivePythonHalfSiteScanEngine(self)
        elif backend == ZFNSearchBackend.PYAHOCORASICK:
            engine = _PyAhoCorasickHalfSiteScanEngine(self)
        elif backend == ZFNSearchBackend.FM_INDEX:
            engine = _FMIndexHalfSiteScanEngine(self)
        else:
            raise ValueError(f"Unsupported ZFN search backend: {backend}")

        self._scan_engines[backend] = engine
        return engine

    def _recommended_worker_cap(self, chrom_sequences: dict[str, str], params: ZFNDesignParameters) -> int:
        """Return a conservative worker cap based on requested workers and live memory.

        The search path should follow the requested sharding profile unless the
        host lacks enough free memory to sustain that concurrency without risking
        an OOM kill.
        """
        del chrom_sequences
        return min(params.sharding.max_workers, self._memory_based_worker_cap(params.sharding.max_workers))

    @staticmethod
    def _available_memory_gb() -> float | None:
        """Return available system memory in GiB by reading /proc/meminfo, or None."""
        try:
            with Path("/proc/meminfo").open(encoding="ascii") as fh:
                for line in fh:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) / (1024 * 1024)
        except (OSError, ValueError, IndexError):
            pass
        return None

    def _memory_based_worker_cap(self, requested: int) -> int:
        """Cap workers based on live available memory to avoid OOM on small hosts.

        Each concurrent shard worker builds per-shard hit lists in memory.
        On a heavily-loaded 12 GB machine (VS Code, Docker, other experiments)
        running 8 workers on chr3 was enough to invoke the OOM killer.
        By default we reserve 2 GiB for OS + other apps and allow ~0.5 GiB per worker slot.
        """
        avail_gb = self._available_memory_gb()
        if avail_gb is None:
            return requested
        reserve_gb = self._DEFAULT_MEMORY_RESERVE_GB
        per_worker_gb = self._DEFAULT_PER_WORKER_GB
        usable_gb = max(0.0, avail_gb - reserve_gb)
        mem_cap = max(1, int(usable_gb / per_worker_gb))
        if mem_cap < requested:
            logger.warning(
                "ZFN worker count memory-capped from %s to %s "
                "(%.1f GiB available, %.1f GiB usable after %.1f GiB reserve)",
                requested,
                mem_cap,
                avail_gb,
                usable_gb,
                reserve_gb,
            )
        return mem_cap

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Return reverse complement for one uppercase DNA/IUPAC sequence."""
        return _reverse_complement_seq(seq)

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
        scan_engine = self._scan_engine_for(params)
        all_sites: list[ZFNOffTargetSite] = []

        worker_cap = self._recommended_worker_cap(chrom_sequences, params)
        workers = min(worker_cap, len(shard_specs)) if shard_specs else 1
        if workers < params.sharding.max_workers:
            logger.warning(
                "ZFN worker count auto-capped from %s to %s to limit peak memory usage "
                "for a large search space (%s bp, max_mismatches=%s)",
                params.sharding.max_workers,
                workers,
                sum(len(seq) for seq in chrom_sequences.values()),
                params.half_site_constraints.max_mismatches,
            )
        logger.info("Starting ZFN shard search: %s shard(s), workers=%s", len(shard_specs), workers)
        logger.info("Using ZFN scan backend: %s", params.search_backend.value)
        if workers == 1 and len(shard_specs) > 8:
            logger.info(
                "ZFN search is currently using one shard worker; runtime tuning remains internal and is auto-selected"
            )

        started_at = time.perf_counter()
        all_sites.extend(
            self._run_shard_searches(
                shard_specs=shard_specs,
                chrom_sequences=chrom_sequences,
                params=params,
                scan_engine=scan_engine,
                workers=workers,
                started_at=started_at,
            )
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
        ranked = self._truncate_sites(ranked, params.top_n_sites)
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

    def search_region(
        self,
        params: ZFNDesignParameters,
        chrom: str,
        scan_start0: int,
        scan_end0: int,
        core_start0: int | None = None,
        core_end0: int | None = None,
        annotation: GenomicAnnotationConfig | None = None,
        top_n_sites: int | None = None,
    ) -> list[ZFNOffTargetSite]:
        """Search one bounded genomic region using the same core engine as full search.

        The region is scanned across ``scan_start0..scan_end0`` and then filtered to
        the optional core window so overlapping shards can share context without
        double-reporting the same site.
        """
        fasta_path = self._resolve_search_space_fasta(params)
        if annotation is not None:
            self._resolve_annotation_path(annotation)

        chrom_sequences = self._load_fasta(fasta_path)
        if chrom not in chrom_sequences:
            raise ValueError(f"Chromosome {chrom!r} not present in resolved search space")

        scan_start = max(0, scan_start0)
        scan_end = min(len(chrom_sequences[chrom]), scan_end0)
        if scan_start >= scan_end:
            return []

        core_start = scan_start if core_start0 is None else max(scan_start, core_start0)
        core_end = scan_end if core_end0 is None else min(scan_end, core_end0)
        if core_start >= core_end:
            return []

        shard = _ShardSpec(
            shard_id=f"{chrom}:{core_start + 1}-{core_end}",
            chrom=chrom,
            core_start0=core_start,
            core_end0=core_end,
            scan_start0=scan_start,
            scan_end0=scan_end,
        )
        scan_engine = self._scan_engine_for(params)
        shard_sites = self._search_shard(
            shard=shard,
            chrom_sequences=chrom_sequences,
            params=params,
            scan_engine=scan_engine,
        )
        shard_sites = self._filter_sites_to_core_window(
            shard_sites, chrom=chrom, core_start0=core_start, core_end0=core_end
        )
        return self._postprocess_sites(
            shard_sites,
            params=params,
            annotation=annotation,
            top_n_sites=top_n_sites,
        )

    def _build_shard_specs(self, chrom_sequences: dict[str, str], params: ZFNDesignParameters) -> list[_ShardSpec]:
        """Build chromosome/chunk shard specs with overlap safety guarantees."""
        contig_lengths = {chrom: len(sequence) for chrom, sequence in chrom_sequences.items()}
        return build_zfn_shard_specs(contig_lengths, params)

    def _filter_sites_to_core_window(
        self,
        sites: list[ZFNOffTargetSite],
        chrom: str,
        core_start0: int,
        core_end0: int,
    ) -> list[ZFNOffTargetSite]:
        """Keep only sites fully contained in the shard core window."""
        core_start_1 = core_start0 + 1
        return [
            site
            for site in sites
            if site.chrom == chrom and site.start_1based >= core_start_1 and site.end_1based <= core_end0
        ]

    def _truncate_sites(self, sites: list[ZFNOffTargetSite], top_n_sites: int | None) -> list[ZFNOffTargetSite]:
        """Optionally truncate a ranked site list."""
        if top_n_sites is None or len(sites) <= top_n_sites:
            return sites
        return sites[:top_n_sites]

    def _postprocess_sites(
        self,
        sites: list[ZFNOffTargetSite],
        params: ZFNDesignParameters,
        annotation: GenomicAnnotationConfig | None = None,
        top_n_sites: int | None = None,
    ) -> list[ZFNOffTargetSite]:
        """Apply dedupe, optional annotation, ranking, and truncation."""
        deduped = self._dedupe_sites(sites)
        if annotation and self.annotation_provider:
            deduped = self.annotation_provider.annotate(deduped, annotation)
        ranked = rank_sites(deduped, params)
        return self._truncate_sites(ranked, top_n_sites)

    def _run_shard_searches(
        self,
        shard_specs: list[_ShardSpec],
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        scan_engine: _HalfSiteScanEngine,
        workers: int,
        started_at: float,
    ) -> list[ZFNOffTargetSite]:
        """Execute shard searches with either parallel or serial scheduling."""
        if workers > 1 and len(shard_specs) > 1:
            return self._run_parallel_shard_searches(
                shard_specs=shard_specs,
                chrom_sequences=chrom_sequences,
                params=params,
                scan_engine=scan_engine,
                workers=workers,
                started_at=started_at,
            )
        return self._run_serial_shard_searches(
            shard_specs=shard_specs,
            chrom_sequences=chrom_sequences,
            params=params,
            scan_engine=scan_engine,
            started_at=started_at,
        )

    def _run_parallel_shard_searches(
        self,
        shard_specs: list[_ShardSpec],
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        scan_engine: _HalfSiteScanEngine,
        workers: int,
        started_at: float,
    ) -> list[ZFNOffTargetSite]:
        """Execute shard searches using a thread pool."""
        all_sites: list[ZFNOffTargetSite] = []
        progress_interval = max(1, len(shard_specs) // 20)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    self._search_shard,
                    shard=shard,
                    chrom_sequences=chrom_sequences,
                    params=params,
                    scan_engine=scan_engine,
                ): shard.shard_id
                for shard in shard_specs
            }
            completed_count = 0
            pending = set(futures)
            heartbeat_interval_s = 60.0
            last_heartbeat_at = started_at
            while pending:
                done, pending = wait(pending, timeout=heartbeat_interval_s, return_when=FIRST_COMPLETED)
                if not done:
                    last_heartbeat_at = self._maybe_log_shard_heartbeat(
                        completed_count=completed_count,
                        total_shards=len(shard_specs),
                        pending_count=len(pending),
                        started_at=started_at,
                        last_heartbeat_at=last_heartbeat_at,
                        heartbeat_interval_s=heartbeat_interval_s,
                    )
                    continue

                for future in done:
                    all_sites.extend(future.result())
                    completed_count += 1

                self._log_shard_progress_if_needed(
                    completed_count=completed_count,
                    total_shards=len(shard_specs),
                    progress_interval=progress_interval,
                    started_at=started_at,
                )
        return all_sites

    def _run_serial_shard_searches(
        self,
        shard_specs: list[_ShardSpec],
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        scan_engine: _HalfSiteScanEngine,
        started_at: float,
    ) -> list[ZFNOffTargetSite]:
        """Execute shard searches serially."""
        all_sites: list[ZFNOffTargetSite] = []
        progress_interval = max(1, len(shard_specs) // 20) if shard_specs else 1
        for completed_count, shard in enumerate(shard_specs, start=1):
            all_sites.extend(
                self._search_shard(
                    shard=shard,
                    chrom_sequences=chrom_sequences,
                    params=params,
                    scan_engine=scan_engine,
                )
            )
            self._log_shard_progress_if_needed(
                completed_count=completed_count,
                total_shards=len(shard_specs),
                progress_interval=progress_interval,
                started_at=started_at,
            )
        return all_sites

    def _maybe_log_shard_heartbeat(
        self,
        completed_count: int,
        total_shards: int,
        pending_count: int,
        started_at: float,
        last_heartbeat_at: float,
        heartbeat_interval_s: float,
    ) -> float:
        """Emit a heartbeat when no shard has completed recently."""
        now = time.perf_counter()
        if now - last_heartbeat_at < heartbeat_interval_s:
            return last_heartbeat_at

        elapsed = now - started_at
        logger.info(
            "ZFN shard progress: %s/%s complete (%.1f%%) after %.1fs; %s shard(s) still running",
            completed_count,
            total_shards,
            (completed_count * 100.0) / max(1, total_shards),
            elapsed,
            pending_count,
        )
        return now

    def _log_shard_progress_if_needed(
        self,
        completed_count: int,
        total_shards: int,
        progress_interval: int,
        started_at: float,
    ) -> None:
        """Emit normal progress logs for shard completion milestones."""
        if completed_count != total_shards and completed_count % progress_interval != 0:
            return

        elapsed = time.perf_counter() - started_at
        logger.info(
            "ZFN shard progress: %s/%s complete (%.1f%%) after %.1fs",
            completed_count,
            total_shards,
            (completed_count * 100.0) / max(1, total_shards),
            elapsed,
        )

    def _resolve_target_contigs(self, chrom_sequences: dict[str, str], requested: list[str]) -> list[str]:
        """Resolve requested chromosome filters against loaded FASTA contigs."""
        return resolve_target_contigs(list(chrom_sequences.keys()), requested)

    def _chrom_matches_token(self, chrom: str, normalized: str, token: str) -> bool:
        """Return whether a contig matches one filter token."""
        del normalized
        return _chrom_matches_filter(chrom, token)

    def _normalize_chrom(self, chrom: str) -> str:
        """Normalize chromosome labels for robust alias matching (chr3 == 3)."""
        return _normalize_chrom_token(chrom)

    def _search_shard(
        self,
        shard: _ShardSpec,
        chrom_sequences: dict[str, str],
        params: ZFNDesignParameters,
        scan_engine: _HalfSiteScanEngine,
    ) -> list[ZFNOffTargetSite]:
        """Run one shard search and return local paired sites."""
        left_hits = scan_engine.scan_half_site(
            kind="L",
            query=params.left_half_site,
            chrom_sequences=chrom_sequences,
            params=params,
            target_chrom=shard.chrom,
            region_start0=shard.scan_start0,
            region_end0=shard.scan_end0,
        )
        right_hits = scan_engine.scan_half_site(
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
        """Compatibility shim for the legacy built-in Python half-site scan path."""
        exhaustive_engine = self._scan_engines.get(ZFNSearchBackend.EXHAUSTIVE_PYTHON)
        if exhaustive_engine is None:
            exhaustive_engine = _ExhaustivePythonHalfSiteScanEngine(self)
            self._scan_engines[ZFNSearchBackend.EXHAUSTIVE_PYTHON] = exhaustive_engine

        return exhaustive_engine.scan_half_site(
            kind=kind,
            query=query,
            chrom_sequences=chrom_sequences,
            params=params,
            target_chrom=target_chrom,
            region_start0=region_start0,
            region_end0=region_end0,
        )

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
        return _evaluate_window_match(
            kind=kind,
            query=query,
            observed=observed,
            mode=mode,
            seed_len=seed_len,
            max_mm=max_mm,
            max_seed_mm=max_seed_mm,
            expanded_query=expanded_query,
            seed_positions=seed_positions,
        )

    def _seed_positions(self, kind: str, seq_len: int, seed_len: int | None) -> set[int]:
        """Return seed positions nearest FokI according to half-site side."""
        return _seed_positions_for(kind, seq_len, seed_len)

    def _base_match(self, query_base: str, observed_base: str, mode: IUPACMode) -> bool:
        """Return whether one query base matches one observed base under configured IUPAC mode."""
        return _base_match_for(query_base, observed_base, mode)

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
