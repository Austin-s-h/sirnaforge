"""High-level ZFN pair evaluation workflow."""

from __future__ import annotations

import sys
import time

from sirnaforge import __version__
from sirnaforge.models.zfn import (
    GenomicAnnotationConfig,
    ZFNCandidate,
    ZFNDesignParameters,
    ZFNDesignResult,
    ZFNMutationConstraints,
    ZFNOffTargetFilterCriteria,
    ZFNOffTargetSite,
)
from sirnaforge.utils.logging_utils import get_logger

from .experimental import emit_zfn_experimental_warning
from .interfaces import ZFNOffTargetSearcher
from .search import ExhaustiveZFNOffTargetSearcher

logger = get_logger(__name__)


class ZFNDesigner:
    """Evaluate one provided ZFN pair with exhaustive off-target search."""

    def __init__(self, searcher: ZFNOffTargetSearcher | None = None) -> None:
        """Initialize with a searcher implementation."""
        self.searcher = searcher or ExhaustiveZFNOffTargetSearcher()

    def evaluate_pair(
        self,
        params: ZFNDesignParameters,
        annotation: GenomicAnnotationConfig | None = None,
    ) -> ZFNDesignResult:
        """Run end-to-end pair evaluation and return typed results."""
        # This is the documented Python API entry point, so it has to carry the
        # experimental notice itself -- a caller who never touches the CLI or
        # SiRNAWorkflow would otherwise get ZFN artifacts with no indication of status.
        emit_zfn_experimental_warning()

        start = time.perf_counter()
        sites = self.searcher.search(params=params, annotation=annotation)
        candidate = self._build_candidate(params, sites)
        elapsed = max(0.0, time.perf_counter() - start)

        return ZFNDesignResult(
            parameters=params,
            annotation=annotation,
            candidates=[candidate],
            off_target_sites=sites,
            processing_time_s=elapsed,
            tool_versions=self._tool_versions(),
        )

    def _build_candidate(self, params: ZFNDesignParameters, sites: list[ZFNOffTargetSite]) -> ZFNCandidate:
        """Aggregate site-level data into one candidate summary."""
        filters = params.off_target_filters
        concerning = [site for site in sites if site.score >= filters.min_site_score_to_count]
        exonic = [site for site in concerning if site.region == "exon"]
        promoter = [site for site in concerning if site.region == "promoter"]

        pass_status: bool | str = self._apply_filters(filters, concerning, exonic, promoter)

        if concerning:
            worst_score = min(site.score for site in concerning)
            best_score = max(site.score for site in concerning)
        else:
            worst_score = None
            best_score = None

        on_target_quality = self._score_on_target(params)
        off_target_specificity = self._score_specificity(concerning)
        manufacturability = self._score_manufacturability(params)
        weights = params.scoring
        composite = (
            weights.on_target_quality * on_target_quality
            + weights.off_target_specificity * off_target_specificity
            + weights.manufacturability * manufacturability
        )

        return ZFNCandidate(
            id="ZFN_PAIR_1",
            left_half_site=params.left_half_site,
            right_half_site=params.right_half_site,
            allowed_spacers=list(params.spacer_constraints.allowed_spacer_lengths),
            predicted_sites_total=len(concerning),
            predicted_sites_exonic=len(exonic),
            predicted_sites_promoter=len(promoter),
            worst_site_score=worst_score,
            best_offtarget_score=best_score,
            passes_offtarget_filters=pass_status,
            component_scores={
                "on_target_quality": on_target_quality,
                "off_target_specificity": off_target_specificity,
                "manufacturability": manufacturability,
            },
            composite_score=max(0.0, min(100.0, composite)),
            top_offtargets=sorted(sites, key=lambda s: s.score, reverse=True)[: params.report_n_sites],
        )

    def build_candidate(self, params: ZFNDesignParameters, sites: list[ZFNOffTargetSite]) -> ZFNCandidate:
        """Public wrapper for candidate-summary construction."""
        return self._build_candidate(params, sites)

    def _apply_filters(
        self,
        filters: ZFNOffTargetFilterCriteria,
        concerning: list[ZFNOffTargetSite],
        exonic: list[ZFNOffTargetSite],
        promoter: list[ZFNOffTargetSite],
    ) -> bool | str:
        """Apply post-ranking filter policy."""
        if filters.max_total_sites is not None and len(concerning) > filters.max_total_sites:
            return f"TOTAL_SITES_EXCEEDED ({len(concerning)} > {filters.max_total_sites})"

        if filters.max_exonic_sites is not None and len(exonic) > filters.max_exonic_sites:
            return f"EXONIC_SITES_EXCEEDED ({len(exonic)} > {filters.max_exonic_sites})"

        if filters.max_promoter_sites is not None and len(promoter) > filters.max_promoter_sites:
            return f"PROMOTER_SITES_EXCEEDED ({len(promoter)} > {filters.max_promoter_sites})"

        return True

    def _score_on_target(self, params: ZFNDesignParameters) -> float:
        """Simple on-target construct quality score in [0,100]."""
        half = params.half_site_constraints
        left_len_ok = half.min_len <= len(params.left_half_site) <= half.max_len
        right_len_ok = half.min_len <= len(params.right_half_site) <= half.max_len
        length_score = 100.0 if (left_len_ok and right_len_ok) else 60.0

        spacer_count = len(set(params.spacer_constraints.allowed_spacer_lengths))
        spacer_score = min(100.0, 80.0 + (spacer_count * 5.0))
        return (0.7 * length_score) + (0.3 * spacer_score)

    def _score_specificity(self, concerning: list[ZFNOffTargetSite]) -> float:
        """Specificity score in [0,100] where fewer/high-quality off-targets is better."""
        if not concerning:
            return 100.0

        burden = sum(site.score for site in concerning) / max(1, len(concerning))
        volume_penalty = min(60.0, float(len(concerning)) * 3.0)
        similarity_penalty = burden * 0.4
        return max(0.0, 100.0 - volume_penalty - similarity_penalty)

    def _score_manufacturability(self, params: ZFNDesignParameters) -> float:
        """Manufacturability heuristic in [0,100].

        Considers:
        - IUPAC ambiguity codes (each non-ACGT base adds complexity)
        - Homopolymer runs (impede synthesis)
        - Mutation constraint budgets (per-sub-finger and global, if provided)

        Each zinc finger recognises a 3-bp triplet.  IUPAC ambiguity codes
        within a triplet represent tolerated mutations; the constraint
        budgets set an upper limit on how many such positions are acceptable
        per finger and across the whole half-site.
        """
        left = params.left_half_site
        right = params.right_half_site
        combined = left + right

        # --- Base complexity penalty (IUPAC ambiguity) ---
        iupac_count = sum(1 for base in combined if base not in {"A", "C", "G", "T"})
        complexity_score = max(0.0, 100.0 - (iupac_count * 4.0))

        # --- Homopolymer repeat penalty ---
        repeat_penalty = 0.0
        if "AAAA" in combined or "TTTT" in combined or "CCCC" in combined or "GGGG" in combined:
            repeat_penalty = 15.0

        base_score = max(0.0, complexity_score - repeat_penalty)

        # --- Mutation constraint penalties ---
        constraint_penalty = self._compute_mutation_constraint_penalty(params)

        return max(0.0, base_score - constraint_penalty)

    def score_manufacturability(self, params: ZFNDesignParameters) -> float:
        """Public wrapper for manufacturability scoring logic."""
        return self._score_manufacturability(params)

    # ------------------------------------------------------------------
    # Mutation-constraint helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _count_ambiguous_per_triplet(seq: str) -> list[int]:
        """Split *seq* into 3-bp triplets and count non-ACGT bases per triplet.

        Trailing bases that do not form a complete triplet are grouped as
        one shorter "finger" so every base is accounted for.
        """
        canonical = {"A", "C", "G", "T"}
        triplets: list[int] = []
        for i in range(0, len(seq), 3):
            chunk = seq[i : i + 3]
            triplets.append(sum(1 for b in chunk if b not in canonical))
        return triplets

    @staticmethod
    def count_ambiguous_per_triplet(seq: str) -> list[int]:
        """Public wrapper for triplet ambiguity counting helper."""
        return ZFNDesigner._count_ambiguous_per_triplet(seq)

    def _compute_mutation_constraint_penalty(self, params: ZFNDesignParameters) -> float:
        """Return a [0, 40] penalty based on mutation constraint violations.

        The penalty is 0 when no constraints are set or all budgets are
        satisfied.  Each per-finger violation adds up to 10 pts, and
        overall-budget violations add up to 20 pts, capped at 40 total.
        """
        mc: ZFNMutationConstraints | None = params.mutation_constraints
        if mc is None:
            return 0.0

        left_trips = self._count_ambiguous_per_triplet(params.left_half_site)
        right_trips = self._count_ambiguous_per_triplet(params.right_half_site)
        all_trips = left_trips + right_trips

        penalty = 0.0

        # Per-sub-finger checks (explicit entries override the default)
        explicit_indices: set[int] = set()
        for sf in mc.subfinger_mutations:
            explicit_indices.add(sf.subfinger_index)
            idx = sf.subfinger_index - 1  # 1-based → 0-based
            if 0 <= idx < len(all_trips) and all_trips[idx] > sf.max_mutations:
                excess = all_trips[idx] - sf.max_mutations
                penalty += min(10.0, excess * 5.0)

        # Default mutation allowance for fingers not explicitly constrained
        if mc.default_subfinger_mutation is not None:
            default_max = mc.default_subfinger_mutation.max_mutations
            for idx, count in enumerate(all_trips):
                if (idx + 1) in explicit_indices:
                    continue
                if count > default_max:
                    excess = count - default_max
                    penalty += min(10.0, excess * 5.0)

        # Global / overall mutation budgets
        total_ambiguous = sum(all_trips)
        for ov in mc.overall_mutations:
            if total_ambiguous > ov.max_mutations:
                excess = total_ambiguous - ov.max_mutations
                penalty += min(20.0, excess * 4.0)

        return min(40.0, penalty)

    def _tool_versions(self) -> dict[str, str]:
        """Collect tool versions used in evaluation."""
        return {
            "sirnaforge": __version__,
            "python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        }

    def tool_versions(self) -> dict[str, str]:
        """Public wrapper for tool version metadata."""
        return self._tool_versions()
