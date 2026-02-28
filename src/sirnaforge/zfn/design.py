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
    ZFNOffTargetFilterCriteria,
    ZFNOffTargetSite,
)
from sirnaforge.utils.logging_utils import get_logger

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
        """Manufacturability heuristic in [0,100]."""
        left = params.left_half_site
        right = params.right_half_site
        combined = left + right

        iupac_count = sum(1 for base in combined if base not in {"A", "C", "G", "T"})
        complexity_score = max(0.0, 100.0 - (iupac_count * 4.0))
        repeat_penalty = 0.0
        if "AAAA" in combined or "TTTT" in combined or "CCCC" in combined or "GGGG" in combined:
            repeat_penalty = 15.0

        return max(0.0, complexity_score - repeat_penalty)

    def _tool_versions(self) -> dict[str, str]:
        """Collect tool versions used in evaluation."""
        return {
            "sirnaforge": __version__,
            "python": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        }
