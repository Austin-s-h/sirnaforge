"""ZFN evaluation and off-target search module."""

from .annotation import GTFZFNAnnotationProvider
from .design import ZFNDesigner
from .interfaces import ZFNAnnotationProvider, ZFNOffTargetSearcher
from .rank import rank_sites
from .search import ExhaustiveZFNOffTargetSearcher

# The ZFN arm ships EXPERIMENTAL for 0.6.0: known defects in half-site orientation
# handling, FokI seed-region weighting and off-target region classification are deferred
# to a tracking issue rather than fixed. Every entry point repeats this text verbatim so a
# user cannot reach ZFN output without seeing it.
ZFN_EXPERIMENTAL_WARNING = (
    "The ZFN module is EXPERIMENTAL and has known unfixed defects in half-site orientation "
    "handling, FokI seed-region weighting and off-target region classification. "
    "Do not use ZFN results for decisions without independent validation."
)

__all__ = [
    "ZFN_EXPERIMENTAL_WARNING",
    "ZFNOffTargetSearcher",
    "ZFNAnnotationProvider",
    "GTFZFNAnnotationProvider",
    "ExhaustiveZFNOffTargetSearcher",
    "ZFNDesigner",
    "rank_sites",
]
