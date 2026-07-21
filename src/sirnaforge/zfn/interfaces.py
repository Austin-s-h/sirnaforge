"""Interfaces for ZFN off-target search and annotation engines."""

from __future__ import annotations

from typing import Protocol

from sirnaforge.models.zfn import GenomicAnnotationConfig, ZFNDesignParameters, ZFNOffTargetSite


class ZFNOffTargetSearcher(Protocol):
    """Protocol for ZFN off-target site searchers."""

    def search(
        self,
        params: ZFNDesignParameters,
        annotation: GenomicAnnotationConfig | None = None,
    ) -> list[ZFNOffTargetSite]:
        """Return predicted off-target cut sites for a provided ZFN pair."""


class ZFNAnnotationProvider(Protocol):
    """Optional protocol for site-level region annotation."""

    def annotate(
        self,
        sites: list[ZFNOffTargetSite],
        config: GenomicAnnotationConfig,
    ) -> list[ZFNOffTargetSite]:
        """Attach region and nearest-gene annotation to predicted sites."""
