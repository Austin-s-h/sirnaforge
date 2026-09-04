# ZFN Off-Target Ranking

> **⚠️ EXPERIMENTAL.** The ZFN arm has known unfixed defects tracked in the ZFN experimental-status
> issue; do not use ZFN output for any decision without independent validation. See
> [ZFN Module Guide](zfn_module.md).

`siRNAforge` provides deterministic ZFN off-target ordering via `sirnaforge.zfn.rank.rank_sites`.

## Policy

When two sites have equal score, tie-breaks are applied in this order:

1. Region priority: `exon` > `promoter` > `intron` > `intergenic` > `unknown`
2. Chromosomal location (chromosome, then `start_1based`, then `end_1based`)
3. Stable fallback by `site_id`

This policy is designed to match PROGNOS-style ranking behavior used by the ZFN regression tests.

## API

```python
from sirnaforge.zfn.rank import rank_sites

ranked = rank_sites(sites, params)
```

- `sites`: sequence of `ZFNOffTargetSite`
- `params`: optional `ZFNDesignParameters` (reserved for algorithm-specific extension)
- returns: new `list[ZFNOffTargetSite]` sorted by rank
