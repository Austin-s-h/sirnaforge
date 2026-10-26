# Separately recovered efficacy panels

These source panels are kept outside `source_extracts/` and `built/` so they do
not silently become part of the existing model-review training collection.

## DSIR 2012 endogenous knockdown panel

Source: Filhol et al., “DSIR: assessing the design of highly potent siRNA by
testing a set of cancer-relevant target genes,” *PLOS ONE* 7(10), e48057 (2012),
[doi:10.1371/journal.pone.0048057](https://doi.org/10.1371/journal.pone.0048057).
The article's license statement permits reuse under CC BY; it does not specify a
version, so this record does not infer one.

The compact `dsir_2012.jsonl` extract joins the guide/passenger strands in
Supplementary Table S1 to the measured residual-expression and reported feature
rows in Supplementary Table S3 by `siRNA_id`. The extraction checked that target,
position and DSIR score agree between the tables. It retains the source sequences,
residual mRNA endpoint, reported DSIR score, and source-computed target-accessibility
and off-target summaries. The measured endpoint is percent residual non-cleaved
target mRNA relative to control; lower is more knockdown. The paper used 20 nM
siRNA in HeLa cells, with a 72-hour endpoint except BCL2L1 at 24 hours. It
reports two design rounds; second-round candidates used findings from the first
round. The article reports a 128-sequence panel in prose, while the S1/S3 join
has 129 unique guide
rows. The table-level discrepancy is preserved rather than silently dropping a
record.

The panel does not contain broad measured unintended repression. The paper's
CSNK2A1/CSNK2A2 cross-reaction check measured only three candidate guides against
one possible off-target gene and found no transcript effect. It is insufficient
for validating the genome-wide liability score. The article also analyzed
several same-panel efficacy features, and its DSIR selection restricts the score
range; therefore the new frozen siRNAforge feature check is exploratory, not
independent validation of a complete profile.

The source files are available from PLOS:

| File | URL suffix | SHA-256 |
| --- | --- | --- |
| Supplementary Table S1, sequences | `journal.pone.0048057.s003` | `4e44096679123c5e5355fcbb4e80cb508d52a7ba3fe72ac7579d7aa55d9703be` |
| Supplementary Table S3, measurements/features | `journal.pone.0048057.s005` | `6602f45d52a7f70f3dfd0a5dea7b2988ec058090fea63de0507871205ea0873b` |
| Derived JSONL extract | `dsir_2012.jsonl` | `9ba410513ed88765ca97b39f0199ec6f2818b0a67712b37daaaba43bbd8b1d68` |

The extraction produced 129 matched source IDs, 129 distinct full guides and
10 target groups. Exact identity checks found no 19-nt guide-core or full 21-nt
guide match in the tracked review source extracts or the OligoGym sequence files
audited by `scripts/validate_dsir2012_independent.py`. This local audit does not
prove separation from every external training corpus.

From the repository root, run:

```bash
uv run python scripts/validate_dsir2012_independent.py
```

The script reports per-target Spearman correlations for asymmetry calculated on
the verified 19-nt paired core (an applicability-limited proxy) and the shipped
GC score calculated on the full 21-nt guide against residual expression,
equal-target means and a target-cluster bootstrap interval. The DSIR duplexes
have two-base 3' overhangs, unlike the shipped fully complementary 21-nt duplex
geometry, so the asymmetry result is not direct validation of that geometry.
The analysis is fixed and does not fit weights. It records
panel and code hashes, the current profile identity, the exact sequence-overlap
scope, and the limits above.
