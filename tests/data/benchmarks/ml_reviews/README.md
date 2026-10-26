# Measured datasets from the 2024 and 2026 siRNA ML reviews

This is a **review proposal**, not a complete recovery of all compatible examples
mentioned in either review and not a new validated efficacy benchmark. Both
`coverage_2024.json` and `coverage_2026.json` explicitly record
`all_compatible_examples_recovered: false`. Their entries include inaccessible,
restricted, incompatible, predicted-only and incompletely recovered sources.

Reviews:

- *Machine learning for siRNA efficiency prediction: A systematic review* (2024),
  <https://doi.org/10.1016/j.hsr.2024.100157>.
- *From rules to foundation models: a comprehensive review of machine learning
  approaches for siRNA design* (2026),
  <https://doi.org/10.1093/nargab/lqag099>.

## What is included

`collections.json` names the source extracts, original-file checksums, source
URLs, attributed licenses and interpretation limits. `source_extracts/*.jsonl`
retain individual source records and measurement strings. `built/` contains the
normalized audit records, compatible sequence records and their checksum manifest.
The manifest is the authority for exact per-collection counts.

The current build has 8,287 source/compatible assay records in eight collections
and 4,091 distinct paired-guide sequence clusters. These are **not 8,287 new,
independent measurements**. In particular, Huesken observations recur through the
Ichihara collection and siRNAEfficacyDB; database aggregations and categorical
reuses also overlap original studies. The manifest flags 6,158 accepted rows
whose paired sequences overlap the existing OligoGym table. Neither sequence
identity nor a different collection ID establishes an independent assay.

Multiple doses/replicates remain separate source observations. Sequence hashes
are comparison/grouping aids, not validated development/held-out partitions. No
validation split is assigned, and the old regression fixtures are unchanged.

One source-specific exception is the 32-record `shabalina_independent` subset. The
Shabalina et al. paper describes these as an additional set of 32 siRNAs that was
not included in or analyzed with its 653-guide training set, and reports them as
two low-concentration experiments. The checked-in sequence hashes confirm zero
exact-sequence overlap with the 653 extracted training records: 32 observations,
31 distinct guide sequences, and 10 source accession groups. The paper's results
for its own models (R² 0.258 and 0.281) are not validation results for
siRNAforge.

This small set can support a predeclared, exploratory external check of applicable
sequence-only efficacy features. It does not validate the complete default score
or its gates: all entries are paired cores without full duplex geometry, chemistry
is unstated, and target identity remains unestablished. It contains no measured
unintended-repression endpoint. Any use must keep residual expression as its own
endpoint, group uncertainty by source accession, and report these limitations.
See the [primary paper](https://doi.org/10.1186/1471-2105-7-65), Results section.

A fixed, exploratory check of the shipped GC feature is reproducible with
`uv run python scripts/validate_shabalina_independent.py`. At scoring version
4.1.0, the default design GC score has a pooled Spearman ρ of −0.215 against
residual expression (the expected direction is negative), but its accession
fixed-effect slope is +12.47 residual-expression percentage points per score
unit, with within-accession Pearson r = +0.159. These are point estimates over
32 observations and 10 unverified accession groups, not inferential validation;
the opposing pooled and within-accession directions do not support promoting the
GC weight. The script emits the source checksums, profile identity, group counts,
and limitations with the metrics.

## Interpretation boundaries

- `compatible_core_only` means a reported 19–23 nt core can be represented for
  sequence-level regression. It does **not** establish full strands or overhangs.
- Full strands and overhangs are retained only when the extract supplies them.
  Core-only sense sequences are reverse-complemented into guide orientation;
  literal source fields remain in the JSONL audit records.
- The canonical comparison alphabet is DNA A/C/G/T. Original notation and
  chemistry are retained separately; a chemically modified record is not thereby
  eligible for validating an unmodified-siRNA potency profile.
- Numeric labels, censored measurements, categorical observations, endpoint,
  direction and units remain distinct. IC50 is not percent knockdown; remaining
  expression is not inhibition. No universal label rescaling is performed.
- Every normalized row has `target_identity_status: not_established`. Source
  accessions/coordinates are annotations, not independently verified native
  transcript identity or accessibility evidence.
- `records.csv` is a review fixture table, **not** a registered
  `sirnaforge benchmark prepare` panel. It must not be fed to that interface as
  though its heterogeneous geometry and endpoints were a single existing panel.
- The coverage inventories audit review entries, not distinct experimentally
  independent datasets. Their common partial/gap status is deliberately
  conservative and is not a claim that every downloaded collection is incomplete.

## Rebuild and verification

From the source checkout, using its supported Python environment:

```bash
uv run python scripts/build_review_benchmarks.py \
  --source-dir tests/data/benchmarks/ml_reviews \
  --existing-records-csv tests/data/benchmarks/oligogym/records.csv \
  --out-dir /path/to/rebuilt-review-data

make test-dev
```

The builder is offline and deterministic. It neither downloads source data nor
runs RNA folding, alignment, model fitting or screening. Tests rebuild to an
isolated directory and compare compatible CSV records with the checked-in build.
Original upstream extraction scripts and raw supplements are not all vendored;
the recorded URLs/checksums are provenance handles, not a claim of end-to-end
upstream-download regeneration.

## Third-party attribution and licenses

The repository's software license does **not** relicense third-party datasets.
Retain the collection-level provenance and license notices when reusing records;
several sources have noncommercial restrictions. Assess suitability and any
additional permissions before commercial redistribution or use. Consult the
linked original sources for the governing terms.

| Collection | Primary/source link | License recorded in catalog |
| --- | --- | --- |
| Ichihara 2007 (includes Huesken reuse) | <https://doi.org/10.1093/nar/gkm699> | CC BY-NC 2.0 UK |
| Katoh 2007 | <https://doi.org/10.1093/nar/gkl1120> | CC BY-NC 2.0 UK |
| Shabalina 2006 and independent set | <https://doi.org/10.1186/1471-2105-7-65> | CC BY 2.0 |
| Wang 2009 validation | <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2794195/supplementaryFiles> | CC BY-NC 2.5 |
| siRNAEfficacyDB aggregation | <https://doi.org/10.1049/syb2.12102> | CC BY-NC 2.0 UK (as recorded; constituent-source terms also matter) |
| Kliuchnikov 2024 | <https://doi.org/10.1016/j.omtn.2024.102415> | CC BY 4.0 |
| Takasaki measured categorical set | <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1764479/supplementaryFiles> | CC BY 2.0 |

Do not describe a record as independent external validation merely because it
arrives through another model repository, review or database. Independence must
be established from primary-source study design and checked against the actual
training sequences; the `shabalina_independent` subset above is documented on
that basis, with its limits stated explicitly.
