# Target-accessibility regression fixture

Provenance, citations and regeneration notes for the files backing
`tests/unit/test_target_accessibility.py`. Issue #95. Issue #102 added the dataset checksums and the
predeclared splits below. The offline orthologue mapping is documented at the bottom.

## Citations

**Measured siRNA efficacy** — every efficacy value in `sirna_efficacy_subset.csv` comes from:

> Huesken D, Lange J, Mickanin C, Weiler J, Asselbergs F, Warner J, Meloon B, Engel S,
> Rosenberg A, Cohen D, Labow M, Reinhardt M, Natt F, Hall J.
> "Design of a genome-wide siRNA library using an artificial neural network."
> Nat Biotechnol. 2005 Aug;23(8):995-1001. doi:10.1038/nbt1118. PMID: 16025102.

Obtained not from the paper but from the third-party redistribution
<https://github.com/apkrfi/unMod-siRNA-Pred>, file `raw_data/All_Dataset.csv`.

⚠️ **Limitation, stated so it is not discovered later.** The Huesken paper is not open access, so
the efficacy values could **not** be verified against the primary source. This fixture reproduces a
third-party redistribution of them. That bounds what the regression test proves: it pins the
scoring geometry, normalisation and W/L defaults against a _specific published-derived table_, not
against independently verified biology.

The redistributed dataset also carries 385 rows labelled `sources = Other`, whose provenance could
not be established and which therefore cannot be attributed. **They are excluded from this
fixture**, so every vendored row traces to the one citation above. Sources split cleanly by
accession in that dataset, so the exclusion drops whole transcripts rather than thinning any of
them; it cost 0.03 of Spearman rho (+0.275 mixed, +0.2437 Huesken-only) and moved both correlations
_closer_ to the full-set values, so nothing material was lost. The untracked full-dataset script
(`scripts/validate_target_accessibility.py`) still uses all 2,816 rows, because it is a calibration
record rather than a shipped artifact.

**Transcript sequences** — not part of the Huesken paper. They are NCBI Nucleotide (`nuccore`)
records, fetched full-length:

```bash
# accession list below, joined with commas
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=<ACCESSIONS>"
```

Full-length records, deliberately: the ~59-nt local context shipped with the efficacy dataset must
**not** be folded, because local opening probability depends on context well outside the site.

**TP53-201** (`tp53_201.fa`) — `ENST00000269305.9` cDNA from Ensembl, 2,512 nt. TP53 is this
repository's public example gene throughout.

## `sirna_efficacy_subset.csv` + `sirna_efficacy_subset_transcripts.fa`

180 siRNAs with measured knockdown across 8 transcripts (9,051 nt); 17 KB together.

- `efficacy` is **inhibition**: higher means more knockdown.
- `guide_sequence` is the 21 nt guide as published. Its reverse complement is the target site.
- `site_start` is the 0-based offset of that target site in the transcript, precomputed so the test
  does not re-derive coordinates. The test asserts the offsets still describe the sequences.
- **Accessions retained** (8 of the full set's 41): `NM_003969`, `NM_004223`, `NM_006357`,
  `NM_007019`, `NM_012864`, `NM_014176`, `NM_016406`, `XM_214061`.

### Regeneration

```bash
uv run python scripts/build_accessibility_test_fixture.py \
    --benchmark-csv work/sirna_bench.csv \
    --transcripts work/bench_tx.fa \
    --out-dir tests/unit/data
```

That script is the definition of the subset; its docstring states the selection rule. In summary:
keep `sources == Huesken` rows only; keep rows whose target site occurs exactly once in its
transcript; take accessions in descending siRNA-per-kb order until the 10,000 nt transcript budget
is spent (density, because the fixture's cost is transcript length while its power is row count);
then order by (accession, `site_start`) and keep every 4th row, capped at 180. Fully deterministic,
no randomness, and **not** tuned against the resulting correlation.

### What it is asserted to reproduce

| Quantity (RNAplfold W=150, L=100)      | Full set (n=2,779) | This subset (n=180) |
| -------------------------------------- | ------------------ | ------------------- |
| Spearman rho, scored seed-end 8-mer    | +0.267             | +0.2437             |
| Spearman rho, non-seed 8-mer (control) | +0.067             | +0.0639             |

The subset is smaller and noisier, so the test asserts floors and ceilings with headroom rather
than these values. The full set is checked by `scripts/validate_target_accessibility.py`, which
reads the untracked full data from a path.

The pSUPER/pSuper-Retro anti-p53 target site `GACTCCAGTGGTAATCTAC` (0-based offset 916, unique in
TP53-201) must not land in the bottom decile of the scored term across the transcript's 19-mer
windows. It sits at the 93rd percentile.

## Checksums (issue #102)

SHA-256 of every file any scoring measurement was taken from, so a number quoted in a docstring, an
issue or the CHANGELOG can be traced to the exact bytes that produced it.

| File                                              |  Bytes | SHA-256                                                            |
| ------------------------------------------------- | -----: | ------------------------------------------------------------------ |
| `sirna_efficacy_subset.csv` (tracked)             |  7,586 | `6e0a2efea505a1df7110136a8818e5bbbb28243dc12becdfbc73402fa8a5e10d` |
| `sirna_efficacy_subset_transcripts.fa` (tracked)  |  9,293 | `261ea9615878b6138e29f34ab099aa070a33db8dbfa36f8a5b4fe8b95a8588e8` |
| `work/sirna_bench.csv` (**untracked** full panel) | 334,475 | `a53caa3b1deb8fd9d2c9496a68c1e709af5a836f264e540361fb9dbdb1ff1523` |
| `work/bench_tx.fa` (**untracked** transcripts)    | 127,390 | `0f63be5fe04086da7c595996442cf89a1f70d8adc05d70ff3fd4633c978ba5a2` |

`scripts/validate_scoring_profiles.py` prints the panel's SHA-256 on every run, so a report that does
not match the row above was measured on different bytes.

## Predeclared splits (issues #102, #97)

Every efficacy figure attributed to "the panel" in `models/scoring_profile.py`, `docs/scoring.md` and
the CHANGELOG comes from `work/sirna_bench.csv`: **2,816 guides across 41 transcript accessions.**

**The split is by transcript accession, never by row**, because a guide's efficacy is correlated with
its transcript's and a row-wise split would put near-neighbours on both sides. The rule is
`sha256(accession)[0] % 2 == 0 -> development`, implemented in
`scripts/validate_scoring_profiles.py::split_of`. It is deterministic, computable from the accession
alone, independent of every efficacy value and every feature, and was fixed before any number was
computed.

| Split         | Transcripts | Guides |
| ------------- | ----------: | -----: |
| `development` |          21 |  1,605 |
| `held_out`    |          20 |  1,211 |

- `development`: `NM_002559` `NM_003342` `NM_003344` `NM_004223` `NM_004359` `NM_005339` `NM_007019`
  `NM_012864` `NM_014501` `NM_016021` `NM_017346` `NM_020548` `NM_021988` `NM_022005` `NM_025237`
  `NM_031313` `NM_053656` `U47298` `U92436` `XM_214061` `XM_371822`
- `held_out`: `AH001498` `BD135193` `J03132` `M60857` `NM_001001481` `NM_001632` `NM_002046`
  `NM_003337` `NM_003340` `NM_003345` `NM_003347` `NM_003348` `NM_003969` `NM_005450` `NM_006357`
  `NM_014176` `NM_015213` `NM_016406` `U47296` `X75932`

**Overlapping-study de-duplication.** The redistributed dataset carries two source labels, `Huesken`
(2,431 rows / 30 accessions) and `Other` (385 rows / 11 accessions, provenance unestablished). They
partition cleanly **by accession**, so no accession — and therefore no sequence — appears under both
labels or in both splits. The split therefore never places two rows from the same study *and* the
same transcript on opposite sides. `--huesken-only` drops the `Other` rows if a measurement should
cite one paper.

⚠️ **What these splits do and do not buy.** They control transcript-level overfitting **by this
branch**: a number reported on `held_out` was not available to the term-set comparison run here. They
do **not** insulate the choice of term or window — A/U(1-5) and its window were pre-declared from
published analyses of the *whole* 2,816-guide panel (ρ +0.378 at n = 2,816), which includes all 20
held-out accessions, so the held-out rows were fully available to whoever chose them. They do **not**
constitute validation. Both halves are the **same study, the same assay, the same laboratory and the
same chemistry**, so a held-out transcript says nothing about cross-lab replication — issue #97's
open question 4, which 0.7.1 does not address. No term is promoted past `experimental` on the
strength of this split, and no weight in siRNAforge was fitted on either half.

⚠️ **The tracked fixture is not the panel and its numbers must not be quoted as the panel's.**
`sirna_efficacy_subset.csv` is 180 rows across 8 accessions, which straddle both splits (development:
`NM_004223` `NM_007019` `NM_012864` `XM_214061`; held out: `NM_003969` `NM_006357` `NM_014176`
`NM_016406`). A clustered test run on it has **8** clusters against the panel's 41, and its
`target_accessibility` ρ is +0.2437 against the full set's +0.267.

# Offline orthologue mapping — `ortholog_mapping_synthetic.json`

The offline orthologue mapping for `tests/unit/test_hit_class_persistence.py`, issue #101. Fully
synthetic and unrelated to any real gene: `ENSG00000000001` → mouse `ENSMUSG00000000002` matches the
synthetic cDNA headers that test builds its transcript index from. **Rat is deliberately absent**, so
one run still covers both branches — a species the mapping resolves, and a species it states nothing
about, whose hits must fall back to the labelled symbol heuristic. Nothing here is generated; edit it
by hand alongside the fixture headers it mirrors.

It exists because cross-species classification otherwise calls Ensembl Compara: these tests each
spent ~25s failing that call behind a TLS-intercepting proxy, three attempts with 2s + 4s backoff per
route. `WorkflowConfig(ortholog_mapping_file=...)` is the offline path #101 requires, and this file is
the fixture that uses it.
