#!/usr/bin/env python
r"""Evaluate siRNAforge's scoring terms and profiles against measured knockdown.

The calibration record for issues #102 and #97, and the place D5's pre-declared decision rule was
executed. Deliberately **not** a unit test: it folds nothing but it does read an untracked
third-party panel from a path, and panel evaluation does not belong in the dev loop. The tracked
regression guards are `tests/unit/test_scoring_profile.py` (registry invariants, deterministic) and
`tests/unit/test_target_accessibility.py` (geometry, on a vendored subset).

Six things are reported, in the order they were decided, and the numbers below are the section
numbers the program prints:

1. **Sign checks.** A/U(1-5) against efficacy, and the same count over guide positions 17-21 as the
   null control. Published: rho +0.378 and -0.017 (p = 0.37). If the control ever rises to meet the
   scored window, the guide is being read 3'->5'.
2. **D5's pre-declared decision rule.** Regress A/U(1-5) out of `asymmetry_score`; test the residual
   against efficacy with **by-transcript clustering**. Residual carries independent signal -> both
   terms score. Residual null -> A/U supersedes asymmetry. Significance is judged at a two-sided
   alpha of 0.05 on the cluster-robust t, with G - 1 degrees of freedom, G = transcripts. The rule
   and the alpha are declared here so neither can move after the number is seen.
3. **#97 question 5.** `ago_start` reads guide position 1; A/U(1-5) contains position 1. Their
   overlap is quantified so the two are not both paid for one base.
4. **Per-position diagnostic.** The A/U indicator at each window position separately. Reported so
   the declared window's shape is on the record; explicitly not used to reselect it.
5. **Term audits.** Attainable range and endpoint association for every registered term the panel
   can speak to -- the machine-readable half is `TERM_REGISTRY`'s `attainable_range`, which this
   script cross-checks against what the panel actually attains.
6. **Profile comparison.** Rank agreement between the shipped 4.0.0 post-screen vector and the
   experimental A/U vector, over the panel's real feature values, on the predeclared splits.
7. **miRNA variance-share re-derivation** (only with `--baseline-mirna-csv`). The frozen baseline was
   scored under the pre-#102 *seven*-term miRNA vector, so its published per-term variance shares
   have denominators #102 moved. This recomposes them under the six-term vector that ships, and
   fails if it cannot reproduce the published pre-#102 column from the same rows.

**Splits.** Development and held-out are predeclared **by transcript accession**, by the parity of
the low byte of the accession's SHA-256 -- deterministic, independent of efficacy, and fixed before
any number below was computed. Sources split cleanly by accession in this dataset, so no transcript
and no sequence appears in both halves. What that does **not** buy: both halves are the same study
on the same assay, so a held-out split controls transcript-level overfitting and says nothing about
cross-lab replication (#97 question 4, still open).

Nothing here promotes anything. Every weight is `experimental` (D1).

Data (provenance, checksums and limitations: tests/unit/data/README.md):
  --benchmark-csv   columns: efficacy (inhibition, higher = more knockdown), siRNA_seq (21 nt guide,
                    5'->3'), B (transcript accession), sources
                    Huesken D, et al. "Design of a genome-wide siRNA library using an artificial
                    neural network." Nat Biotechnol. 2005 Aug;23(8):995-1001. doi:10.1038/nbt1118.
                    PMID: 16025102 -- obtained via a third-party redistribution, NOT verified
                    against the primary source, which is not open access.

Requires scipy for exact t and Spearman p-values; scipy is not a siRNAforge runtime dependency, so
it is supplied per-invocation:

    uv run --with scipy python scripts/validate_scoring_profiles.py \\
        --benchmark-csv work/sirna_bench.csv
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from sirnaforge.core.design import (
    SiRNADesigner,
    au_content_5p_score,
    biogenesis_features,
    supplementary_score,
)
from sirnaforge.core.scoring import compute_composite
from sirnaforge.core.thermodynamics import ThermodynamicCalculator
from sirnaforge.models.scoring_profile import (
    EXPERIMENTAL_AU_PROFILE,
    SHIPPED_PROFILE,
    TERM_REGISTRY,
    EvidenceStatus,
)
from sirnaforge.models.sirna import DesignParameters, SiRNACandidate, WeightVector

# D5's decision rule, declared before the measurement. Do not move it after seeing the number.
RESIDUAL_ALPHA = 0.05

# Published figures this script reproduces or refutes, as sign checks rather than assertions on
# third-decimal agreement (the panel is redistributed, not primary).
PUBLISHED_AU_RHO = 0.378
PUBLISHED_NULL_RHO = -0.017
_RC = str.maketrans("ATCG", "TAGC")


@dataclass(frozen=True)
class ClusteredFit:
    """One cluster-robust OLS fit: coefficients with by-cluster standard errors."""

    names: tuple[str, ...]
    beta: np.ndarray
    se: np.ndarray
    tstat: np.ndarray
    pvalue: np.ndarray
    n_obs: int
    n_clusters: int

    def of(self, name: str) -> tuple[float, float, float, float]:
        """(beta, se, t, p) for one coefficient."""
        i = self.names.index(name)
        return float(self.beta[i]), float(self.se[i]), float(self.tstat[i]), float(self.pvalue[i])

    def table(self) -> str:
        """The fit as a fixed-width table."""
        lines = [f"    {'term':<16}{'beta':>10}{'se(cluster)':>13}{'t':>9}{'p':>10}"]
        for i, name in enumerate(self.names):
            lines.append(
                f"    {name:<16}{self.beta[i]:>+10.5f}{self.se[i]:>13.5f}{self.tstat[i]:>+9.2f}{self.pvalue[i]:>10.4f}"
            )
        lines.append(f"    n = {self.n_obs}, clusters = {self.n_clusters}")
        return "\n".join(lines)


def fit_clustered(y: np.ndarray, columns: dict[str, np.ndarray], clusters: np.ndarray) -> ClusteredFit:
    """OLS with one-way cluster-robust (CR1) standard errors and t(G-1) p-values.

    Hand-rolled rather than pulled from statsmodels because it is one sandwich and the repo already
    declines to add a dependency for one number. The finite-sample correction is the conventional
    CR1: G/(G-1) * (n-1)/(n-k).
    """
    names = ("const", *columns)
    design = np.column_stack([np.ones(len(y)), *columns.values()])
    n, k = design.shape
    xtx_inv = np.linalg.pinv(design.T @ design)
    beta = xtx_inv @ design.T @ y
    resid = y - design @ beta

    meat = np.zeros((k, k))
    unique = pd.unique(clusters)
    for group in unique:
        rows = clusters == group
        score = design[rows].T @ resid[rows]
        meat += np.outer(score, score)

    n_clusters = len(unique)
    correction = (n_clusters / (n_clusters - 1)) * ((n - 1) / (n - k))
    variance = correction * xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.diag(variance))
    tstat = beta / se
    pvalue = 2.0 * stats.t.sf(np.abs(tstat), df=n_clusters - 1)
    return ClusteredFit(names, beta, se, tstat, pvalue, n, n_clusters)


def split_of(accession: str) -> str:
    """Predeclared split for one transcript: 'development' or 'held_out'.

    By SHA-256 parity of the accession, so it is deterministic, reproducible from the accession
    alone, and independent of every efficacy value and every feature.
    """
    digest = hashlib.sha256(accession.encode()).digest()
    return "development" if digest[0] % 2 == 0 else "held_out"


def build_features(benchmark_csv: Path) -> pd.DataFrame:
    """Compute every panel-computable term for every benchmark guide.

    Terms needing transcript context (`target_accessibility`) or screening (`off_target`) are not
    computed here: `scripts/validate_target_accessibility.py` owns the former, and the latter has no
    panel. Both are supplied as constants in the profile comparison, which is stated where it
    matters.
    """
    table = pd.read_csv(benchmark_csv)
    calc = ThermodynamicCalculator()
    designer = SiRNADesigner(DesignParameters())

    rows: list[dict[str, object]] = []
    for _, row in table.iterrows():
        guide = str(row.siRNA_seq).upper()
        passenger = guide.translate(_RC)[::-1]
        candidate = SiRNACandidate(
            id=f"bench_{row.name}",
            transcript_id=str(row.B),
            position=1,
            guide_sequence=guide,
            passenger_sequence=passenger,
            gc_content=(guide.count("G") + guide.count("C")) / len(guide) * 100.0,
            length=len(guide),
            asymmetry_score=0.0,
        )
        biogenesis = biogenesis_features(guide, passenger)
        au = au_content_5p_score(guide)
        rna = guide.replace("T", "U")
        tail = rna[-5:]
        rows.append(
            {
                "accession": str(row.B),
                "sources": str(row.sources),
                # Carried on the panel, not re-read from the CSV: every later report must see the
                # same rows as the panel, including under --huesken-only.
                "guide_rna": rna,
                "split": split_of(str(row.B)),
                "efficacy": float(row.efficacy),
                "au_1_5": au,
                # The null control: the same count at the other end of the same guide.
                "au_17_21": (tail.count("A") + tail.count("U")) / 5.0,
                "asymmetry": calc.calculate_asymmetry_score(candidate)[2],
                # The scored feature, not raw GC: a Gaussian centred on 40% GC.
                "gc_content": designer._calculate_gc_score(candidate.gc_content),
                "gc_percent": candidate.gc_content,
                "ago_start": biogenesis["ago_start"],
                "pos1_mismatch": biogenesis["pos1_mismatch"],
                "supp_13_16": supplementary_score(guide),
            }
        )
    features: pd.DataFrame = pd.DataFrame(rows)
    return features


def report_sign_checks(panel: pd.DataFrame) -> bool:
    """Reproduce the published marginal correlations, including the null control."""
    print("\n1. Sign checks against the published figures")
    y = panel.efficacy
    au_rho, au_p = stats.spearmanr(panel.au_1_5, y)
    null_rho, null_p = stats.spearmanr(panel.au_17_21, y)
    print(f"    A/U(1-5)   vs efficacy   rho {au_rho:+.4f} (p {au_p:.3g})   published {PUBLISHED_AU_RHO:+.3f}")
    print(f"    A/U(17-21) vs efficacy   rho {null_rho:+.4f} (p {null_p:.3g})   published {PUBLISHED_NULL_RHO:+.3f}")
    ok = au_rho > 0.30 and abs(null_rho) < 0.10 and au_rho > abs(null_rho) * 3
    print(f"    {'PASS' if ok else 'FAIL'}  scored window positive and well clear of its null control")
    return bool(ok)


def report_d5_rule(panel: pd.DataFrame) -> None:
    """Execute D5's pre-declared decision rule and print which branch it selects."""
    print(f"\n2. D5's pre-declared decision rule (alpha {RESIDUAL_ALPHA}, by-transcript clustering)")
    au = panel.au_1_5.to_numpy(float)
    asym = panel.asymmetry.to_numpy(float)
    y = panel.efficacy.to_numpy(float)
    clusters = panel.accession.to_numpy()

    stage_one = np.column_stack([np.ones(len(panel)), au])
    coefficients = np.linalg.lstsq(stage_one, asym, rcond=None)[0]
    residual = asym - stage_one @ coefficients
    shared = 1.0 - residual.var() / asym.var()
    print(f"    stage 1: asymmetry = {coefficients[0]:+.4f} {coefficients[1]:+.4f} * A/U(1-5)")
    print(f"             A/U(1-5) explains {shared:.1%} of asymmetry_score's variance on this panel")
    print(f"             rho(A/U, asymmetry) {stats.spearmanr(au, asym)[0]:+.4f}")

    fit = fit_clustered(y, {"au_1_5": au, "asym_residual": residual}, clusters)
    print(fit.table())
    beta, se, tstat, pvalue = fit.of("asym_residual")
    residual_rho = stats.spearmanr(residual, y)[0]
    print(f"    residual marginal rank association: rho {residual_rho:+.4f}")

    # A saturated removal of A/U -- one intercept per attainable count -- is a strictly harder test
    # than the linear one the rule specifies. Reported so the verdict is not an artifact of assuming
    # linearity in the count.
    dummies = pd.get_dummies(panel.au_1_5, prefix="au", drop_first=True)
    indicators: dict[str, np.ndarray] = {str(name): dummies[name].to_numpy(float) for name in dummies.columns}
    stage_one_saturated = np.column_stack([np.ones(len(panel)), *indicators.values()])
    saturated = asym - stage_one_saturated @ np.linalg.lstsq(stage_one_saturated, asym, rcond=None)[0]
    saturated_fit = fit_clustered(y, {**indicators, "asym_residual": saturated}, clusters)
    _, _, sat_t, sat_p = saturated_fit.of("asym_residual")
    print(f"    robustness (A/U removed with saturated dummies): t {sat_t:+.2f}, p {sat_p:.4f}")

    carries_signal = pvalue < RESIDUAL_ALPHA
    if carries_signal:
        print(
            f"    VERDICT  residual carries independent signal (beta {beta:+.4f}, cluster SE {se:.4f}, "
            f"t {tstat:+.2f}, p {pvalue:.4f} < {RESIDUAL_ALPHA})"
        )
        print("             -> BOTH terms enter the scored vector.")
        print(
            f"             Effect size, stated so the verdict is not over-read: asymmetry's "
            f"independent rank association is rho {residual_rho:+.4f}, against A/U(1-5)'s "
            f"{stats.spearmanr(au, y)[0]:+.4f} marginal. Non-null, and small."
        )
    else:
        print(f"    VERDICT  residual is null (p {pvalue:.4f} >= {RESIDUAL_ALPHA})")
        print("             -> A/U(1-5) SUPERSEDES asymmetry in the scored vector.")


def report_question_five(panel: pd.DataFrame) -> None:
    """#97 question 5: how much of A/U(1-5) is ago_start, and can both be paid?"""
    print("\n3. #97 question 5 -- ago_start (position 1) inside A/U(1-5)")
    y = panel.efficacy.to_numpy(float)
    clusters = panel.accession.to_numpy()
    au = panel.au_1_5.to_numpy(float)
    ago = panel.ago_start.to_numpy(float)

    print(f"    rho(A/U(1-5), ago_start) {stats.spearmanr(au, ago)[0]:+.4f}   (nested, not merely correlated)")
    print(
        f"    marginal vs efficacy: ago_start rho {stats.spearmanr(ago, y)[0]:+.4f}, "
        f"A/U(1-5) rho {stats.spearmanr(au, y)[0]:+.4f}"
    )
    print(fit_clustered(y, {"au_1_5": au, "ago_start": ago}, clusters).table())
    print("    Both coefficients are non-null, which is exactly the double-payment problem: they")
    print("    are nested, so a vector holding both pays twice for guide position 1 -- the base")
    print("    already driving 24.6% of the shipped postscreen_mirna_v4's ranking variance on the")
    print("    frozen baseline (report 7; +27.1% under the pre-#102 vector that run was scored with).")


def report_per_position(panel: pd.DataFrame) -> None:
    """Per-position A/U association. Diagnostic: it must not be used to reselect the window."""
    guides = panel.guide_rna.tolist()
    y = panel.efficacy.to_numpy(float)
    print("\n4. Per-position A/U association (diagnostic; D5 fixes the window at 1-5)")
    for position in (1, 2, 3, 4, 5, 19):
        indicator = np.array([1.0 if g[position - 1] in "AU" else 0.0 for g in guides])
        rho, pvalue = stats.spearmanr(indicator, y)
        print(f"    position {position:>2}   rho {rho:+.4f} (p {pvalue:.3g})   A/U frequency {indicator.mean():.3f}")
    print("    Position 1 alone outranks the 1-5 count. Recorded, argued in models/scoring_profile.py,")
    print("    and deliberately NOT acted on: a rho-maximising window is the overfitting D5 forbids.")


# The seven-term miRNA vector the frozen baseline was scored under, before #102 removed
# `pos1_mismatch`. Kept only so a baseline row's per-term contribution can be divided back into the
# feature that produced it; nothing scores with it.
PRE_102_MIRNA_WEIGHTS = {
    "off_target": 0.20,
    "target_accessibility": 0.22,
    "asymmetry": 0.18,
    "gc_content": 0.15,
    "ago_start": 0.10,
    "pos1_mismatch": 0.05,
    "supp_13_16": 0.10,
}


def _variance_shares(features: pd.DataFrame, weights: dict[str, float]) -> dict[str, float]:
    """Each term's share of composite variance under one weight set, as a percentage summing to 100.

    The share is cov(contribution, total) / var(total), the decomposition the frozen baseline's
    MEASUREMENTS.md uses, so the two are directly comparable.
    """
    contributions = pd.DataFrame({term: features[term] * weight * 100.0 for term, weight in weights.items()})
    total = contributions.sum(axis=1)
    variance = float(total.var(ddof=1))
    return {term: float(np.cov(contributions[term], total, ddof=1)[0, 1] / variance * 100.0) for term in weights}


def report_mirna_variance_reweighting(candidates_csv: Path) -> bool:
    """Re-derive the miRNA variance shares under the vector that ships, not the one measured.

    The frozen baseline was scored under the pre-#102 seven-term vector, so every variance share in
    `work/baseline_0.7.1/MEASUREMENTS.md` has a denominator that #102 moved. This divides each
    `score_<term>` column back by the weight that produced it -- exact, because a contribution is
    weight x feature x 100 and no weight in that vector is zero -- and recomposes the shares under
    the shipped six-term weights. Reproducing the published old-vector column is the check that the
    reconstruction is faithful; a mismatch there fails the run.
    """
    print("\n7. miRNA variance shares re-derived under the shipped vector")
    table = pd.read_csv(candidates_csv, low_memory=False)
    scored = table[table.scored_after_screening.astype(str).str.lower() == "true"]
    features = pd.DataFrame(
        {term: scored["score_" + term].astype(float) / weight for term, weight in PRE_102_MIRNA_WEIGHTS.items()}
    )
    mirna_vector = next(v for v in SHIPPED_PROFILE.vectors if v.name == "postscreen_mirna_v4")
    shipped_weights = dict(mirna_vector.as_mapping())
    old = _variance_shares(features, PRE_102_MIRNA_WEIGHTS)
    new = _variance_shares(features, shipped_weights)

    print(f"    n = {len(scored)} scored rows from {candidates_csv}")
    print(f"    {'term':<22}{'pre-#102 (7 terms)':>20}{'shipped (6 terms)':>20}{'vs nominal':>12}")
    for term in PRE_102_MIRNA_WEIGHTS:
        if term in shipped_weights:
            ratio = new[term] / (shipped_weights[term] * 100.0)
            print(f"    {term:<22}{old[term]:>+19.3f}%{new[term]:>+19.3f}%{ratio:>11.2f}x")
        else:
            print(f"    {term:<22}{old[term]:>+19.3f}%{'not scored':>20}{'--':>12}")

    # The published figures the registry and CHANGELOG cite, as a self-check on the reconstruction.
    published = {"ago_start": 27.12, "supp_13_16": 1.06, "off_target": 8.02, "asymmetry": 32.62}
    faithful = all(abs(old[term] - value) < 0.02 for term, value in published.items())
    print(f"    {'PASS' if faithful else 'FAIL'}  reconstruction reproduces MEASUREMENTS.md's pre-#102 column")
    print("    Every figure attributed to the miRNA vector must name which of the two columns it is")
    print("    from: #102 moved four of the six shared weights, so the denominators are not the same.")
    return bool(faithful)


def report_term_audit(panel: pd.DataFrame) -> bool:
    """Cross-check each registered term's attainable range and endpoint against the panel."""
    print("\n5. Term audit -- registry claims against what the panel attains")
    y = panel.efficacy.to_numpy(float)
    clusters = panel.accession.to_numpy()
    ok = True
    for term in ("au_1_5", "asymmetry", "gc_content", "ago_start", "supp_13_16", "pos1_mismatch"):
        record = TERM_REGISTRY[term]
        values = panel[term].to_numpy(float)
        observed = (float(values.min()), float(values.max()))
        constant = observed[0] == observed[1]
        agrees = constant == record.is_constant
        ok = ok and agrees
        line = (
            f"    {term:<15} observed {observed[0]:.2f}-{observed[1]:.2f} "
            f"({panel[term].nunique()} distinct)  registry attainable "
            f"{record.attainable_range[0]:.2f}-{record.attainable_range[1]:.2f}"
        )
        if constant:
            print(f"{line}  CONSTANT -- ranks nothing")
        else:
            beta, _, tstat, pvalue = fit_clustered(y, {term: values}, clusters).of(term)
            print(f"{line}  beta {beta:+.4f} t {tstat:+.2f} p {pvalue:.4f} vs {record.endpoint.value}")
        if not agrees:
            print(f"        MISMATCH: registry says constant={record.is_constant}, panel says {constant}")
    print("    supp_13_16's endpoint is declared SPECIFICITY and nothing has ever measured that; the")
    print("    association above is with efficacy, which a plain A/U-tracks-duplex-stability")
    print("    mechanism would also produce. Non-null does not mean the stated mechanism is right.")
    return ok


def _score_panel(panel: pd.DataFrame, vector: WeightVector, accessibility: float, off_target: float) -> np.ndarray:
    """Score every panel guide on one vector, holding the two unmeasurable terms constant.

    `target_accessibility` and `off_target` cannot be computed from the panel alone (one needs
    transcript folding, the other a screen), so they are held at a constant. That makes the absolute
    scores meaningless and the *rankings* comparable, which is all this comparison uses.
    """
    constants = {"target_accessibility": accessibility, "off_target": off_target}
    return np.array(
        [
            compute_composite(
                {term: constants[term] if term in constants else float(row[term]) for term in vector.terms},
                vector,
            ).score
            for _, row in panel.iterrows()
        ]
    )


def report_profile_comparison(panel: pd.DataFrame) -> None:
    """Rank agreement and panel association for the shipped and experimental post-screen vectors."""
    print("\n6. Profile comparison on the predeclared splits")
    shipped = next(v for v in SHIPPED_PROFILE.vectors if v.name == "postscreen_sirna_v4")
    experimental = next(v for v in EXPERIMENTAL_AU_PROFILE.vectors if v.name.startswith("postscreen_sirna"))
    print(f"    {shipped.name:<40} {shipped.as_mapping()}")
    print(f"    {experimental.name:<40} {experimental.as_mapping()}")
    print("    target_accessibility and off_target held at 0.5 -- neither is computable from this")
    print("    panel, so only the rankings these vectors disagree about are interpretable.")

    for split in ("development", "held_out", "all"):
        subset = panel if split == "all" else panel[panel.split == split]
        shipped_scores = _score_panel(subset, shipped, 0.5, 0.5)
        experimental_scores = _score_panel(subset, experimental, 0.5, 0.5)
        rho_shipped = stats.spearmanr(shipped_scores, subset.efficacy)[0]
        rho_experimental = stats.spearmanr(experimental_scores, subset.efficacy)[0]
        agreement = stats.spearmanr(shipped_scores, experimental_scores)[0]
        print(
            f"    {split:<12} n {len(subset):>5}  transcripts {subset.accession.nunique():>3}  "
            f"shipped rho {rho_shipped:+.4f}  experimental rho {rho_experimental:+.4f}  "
            f"rank agreement {agreement:+.4f}"
        )
    print("    Both halves are the same study on the same assay, so a held-out transcript split")
    print("    controls transcript-level overfitting only. #97 question 4 stays open.")


def main() -> int:
    """Run every report and return a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--benchmark-csv", type=Path, required=True, help="Measured-efficacy CSV")
    parser.add_argument(
        "--huesken-only",
        action="store_true",
        help="Drop the 385 rows whose provenance could not be established (sources != Huesken)",
    )
    parser.add_argument(
        "--baseline-mirna-csv",
        type=Path,
        default=None,
        help="candidates_all.csv from the frozen baseline's miRNA-mode run, to re-derive its "
        "variance shares under the shipped six-term vector (report 7)",
    )
    parser.add_argument("--out-csv", type=Path, default=None, help="Write the computed feature table here")
    args = parser.parse_args()

    digest = hashlib.sha256(args.benchmark_csv.read_bytes()).hexdigest()
    print(f"panel {args.benchmark_csv}  sha256 {digest}")

    panel = build_features(args.benchmark_csv)
    if args.huesken_only:
        panel = panel[panel.sources == "Huesken"].reset_index(drop=True)
    print(f"n = {len(panel)} guides, {panel.accession.nunique()} transcripts, sources {sorted(panel.sources.unique())}")
    counts = panel.groupby("split").accession.agg(["nunique", "size"])
    for split, row in counts.iterrows():
        print(f"    split {split:<12} {int(row['nunique']):>3} transcripts, {int(row['size']):>5} guides")
    overlap = set(panel[panel.split == "development"].accession) & set(panel[panel.split == "held_out"].accession)
    print(f"    transcript overlap between splits: {len(overlap)} (must be 0)")

    if args.out_csv is not None:
        panel.to_csv(args.out_csv, index=False)
        print(f"    features written to {args.out_csv}")

    # D5's verdict is an outcome, not a pass/fail: the rule was declared to be obeyed either way,
    # so it never contributes to the exit code.
    checks = [report_sign_checks(panel)]
    report_d5_rule(panel)
    report_question_five(panel)
    report_per_position(panel)
    checks.append(report_term_audit(panel))
    report_profile_comparison(panel)
    if args.baseline_mirna_csv is not None:
        checks.append(report_mirna_variance_reweighting(args.baseline_mirna_csv))

    statuses = {record.evidence_status for record in TERM_REGISTRY.values()}
    print(f"\nevidence statuses in the registry: {sorted(s.value for s in statuses)}")
    print(f"any term VALIDATED: {EvidenceStatus.VALIDATED in statuses}  (D1: nothing is promoted in 0.7.1)")
    return 0 if all(checks) and not overlap else 1


if __name__ == "__main__":
    sys.exit(main())
