"""The scoring term registry: what each term reads, and what evidence stands behind it.

Issue #102. A weight vector says how much a term is worth; it does not say what the term *is*.
Before this module the answer lived in a `Field(description=...)` one-liner, a docstring paragraph
and an issue comment, and the three did not always agree -- ``pos1_mismatch`` held 0.05 of
``postscreen_mirna_v4`` while being exactly constant on 13,415 scored rows, and no artifact
recorded that its range was unattainable.

So each term gets a :class:`TermRecord` naming the molecule and strand it reads, the positions, the
endpoint it claims to predict, when it is applicable at all, its formula and transform, what happens
when it cannot be computed, its **declared** and **attainable** feature ranges, its evidence source
and its evidence status. :class:`ScoringProfile` then bundles a named vector set with the records
behind it, and refuses at construction to score a term whose attainable range is a single point --
which is the machine-checked form of "a constant term must not consume weight".

**Nothing here is promoted (D1).** Every weight and threshold in siRNAforge 0.7.1 is
``EvidenceStatus.EXPERIMENTAL``. This registry is documentation and validation of a *ranking*
utility: it is not calibrated potency, it is not a probability of any safety event, and no entry
below claims #98's 1.0.0 exit criteria. `validated` exists in the enum because the enum has to be
able to express a promotion that has not happened, not because anything reaches it.

This module resolves nothing and scores nothing. ``core/scoring.py::compute_composite`` remains the
only scorer and stays a pure weighted sum over a validated vector;
``ScoringWeights.vector_for`` remains the only vector selector. A profile here is a *description* of
a vector set plus its evidence, and -- for the experimental profile -- a place to write down a
candidate vector so its numbers can be reviewed and compared offline
(``scripts/validate_scoring_profiles.py``) without becoming a second scoring engine.
"""

from enum import Enum
from typing import ClassVar

from pydantic import BaseModel, ConfigDict, Field, model_validator

from sirnaforge.models.sirna import (
    DesignWeights,
    PostScreenMiRNAWeights,
    PostScreenSiRNAWeights,
    WeightVector,
)

# The panel every efficacy figure in this module was measured on. One study, one assay: it bounds
# every "measured" claim below, and cross-lab replication (#97 question 4) is not addressed by it.
HUESKEN_PANEL = (
    "Huesken et al., Nat Biotechnol 2005;23(8):995-1001, doi:10.1038/nbt1118, PMID 16025102, "
    "n = 2,816 guides / 41 transcripts, via the third-party redistribution recorded in "
    "tests/unit/data/README.md"
)

# The frozen public TP53 screen. Establishes term *behaviour* (range, variance share, saturation)
# and, by its own MEASUREMENTS.md F.8, nothing whatever about knockdown.
TP53_BASELINE = "frozen public TP53 baseline, work/baseline_0.7.1/MEASUREMENTS.md"


class EvidenceStatus(str, Enum):
    """How much a term's weight has actually earned.

    Attributes:
        EXPERIMENTAL: Declared expert prior, or an association measured on one panel. The status of
            every weight and threshold that ships in 0.7.1 (D1).
        VALIDATED: Held-out and independently replicated against a predeclared endpoint. Nothing in
            siRNAforge holds this status; the member exists so a future promotion has a spelling.
        DEPRECATED: Retained as a reported diagnostic, deliberately absent from every default
            vector. A profile may not score a deprecated term.
    """

    EXPERIMENTAL = "experimental"
    VALIDATED = "validated"
    DEPRECATED = "deprecated"


class Molecule(str, Enum):
    """Which molecule a term reads."""

    GUIDE = "guide"
    PASSENGER = "passenger"
    DUPLEX = "duplex"
    TARGET_MRNA = "target_mrna"
    TRANSCRIPTOME = "transcriptome"


class Endpoint(str, Enum):
    """What the term claims to predict, kept separate from what has been measured about it.

    ``supp_13_16`` is the reason this is its own field: it is documented as a *specificity*
    heuristic and the only thing ever measured on it is an efficacy association.
    """

    EFFICACY = "efficacy"
    RISC_LOADING = "risc_loading"
    SPECIFICITY = "specificity"
    ISOFORM_COVERAGE = "isoform_coverage"
    NONE_DECLARED = "none_declared"


class TermRecord(BaseModel):
    """Everything about one scoring term except how much it is worth.

    Frozen and ``extra="forbid"``: a field this does not have is added here rather than smuggled
    into ``notes``.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    term: str = Field(description="Feature key, exactly as compute_composite receives it")
    molecule: Molecule = Field(description="Molecule the term reads")
    strand: str = Field(description="Strand or role read, in words (e.g. 'guide 5'->3'')")
    positions: str = Field(description="Positions read, 1-based and inclusive, or 'whole' / 'n/a'")
    endpoint: Endpoint = Field(description="Endpoint the term is *claimed* to predict")
    applicability: str = Field(description="Guide length, chemistry and assay conditions under which it applies")
    formula: str = Field(description="How the raw quantity is computed")
    units: str = Field(description="Units of the raw quantity before the transform")
    transform: str = Field(description="Map from the raw quantity onto the feature scale")
    missing_value_policy: str = Field(description="What happens when the raw quantity cannot be computed")
    declared_range: tuple[float, float] = Field(
        description="Range the feature is declared over; the scorer enforces it"
    )
    attainable_range: tuple[float, float] = Field(
        description=(
            "Range actually reachable under the duplex construction and reference scope siRNAforge "
            "ships. A single point here means the term cannot rank anything, and a profile refuses "
            "to score it."
        )
    )
    evidence_source: str = Field(description="Where the evidence for this term comes from")
    evidence_status: EvidenceStatus = Field(description="How far that evidence goes")
    notes: str = Field(default="", description="Limitations that a reader would otherwise have to rediscover")

    @model_validator(mode="after")
    def ranges_are_ordered_and_nested(self) -> "TermRecord":
        """Reject an inverted range, or an attainable range outside the declared one."""
        for label, (low, high) in (("declared", self.declared_range), ("attainable", self.attainable_range)):
            if low > high:
                raise ValueError(f"{self.term}: {label}_range {low}..{high} is inverted")
        if self.attainable_range[0] < self.declared_range[0] or self.attainable_range[1] > self.declared_range[1]:
            raise ValueError(
                f"{self.term}: attainable_range {self.attainable_range} escapes declared_range {self.declared_range}"
            )
        return self

    @property
    def is_constant(self) -> bool:
        """True when the attainable range is a single point, so the term ranks nothing."""
        return self.attainable_range[0] == self.attainable_range[1]


class ScoringProfile(BaseModel):
    """A named vector set plus the term records standing behind it.

    Selects nothing and scores nothing. Its validator is the point: a profile cannot declare a
    weight on a term the registry does not describe, on a term whose attainable range is a single
    point, or on a term the registry marks deprecated.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    profile_id: str = Field(description="Stable identifier, used in reports and in the CHANGELOG")
    description: str = Field(description="What this profile is for")
    status: EvidenceStatus = Field(description="Strongest status any weight in this profile has earned")
    ships_as_default: bool = Field(description="Whether ScoringWeights' defaults are these vectors")
    vectors: tuple[WeightVector, ...] = Field(description="The vectors this profile declares")

    @model_validator(mode="after")
    def every_scored_term_is_registered_and_can_rank(self) -> "ScoringProfile":
        """Refuse a profile that pays weight to an undescribed, constant or deprecated term."""
        if self.status is EvidenceStatus.VALIDATED:
            raise ValueError(
                f"profile '{self.profile_id}' claims VALIDATED; no siRNAforge weight has been "
                "held out and independently replicated, so no profile may claim it (D1)"
            )
        for vector in self.vectors:
            weights = vector.as_mapping()
            for term, weight in weights.items():
                record = TERM_REGISTRY.get(term)
                if record is None:
                    raise ValueError(
                        f"profile '{self.profile_id}' vector '{vector.name}' scores '{term}', "
                        "which has no TermRecord; register the term or stop paying it weight"
                    )
                if weight <= 0.0:
                    continue
                if record.is_constant:
                    raise ValueError(
                        f"profile '{self.profile_id}' vector '{vector.name}' pays {weight} to "
                        f"'{term}', whose attainable range is the single point "
                        f"{record.attainable_range[0]}; a constant term ranks nothing and must "
                        "not consume weight"
                    )
                if record.evidence_status is EvidenceStatus.DEPRECATED:
                    raise ValueError(
                        f"profile '{self.profile_id}' vector '{vector.name}' pays {weight} to deprecated term '{term}'"
                    )
        return self

    def term_records(self) -> tuple[TermRecord, ...]:
        """The records for every term any of this profile's vectors scores, in registry order."""
        scored = {term for vector in self.vectors for term in vector.terms}
        return tuple(record for term, record in TERM_REGISTRY.items() if term in scored)


# --------------------------------------------------------------------------------------------
# The registry. Order is reporting order: the four shared scored terms, the biogenesis terms,
# the experimental A/U term, then the terms that are computed and deliberately never scored.
# --------------------------------------------------------------------------------------------

_RECORDS: tuple[TermRecord, ...] = (
    TermRecord(
        term="off_target",
        molecule=Molecule.TRANSCRIPTOME,
        strand="guide, aligned against every screened transcriptome",
        positions="whole guide",
        endpoint=Endpoint.SPECIFICITY,
        applicability=(
            "post-screen only; requires at least one screened species whose reference resolved and "
            "indexed. Undefined before screening, which is why design_v4 does not contain it."
        ),
        formula="exp(-genuine_off_target_count / 10.0); on-target, ortholog and repeat hits excluded",
        units="count of genuine off-target sites",
        transform="exponential decay, constant 10.0, clamped to [0, 1]",
        missing_value_policy=(
            "no substitution: an unscreened candidate has no composite_score at all, because weights "
            "are never renormalised over the terms that happened to be available"
        ),
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source=(
            f"{TP53_BASELINE}: +9.68% of composite variance at nominal 0.25 on a human-only screen "
            "(0.39x nominal), rising to +19.31% when one more species is screened"
        ),
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "Its measured influence is a property of the screening scope, not of the term: on a "
            "single-species screen the contribution compresses against its 25.0 ceiling (mean "
            "18.465, 142 distinct values on 34,861 rows). Never quote the share without the scope. "
            "The decay constant 10.0 is a declared prior and has never been fitted."
        ),
    ),
    TermRecord(
        term="target_accessibility",
        molecule=Molecule.TARGET_MRNA,
        strand="target mRNA, sense",
        positions="the 8 target bases pairing guide positions 1-8, i.e. the target site's 3' end",
        endpoint=Endpoint.EFFICACY,
        applicability=(
            "requires full transcript context; the whole transcript is folded once with RNAplfold "
            "(W=150, L=100). Not computable from a bare site, and None when the site sits closer to "
            "the transcript 5' end than the window is long."
        ),
        formula="P(the 8-mer ending at the target site's 3'-most base is unpaired), RNA.pfl_fold_up",
        units="probability",
        transform="(clamp(log10 P, -5.0, 0) + 5.0) / 5.0; the floor is fixed, never per-transcript",
        missing_value_policy=(
            "None, and the term is then absent from the feature map: no design_score is produced "
            "rather than one computed as if the site were open. Exactly 0 probability is a computed "
            "'inaccessible' and maps to 0.0, which is not the same thing."
        ),
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source=(
            f"{HUESKEN_PANEL}: Spearman rho +0.267 (n=2,779 locatable), against +0.067 for the same "
            f"statistic at the non-seed end of the same site. {TP53_BASELINE}: +29.06% of variance, "
            "4,677 distinct contributions."
        ),
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "Geometry is the thing that is easy to invert: guide G[i] pairs target T[L+1-i], so the "
            "guide seed reads the target site's 3' end. The non-seed-end control exists to catch "
            "exactly that inversion. W/L are configurable and change the feature's scale, so a run "
            "that alters them is not comparable with a default run."
        ),
    ),
    TermRecord(
        term="asymmetry",
        molecule=Molecule.DUPLEX,
        strand="both ends of the guide/passenger duplex",
        positions="the terminal 5 bp at each end (END_WINDOW_NT)",
        endpoint=Endpoint.RISC_LOADING,
        applicability=(
            "needs ViennaRNA; falls back to a proxy if the backend is unavailable, and yields NaN if "
            "even that fails. Computed from sequence alone, so available at design time."
        ),
        formula="dG(guide 5' end) - dG(guide 3' end), ViennaRNA duplex end stabilities",
        units="kcal/mol difference",
        transform="(dG_5p - dG_3p + 5.0) / 10.0, clamped to [0, 1]",
        missing_value_policy="NaN is dropped from the feature map, so the candidate gets no score on that vector",
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source=(
            f"{HUESKEN_PANEL}: Spearman rho +0.254 marginal. After A/U(1-5) is regressed out, the "
            "residual keeps a small but non-null association with efficacy -- beta +0.0533 per unit, "
            "by-transcript cluster-robust SE 0.0213, t = 2.50, p = 0.017, G = 41 clusters. "
            f"{TP53_BASELINE}: rank-1 scored term at +31.36% of variance."
        ),
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "It correlates rho +0.544 with A/U(1-5) on the panel and shares 30% of its variance with "
            "it, which is why D5 pre-declared a residual test rather than adding both blind. The "
            "residual survived, so both terms are kept -- and it survives dropping the 385 panel rows "
            "whose provenance could not be established (beta +0.0532, t = 2.48, p = 0.019, 30 "
            "clusters). But asymmetry's *independent* rank association is only rho +0.069, and the "
            "0.40 / 0.25 / 0.20 weights predate that number and were not fitted to it."
        ),
    ),
    TermRecord(
        term="gc_content",
        molecule=Molecule.GUIDE,
        strand="guide, whole",
        positions="whole",
        endpoint=Endpoint.EFFICACY,
        applicability="any guide of any supported length; pure sequence",
        formula="exp(-((GC% - 40) / 10)^2): a Gaussian preference centred on 40% GC",
        units="percent GC",
        transform="the Gaussian above; already in (0, 1]",
        missing_value_policy="not reachable: GC content is always computable",
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source=(
            f"declared expert prior; {TP53_BASELINE}: +29.91% of variance on only 6 distinct values. "
            f"{HUESKEN_PANEL}: beta +0.152 (cluster-robust t = 7.02, G = 41), so the 40%-centred "
            "preference does point the right way on measured knockdown."
        ),
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "Two things an auditor should know. Six distinct contribution values over 34,861 baseline "
            "rows means it behaves as a coarse stratifier rather than a continuous score, and its "
            "1.50x over-delivery against nominal follows from that coarseness. And the centre is a "
            "hard-coded 40%, which is **not** the midpoint of the default GC filter window "
            "(gc_min 35 / gc_max 60, midpoint 47.5): the score's optimum and the gate's optimum are "
            "different numbers, so a candidate at 47% GC sits mid-window and is scored as mildly "
            "off-optimum. Recorded here rather than changed -- filter defaults and bounds are #101's."
        ),
    ),
    TermRecord(
        term="ago_start",
        molecule=Molecule.GUIDE,
        strand="guide, 5'->3'",
        positions="1",
        endpoint=Endpoint.RISC_LOADING,
        applicability="any guide; pure sequence, so computable on every path including dirty controls",
        formula="1.0 if guide position 1 is A or U (read as RNA), else 0.0",
        units="indicator",
        transform="none; already an indicator in {0, 1}",
        missing_value_policy="not reachable for a non-empty guide",
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source=(
            f"{HUESKEN_PANEL}: Spearman rho +0.415 vs efficacy -- the strongest single feature "
            "measured on this panel, stronger than the A/U(1-5) count it is nested inside. "
            f"{TP53_BASELINE}: +27.12% of variance at nominal 0.10 (2.71x), on two distinct values."
        ),
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "Audited for #102: the range is genuinely attainable (A/U at position 1 on 50.8% of panel "
            "guides, 61.0% of baseline candidates), and a binary term over a near-balanced split is "
            "exactly why 0.10 of weight moves 27% of the ranking. It is scored only in "
            "postscreen_mirna_v4; the two siRNA vectors read no base position at all. Overlap with "
            "au_1_5 is total, not partial -- position 1 is inside 1-5 -- so the two must never both "
            "be paid (#97 question 5)."
        ),
    ),
    TermRecord(
        term="supp_13_16",
        molecule=Molecule.GUIDE,
        strand="guide, 5'->3'",
        positions="13-16",
        endpoint=Endpoint.SPECIFICITY,
        applicability=(
            "guides of at least 16 nt. Shorter guides receive a fabricated 0.5, which is unreachable "
            "in practice because the supported design range is 19-23 nt."
        ),
        formula="fraction of guide positions 13-16 that are A or U (read as RNA)",
        units="fraction",
        transform="none; already a fraction in [0, 1] on a 4-base window, so 5 distinct values",
        missing_value_policy="0.5 for a guide shorter than 16 nt -- a substituted midpoint, not a measurement",
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source=(
            f"declared expert prior on 3' supplementary pairing. {HUESKEN_PANEL}: associates with "
            "*efficacy* at rho +0.165, beta +0.148 (cluster-robust t = 8.06, G = 41), and survives "
            f"controlling for A/U(1-5) (beta +0.104, t = 5.43). {TP53_BASELINE}: +1.06% of variance "
            "at nominal 0.10, i.e. 0.11x."
        ),
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "Audited for #102, and the audit found an endpoint mismatch rather than a dead term. Its "
            "declared endpoint is specificity -- less 3' supplementary pairing -- and *nothing has "
            "ever measured that*. What is measured is an efficacy association, which is also what a "
            "plain A/U-content-tracks-duplex-stability mechanism would produce, so the number does "
            "not support the stated mechanism even though it is not null. It under-delivers 9x "
            "against nominal on the baseline while holding the same 0.10 as ago_start."
        ),
    ),
    TermRecord(
        term="pos1_mismatch",
        molecule=Molecule.DUPLEX,
        strand="guide position 1 against passenger position L",
        positions="1",
        endpoint=Endpoint.RISC_LOADING,
        applicability=(
            "**inapplicable to every duplex siRNAforge builds.** The passenger is constructed as the "
            "exact reverse complement of the guide, so guide position 1 always forms a Watson-Crick "
            "pair and the feature is always 0.0. It would become applicable only if the designer "
            "gained deliberately mismatched passengers."
        ),
        formula="1.0 if guide position 1 forms a G:U wobble or a mismatch with the passenger 3' base, else 0.0",
        units="indicator",
        transform="none",
        missing_value_policy="not reachable",
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 0.0),
        evidence_source=(
            f"{TP53_BASELINE}: min = max = mean = sd = 0.0, one distinct value over all 13,415 scored "
            "miRNA-mode candidates, 0.000% of variance, while holding 0.05 of postscreen_mirna_v4. "
            "Independently reproduced on all 2,816 Huesken panel guides with exact-reverse-complement "
            "passengers."
        ),
        evidence_status=EvidenceStatus.DEPRECATED,
        notes=(
            "Removed from postscreen_mirna_v4 in the 4.0.0 revision. It is still computed and still "
            "written to component_scores and the score_pos1_mismatch column, so a run remains "
            "auditable and the term can return to a vector if mismatched passengers are ever "
            "designed. The 0.05 it held was not reassigned by judgement: the four shared terms were "
            "restored to exactly 0.80 x postscreen_sirna_v4, which is the construction rule that "
            "vector's own docstring already declared."
        ),
    ),
    TermRecord(
        term="au_1_5",
        molecule=Molecule.GUIDE,
        strand="guide, 5'->3'",
        positions="1-5",
        endpoint=Endpoint.EFFICACY,
        applicability=(
            "any guide of at least 5 nt; pure sequence, so unlike asymmetry it needs no ViennaRNA and "
            "is computable at design time, on design_from_sequence, and on injected dirty controls."
        ),
        formula="count of A or U among guide positions 1-5 (read as RNA)",
        units="count out of 5",
        transform="count / 5, giving 6 attainable values",
        missing_value_policy="0.0 window bases are never missing for a supported 19-23 nt guide",
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source=(
            f"{HUESKEN_PANEL}: Spearman rho +0.378 vs efficacy, against rho -0.017 (p = 0.37) for the "
            "same count over guide positions 17-21, which is the null control this window's "
            "sign is checked against."
        ),
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "**Computed and reported; scored by no default vector in 0.7.1.** D5 accepted it as a "
            "scored term with the window pre-declared at 1-5, and D5's addendum -- run and reported "
            "in scripts/validate_scoring_profiles.py -- put both it and asymmetry in the vector "
            "because the asymmetry residual survived clustered inference. Wiring it into "
            "postscreen_sirna_v4 additionally requires the post-screen feature assembly in "
            "workflow.py to forward it, which #100 owns; until then the shipped default vectors are "
            "unchanged and the candidate vector lives in the experimental profile below. "
            "Two caveats that must travel with the rho. The first is what it buys: on the "
            "predeclared held-out transcripts, giving au_1_5 0.10 out of asymmetry's 0.25 moves "
            "postscreen_sirna_v4's rank correlation with efficacy from +0.335 to +0.343, and the two "
            "vectors agree on ranking at rho +0.967. The gain on data it was not chosen against is "
            "0.008 of rho -- real, and small enough that shipping it as a default off this one panel "
            "would not be honest. "
            "The second: the 1-5 window is dominated by position 1 "
            "alone (rho +0.415 there, against +0.079 / +0.067 / +0.073 at positions 3, 4 and 5), so "
            "a flat count over 1-5 is measurably *worse* than the single base it contains. The "
            "window was pre-declared to prevent rho-maximising selection and is deliberately left as "
            "declared; #97 question 3 (flat count vs position-weighted) is still open, and this is "
            "the number that makes it worth answering."
        ),
    ),
    TermRecord(
        term="empirical",
        molecule=Molecule.GUIDE,
        strand="guide, 5'->3'",
        positions="19",
        endpoint=Endpoint.EFFICACY,
        applicability="guides of at least 19 nt; pure sequence",
        formula="0.5, +0.1 for A/U at guide position 19, -0.1 for C there (simplified Reynolds)",
        units="rubric points",
        transform="clamped to [0.4, 0.6]; note this is NOT a [0, 1] feature",
        missing_value_policy="the 0.5 base score stands when no rule fires",
        declared_range=(0.4, 0.6),
        attainable_range=(0.4, 0.6),
        evidence_source=(
            f"{HUESKEN_PANEL} contradicts the surviving rule's direction: C at guide position 19 "
            "associates with MORE knockdown and A/U there with LESS. Reynolds' criteria are numbered "
            "on the sense strand, so they land on the guide 3' end, where this panel measures no "
            f"usable signal. {TP53_BASELINE}: the min_empirical_score gate produced 0 failures."
        ),
        evidence_status=EvidenceStatus.DEPRECATED,
        notes=(
            "Gate-only since #96 and in no vector. The gate is inert by default "
            "(min_empirical_score = 0.4 = the attainable minimum) which is the only reason a rule "
            "measured to point the wrong way costs nothing today. Deprecated here as a *scoring* "
            "term; whether the gate should exist at all is a filter question and #101's, not this "
            "branch's."
        ),
    ),
    TermRecord(
        term="isoform_coverage",
        molecule=Molecule.TRANSCRIPTOME,
        strand="guide against the query gene's protein-coding transcripts",
        positions="whole guide",
        endpoint=Endpoint.ISOFORM_COVERAGE,
        applicability=(
            "post-screen, and only when the guide -> source-transcript map exists. Absent on "
            "design_from_sequence and on miRNA paths."
        ),
        formula="protein-coding transcripts of the query gene containing this guide / all of them",
        units="fraction",
        transform="none",
        missing_value_policy=(
            "None when the denominator is 0, so an annotation gap does not read as zero coverage; the "
            "optional min_isoform_coverage gate skips a None"
        ),
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source="declared expert prior; no ceiling has been calibrated against truth data",
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes="Reported plus an opt-in gate; in no vector since #96, because it is legitimately None on some run shapes.",
    ),
    TermRecord(
        term="conservation",
        molecule=Molecule.TRANSCRIPTOME,
        strand="guide against non-query-species orthologues",
        positions="whole guide",
        endpoint=Endpoint.NONE_DECLARED,
        applicability="post-screen, and only when at least one non-query species was requested",
        formula="non-query species with at least one orthologue hit / non-query species requested",
        units="fraction",
        transform="none",
        missing_value_policy="None when no non-query species was requested; a single-species run has no orthologue evidence",
        declared_range=(0.0, 1.0),
        attainable_range=(0.0, 1.0),
        evidence_source="reported only; no evidence links it to any endpoint",
        evidence_status=EvidenceStatus.EXPERIMENTAL,
        notes=(
            "Reported and in no vector since #96. Cross-species conservation is desirable for one "
            "programme and irrelevant to another, so it has no single declared endpoint and is left "
            "for the reader to weigh."
        ),
    ),
)

TERM_REGISTRY: dict[str, TermRecord] = {record.term: record for record in _RECORDS}


class ExperimentalAUPostScreenSiRNAWeights(WeightVector):
    """The ``postscreen_sirna_v4`` candidate D5's decision rule produced -- **not a default**.

    D5's addendum said: regress A/U(1-5) out of ``asymmetry_score`` and test the residual against
    measured efficacy with by-transcript clustering; if the residual carries independent signal,
    both terms enter the vector. It does (beta +0.0533, cluster-robust t = 2.50, p = 0.017, G = 41),
    so this vector scores both.

    The numbers are the smallest change that admits the new term: ``asymmetry`` gives up 0.10 to
    ``au_1_5`` and nothing else moves. That is a declared prior, not a fit -- deliberately, because
    fitting five weights on the one panel that motivated the term is the single-dataset tuning this
    repository has already been burned by twice. It is here to be *compared*, by
    ``scripts/validate_scoring_profiles.py``, not to be shipped.
    """

    VECTOR_NAME: ClassVar[str] = "postscreen_sirna_v4_au_experimental"
    TERM_NAMES: ClassVar[tuple[str, ...]] = ("off_target", "target_accessibility", "asymmetry", "gc_content", "au_1_5")

    off_target: float = Field(default=0.25, ge=0, le=1, description="As postscreen_sirna_v4")
    target_accessibility: float = Field(default=0.30, ge=0, le=1, description="As postscreen_sirna_v4")
    asymmetry: float = Field(default=0.15, ge=0, le=1, description="0.25 in postscreen_sirna_v4, less au_1_5's 0.10")
    gc_content: float = Field(default=0.20, ge=0, le=1, description="As postscreen_sirna_v4")
    au_1_5: float = Field(default=0.10, ge=0, le=1, description="A/U count over guide positions 1-5, window per D5")


class ExperimentalAUPostScreenMiRNAWeights(WeightVector):
    """The ``postscreen_mirna_v4`` candidate -- **not a default**.

    ``ago_start`` is absent, and that is the answer to #97 question 5 rather than an omission.
    ``ago_start`` reads guide position 1; ``au_1_5`` reads positions 1-5, which *contains* it. On
    the panel their rank correlation is +0.476 and each still adds to the other in a joint model,
    because they are nested rather than merely correlated -- so scoring both pays twice for one
    base, and pays it in the term that already drives 27.1% of this vector's ranking variance.
    One of the two has to go; the superset stays, so the window stays as D5 declared it.

    The cost is stated rather than hidden: ``ago_start`` alone tracks efficacy *better* than the
    A/U(1-5) count that replaces it (rho +0.415 vs +0.378). This vector therefore trades measured
    marginal association for a term set in which no base is paid twice, and #97 question 3 --
    whether the 1-5 window should be position-weighted, which would recover most of that -- is left
    open on purpose rather than settled by a sweep.
    """

    VECTOR_NAME: ClassVar[str] = "postscreen_mirna_v4_au_experimental"
    TERM_NAMES: ClassVar[tuple[str, ...]] = (
        "off_target",
        "target_accessibility",
        "asymmetry",
        "gc_content",
        "au_1_5",
        "supp_13_16",
    )

    off_target: float = Field(default=0.20, ge=0, le=1, description="As the shipped 4.0.0 miRNA vector")
    target_accessibility: float = Field(default=0.24, ge=0, le=1, description="As the shipped 4.0.0 miRNA vector")
    asymmetry: float = Field(default=0.10, ge=0, le=1, description="0.20 shipped, less au_1_5's 0.10")
    gc_content: float = Field(default=0.16, ge=0, le=1, description="As the shipped 4.0.0 miRNA vector")
    au_1_5: float = Field(default=0.20, ge=0, le=1, description="ago_start's 0.10 plus 0.10 from asymmetry")
    supp_13_16: float = Field(default=0.10, ge=0, le=1, description="As the shipped 4.0.0 miRNA vector")


SHIPPED_PROFILE = ScoringProfile(
    profile_id="sirnaforge_4_0_0",
    description=(
        "The three vectors siRNAforge 0.7.1 actually scores with. Four terms in siRNA mode, six in "
        "miRNA mode, three at design time; pos1_mismatch is computed and reported but scored "
        "nowhere, because its attainable range is a single point."
    ),
    status=EvidenceStatus.EXPERIMENTAL,
    ships_as_default=True,
    vectors=(DesignWeights(), PostScreenSiRNAWeights(), PostScreenMiRNAWeights()),
)

EXPERIMENTAL_AU_PROFILE = ScoringProfile(
    profile_id="au_1_5_experimental",
    description=(
        "The post-screen vectors D5's decision rule produced, for offline comparison only. Nothing "
        "in the pipeline selects this profile: the post-screen feature assembly does not forward "
        "au_1_5, and making it do so is #100's file, not #102's."
    ),
    status=EvidenceStatus.EXPERIMENTAL,
    ships_as_default=False,
    vectors=(ExperimentalAUPostScreenSiRNAWeights(), ExperimentalAUPostScreenMiRNAWeights()),
)

PROFILES: dict[str, ScoringProfile] = {
    profile.profile_id: profile for profile in (SHIPPED_PROFILE, EXPERIMENTAL_AU_PROFILE)
}
