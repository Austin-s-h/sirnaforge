"""Core siRNA design algorithms and functionality."""

import hashlib
import logging
import math
import sys
import time
from collections.abc import Mapping

import Bio
from Bio import SeqIO
from Bio.Seq import Seq

from sirnaforge import __version__
from sirnaforge.config.run_policy import FILTER_SPEC_BY_ID
from sirnaforge.core.repeat_detection import RepeatObservation, normalize_guide_sequence
from sirnaforge.core.scoring import ScoringError, compute_composite, target_accessibility_sub_score
from sirnaforge.core.thermodynamics import (
    SEED_ANCHORED_WINDOW_NT,
    TargetAccessibilityProfile,
    TargetSiteAccessibility,
    ThermodynamicCalculator,
)
from sirnaforge.models.policy import FilterAction
from sirnaforge.models.sirna import (
    EMPIRICAL_SCORE_MAX,
    EMPIRICAL_SCORE_MIN,
    DesignParameters,
    DesignResult,
    SiRNACandidate,
    ranking_score,
)
from sirnaforge.models.sirna import SiRNACandidate as _ModelCandidate

logger = logging.getLogger(__name__)

# Duplex ΔG scales with duplex length, so it is normalised per nucleotide before
# scoring. Window taken from the ViennaRNA range measured for canonical siRNA
# duplexes at 37 °C: strong -> 1.0, weak -> 0.0.
DUPLEX_DG_PER_NT_STRONG = -2.1
DUPLEX_DG_PER_NT_WEAK = -1.4

# Sanitized transcript ids embedded in candidate ids are truncated to this length. Beyond it,
# the last ID_DIGEST_LEN chars are replaced by a deterministic digest of the FULL sanitized id
# (hashlib, not hash(), which is PYTHONHASHSEED-dependent) so two ids that share a long common
# prefix cannot collide after truncation.
TRANSCRIPT_ID_MAX_LEN = 24
ID_DIGEST_LEN = 8


def _as_rna(sequence: str) -> str:
    """Read a stored (DNA) sequence as RNA so T and U compare equal."""
    return sequence.upper().replace("T", "U")


class SiRNADesigner:
    """Main siRNA design engine following the algorithm specification."""

    def __init__(
        self,
        parameters: DesignParameters,
        *,
        filter_actions: Mapping[str, FilterAction] | None = None,
    ) -> None:
        """Initialize designer with given parameters.

        Args:
            parameters: Thresholds and scoring settings.
            filter_actions: Per-filter actions from the resolved run policy. Keyword-only with a
                default because the declared defaults are the right answer when a caller has no
                policy in hand -- ``DesignParameters`` carries thresholds and says nothing about
                whether exceeding one rejects a candidate, so without this the designer could only
                ever reject and ``FilterAction.WARN`` would be inert.
        """
        self.parameters = parameters
        self._filter_actions: Mapping[str, FilterAction] = filter_actions or {}
        self.last_guide_to_transcripts: dict[str, set[str]] | None = None

    def _action_for(self, filter_id: str) -> FilterAction:
        """The action in force for one filter: the caller's policy, else the declared default."""
        override = self._filter_actions.get(filter_id)
        if override is not None:
            return override
        spec = FILTER_SPEC_BY_ID.get(filter_id)
        return spec.default_action if spec is not None else FilterAction.FAIL

    def design_from_file(self, input_file: str) -> DesignResult:
        """Design siRNAs from input FASTA file."""
        start_time = time.perf_counter()

        # Parse input sequences
        sequences = list(SeqIO.parse(input_file, "fasta"))
        if not sequences:
            raise ValueError(f"No sequences found in {input_file}")

        all_candidates: list[SiRNACandidate] = []
        rejected_pool: list[SiRNACandidate] = []
        # Map guide_sequence -> set of transcript_ids where it appears
        guide_to_transcripts: dict[str, set[str]] = {}

        # Process each sequence
        for seq_record in sequences:
            transcript_id = seq_record.id
            sequence = str(seq_record.seq).upper()

            # Generate candidates for this sequence
            candidates, rejected = self._enumerate_candidates(sequence, transcript_id)
            rejected_pool.extend(rejected)

            # Apply filters
            filtered_candidates = self._apply_filters(candidates)

            # Score candidates. The transcript is passed through so target_accessibility can be
            # folded from the real mRNA context rather than left inactive.
            scored_candidates = self._score_candidates(filtered_candidates, sequence)

            # Track which transcripts each guide appears in
            for c in scored_candidates:
                guide_to_transcripts.setdefault(c.guide_sequence, set()).add(c.transcript_id)

            all_candidates.extend(scored_candidates)

        # Sort by design_score (descending); composite_score does not exist until screening.
        all_candidates.sort(key=ranking_score, reverse=True)

        # Get top candidates only from those passing filters; fallback to all if none pass
        passing = [
            c
            for c in all_candidates
            if (c.passes_filters is True)
            or (hasattr(_ModelCandidate, "FilterStatus") and c.passes_filters == _ModelCandidate.FilterStatus.PASS)
        ]
        top_candidates = (passing or all_candidates)[: self.parameters.top_n]

        processing_time = max(0.0, time.perf_counter() - start_time)  # Ensure non-negative
        # Compute transcript hit metrics for each candidate (how many input transcripts contain the guide)
        total_seqs = len(sequences)
        for c in all_candidates:
            hits = len(guide_to_transcripts.get(c.guide_sequence, {c.transcript_id}))
            c.transcript_hit_count = hits
            c.transcript_hit_fraction = hits / total_seqs if total_seqs > 0 else 0.0

        # Expose guide_to_transcripts mapping for workflow's protein-coding coverage computation
        self.last_guide_to_transcripts = guide_to_transcripts

        return DesignResult(
            input_file=input_file,
            parameters=self.parameters,
            candidates=all_candidates,
            top_candidates=top_candidates,
            total_sequences=len(sequences),
            total_candidates=len(all_candidates),
            filtered_candidates=len(
                [
                    c
                    for c in all_candidates
                    if (c.passes_filters is True)
                    or (
                        hasattr(_ModelCandidate, "FilterStatus")
                        and c.passes_filters == _ModelCandidate.FilterStatus.PASS
                    )
                ]
            ),
            processing_time=processing_time,
            tool_versions=self._get_tool_versions(),
            rejected_candidates=rejected_pool,
        )

    def design_from_sequence(self, sequence: str, transcript_id: str = "seq1") -> DesignResult:
        """Design siRNAs from a single sequence."""
        start_time = time.perf_counter()

        sequence = sequence.upper()

        # Generate candidates
        candidates, rejected = self._enumerate_candidates(sequence, transcript_id)

        # Apply filters
        filtered_candidates = self._apply_filters(candidates)

        # Score candidates, with the transcript in scope for target_accessibility
        scored_candidates = self._score_candidates(filtered_candidates, sequence)

        # Sort by design_score (descending); composite_score does not exist until screening.
        scored_candidates.sort(key=ranking_score, reverse=True)

        # Get top candidates only from those passing filters; fallback to all if none pass
        passing = [
            c
            for c in scored_candidates
            if (c.passes_filters is True)
            or (hasattr(_ModelCandidate, "FilterStatus") and c.passes_filters == _ModelCandidate.FilterStatus.PASS)
        ]
        top_candidates = (passing or scored_candidates)[: self.parameters.top_n]

        processing_time = max(0.0, time.perf_counter() - start_time)
        # For single-sequence runs, transcript hit metrics are trivial (hits=1, fraction=1.0)
        for c in scored_candidates:
            c.transcript_hit_count = 1
            c.transcript_hit_fraction = 1.0
        return DesignResult(
            input_file="<direct_input>",
            parameters=self.parameters,
            candidates=scored_candidates,
            top_candidates=top_candidates,
            total_sequences=1,
            total_candidates=len(scored_candidates),
            filtered_candidates=len(
                [
                    c
                    for c in scored_candidates
                    if (c.passes_filters is True)
                    or (
                        hasattr(_ModelCandidate, "FilterStatus")
                        and c.passes_filters == _ModelCandidate.FilterStatus.PASS
                    )
                ]
            ),
            processing_time=processing_time,
            tool_versions=self._get_tool_versions(),
            rejected_candidates=rejected,
        )

    def _enumerate_candidates(
        self, sequence: str, transcript_id: str
    ) -> tuple[list[SiRNACandidate], list[SiRNACandidate]]:
        """Enumerate all possible siRNA candidates and record those failing early filters."""
        candidates: list[SiRNACandidate] = []
        rejected: list[SiRNACandidate] = []
        sirna_length = self.parameters.sirna_length
        filters = self.parameters.filters

        # Slide window across sequence
        for i in range(len(sequence) - sirna_length + 1):
            target_seq = sequence[i : i + sirna_length]

            # Generate guide (antisense) and passenger (sense) sequences
            guide_seq = str(Seq(target_seq).reverse_complement())
            passenger_seq = target_seq

            # Early filtering for computational efficiency
            gc_content = self._calculate_gc_content(guide_seq)
            poly_run = self._longest_poly_run(guide_seq)
            fail_reason: SiRNACandidate.FilterStatus | None = None
            if not (filters.gc_min <= gc_content <= filters.gc_max):
                fail_reason = SiRNACandidate.FilterStatus.GC_OUT_OF_RANGE
            elif poly_run > filters.max_poly_runs:
                fail_reason = SiRNACandidate.FilterStatus.POLY_RUNS

            # Create candidate ID with project moniker and sanitized transcript id
            # Format: SIRNAF_<TRANSCRIPT>_<start>_<end>
            # Sanitize transcript_id: keep alphanumerics and underscore, replace others with '-'
            safe_tid = "".join([c if (c.isalnum() or c == "_") else "-" for c in transcript_id])
            # Truncate long transcript ids, disambiguating with a digest of the full id so two
            # ids sharing a >=24-char prefix (e.g. de-novo assembly isoforms) cannot collide.
            if len(safe_tid) > TRANSCRIPT_ID_MAX_LEN:
                digest = hashlib.sha256(safe_tid.encode()).hexdigest()[:ID_DIGEST_LEN]
                keep = TRANSCRIPT_ID_MAX_LEN - ID_DIGEST_LEN - 1
                safe_tid = f"{safe_tid[:keep]}-{digest}"
            candidate_id = f"SIRNAF_{safe_tid}_{i + 1}_{i + sirna_length}"

            candidate = SiRNACandidate(
                id=candidate_id,
                transcript_id=transcript_id,
                position=i + 1,  # 1-based
                guide_sequence=guide_seq,
                passenger_sequence=passenger_seq,
                gc_content=gc_content,
                length=sirna_length,
                asymmetry_score=0.0,  # Will be calculated in scoring
            )

            self._record_enumeration_verdicts(candidate, gc_content, poly_run)

            if fail_reason is not None:
                candidate.passes_filters = fail_reason
                issues = list(candidate.quality_issues or [])
                label = fail_reason.value if hasattr(fail_reason, "value") else str(fail_reason)
                issues.append(label)
                candidate.quality_issues = issues
                rejected.append(candidate)
                continue

            candidates.append(candidate)

        return candidates, rejected

    def _record_enumeration_verdicts(self, candidate: SiRNACandidate, gc_content: float, poly_run: int) -> None:
        """Record what the three enumeration-time gates observed on this candidate.

        These gates decide during enumeration, and used to set ``passes_filters`` directly without
        recording a verdict. That left ``gc_content_min``, ``gc_content_max`` and ``max_poly_runs``
        exporting ``not_evaluated`` and an empty observed value on every row of every run -- so a
        consumer could not tell a gate that passed from one that never ran, and no candidate could be
        shown as clean. ``passes_filters`` is still set by the caller, which owns the first-label rule.
        """
        filters = self.parameters.filters
        for filter_id, observed, passed in (
            ("gc_content_min", gc_content, gc_content >= filters.gc_min),
            ("gc_content_max", gc_content, gc_content <= filters.gc_max),
            ("max_poly_runs", float(poly_run), poly_run <= filters.max_poly_runs),
        ):
            candidate.record_filter_verdict(
                filter_id,
                observed=observed,
                passed=passed,
                action=self._action_for(filter_id),
                status=_ModelCandidate.FilterStatus.GC_OUT_OF_RANGE
                if filter_id.startswith("gc_content")
                else _ModelCandidate.FilterStatus.POLY_RUNS,
            )

    def _apply_filters(self, candidates: list[SiRNACandidate]) -> list[SiRNACandidate]:
        """Apply remaining filters (early GC and poly-run filtering already done in enumeration)."""
        filtered = []

        for candidate in candidates:
            issues: list[str] = []
            status: bool | _ModelCandidate.FilterStatus = True

            # Note: GC content and poly-run filtering already done in _enumerate_candidates
            # This is mainly for any additional filters or post-processing

            # Update candidate with filter results
            candidate.passes_filters = status
            candidate.quality_issues = issues

            filtered.append(candidate)

        return filtered

    def _build_accessibility_profile(self, transcript_sequence: str | None) -> TargetAccessibilityProfile | None:
        """Fold the transcript once for target-site accessibility, or return None if impossible.

        One RNAplfold pass per transcript is both correct and cheaper than the per-candidate guide
        folds it replaced. Returning None (rather than raising or substituting a default) is what
        makes the term inactive for callers that score candidates without transcript context.
        """
        if not transcript_sequence:
            return None

        config = self.parameters.target_accessibility
        try:
            return TargetAccessibilityProfile.fold(
                transcript_sequence,
                window_size=config.window_size,
                max_bp_span=config.max_bp_span,
                u_max=max(self.parameters.sirna_length, SEED_ANCHORED_WINDOW_NT),
            )
        except Exception as exc:
            # A failed fold must leave the term inactive, never score as accessible.
            logger.warning(f"RNAplfold failed on the transcript ({exc}); target_accessibility is inactive")
            return None

    def _store_target_accessibility(
        self, candidate: SiRNACandidate, profile: TargetAccessibilityProfile | None
    ) -> float | None:
        """Record a candidate's target-site opening probabilities and return the scored feature.

        ``candidate.position`` is the 1-based start of the *target site* on the transcript (the
        passenger strand's coordinates), which is what the profile is indexed by.
        """
        site: TargetSiteAccessibility = (
            profile.site_accessibility(candidate.position - 1, candidate.length)
            if profile is not None
            else TargetSiteAccessibility(None, None, None)
        )
        candidate.target_accessibility_p = site.seed_end_8mer
        candidate.target_accessibility_p_17mer = site.seed_anchored_17mer
        candidate.target_accessibility_p_site = site.whole_site

        return target_accessibility_sub_score(
            site.seed_end_8mer, log_floor=self.parameters.target_accessibility.log_floor
        )

    def _score_candidates(
        self, candidates: list[SiRNACandidate], transcript_sequence: str | None = None
    ) -> list[SiRNACandidate]:
        """Score candidates using composite scoring algorithm.

        Args:
            candidates: Candidates to score in place.
            transcript_sequence: The transcript the candidates were enumerated from, folded once
                for the target_accessibility term. Omitting it leaves the term uncomputable, and
                since no weight is ever redistributed, the candidate then has no design_score at
                all -- it is never scored as if the site were accessible.
        """
        profile = self._build_accessibility_profile(transcript_sequence)

        for candidate in candidates:
            # Calculate component scores
            # Thermodynamic end stabilities and asymmetry
            dg5: float = float("nan")
            dg3: float = float("nan")
            try:
                calc = ThermodynamicCalculator()
                dg5, dg3, asym_score = calc.calculate_asymmetry_score(candidate)
            except Exception:
                # Fallback if ViennaRNA or calculation not available
                asym_score = self._calculate_asymmetry_score(candidate)
                dg5, dg3 = float("nan"), float("nan")
            # Thermodynamic duplex stability (ΔG) and score normalization
            dg_score, duplex_dg = self._calculate_duplex_score(candidate)
            candidate.duplex_stability = duplex_dg
            gc_score = self._calculate_gc_score(candidate.gc_content)
            # Guide self-structure: reported and the EXCESS_PAIRING gate input, not a scoring term.
            self._calculate_guide_structure(candidate)
            access_score = self._store_target_accessibility(candidate, profile)
            ot_score = self._calculate_off_target_score(candidate)
            empirical_score = self._calculate_empirical_score(candidate)
            self._apply_score_filters(candidate, asym_score, empirical_score)

            # Optional melting temperature estimation (rough)
            tm_c = float("nan")
            try:
                if "calc" in locals():
                    tm_c = calc.calculate_melting_temperature(candidate.guide_sequence, candidate.passenger_sequence)
                else:
                    calc2 = ThermodynamicCalculator()
                    tm_c = calc2.calculate_melting_temperature(candidate.guide_sequence, candidate.passenger_sequence)
            except Exception:
                tm_c = float("nan")

            # Store component scores
            # Combine thermodynamic components: favor asymmetry with contribution from duplex stability
            thermo_combo = 0.7 * asym_score + 0.3 * dg_score

            candidate.component_scores = {
                "asymmetry": asym_score,
                # store as float; if missing, use NaN to satisfy type expectations
                "duplex_stability_dg": float(duplex_dg) if duplex_dg is not None else float("nan"),
                "duplex_stability_score": dg_score,
                "thermo_combo": thermo_combo,
                "dg_5p": float(dg5),
                "dg_3p": float(dg3),
                "delta_dg_end": float(dg5 - dg3) if (not math.isnan(dg5) and not math.isnan(dg3)) else float("nan"),
                "melting_temp_c": float(tm_c),
                "gc_content": gc_score,
                "guide_paired_fraction": candidate.paired_fraction,  # Diagnostic; gates EXCESS_PAIRING
                "design_off_target_proxy": ot_score,  # Diagnostic only, not fed into composite
                "empirical": empirical_score,
            }
            # Omitted entirely when inactive: post-screen rescoring reads component_scores back as
            # its design-time feature set, and a stored NaN or 0.0 would be indistinguishable from
            # a site that is genuinely closed.
            if access_score is not None:
                candidate.component_scores["target_accessibility"] = access_score
            # Reported, not scored (issue #97 / D5). Same omit-when-inactive rule as above.
            au_score = au_content_5p_score(candidate.guide_sequence)
            if au_score is not None:
                candidate.component_scores["au_1_5"] = au_score

            self._apply_design_score(candidate, asym_score, gc_score, access_score)
            candidate.asymmetry_score = asym_score

        return candidates

    def _apply_design_score(
        self,
        candidate: SiRNACandidate,
        asym_score: float,
        gc_score: float,
        access_score: float | None,
    ) -> None:
        """Write ``design_score`` on the ``design_v4`` vector, and its per-term contributions.

        ``composite_score`` stays None: it needs ``off_target``, which does not exist until
        screening has run. A NaN sub-score (a ViennaRNA failure) or a None accessibility means "no
        evidence", and since weights are never renormalised there is then no design_score to
        report -- the field stays None rather than being computed over a smaller term set.
        """
        features: dict[str, float] = {}
        if not math.isnan(asym_score):
            features["asymmetry"] = asym_score
        if not math.isnan(gc_score):
            features["gc_content"] = gc_score
        if access_score is not None:
            features["target_accessibility"] = access_score

        candidate.scored_after_screening = False
        vector = self.parameters.scoring.vector_for(post_screen=False)
        try:
            result = compute_composite(features, vector)
        except ScoringError as e:
            logger.warning(f"Design scoring skipped for candidate {candidate.id}: {e}")
            candidate.design_score = None
            candidate.weight_set_version = ""
            candidate.weight_vector = ""
            return

        candidate.design_score = result.score
        candidate.weight_set_version = result.weight_set_version
        candidate.weight_vector = result.vector_name
        candidate.score_asymmetry = result.contributions.get("asymmetry")
        candidate.score_gc_content = result.contributions.get("gc_content")
        candidate.score_target_accessibility = result.contributions.get("target_accessibility")
        # Terms outside design_v4 stay None at design time
        candidate.score_off_target = None
        candidate.score_ago_start = None
        candidate.score_pos1_mismatch = None
        candidate.score_supp_13_16 = None

    def _calculate_duplex_score(self, candidate: SiRNACandidate) -> tuple[float, float | None]:
        """Compute duplex stability ΔG and a normalized score in [0,1].

        ΔG is normalised per nucleotide so 19-23 nt designs stay comparable:
        DUPLEX_DG_PER_NT_STRONG -> 1.0, DUPLEX_DG_PER_NT_WEAK -> 0.0, clamped outside.
        On failure or missing backend, returns (asymmetry_score, None) as a fallback.
        """
        try:
            calc = ThermodynamicCalculator()
            dg = calc.calculate_duplex_stability(candidate.guide_sequence, candidate.passenger_sequence)
            # Normalize per nucleotide: more negative is better
            dg_per_nt = dg / len(candidate.guide_sequence)
            span = DUPLEX_DG_PER_NT_WEAK - DUPLEX_DG_PER_NT_STRONG
            score = (DUPLEX_DG_PER_NT_WEAK - dg_per_nt) / span
            score = max(0.0, min(1.0, score))
            return score, float(dg)
        except Exception:
            # Fallback: use asymmetry as proxy if duplex calc not available
            try:
                asym = self._calculate_asymmetry_score(candidate)
            except Exception:
                asym = 0.5
            return asym, None

    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content percentage."""
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100

    @staticmethod
    def _longest_poly_run(sequence: str) -> int:
        """Length of the longest run of one nucleotide -- the quantity ``max_poly_runs`` compares.

        Returns the length rather than a bool so the gate can record what it observed. Reporting the
        bool left ``max_poly_run_length`` unexported and the gate permanently ``unknown``.
        """
        longest = current = 1
        for previous, base in zip(sequence, sequence[1:], strict=False):
            current = current + 1 if base == previous else 1
            longest = max(longest, current)
        return longest if sequence else 0

    def _has_poly_runs(self, sequence: str, max_runs: int) -> bool:
        """Check for runs of identical nucleotides exceeding threshold."""
        return self._longest_poly_run(sequence) > max_runs

    def _calculate_asymmetry_score(self, candidate: SiRNACandidate) -> float:
        """Calculate thermodynamic asymmetry score via ViennaRNA."""
        calc = ThermodynamicCalculator()
        _, _, asymmetry_score = calc.calculate_asymmetry_score(candidate)
        return asymmetry_score

    def _calculate_gc_score(self, gc_content: float) -> float:
        """Calculate GC content score with Gaussian penalty around 40%."""
        # GC_score = exp(-((GC-40)/10)^2)
        return math.exp(-(((gc_content - 40) / 10) ** 2))

    def _calculate_guide_structure(self, candidate: SiRNACandidate) -> None:
        """Fold the guide against itself and record its own structure.

        This is guide self-structure, not target accessibility: it says nothing about whether the
        mRNA site is open (that is `target_accessibility`, folded from the transcript). It is kept
        because it is reported and because it is the EXCESS_PAIRING gate input.

        `mfe == 0.0` with an all-dots structure is the open chain, the physical floor of the MFE,
        and is a correct answer for a short unstructured guide -- not a failed fold.
        """
        calc = ThermodynamicCalculator()
        # Guides are stored as DNA; ViennaRNA needs RNA, as calculate_melting_temperature already does.
        structure, mfe, paired_fraction = calc.calculate_secondary_structure(_as_rna(candidate.guide_sequence))

        candidate.structure = structure
        candidate.mfe = mfe
        candidate.paired_fraction = paired_fraction
        self._flag_excess_pairing(candidate, paired_fraction)

    def _flag_excess_pairing(self, candidate: SiRNACandidate, paired_fraction: float) -> None:
        """Flag a candidate whose guide is too structured to be accessible."""
        ceiling = self.parameters.filters.max_paired_fraction
        candidate.record_filter_verdict(
            "max_paired_fraction",
            observed=paired_fraction,
            passed=paired_fraction <= ceiling,
            action=self._action_for("max_paired_fraction"),
            status=_ModelCandidate.FilterStatus.EXCESS_PAIRING,
        )

    def _calculate_off_target_score(self, candidate: SiRNACandidate) -> float:
        """Score internal sequence repetitiveness as a design-time off-target proxy.

        This looks only at repeated 7-mers *within* the guide. It is not informed by
        transcriptome or miRNA screening: those run after design and are applied as
        a pass/fail gate (see SiRNAWorkflow._integrate_offtarget_results), never fed
        back into composite_score. `off_target_screened` records whether a candidate
        reached that stage.
        """
        # Comprehensive off-target analysis would require external databases and
        # more complex alignment tools
        guide = candidate.guide_sequence

        # Simple penalty for repetitive sequences
        penalty = 0
        for i in range(len(guide) - 6):
            seed = guide[i : i + 7]
            # Count occurrences of this 7-mer in the sequence
            if guide.count(seed) > 1:
                penalty += 10

        # Transform penalty to score: OT_score = exp(-penalty/50)
        candidate.off_target_penalty = penalty
        return math.exp(-penalty / 50)

    def _calculate_empirical_score(self, candidate: SiRNACandidate) -> float:
        """Calculate empirical score using Reynolds et al. rules (simplified).

        Gate only since issue #96: this score is reported and read by `min_empirical_score`, and is
        not a term in any weight vector. Guides are stored as DNA, so the sequence is read as RNA
        (T is U) before the position-19 test; otherwise a T there never earned the A/U bonus. The
        attainable range is EMPIRICAL_SCORE_MIN..EMPIRICAL_SCORE_MAX (0.4-0.6), not 0..1.
        """
        guide = candidate.guide_sequence.upper().replace("T", "U")
        score = 0.5  # Base score

        # Some simplified Reynolds rules
        # Prefer A/U at position 19 (3' end of guide)
        if len(guide) >= 19 and guide[18] in ("A", "U"):
            score += 0.1

        # Avoid C at position 19
        if len(guide) >= 19 and guide[18] == "C":
            score -= 0.1

        # No rule here judges guide position 1. There used to be a +0.1 for G/C there, which
        # contradicted the biogenesis rule rewarding A/U at the same base (see biogenesis_features):
        # G/C gained +1.6 empirical points and lost 7.9 to the miRNA adjustment, so a 0.15-weight
        # declared term was overridden ~5x by an undeclared one and `empirical` ended up with a
        # NEGATIVE variance share. A/U wins; the clause is gone. Do not reinstate it.
        return max(EMPIRICAL_SCORE_MIN, min(EMPIRICAL_SCORE_MAX, score))

    def _apply_score_filters(self, candidate: SiRNACandidate, asymmetry_score: float, empirical_score: float) -> None:
        """Record the asymmetry and empirical-rule verdicts, rejecting only where the action says to.

        Each threshold gates the quantity it is named after: min_asymmetry_score
        gates the thermodynamic asymmetry score, min_empirical_score gates the
        empirical design-rule score.

        Both verdicts are recorded unconditionally. The old ``passes_filters is not True`` early
        return meant a candidate already rejected by GC or pairing was never *measured* against these
        thresholds, so the two gates' reported counts were a function of gate ordering rather than of
        the candidates -- which is how LOW_ASYMMETRY came to report 6,464 rejections on a run where
        it independently rejected 26,431. ``passes_filters`` still keeps the first label.
        """
        filters = self.parameters.filters
        candidate.record_filter_verdict(
            "min_asymmetry_score",
            observed=asymmetry_score,
            passed=ThermodynamicCalculator.meets_asymmetry_threshold(asymmetry_score, filters.min_asymmetry_score),
            action=self._action_for("min_asymmetry_score"),
            status=_ModelCandidate.FilterStatus.LOW_ASYMMETRY,
        )
        candidate.record_filter_verdict(
            "min_empirical_score",
            observed=empirical_score,
            passed=empirical_score >= filters.min_empirical_score,
            action=self._action_for("min_empirical_score"),
            status=_ModelCandidate.FilterStatus.LOW_EMPIRICAL_SCORE,
        )

    @staticmethod
    def stamp_repeat_verdict(candidate: SiRNACandidate, observations: dict[str, RepeatObservation]) -> None:
        """Stamp repeat metadata and verdict on a single candidate if its guide is flagged.

        The REPEAT_ELEMENT verdict is applied only if the candidate is currently passing
        (passes_filters is True or PASS). A candidate that already failed for another
        reason (GC, asymmetry, etc.) retains its earlier verdict — precedence is:
        existing failure > REPEAT_ELEMENT > PASS.

        Args:
            candidate: Candidate to potentially flag.
            observations: Mapping from normalized guide sequence to RepeatObservation.
        """
        norm_guide = normalize_guide_sequence(candidate.guide_sequence)
        obs = observations.get(norm_guide)
        if obs is None:
            return

        # Write repeat metadata regardless of verdict
        candidate.repeat_flagged = obs.is_repeat
        candidate.repeat_transcript_fraction = obs.transcript_fraction

        # Apply REPEAT_ELEMENT verdict only if currently passing
        if obs.is_repeat and (
            candidate.passes_filters is True or candidate.passes_filters == _ModelCandidate.FilterStatus.PASS
        ):
            candidate.passes_filters = _ModelCandidate.FilterStatus.REPEAT_ELEMENT

    def _get_tool_versions(self) -> dict[str, str]:
        """Get versions of tools used in the analysis."""
        try:
            biopython_version = Bio.__version__
        except AttributeError:
            biopython_version = "unknown"

        python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

        return {
            "python": python_version,
            "biopython": biopython_version,
            "sirnaforge": __version__,
        }


# The three miRNA biogenesis features `biogenesis_features` produces. All three are pure functions
# of the guide and passenger sequences, so they are computable for every candidate -- including
# dirty controls, which never pass through MiRNADesigner. Two of them are scored terms of
# postscreen_mirna_v4; `pos1_mismatch` is computed and reported only, because issue #102 measured it
# exactly constant. Read a vector's own TERM_NAMES for what is scored -- this tuple is what is
# *computed*, and the two are deliberately not the same list.
MIRNA_TERM_NAMES = ("ago_start", "pos1_mismatch", "supp_13_16")

# Guide positions 13-16 (1-based), the 3' supplementary pairing region.
SUPP_REGION_SLICE = slice(12, 16)
SUPP_REGION_MIN_LEN = 16

# Guide positions 1-5 (1-based), the A/U window of issue #97. The window is **pre-declared** at 1-5
# by decision D5 and must not be widened, narrowed or shifted to raise a correlation: single-dataset
# tuning has already misfired twice here (the off-target cap of 3, the 0.65 asymmetry floor).
AU_5P_WINDOW_SLICE = slice(0, 5)
AU_5P_WINDOW_LEN = 5

# Watson-Crick pairs and the G:U wobble, read as RNA.
_PERFECT_PAIRS = frozenset({("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")})
_WOBBLE_PAIRS = frozenset({("G", "U"), ("U", "G")})


def classify_pos1_pairing(guide_base: str, passenger_base: str) -> str:
    """Classify the pairing state at guide position 1: perfect, wobble or mismatch.

    Bases are stored as DNA, so they are read as RNA before lookup: otherwise A:T is not found in
    the Watson-Crick set and every A:U pair is called a mismatch.
    """
    pair = (_as_rna(guide_base), _as_rna(passenger_base))
    if pair in _PERFECT_PAIRS:
        return "perfect"
    if pair in _WOBBLE_PAIRS:
        return "wobble"
    return "mismatch"


def supplementary_score(guide: str) -> float:
    """3' supplementary pairing sub-score from guide positions 13-16, in [0, 1].

    High A/U content there means low pairing stability, which is the desirable direction (less 3'
    supplementary pairing, better specificity), so the sub-score rises with A/U as every other
    term rises with the thing it wants.
    """
    if len(guide) < SUPP_REGION_MIN_LEN:
        return 0.5  # Default for short sequences

    supp_region = _as_rna(guide[SUPP_REGION_SLICE])
    au_count = supp_region.count("A") + supp_region.count("U")
    return au_count / len(supp_region) if supp_region else 0.5


def au_content_5p_score(guide: str) -> float | None:
    """A/U content over guide positions 1-5, as a fraction in [0, 1]. Issue #97, window per D5.

    Reported on every candidate and scored by no default vector in 0.7.1: promoting it needs the
    post-screen feature assembly to forward it, which is not this module's file. Returns None for a
    guide shorter than the window rather than a substituted midpoint -- a missing input must not
    score as a good one.

    Pure sequence, so it is computable wherever a guide exists, including on paths that have no
    transcript context and no ViennaRNA.
    """
    if len(guide) < AU_5P_WINDOW_LEN:
        return None
    window = _as_rna(guide[AU_5P_WINDOW_SLICE])
    return (window.count("A") + window.count("U")) / AU_5P_WINDOW_LEN


def biogenesis_features(guide: str, passenger: str) -> dict[str, float]:
    """The three miRNA biogenesis sub-scores, from sequence alone.

    Computed here rather than read back off the candidate so post-screen scoring cannot be handed a
    partially-populated row: every term in postscreen_mirna_v4 must be present, and re-deriving them
    from the sequences means they always are. Each is in [0, 1] like every other term -- they are
    features now, not bonuses added to a finished score.
    """
    pos1_state = classify_pos1_pairing(guide[0] if guide else "", passenger[-1] if passenger else "")
    return {
        # Argonaute loading prefers A/U at guide position 1.
        "ago_start": 1.0 if _as_rna(guide[:1]) in ("A", "U") else 0.0,
        # A G:U wobble or mismatch at position 1 is preferred over a perfect pair.
        "pos1_mismatch": 1.0 if pos1_state in ("wobble", "mismatch") else 0.0,
        "supp_13_16": supplementary_score(guide),
    }


class MiRNADesigner(SiRNADesigner):
    """miRNA-biogenesis-aware siRNA designer with specialized scoring.

    Extends SiRNADesigner with scoring rules optimized for miRNA-like processing:
    - Argonaute selection preferences (pos1 A/U scored; pos1 pairing state reported only)
    - 3' supplementary pairing analysis (positions 13-16)
    - Conservative thermodynamic thresholds
    - Seed region quality assessment
    """

    def __init__(
        self,
        parameters: DesignParameters,
        *,
        filter_actions: Mapping[str, FilterAction] | None = None,
    ) -> None:
        """Initialize miRNA designer with miRNA-specific config validation."""
        super().__init__(parameters, filter_actions=filter_actions)
        # Optionally apply miRNA-specific filter adjustments based on MiRNADesignConfig
        # For now, we rely on the caller to set appropriate filters for miRNA mode

    def _score_candidates(
        self, candidates: list[SiRNACandidate], transcript_sequence: str | None = None
    ) -> list[SiRNACandidate]:
        """Record the miRNA biogenesis evidence, then score exactly as siRNA mode does.

        The design stage uses one vector (``design_v4``) in both modes, so this method no longer
        touches any score: two of the three biogenesis quantities -- ``ago_start`` and
        ``supp_13_16`` -- are terms of ``postscreen_mirna_v4`` and only enter once ``off_target``
        exists, and ``pos1_mismatch`` is scored by no vector at all since #102. What this method does
        do is record all three -- as reported fields and as ``component_scores`` entries -- so the CSV
        shows why a miRNA run ranks as it does.

        Before issue #96 this method folded the bonuses into ``composite_score`` and divided the
        result by ``1 + max_bonus``, which scaled every declared weight by 0.80 in miRNA mode.

        Args:
            candidates: Candidates to score in place.
            transcript_sequence: Forwarded to the base scorer for target_accessibility.
        """
        candidates = super()._score_candidates(candidates, transcript_sequence)

        for candidate in candidates:
            guide = candidate.guide_sequence
            passenger = candidate.passenger_sequence

            # Reported fields keep the stored (DNA) spelling of the base; the tests do not.
            candidate.guide_pos1_base = guide[0] if guide else ""
            candidate.pos1_pairing_state = classify_pos1_pairing(
                candidate.guide_pos1_base, passenger[-1] if passenger else ""
            )
            candidate.supp_13_16_score = supplementary_score(guide)
            candidate.seed_class = self._classify_seed_region(guide)

            # Diagnostics on the row; post-screen scoring re-derives them from the sequences so a
            # candidate that never reached this method is still fully scorable.
            candidate.component_scores.update(biogenesis_features(guide, passenger))

        return candidates

    def _classify_pos1_pairing(self, guide_base: str, passenger_base: str) -> str:
        """Classify pairing state at guide position 1 (delegates to `classify_pos1_pairing`)."""
        return classify_pos1_pairing(guide_base, passenger_base)

    def _calculate_supplementary_score(self, guide: str) -> float:
        """3' supplementary pairing sub-score (delegates to `supplementary_score`)."""
        return supplementary_score(guide)

    def _classify_seed_region(self, guide: str) -> str:
        """Classify seed match class based on guide positions 2-8.

        Args:
            guide: Guide strand sequence

        Returns:
            Seed class: "6mer", "7mer-m8", "7mer-a1", or "8mer"
        """
        # This is a simplified classification based on seed length
        # In practice, seed class depends on target matching, which happens during off-target analysis
        # For now, we just categorize based on sequence properties
        if len(guide) < 8:
            return "6mer"
        # This is a placeholder - actual seed class determined during off-target matching
        return "8mer"
