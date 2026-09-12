"""Build the report payload from a finished run directory.

Issue #103. This module reads ``candidates_all.csv``, the aggregated hit tables and ``manifest.json``
and returns a typed payload. It performs **no** alignment, folding, network access or classification:
``hit_class`` is read from the table the run published, so the report can never disagree with the run.

The report unit is the **guide sequence**, not the candidate id. A guide is enumerated once per
transcript it was found on -- median 6 and up to 34 rows on the public 0.7.1 baseline -- while
off-target screening runs once per distinct guide. Joining evidence on ``id`` therefore leaves most
candidates falsely reading as having zero off-targets, so the guide sequence is the join key and the
per-transcript rows collapse into an isoform sub-table.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from collections.abc import Container, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from sirnaforge.config.run_policy import (
    EntryPoint,
    ResolvedRunPolicy,
    RunPolicyError,
    filters_from_manifest,
    resolve_run_policy,
)
from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS, hit_class_of, is_annotated
from sirnaforge.core.hit_classification import HitClass
from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation
from sirnaforge.reporting.structure import layouts_for
from sirnaforge.reporting.tracks import transcript_regions, uncovered_stretches

#: Bump when the payload's shape changes, so a report and the run it describes can never be
#: silently mismatched.
PAYLOAD_SCHEMA_VERSION = "1.0.0"

#: Hit rows embedded per guide. Everything outside this scope is carried as counts only, which is a
#: deliberate scope decision (#103), not a limitation -- and a guide with hits only outside it must
#: render its counts rather than an empty panel.
EMBED_SPECIES = "human"
EMBED_MAX_NM = 2

#: Classes that are a liability. On-target isoforms, orthologues and repeat-mediated hits are
#: displayed, but never counted as off-targets.
LIABILITY = frozenset({HitClass.OFF_TARGET, HitClass.UNDETERMINED})

_UNKNOWN_SYMBOL = "unknown"


class ReportInputError(ValueError):
    """A run directory that cannot produce a truthful report."""


@dataclass
class GuideEntry:
    """Everything the report shows for one guide sequence."""

    guide: str
    passenger: str | None
    overhang: str | None
    modifications: str | None
    n_rows: int
    design_score: float | None
    composite_score: float | None
    weight_vector: str | None
    metrics: dict[str, Any]
    gates: list[list[Any]]
    isoforms: list[dict[str, Any]]
    offtarget_by_symbol: list[dict[str, Any]]
    offtarget_matrix: list[dict[str, Any]]
    offtarget_rows: list[dict[str, Any]]
    offtarget_embedded_scope_empty_but_counts_exist: bool
    liability_count: int
    mirna: list[dict[str, Any]]
    #: species -> {nm, seed_mismatches} for the best ortholog alignment. Empty when none was found.
    ortholog: dict[str, dict[str, Any]] = field(default_factory=dict)
    run_verdict: str | None = None
    structure: str | None = None
    transcript_hits: int | None = None

    @property
    def undeclared_run_rejection(self) -> bool:
        """The run rejected this guide everywhere, and no declared gate accounts for it.

        ``REPEAT_ELEMENT`` is the live case: the pipeline stamps it, and the 16-filter registry
        declares no repeat gate, so the report has no descriptor that can re-derive the rejection. It
        must not therefore call the guide clean -- on one MSH3 run that would have published 185
        guides as passing that the run threw out.
        """
        return self.run_verdict not in (None, "PASS") and not (self.n_gates_failed or self.n_gates_unknown)

    @property
    def n_gates_failed(self) -> int:
        """Gates this guide independently fails -- all of them, not only the first to fire."""
        return sum(1 for g in self.gates if g[1] == _VERDICT_CODE[FilterEvaluation.FAIL.value])

    @property
    def n_gates_unknown(self) -> int:
        """Gates that could not be evaluated. An unknown is never folded into a pass."""
        return sum(1 for g in self.gates if g[1] == _VERDICT_CODE[FilterEvaluation.UNKNOWN.value])

    @property
    def n_gates_warned(self) -> int:
        """Gates this guide exceeds whose action is ``warn``: reported, and not a rejection."""
        return sum(1 for g in self.gates if g[1] == VERDICT_WARN)

    @property
    def status(self) -> str:
        """``fail`` / ``unknown`` / ``warn`` / ``pass`` -- four states, because fewer would lie.

        A guide failing no gate is only clean when every gate could actually be evaluated. On the
        public 0.7.1 baseline 450 guides fail nothing here while the run failed them, because the
        gates that failed them read counters the run does not export. Collapsing that into "pass"
        is the fabricated-evidence mistake this whole report exists to make visible.

        ``warn`` is separate from ``fail`` for the mirror-image reason. A warn-action gate the guide
        exceeds is a real finding, but the run did not reject the guide for it, so calling it ``fail``
        would make the report contradict a run PASS it actually agrees with.

        A rejection no declared gate can express is also ``unknown``: see
        :attr:`undeclared_run_rejection`. The report may report less than the run. It may not report
        more.
        """
        if self.n_gates_failed:
            return "fail"
        if self.n_gates_unknown or self.undeclared_run_rejection:
            return "unknown"
        return "warn" if self.n_gates_warned else "pass"


@dataclass
class ReportPayload:
    """The whole report, ready to render."""

    schema_version: str
    run: dict[str, Any]
    filters: list[dict[str, Any]]
    guides: list[GuideEntry]
    provenance: dict[str, Any]
    caveats: list[str] = field(default_factory=list)


def _normalise_guide(seq: object) -> str:
    """Uppercase RNA spelling, so a U/T difference cannot split one guide into two."""
    return str(seq or "").strip().upper().replace("T", "U")


def _num(value: object) -> float | int | None:
    """A number, or None for NaN/blank -- never a fabricated zero.

    Coerces through ``float`` rather than testing ``isinstance(value, (int, float))``, because
    ``numpy.int64`` is **not** a Python ``int`` while ``numpy.float64`` *is* a Python ``float``, and a
    row taken with ``DataFrame.iloc`` hands back numpy scalars. The isinstance test therefore nulled
    every integer-dtype column and let float columns through, which made ``max_off_target_count``
    read ``unknown`` for every guide on a real run while its threshold was rejecting 65% of them.
    """
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, str) and not value.strip():
        return None
    try:
        number = float(value)  # type: ignore[arg-type]  # non-numeric raises, which is the answer
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):  # NaN or infinite is not a value
        return None
    return int(number) if number.is_integer() else round(number, 6)


_COMPARE = {
    FilterComparator.LE: lambda v, t: v <= t,
    FilterComparator.GE: lambda v, t: v >= t,
    FilterComparator.LT: lambda v, t: v < t,
    FilterComparator.GT: lambda v, t: v > t,
}


#: Why a gate reached its verdict, as a code. The descriptor is emitted once for the whole report, so
#: a per-guide outcome carries only what varies -- value, verdict and this. Repeating filter_id,
#: column, comparator, threshold and scope on every guide cost 9 MB of the first render's 38 MB.
REASON_OK = 0
REASON_MISSING_COLUMN = 1
REASON_EMPTY_VALUE = 2
REASON_FILTER_OFF = 3
REASON_NO_THRESHOLD = 4
#: The run itself recorded UNKNOWN for this gate. A distinct reason from an empty or absent column,
#: because the run measured nothing on purpose and said so, rather than the report failing to find it.
REASON_RUN_UNKNOWN = 5

_VERDICT_CODE = {
    FilterEvaluation.PASS.value: 0,
    FilterEvaluation.FAIL.value: 1,
    FilterEvaluation.UNKNOWN.value: 2,
    FilterEvaluation.NOT_EVALUATED.value: 3,
}

#: A gate whose threshold was exceeded but whose action is ``warn``: a real finding the reader should
#: see, and not a rejection. It needs its own code because folding it into ``fail`` would make the
#: report contradict the run on every warn-flagged guide -- inverting the one metric
#: (``contradicted_run_pass``) that exists to prove the report and the pipeline agree.
VERDICT_WARN = 4

#: Guide statuses, worst first. ``warn`` belongs in the header tally like the rest: counting only
#: pass/unknown/fail dropped every warn-status guide out of a total that claims to be all of them.
STATUSES = ("fail", "unknown", "warn", "pass")


def observed_column(descriptor: Any, populated: Container[str]) -> str | None:
    """The column carrying the value this gate compared, or None if the run exports neither.

    ``<filter_id>_observed`` is preferred over the descriptor's own ``column``, because it is the
    number the pipeline itself compared. It is the answer for the gates that read human-stratified
    counters: 0.7.1 does not export ``transcriptome_hits_1mm_human`` under that name (#101), but it
    does export ``max_transcriptome_hits_1mm_observed``, and on a four-species MSH3 run the observed
    column reproduces each gate's own verdict on 100% of 40,081 rows while the same-named all-species
    column disagrees -- 17,600 hits against 63,801. Reading the descriptor's column instead would let
    the report contradict the run using a counter with a wider scope than the gate's.

    ``populated`` must hold only columns that carry at least one value. A column present but empty
    for every row is the shape a gate takes when the run never recorded its verdict, and preferring
    it would turn a gate the descriptor column *can* answer into an unknown.
    """
    for candidate in (f"{descriptor.filter_id}_observed", descriptor.column):
        if candidate in populated:
            return candidate
    return None


def _recorded_verdict(row: pd.Series, filter_id: str) -> str | None:
    """The verdict the RUN recorded for one gate on one row, from its ``<filter_id>_verdict`` column.

    Distinct from :func:`_run_verdict`, which reads the single ``passes_filters`` label over a guide's
    rows. Only this per-gate record can say "in force, evidence unavailable" -- a state no threshold
    comparison reconstructs.
    """
    value = row.get(f"{filter_id}_verdict")
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    return str(value)


def _evaluate(descriptor: Any, row: pd.Series, column: str | None) -> tuple[float | int | None, int, int]:
    """Evaluate one descriptor against one candidate row, independently of every other gate.

    Returns ``unknown`` rather than a pass when the run exports no column the threshold can read.
    Reporting those as passes would be the fabricated-evidence failure this report exists to make
    visible.

    Args:
        descriptor: The gate as configured for this run.
        row: One candidate row.
        column: Column to read, from :func:`observed_column`; None when the run exports neither.

    Returns:
        ``(value, verdict_code, reason_code)``. The descriptor itself is emitted once per report.
    """
    # Read the value first, and report it whatever the verdict turns out to be. A gate that is off or
    # has no threshold still *measured* something the reader wants: three of the four off gates on the
    # reference run carry a real per-guide number, and `max_mirna_1mm_seed`'s own policy definition says
    # "the number is reported and nothing acts on it" -- which returning None made untrue.
    value = _num(row.get(column)) if column is not None else None
    not_evaluated = _VERDICT_CODE[FilterEvaluation.NOT_EVALUATED.value]
    if descriptor.action.value == FilterAction.OFF.value:
        return value, not_evaluated, REASON_FILTER_OFF
    if descriptor.threshold is None:
        return value, not_evaluated, REASON_NO_THRESHOLD
    # The run's own UNKNOWN wins, before any comparison: re-thresholding moves a ceiling, it cannot
    # conjure the measurement. Load-bearing, not defensive -- the gate writes UNKNOWN with an empty
    # observed value so nothing re-derives a pass, but `observed_column` falls back to the descriptor's
    # own column when the observed one is empty for every row, and three of those fallbacks are exported
    # and default to 0. Without this an unscreened guide's `max_off_target_count` read as a pass at 0.
    if _recorded_verdict(row, descriptor.filter_id) == FilterEvaluation.UNKNOWN.value:
        return None, _VERDICT_CODE[FilterEvaluation.UNKNOWN.value], REASON_RUN_UNKNOWN
    if column is None:
        return None, _VERDICT_CODE[FilterEvaluation.UNKNOWN.value], REASON_MISSING_COLUMN
    if value is None:
        return None, _VERDICT_CODE[FilterEvaluation.UNKNOWN.value], REASON_EMPTY_VALUE

    passed = _COMPARE[descriptor.comparator](value, descriptor.threshold)
    if passed:
        code = _VERDICT_CODE[FilterEvaluation.PASS.value]
    elif descriptor.action.value == FilterAction.WARN.value:
        code = VERDICT_WARN
    else:
        code = _VERDICT_CODE[FilterEvaluation.FAIL.value]
    return value, code, REASON_OK


def _transcript_hits(rows: pd.DataFrame) -> int | None:
    """How many distinct transcripts carry this guide -- the numerator of isoform coverage.

    Not the row count. A guide's site can occur twice in one transcript, so enumerations exceed
    isoforms: on one MSH3 run 14 guides have more rows than transcripts and one has 19 rows over 10.
    Prefers the run's own ``transcript_hit_count`` (which equals the distinct count on every row of
    that run) and falls back to counting, so the report agrees with the column when it exists.
    """
    stated = _num(rows["transcript_hit_count"].iloc[0]) if "transcript_hit_count" in rows.columns else None
    if stated is not None:
        return int(stated)
    return int(rows["transcript_id"].nunique()) if "transcript_id" in rows.columns else None


def _run_verdict(rows: pd.DataFrame) -> str | None:
    """The run's own verdict for a guide: PASS if any of its rows passed, else its first label.

    Guide-level because screening is: one passing enumeration is enough for the run to keep the
    guide, so the report must read the group rather than the best row.
    """
    if "passes_filters" not in rows.columns:
        return None
    labels = [str(v).split(" (")[0] for v in rows["passes_filters"].dropna()]
    if not labels:
        return None
    return "PASS" if "PASS" in labels else labels[0]


def _scope_label(descriptor: Any) -> str:
    """Human-readable scope, emitted once per filter rather than once per guide."""
    label = ",".join(sorted(descriptor.scope.species)) or "all species"
    if descriptor.scope.max_mismatches is not None:
        label += f", nm<={descriptor.scope.max_mismatches}"
    return label


def _gene_query(manifest: Mapping[str, Any], candidates: pd.DataFrame) -> str:
    """The run's target: the manifest's own ``gene_query``, else a gene column, else unknown.

    The manifest is the authority because ``candidates_all.csv`` carries no gene column at all -- an
    earlier version read this from a loop variable that had leaked past its loop and reported a
    transcript accession as the gene.
    """
    stated = str(manifest.get("gene_query") or "").strip()
    if stated:
        return stated
    for column in ("gene_name", "gene_symbol", "gene_id"):
        if column in candidates.columns:
            values = candidates[column].dropna().astype(str)
            values = values[values.str.strip().ne("")]
            if not values.empty:
                names = sorted(values.unique())
                return names[0] if len(names) == 1 else f"{names[0]} +{len(names) - 1} more"
    return "unknown"


def _verdict_agreement(candidates: pd.DataFrame, guides: list[GuideEntry]) -> dict[str, Any]:
    """Compare the report's independent verdicts against the run's own ``passes_filters`` column.

    Published in the report header, because a reader has no other way to tell whether the gate panel
    agrees with the pipeline. The asymmetry is the point: the report must never contradict a run PASS,
    while a run failure it cannot re-derive is reported ``unknown`` rather than quietly flipped.
    """
    if "passes_filters" not in candidates.columns:
        return {"comparable": False, "reason": "the run exports no passes_filters column"}
    run_pass = set(candidates.loc[candidates["passes_filters"] == "PASS", "_guide"])
    by_status = {g.guide: g.status for g in guides}
    contradicted = sorted(g for g in run_pass if by_status.get(g) == "fail")
    not_rederivable = sorted(g for g, s in by_status.items() if s == "unknown" and g not in run_pass)
    overruled = sorted(g.guide for g in guides if g.status in ("pass", "warn") and g.guide not in run_pass)
    return {
        "comparable": True,
        "run_pass_guides": len(run_pass),
        "contradicted_run_pass": len(contradicted),
        "run_failed_not_rederivable": len(not_rederivable),
        # Must stay 0. It counts guides the report calls clean while the run rejected them -- the
        # fabricated-evidence direction, which no earlier metric could see.
        "overruled_run_fail": len(overruled),
        "undeclared_run_rejections": sum(1 for g in guides if g.undeclared_run_rejection),
    }


def _read_hits(path: Path, caveats: list[str]) -> pd.DataFrame:
    """Read a published hit table, requiring the classification columns the report displays."""
    empty: pd.DataFrame = pd.DataFrame()
    if not path.exists():
        caveats.append(f"no hit table at {path.name}; every off-target panel is empty for want of input")
        return empty
    hits: pd.DataFrame = pd.read_csv(path, sep="\t", low_memory=False)
    missing = [c for c in CLASSIFICATION_COLUMNS if c not in hits.columns]
    if missing:
        raise ReportInputError(
            f"{path.name} is missing {missing}. The report will not classify hits itself -- it would then "
            "be able to disagree with the run. Re-publish the table through the workflow, which writes "
            "these columns, and report again."
        )
        # Note for reviewers: this is deliberately fatal. A table without hit_class would render
        # on-target isoform alignments as off-targets, which on the public baseline is 50.1% of rows.
    unannotated = int((~hits.apply(is_annotated, axis=1)).sum()) if len(hits) else 0
    if unannotated:
        caveats.append(f"{unannotated} hit rows carry no classification and are shown as undetermined")
    return hits


@dataclass(frozen=True)
class _GatePanel:
    """The gates the report will apply, and where they came from."""

    filters: tuple[Any, ...]
    profile_name: str
    run_mode: str
    source: str


def _gate_panel(policy: ResolvedRunPolicy | None, manifest: Mapping[str, Any], caveats: list[str]) -> _GatePanel:
    """The run's own gates: the live policy, else the manifest, else library defaults with a caveat.

    Library defaults are the last resort and must be declared, because they are not this run's gates.
    Defaulting silently published a GC ceiling of 60 against a run that set 65, which failed 227
    guides on a threshold the run never applied and made the report contradict 41 of its PASSes.
    """
    if policy is not None:
        return _GatePanel(tuple(policy.filters), policy.profile.name, policy.run_mode.value, "the live run policy")

    block = manifest.get("run_policy")
    if isinstance(block, Mapping):
        try:
            filters = filters_from_manifest(block)
        except RunPolicyError as exc:
            caveats.append(f"the manifest's gate registry is unusable ({exc}); gates below are library defaults")
        else:
            if filters:
                profile = block.get("profile") or {}
                return _GatePanel(
                    filters,
                    str(profile.get("name") or "unknown") if isinstance(profile, Mapping) else "unknown",
                    str(block.get("run_mode") or "unknown"),
                    "the run's own manifest",
                )
            caveats.append("the manifest declares no gates; gates below are library defaults")
    else:
        caveats.append("this run published no policy in its manifest; gates below are library defaults")

    fallback = resolve_run_policy(
        entry_point=EntryPoint.SCREENING_WORKFLOW, query_species="human", screen_species=["human"]
    )
    return _GatePanel(tuple(fallback.filters), fallback.profile.name, fallback.run_mode.value, "library defaults")


def build_payload(run_dir: Path | str, *, policy: ResolvedRunPolicy | None = None) -> ReportPayload:
    """Build the payload for a finished run directory.

    Args:
        run_dir: A completed output directory. No live pipeline state is required.
        policy: Resolved policy supplying the gate descriptors. Read from the run's own
            ``manifest.json`` when omitted.

    Returns:
        A :class:`ReportPayload`.

    Raises:
        ReportInputError: The directory has no candidate table, or a hit table with no
            classification columns.
    """
    run_dir = Path(run_dir)
    caveats: list[str] = []

    candidates_csv = next(iter(sorted(run_dir.glob("**/candidates_all.csv"))), None)
    if candidates_csv is None:
        raise ReportInputError(f"no candidates_all.csv under {run_dir}")
    candidates = pd.read_csv(candidates_csv, low_memory=False)

    agg = run_dir / "off_target" / "results" / "aggregated"
    hits = _read_hits(agg / "combined_offtargets.tsv", caveats)
    mirna_path = agg / "combined_mirna_hits.tsv"
    mirna = pd.read_csv(mirna_path, sep="\t", low_memory=False) if mirna_path.exists() else pd.DataFrame()

    manifest_path = next(iter(sorted(run_dir.glob("**/manifest.json"))), None)
    manifest: dict[str, Any] = {}
    if manifest_path is not None:
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:  # a bad manifest costs provenance, not the report
            caveats.append(f"manifest.json could not be read ({exc}); the header shows less provenance")

    panel = _gate_panel(policy, manifest, caveats)

    candidates["_guide"] = candidates["guide_sequence"].map(_normalise_guide)
    if not hits.empty:
        hits["_guide"] = hits["qseq"].map(_normalise_guide)
        hits["_nm"] = pd.to_numeric(hits["nm"], errors="coerce")
        hits["_class"] = hits.apply(lambda r: hit_class_of(r).value, axis=1)
    if not mirna.empty and "qseq" in mirna.columns:
        mirna["_guide"] = mirna["qseq"].map(_normalise_guide)

    hits_by_guide = dict(tuple(hits.groupby("_guide"))) if not hits.empty else {}
    mirna_by_guide = dict(tuple(mirna.groupby("_guide"))) if not mirna.empty and "_guide" in mirna else {}

    descriptors = [f.descriptor for f in panel.filters]
    populated = {c for c in candidates.columns if candidates[c].notna().any()}
    gate_columns = [observed_column(d, populated) for d in descriptors]
    register = _register_index(candidates)
    guides: list[GuideEntry] = []

    # Sorted so `best` is the guide's best-scoring enumeration rather than whichever row the CSV
    # happens to list first; the gates read constant-per-guide columns either way.
    by = [c for c in ("composite_score", "design_score") if c in candidates.columns]
    ranked = candidates.sort_values(by, ascending=False, na_position="last") if by else candidates
    for guide, rows in ranked.groupby("_guide", sort=False):
        guides.append(
            _build_guide(
                guide=str(guide),
                rows=rows,
                best=rows.iloc[0],
                descriptors=descriptors,
                gate_columns=gate_columns,
                hits=hits_by_guide.get(guide, pd.DataFrame()),
                mirna=mirna_by_guide.get(guide, pd.DataFrame()),
                register=register,
            )
        )

    guides.sort(key=lambda g: (-(g.composite_score or g.design_score or -1), g.guide))

    undeclared = sorted({g.run_verdict for g in guides if g.undeclared_run_rejection if g.run_verdict})
    if undeclared:
        caveats.append(
            f"{sum(1 for g in guides if g.undeclared_run_rejection)} guides were rejected by the run as "
            f"{', '.join(undeclared)}, which no declared gate expresses; they are reported not "
            "established rather than clean"
        )

    n_liab = sum(g.liability_count for g in guides)
    agreement = _verdict_agreement(candidates, guides)
    run = {
        "candidate_rows": int(len(candidates)),
        "guides": len(guides),
        "hit_rows": int(len(hits)),
        "liability_rows": n_liab,
        "non_liability_rows": int(len(hits)) - n_liab,
        "mirna_rows": int(len(mirna)),
        "embed_scope": f"{EMBED_SPECIES}, nm<={EMBED_MAX_NM}",
        "gene_query": _gene_query(manifest, candidates),
        "status_counts": {s: sum(1 for g in guides if g.status == s) for s in STATUSES},
        "agreement": agreement,
        # The isoform-coverage denominator. Every transcript the design step enumerated on, which is
        # not the same as the transcripts the map can plot -- a transcript can be present with no
        # scoreable row and still be an isoform this guide either does or does not reach.
        "transcript_ids": sorted(candidates["transcript_id"].dropna().astype(str).unique())
        if "transcript_id" in candidates.columns
        else [],
        "transcripts": _transcript_maps(candidates, guides, run_dir, caveats),
        # Keyed by dot-bracket and computed once per distinct structure, which is what makes the
        # layouts small enough to embed: 40,079 candidates carry 1,333 distinct structures.
        "structure_layouts": layouts_for(g.structure for g in guides),
    }
    return ReportPayload(
        schema_version=PAYLOAD_SCHEMA_VERSION,
        run=run,
        filters=[
            {
                **f.descriptor.model_dump(mode="json"),
                "setting_key": f.setting_key,
                "definition": f.definition,
                "scope_label": _scope_label(f.descriptor),
                "read_column": column,
            }
            for f, column in zip(panel.filters, gate_columns, strict=True)
        ],
        guides=guides,
        provenance={
            "run_dir": str(run_dir),
            "candidates_csv": str(candidates_csv.relative_to(run_dir)),
            "manifest_present": manifest_path is not None,
            "tool_version": str(manifest.get("tool_version") or "unknown"),
            "run_timestamp": str(manifest.get("run_timestamp") or "unknown"),
            "policy_profile": panel.profile_name,
            "run_mode": panel.run_mode,
            "policy_source": panel.source,
            "payload_schema_version": PAYLOAD_SCHEMA_VERSION,
        },
        caveats=caveats,
    )


#: Transcript stretches shorter than this are not reported as carrying no candidate: a 23-mer window
#: cannot start in the last 22 nt of a transcript, so short tails are arithmetic, not a design choice.
MIN_UNCOVERED_NT = 40

#: Value plotted on the design map, in order of preference. ``composite_score`` does not exist until
#: post-screen scoring, so a design-only run or one whose scoring failed has only ``design_score``.
#: Without the fallback those runs lose the whole map to an all-NaN column.
_MAP_VALUE_COLUMNS = ("composite_score", "design_score")


def _transcript_maps(
    candidates: pd.DataFrame, guides: list[GuideEntry], run_dir: Path, caveats: list[str]
) -> list[dict[str, Any]]:
    """Per-transcript position series for the design map, most-enumerated transcript first.

    A row the run rejected is plotted as a rejection; every other row carries its guide's report
    status, so a window the report cannot establish is not drawn as passing. Regions come from the
    run's ORF report; a transcript missing from it still gets a map, without a region bar.
    """
    required = {"position", "transcript_id", "passes_filters"}
    missing = sorted(required - set(candidates.columns))
    if missing:
        caveats.append(f"no design map: this run exports no {', '.join(missing)}")
        return []
    value_column = next((c for c in _MAP_VALUE_COLUMNS if c in candidates.columns), None)
    if value_column is None:
        caveats.append(f"no design map: this run exports none of {', '.join(_MAP_VALUE_COLUMNS)}")
        return []

    regions = transcript_regions(run_dir)
    if not regions:
        caveats.append("no ORF report in this run, so the design map cannot label CDS and UTR")

    status_by_guide = {g.guide: g.status for g in guides}
    ids = candidates["id"].astype(str) if "id" in candidates.columns else pd.Series("", index=candidates.index)
    rows = candidates[~ids.str.contains("DIRTY", na=False)].copy()
    rows["_pos"] = pd.to_numeric(rows["position"], errors="coerce")
    rows["_val"] = pd.to_numeric(rows[value_column], errors="coerce")
    rows = rows.dropna(subset=["_pos", "_val"])
    run_pass = rows["passes_filters"].astype(str).str.split(" (", regex=False).str[0].eq("PASS")
    rows["_class"] = [
        "fail" if not passed else status_by_guide.get(guide, "unknown")
        for guide, passed in zip(rows["_guide"], run_pass, strict=True)
    ]

    # The run exports no window-length column, so take it from the guides themselves.
    window = int(rows["guide_sequence"].astype(str).str.len().max() or 0) if len(rows) else 0
    out: list[dict[str, Any]] = []
    for tid, group in rows.groupby("transcript_id", sort=False):
        region = regions.get(str(tid))
        length = region.length if region else int(group["_pos"].max()) + max(window - 1, 0)
        series = {
            name: [[int(p), round(float(v), 3)] for p, v in zip(part["_pos"], part["_val"], strict=True)]
            for name, part in group.groupby("_class", sort=False)
        }
        out.append(
            {
                "transcript_id": str(tid),
                "length": int(length),
                "cds_start": region.cds_start if region else None,
                "cds_end": region.cds_end if region else None,
                "windows": int(len(group)),
                "value_column": value_column,
                "series": series,
                "gaps": [
                    [start, end]
                    for start, end in uncovered_stretches(
                        (int(p) for p in group["_pos"]), int(length), window=max(window, 1), min_nt=MIN_UNCOVERED_NT
                    )
                ]
                if window
                else [],
            }
        )
    out.sort(key=lambda d: (-d["windows"], d["transcript_id"]))
    return out


def _build_guide(
    *,
    guide: str,
    rows: pd.DataFrame,
    best: pd.Series,
    descriptors: list[Any],
    gate_columns: list[str | None],
    hits: pd.DataFrame,
    mirna: pd.DataFrame,
    register: dict[str, list[int]],
) -> GuideEntry:
    gates = [list(_evaluate(d, best, c)) for d, c in zip(descriptors, gate_columns, strict=True)]

    isoforms = _isoform_table(rows, register)
    by_symbol, matrix, embedded, liability = _offtarget_views(hits)

    counts_exist = bool(len(hits)) and not embedded
    return GuideEntry(
        run_verdict=_run_verdict(rows),
        transcript_hits=_transcript_hits(rows),
        structure=(str(best.get("structure")) if pd.notna(best.get("structure")) else None),
        guide=guide,
        passenger=(str(best.get("passenger_sequence")) if pd.notna(best.get("passenger_sequence")) else None),
        overhang=(str(best.get("passenger_overhang")) if pd.notna(best.get("passenger_overhang")) else None),
        modifications=(
            str(best.get("passenger_modifications")) if pd.notna(best.get("passenger_modifications")) else None
        ),
        n_rows=int(len(rows)),
        design_score=_num(best.get("design_score")),
        composite_score=_num(best.get("composite_score")),
        weight_vector=(str(best.get("weight_vector")) if pd.notna(best.get("weight_vector")) else None),
        metrics={
            k: _num(best.get(k))
            for k in (
                "gc_content",
                "asymmetry_score",
                "paired_fraction",
                "mfe",
                "empirical_score",
                "target_accessibility_p",
                "melting_temperature",
                "off_target_count",
                "undetermined_hits",
                "transcriptome_hits_total",
                "mirna_hits_0mm_seed",
            )
            if k in best.index
        },
        gates=gates,
        isoforms=isoforms,
        offtarget_by_symbol=by_symbol,
        offtarget_matrix=matrix,
        offtarget_rows=embedded,
        offtarget_embedded_scope_empty_but_counts_exist=counts_exist,
        liability_count=liability,
        mirna=_mirna_table(mirna),
        ortholog=_ortholog_conservation(hits),
    )


#: Positions this close on one transcript share a window, so their designs can score identically.
REGISTER_NEIGHBOUR_NT = 2


def _register_index(candidates: pd.DataFrame) -> dict[str, list[int]]:
    """Transcript -> every enumerated position, across **all** guides.

    Register neighbours are cross-guide by nature: the point is that two *different* designs one
    nucleotide apart share a window. Built once from the whole candidate table, because a per-guide
    map can only ever see the guide's own rows and would report no neighbours at all.
    """
    index: dict[str, list[int]] = defaultdict(list)
    for tx, pos in zip(candidates.get("transcript_id", []), candidates.get("position", []), strict=False):
        p = _num(pos)
        if p is not None:
            index[str(tx or "")].append(int(p))
    return index


def _isoform_table(rows: pd.DataFrame, register: dict[str, list[int]]) -> list[dict[str, Any]]:
    """One entry per transcript the guide was enumerated on, flagging register neighbours.

    A shorter guide starting at *s* sits inside the longer window at *s-1*, so two designs one
    nucleotide apart can share a window and receive an identical score. That much is a property of the
    enumeration and holds for any target.

    An earlier version of this report told the reader "two such designs differed 1.4x in measured
    knockdown" in the rendered card. That claim is now confined to this docstring, for two reasons.
    It is traceable but **mislabelled**: 1.4x is the ratio of *fraction remaining* between AZ's HD-001
    (0.49) and HD-002 (0.35) at transcript positions 1982/1983, and the ratio of *knockdown* for the
    same pair is 1.27x. And it is not the strongest case on that panel -- the pair at 2733/2735
    differs 1.97x in fraction remaining. Neither number belongs in a card that ships with the tool
    and is read against targets that panel says nothing about.
    """
    out: list[dict[str, Any]] = []
    for _, r in rows.iterrows():
        tx = str(r.get("transcript_id") or "")
        pos = _num(r.get("position"))
        neighbours = [
            p
            for p in register.get(tx, [])
            if pos is not None and p != int(pos) and abs(p - int(pos)) <= REGISTER_NEIGHBOUR_NT
        ]
        out.append(
            {
                "candidate_id": str(r.get("id") or ""),
                "transcript": tx,
                "position": pos,
                "register_neighbours": sorted(neighbours),
            }
        )
    out.sort(key=lambda d: (d["transcript"], d["position"] if d["position"] is not None else -1))
    return out


def _offtarget_views(
    hits: pd.DataFrame,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], int]:
    """Group hits by gene symbol, build the complete count matrix, and embed the in-scope rows."""
    if hits.empty:
        return [], [], [], 0

    liability = int(hits["_class"].isin({c.value for c in LIABILITY}).sum())

    grouped: Counter[tuple[str, str, str]] = Counter()
    min_nm: dict[tuple[str, str, str], float] = {}
    for _, r in hits.iterrows():
        symbol = str(r.get("hit_symbol") or r.get("matched_symbol") or "").strip() or _UNKNOWN_SYMBOL
        key = (symbol, str(r["_class"]), str(r.get("species") or "").strip() or "query")
        grouped[key] += 1
        nm = r["_nm"]
        if pd.notna(nm):
            min_nm[key] = min(min_nm.get(key, float(nm)), float(nm))
    by_symbol = [
        {
            "symbol": s,
            "hit_class": c,
            "species": sp,
            "n": n,
            "min_nm": min_nm.get((s, c, sp)),
            "is_liability": c in {x.value for x in LIABILITY},
        }
        for (s, c, sp), n in grouped.most_common()
    ]

    matrix_counter: Counter[tuple[str, str, str]] = Counter()
    for _, r in hits.iterrows():
        nm = r["_nm"]
        band = "unknown" if pd.isna(nm) else ("0" if nm == 0 else "1" if nm == 1 else "2" if nm == 2 else ">=3")
        matrix_counter[(str(r.get("species") or "query"), band, str(r["_class"]))] += 1
    matrix = [
        {"species": sp, "nm": band, "hit_class": cls, "n": n} for (sp, band, cls), n in sorted(matrix_counter.items())
    ]

    # Row-level drill-down exists to inspect liabilities, so only liability rows in the embedded
    # scope carry per-row detail. On-target isoform, ortholog and repeat alignments are still fully
    # displayed -- by gene in offtarget_by_symbol and completely in offtarget_matrix -- they just do
    # not each need a row. Embedding every class cost 14.7 MB of the first render's 38 MB.
    in_scope = hits[
        (hits["species"].astype(str).str.strip().isin({EMBED_SPECIES, ""}))
        & (hits["_nm"].notna())
        & (hits["_nm"] <= EMBED_MAX_NM)
        & (hits["_class"].isin({c.value for c in LIABILITY}))
    ]
    embedded = [
        {
            "t": str(r.get("rname") or ""),
            "s": str(r.get("hit_symbol") or "").strip() or _UNKNOWN_SYMBOL,
            "c": str(r["_class"]),
            "nm": _num(r.get("nm")),
            "sm": _num(r.get("seed_mismatches")),
            "cig": str(r.get("cigar") or ""),
            "st": str(r.get("strand") or ""),
        }
        for _, r in in_scope.iterrows()
    ]
    embedded.sort(key=lambda d: (d["nm"] if d["nm"] is not None else 99, d["s"]))
    return by_symbol, matrix, embedded, liability


#: Columns the miRNA name may arrive under. ``mirna_id`` is what the aggregate actually publishes
#: (``Hsa-Mir-24-P2_3p``); the table has no ``rname``, so reading only that returned an empty string
#: for every row and the panel listed 18,078 anonymous seed matches. Which miRNA is mimicked is the
#: whole question -- a perfect seed match to a cardiac or neuronal family is a different finding from
#: one to an unexpressed paralogue, and a retraction in this programme turned on exactly that.
_MIRNA_NAME_COLUMNS = ("mirna_id", "rname", "mirna", "name")


def _ortholog_conservation(hits: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Best ortholog alignment per species: how conserved this guide's site is, per species.

    Cross-species conservation is already in the screen and is thrown away twice over. The screen
    aligns each guide against every requested species' transcriptome and classifies a hit on the
    orthologous gene as ``ortholog`` -- explicitly *not* a liability -- and the candidate table then
    summarises those hits as ``conservation_score``, which is ``(species hit) / 3``.

    That summary cannot answer the question a cross-species programme asks, for two reasons. It is
    **mismatch-blind**: on one MSH3 run it counted a species as conserved on alignments up to 8
    mismatches, and mouse ortholog hits ran 1,413 at nm=0 against 1,493 at nm>=3. And it is
    **seed-blind**: mouse nm=1 split 96 seed-intact against 79 seed-hit, and one mismatch outside
    positions 2-8 is a different molecule from one inside them -- allowing it took the
    mouse-and-macaque pool from 49 guides to 113.

    So the per-species best ``(nm, seed_mismatches)`` is published instead and the thresholding is
    left to the reader. Requiring "perfect in mouse and macaque" is a programme decision, not a
    property of the target, and hard-coding it as a gate would put one programme's requirement in
    every run.
    """
    if hits.empty or "_class" not in hits.columns:
        return {}
    ortholog = hits[hits["_class"].eq(HitClass.ORTHOLOG.value)]
    if ortholog.empty:
        return {}
    out: dict[str, dict[str, Any]] = {}
    for species, group in ortholog.groupby(ortholog["species"].astype(str).str.strip()):
        nm = pd.to_numeric(group["nm"], errors="coerce")
        seed = pd.to_numeric(group.get("seed_mismatches"), errors="coerce")
        order = nm.fillna(99) * 100 + seed.fillna(99)
        best = order.idxmin()
        out[species or "query"] = {"nm": _num(nm.get(best)), "seed_mismatches": _num(seed.get(best))}
    return out


def _mirna_table(mirna: pd.DataFrame) -> list[dict[str, Any]]:
    """Named miRNA seed hits at 0 and 1 seed mismatch, with the database and match coordinate."""
    if mirna.empty:
        return []
    name_column = next((c for c in _MIRNA_NAME_COLUMNS if c in mirna.columns), None)
    out: list[dict[str, Any]] = []
    for _, r in mirna.iterrows():
        sm = _num(r.get("seed_mismatches"))
        if sm is not None and sm > 1:
            continue
        out.append(
            {
                "mirna": str(r.get(name_column) or "") if name_column else "",
                "source": str(r.get("species") or r.get("source") or ""),
                "database": str(r.get("database") or ""),
                "seed_mismatches": sm,
                "nm": _num(r.get("nm")),
                # The seed offset. A v0.5.1 defect counted a guide-seed motif matching anywhere on a
                # miRNA as a perfect seed hit, and coord != 1 is what distinguished the artifact from
                # the real thing, so the report shows it rather than asking a reader to trust the fix.
                "coord": _num(r.get("coord")),
            }
        )
    out.sort(key=lambda d: (d["seed_mismatches"] if d["seed_mismatches"] is not None else 99, d["mirna"]))
    return out
