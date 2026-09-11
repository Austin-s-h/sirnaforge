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
from collections import Counter, defaultdict
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from sirnaforge.config.run_policy import EntryPoint, ResolvedRunPolicy, resolve_run_policy
from sirnaforge.core.hit_annotation import CLASSIFICATION_COLUMNS, hit_class_of, is_annotated
from sirnaforge.core.hit_classification import HitClass
from sirnaforge.models.policy import FilterAction, FilterComparator, FilterEvaluation

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
        """
        if self.n_gates_failed:
            return "fail"
        if self.n_gates_unknown:
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
    """A number, or None for NaN/blank -- never a fabricated zero."""
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, str):
        value = value.strip()
        if not value:
            return None
        try:
            value = float(value)
        except ValueError:
            return None
    if not isinstance(value, (int, float)) or pd.isna(value):
        return None
    return int(value) if float(value).is_integer() else round(float(value), 6)


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


def _evaluate(descriptor: Any, row: pd.Series) -> tuple[float | int | None, int, int]:
    """Evaluate one descriptor against one candidate row, independently of every other gate.

    Returns ``unknown`` rather than a pass when the run does not export the column the threshold
    reads. Six of the twelve active gates read human-stratified counters that 0.7.1 does not export
    yet (#101), and reporting those as passes would be the fabricated-evidence failure this report
    exists to make visible.

    Returns:
        ``(value, verdict_code, reason_code)``. The descriptor itself is emitted once per report.
    """
    if descriptor.action.value == FilterAction.OFF.value:
        return None, _VERDICT_CODE[FilterEvaluation.NOT_EVALUATED.value], REASON_FILTER_OFF
    if descriptor.threshold is None:
        return None, _VERDICT_CODE[FilterEvaluation.NOT_EVALUATED.value], REASON_NO_THRESHOLD
    if descriptor.column not in row.index:
        return None, _VERDICT_CODE[FilterEvaluation.UNKNOWN.value], REASON_MISSING_COLUMN

    value = _num(row.get(descriptor.column))
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
    return {
        "comparable": True,
        "run_pass_guides": len(run_pass),
        "contradicted_run_pass": len(contradicted),
        "run_failed_not_rederivable": len(not_rederivable),
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


def build_payload(run_dir: Path | str, *, policy: ResolvedRunPolicy | None = None) -> ReportPayload:
    """Build the payload for a finished run directory.

    Args:
        run_dir: A completed output directory. No live pipeline state is required.
        policy: Resolved policy supplying the gate descriptors. Resolved from the run's own
            parameters when omitted.

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

    if policy is None:
        policy = resolve_run_policy(
            entry_point=EntryPoint.SCREENING_WORKFLOW, query_species="human", screen_species=["human"]
        )

    candidates["_guide"] = candidates["guide_sequence"].map(_normalise_guide)
    if not hits.empty:
        hits["_guide"] = hits["qseq"].map(_normalise_guide)
        hits["_nm"] = pd.to_numeric(hits["nm"], errors="coerce")
        hits["_class"] = hits.apply(lambda r: hit_class_of(r).value, axis=1)
    if not mirna.empty and "qseq" in mirna.columns:
        mirna["_guide"] = mirna["qseq"].map(_normalise_guide)

    hits_by_guide = dict(tuple(hits.groupby("_guide"))) if not hits.empty else {}
    mirna_by_guide = dict(tuple(mirna.groupby("_guide"))) if not mirna.empty and "_guide" in mirna else {}

    descriptors = [f.descriptor for f in policy.filters]
    register = _register_index(candidates)
    guides: list[GuideEntry] = []

    for guide, rows in candidates.groupby("_guide"):
        best = rows.iloc[0]
        guides.append(
            _build_guide(
                guide=str(guide),
                rows=rows,
                best=best,
                descriptors=descriptors,
                hits=hits_by_guide.get(guide, pd.DataFrame()),
                mirna=mirna_by_guide.get(guide, pd.DataFrame()),
                register=register,
            )
        )

    guides.sort(key=lambda g: (-(g.composite_score or g.design_score or -1), g.guide))

    manifest_path = next(iter(sorted(run_dir.glob("**/manifest.json"))), None)
    manifest: dict[str, Any] = {}
    if manifest_path is not None:
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:  # a bad manifest costs provenance, not the report
            caveats.append(f"manifest.json could not be read ({exc}); the header shows less provenance")
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
        "status_counts": {s: sum(1 for g in guides if g.status == s) for s in ("pass", "unknown", "fail")},
        "agreement": agreement,
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
            }
            for f in policy.filters
        ],
        guides=guides,
        provenance={
            "run_dir": str(run_dir),
            "candidates_csv": str(candidates_csv.relative_to(run_dir)),
            "manifest_present": manifest_path is not None,
            "tool_version": str(manifest.get("tool_version") or "unknown"),
            "run_timestamp": str(manifest.get("run_timestamp") or "unknown"),
            "policy_profile": policy.profile.name,
            "run_mode": policy.run_mode.value,
            "payload_schema_version": PAYLOAD_SCHEMA_VERSION,
        },
        caveats=caveats,
    )


def _build_guide(
    *,
    guide: str,
    rows: pd.DataFrame,
    best: pd.Series,
    descriptors: list[Any],
    hits: pd.DataFrame,
    mirna: pd.DataFrame,
    register: dict[str, list[int]],
) -> GuideEntry:
    gates = [list(_evaluate(d, best)) for d in descriptors]

    isoforms = _isoform_table(rows, register)
    by_symbol, matrix, embedded, liability = _offtarget_views(hits)

    counts_exist = bool(len(hits)) and not embedded
    return GuideEntry(
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
    nucleotide apart can share a window and receive an identical score. On the reference panel two
    such designs differed 1.4x in measured knockdown, so the report says so rather than letting them
    look like duplicates.
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


def _mirna_table(mirna: pd.DataFrame) -> list[dict[str, Any]]:
    """Named miRNA seed hits at 0 and 1 seed mismatch, by source database."""
    if mirna.empty:
        return []
    out: list[dict[str, Any]] = []
    for _, r in mirna.iterrows():
        sm = _num(r.get("seed_mismatches"))
        if sm is not None and sm > 1:
            continue
        out.append(
            {
                "mirna": str(r.get("rname") or r.get("mirna") or ""),
                "source": str(r.get("species") or r.get("source") or ""),
                "seed_mismatches": sm,
                "nm": _num(r.get("nm")),
            }
        )
    out.sort(key=lambda d: (d["seed_mismatches"] if d["seed_mismatches"] is not None else 99, d["mirna"]))
    return out
