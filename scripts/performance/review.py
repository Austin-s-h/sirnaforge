#!/usr/bin/env python3
"""Isolated performance investigations; no installation or default configuration changes.

Run with Python 3.11+ (or Python 3.10 with tomli) and --help. Results record failures,
commands and checksums. Warm/cold labels describe user preparation, not OS cache eviction.
"""

# ruff: noqa: D103

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import signal
import statistics
import subprocess
import tempfile
import time
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

import ahocorasick

try:
    import tomllib
except ImportError:
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[2]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
TOOLS = {
    "pixi": ["--version"],
    "uv": ["--version"],
    "micromamba": ["--version"],
    "nextflow": ["-version"],
    "java": ["-version"],
    "bwa-mem2": ["version"],
    "bwa-mem3": ["version"],
    "samtools": ["--version"],
    "RNAfold": ["--version"],
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def machine() -> dict:
    return {
        "platform": platform.platform(),
        "architecture": platform.machine(),
        "cpu_count": os.cpu_count(),
        "python": platform.python_version(),
        "git_revision": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, capture_output=True, check=False
        ).stdout.strip(),
    }


def bwa_comparison_succeeded(result: dict, allow_sam_differences: bool = False) -> bool:
    """Require complete execution and matching references; optionally report SAM differences."""
    parity = result.get("parity")
    backends = result.get("backends", {})
    return (
        bool(backends)
        and all(backend.get("status") == "ok" for backend in backends.values())
        and parity is not None
        and parity.get("reference_headers_equal") is True
        and (allow_sam_differences or parity.get("equal") is True)
    )


def run(command: list[str], timeout: float, stdout_path: Path | None = None) -> dict:
    """Measure one subprocess, with GNU time resource metrics when installed.

    stdout is bounded in metadata; alignment SAM goes directly to a file. A timed-out
    process group is killed so workers do not silently survive the benchmark.
    """
    result = {"command": command, "timeout_seconds": timeout}
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="sirnaforge-perf-") as temp:
        metrics = Path(temp) / "resources.txt"
        timer = Path("/usr/bin/time")
        # GNU time is available on Linux; BSD time uses incompatible flags.
        timed = platform.system() == "Linux" and timer.exists()
        argv = [str(timer), "-f", "%U %S %M", "-o", str(metrics), *command] if timed else command
        output_file = stdout_path.open("w") if stdout_path else None
        output = output_file if output_file is not None else subprocess.PIPE
        try:
            try:
                process = subprocess.Popen(
                    argv, stdout=output, stderr=subprocess.PIPE, text=True, start_new_session=True
                )
                try:
                    stdout, stderr = process.communicate(timeout=timeout)
                    result["status"] = "ok" if process.returncode == 0 else "failed"
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    stdout, stderr = process.communicate()
                    result["status"] = "timeout"
                result.update(returncode=process.returncode, stdout=(stdout or "")[-4000:], stderr=stderr[-4000:])
            except OSError as error:
                result.update(status="unavailable", error=str(error))
        finally:
            if output_file is not None:
                output_file.close()
        if metrics.exists():
            lines = metrics.read_text().splitlines()
            fields = lines[-1].split() if lines else []
            if len(fields) == 3:
                result.update(
                    cpu_user_seconds=float(fields[0]), cpu_system_seconds=float(fields[1]), peak_rss_kib=int(fields[2])
                )
    result["wall_seconds"] = time.perf_counter() - started
    return result


def inventory(timeout: float) -> dict:
    tools = {}
    for name, args in TOOLS.items():
        executable = shutil.which(name)
        tools[name] = {"path": executable, "status": "unavailable"}
        if executable:
            tools[name] = {"path": executable, **run([executable, *args], timeout)}
    paths = ["pyproject.toml", "uv.lock", "environment-dev.yml", "docker/environment-nextflow.yml", "docker/Dockerfile"]
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    return {
        "kind": "inventory",
        "machine": machine(),
        "tools": tools,
        "files": {name: {"sha256": digest(ROOT / name)} for name in paths},
        "python_dependencies": project["project"]["dependencies"],
        "dev_dependencies": project["dependency-groups"]["dev"],
        "conda_environment_sources": {name: (ROOT / name).read_text() for name in paths if name.endswith(".yml")},
        "caveat": "Tool version probes are PATH observations, not successful workflow execution. Nextflow may try network access; probe is bounded by timeout.",
    }


def sam_records(path: Path) -> Counter:
    """Compare SAM semantics as a multiset, ignoring record/tag order and PG identifiers.

    Coordinates, flags, MAPQ, CIGAR, query sequence and alignment tags are retained.
    PG is program provenance, recorded independently through binary checksum/version.
    Header reference names/lengths are compared separately by sam_headers().
    """
    records: Counter[tuple[str, ...]] = Counter()
    for line in path.read_text().splitlines():
        if line.startswith("@") or not line:
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            raise ValueError(f"Malformed SAM record in {path}: {line[:120]}")
        tags = sorted(field for field in fields[11:] if not field.startswith("PG:"))
        records[tuple(fields[:11] + tags)] += 1
    return records


def sam_headers(path: Path) -> list:
    return sorted(sorted(line.split("\t")[1:]) for line in path.read_text().splitlines() if line.startswith("@SQ\t"))


def fasta_records(path: Path) -> dict[str, str]:
    records = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError(f"Duplicate FASTA ID: {name}")
            records[name] = ""
        elif line.strip():
            if name is None:
                raise ValueError("FASTA sequence before header")
            records[name] += line.strip().upper()
    return records


def iter_fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    """Yield FASTA records without retaining the complete reference in memory."""
    name = None
    chunks: list[str] = []
    seen: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                if name in seen:
                    raise ValueError(f"Duplicate FASTA ID: {name}")
                seen.add(name)
                chunks = []
            elif line.strip():
                if name is None:
                    raise ValueError("FASTA sequence before header")
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def exact_oracle(reference: Path, queries: Path, sam: Path) -> dict:
    """Independently enumerate full-length exact matches, including reverse/repeated loci."""
    guides = fasta_records(queries)
    patterns: dict[str, list[tuple[str, str]]] = {}
    complement = str.maketrans("ACGT", "TGCA")
    for qname, query in guides.items():
        if not set(query) <= set("ACGT"):
            continue
        patterns.setdefault(query, []).append((qname, "+"))
        patterns.setdefault(query.translate(complement)[::-1], []).append((qname, "-"))
    automaton = ahocorasick.Automaton()
    for pattern, targets in patterns.items():
        automaton.add_word(pattern, (pattern, targets))
    automaton.make_automaton()
    expected: set[tuple[str, str, int, str]] = set()
    for rname, sequence in iter_fasta_records(reference):
        for end, (pattern, targets) in automaton.iter(sequence):
            position = end - len(pattern) + 2
            expected.update((qname, rname, position, strand) for qname, strand in targets)
    observed = set()
    for record in sam_records(sam):
        qname, flag, rname, position, _, cigar = record[:6]
        if not int(flag) & 4 and cigar == f"{len(guides[qname])}M" and "NM:i:0" in record[11:]:
            observed.add((qname, rname, int(position), "-" if int(flag) & 16 else "+"))
    return {
        "equal": expected == observed,
        "expected_count": len(expected),
        "observed_count": len(observed),
        "missing": sorted(expected - observed),
        "unexpected": sorted(observed - expected),
        "scope": "Full-length ungapped exact matches only; other mismatch/gap placements need reviewed truth.",
    }


def differences(left: Counter, right: Counter) -> dict:
    def sample(records: Counter) -> list:
        return [{"record": list(key), "count": count} for key, count in sorted(records.items())[:20]]

    return {
        "equal": left == right,
        "left_records": sum(left.values()),
        "right_records": sum(right.values()),
        "left_only_count": sum((left - right).values()),
        "right_only_count": sum((right - left).values()),
        "left_only_sample": sample(left - right),
        "right_only_sample": sample(right - left),
    }


def bwa(args: argparse.Namespace) -> dict:  # noqa: PLR0912 - explicit backend-specific acquisition and parity stages
    args.work_dir.mkdir(parents=True, exist_ok=True)
    query_records = fasta_records(args.queries)
    # Prevent stale outputs or a prior backend's index from contaminating a run.
    session = Path(tempfile.mkdtemp(prefix="bwa-review-", dir=args.work_dir))
    result = {
        "kind": "bwa_comparison",
        "machine": machine(),
        "session": str(session),
        "threads": args.threads,
        "mem3_mode": args.mem3_mode,
        "repetitions": args.repetitions,
        "query_count": len(query_records),
        "query_length_counts": dict(sorted(Counter(map(len, query_records.values())).items())),
        "cache_state": "fresh per-backend index unless MEM2 reuse is requested; OS page cache uncontrolled; alignment repeats reuse index",
        "inputs": {str(path): digest(path) for path in (args.reference, args.queries)},
        "backends": {},
        "decision": "defer",
        "parity": None,
        "limitations": [
            "SAM parity alone does not prove exhaustive biological correctness or downstream gate parity.",
            "Alignment repetitions reuse the index built during this invocation; run separate invocations for index-build repetitions.",
            "Throughput is query count divided by alignment wall time and includes index loading; OS page cache is uncontrolled.",
            "No cross-backend index reuse is attempted; format compatibility remains unproven.",
        ],
    }
    if args.mem2_index:
        if not args.mem2_index.is_file() or digest(args.mem2_index) != digest(args.reference):
            raise ValueError("--mem2-index must be a reference FASTA with bytes identical to --reference")
        for suffix in (".0123", ".bwt.2bit.64", ".amb", ".ann", ".pac"):
            member = Path(f"{args.mem2_index}{suffix}")
            if not member.is_file() or member.stat().st_size == 0:
                raise ValueError(f"--mem2-index lacks a nonempty BWA-MEM2 index member: {member}")
        result["mem2_reused_index"] = str(args.mem2_index.resolve())
    sams = {}
    for label, requested in (("mem2", args.mem2), ("mem3", args.mem3)):
        executable = shutil.which(requested)
        if executable is None:
            result["backends"][label] = {"status": "unavailable", "requested": requested}
            continue
        executable = str(Path(executable).resolve())
        folder = session / label
        folder.mkdir()
        reuse_mem2 = label == "mem2" and args.mem2_index is not None
        reference = args.mem2_index.resolve() if reuse_mem2 else folder / "reference.fasta"
        if not reuse_mem2:
            shutil.copyfile(args.reference, reference)
        backend: dict = {
            "path": executable,
            "binary_sha256": digest(Path(executable)),
            "version": run([executable, "version"], args.timeout),
        }
        result["backends"][label] = backend
        index_command = [executable, "index"]
        if label == "mem3":
            index_command += ["-t", str(args.threads)]
        index_command.append(str(reference))
        backend["index"] = (
            {"status": "ok", "reused": True, "reference_sha256": digest(reference)}
            if reuse_mem2
            else run(["env", f"OMP_NUM_THREADS={args.threads}", *index_command], args.timeout)
        )
        if backend["index"]["status"] != "ok":
            backend["status"] = "index_failed"
            continue
        index_folder = reference.parent if reuse_mem2 else folder
        index_suffixes = {".0123", ".bwt.2bit.64", ".amb", ".ann", ".pac"}
        backend["index_members"] = {
            path.name: {"bytes": path.stat().st_size, "sha256": digest(path)}
            for path in sorted(index_folder.iterdir())
            if path.name in {f"{reference.name}{suffix}" for suffix in index_suffixes}
        }
        backend["index_size_bytes"] = sum(member["bytes"] for member in backend["index_members"].values())
        trials = []
        baseline = None
        for repetition in range(args.repetitions):
            sam = folder / f"alignment-{repetition}.sam"
            command = [
                executable,
                "mem",
                "-t",
                str(args.threads),
                "-a",
                "-k",
                "12",
                "-T",
                "15",
                "-w",
                "100",
                str(reference),
                str(args.queries.resolve()),
            ]
            if label == "mem3" and args.mem3_mode == "compat":
                command.insert(2, "--compat=bwa-mem2")
            trial = run(command, args.timeout, sam)
            trials.append(trial)
            if trial["status"] == "ok":
                current = sam_records(sam)
                if baseline is None:
                    baseline = current
                    sams[label] = sam
                elif current != baseline:
                    trial["status"] = "nondeterministic_output"
        backend["alignments"] = trials
        if label in sams:
            expected_queries = set(query_records)
            observed_queries = {record[0] for record in sam_records(sams[label])}
            backend["query_membership"] = {
                "complete": expected_queries == observed_queries,
                "missing": sorted(expected_queries - observed_queries),
                "unexpected": sorted(observed_queries - expected_queries),
            }
            backend["exact_match_oracle"] = exact_oracle(args.reference, args.queries, sams[label])
        backend["status"] = "ok" if all(t["status"] == "ok" for t in trials) else "alignment_failed"
        if label in sams and not backend["query_membership"]["complete"]:
            backend["status"] = "incomplete_query_membership"
        successful = [trial["wall_seconds"] for trial in trials if trial["status"] == "ok"]
        if successful:
            backend["wall_summary"] = {
                "median": statistics.median(successful),
                "min": min(successful),
                "max": max(successful),
                "samples": len(successful),
            }
            backend["median_query_throughput_per_second"] = len(query_records) / backend["wall_summary"]["median"]
    if len(sams) == 2 and all(backend["status"] == "ok" for backend in result["backends"].values()):
        result["parity"] = differences(sam_records(sams["mem2"]), sam_records(sams["mem3"]))
        result["parity"]["reference_headers_equal"] = sam_headers(sams["mem2"]) == sam_headers(sams["mem3"])
    result["reason"] = (
        "Adoption requires representative benchmarks, downstream evidence/gate parity, and integration checks."
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=30)
    commands = parser.add_subparsers(dest="kind", required=True)
    commands.add_parser("inventory")
    measure = commands.add_parser("measure", help="Measure an explicit argv; use -- before the command")
    measure.add_argument("--cache-state", required=True)
    measure.add_argument("--repetitions", type=int, default=3)
    measure.add_argument("command", nargs=argparse.REMAINDER)
    align = commands.add_parser("bwa")
    align.add_argument("--mem2", default="bwa-mem2")
    align.add_argument("--mem3", default="bwa-mem3")
    align.add_argument("--mem3-mode", choices=["plain", "compat"], default="plain")
    align.add_argument("--reference", type=Path, default=FIXTURES / "reference.fasta")
    align.add_argument("--queries", type=Path, default=FIXTURES / "guides.fasta")
    align.add_argument(
        "--mem2-index", type=Path, help="Reuse a previously built MEM2 index whose FASTA matches --reference"
    )
    align.add_argument("--work-dir", type=Path, required=True)
    align.add_argument("--threads", type=int, default=1)
    align.add_argument("--repetitions", type=int, default=3)
    align.add_argument(
        "--allow-sam-differences",
        action="store_true",
        help="Report SAM differences without failing if both backends complete and reference headers match.",
    )
    args = parser.parse_args()
    if args.timeout <= 0 or getattr(args, "repetitions", 1) < 1 or getattr(args, "threads", 1) < 1:
        parser.error("timeout, repetitions and threads must be positive")
    if args.kind == "inventory":
        result = inventory(args.timeout)
    elif args.kind == "bwa":
        result = bwa(args)
    else:
        command = args.command[1:] if args.command[:1] == ["--"] else args.command
        if not command:
            parser.error("measure requires a command")
        trials = [run(command, args.timeout) for _ in range(args.repetitions)]
        success = [trial["wall_seconds"] for trial in trials if trial["status"] == "ok"]
        result = {
            "kind": "command_measurement",
            "machine": machine(),
            "cache_state": args.cache_state,
            "trials": trials,
            "wall_median_seconds": statistics.median(success) if success else None,
            "complete": len(success) == len(trials),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {args.output}")
    if args.kind == "bwa":
        return 0 if bwa_comparison_succeeded(result, args.allow_sam_differences) else 1
    if args.kind == "measure":
        return 0 if result["complete"] else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
