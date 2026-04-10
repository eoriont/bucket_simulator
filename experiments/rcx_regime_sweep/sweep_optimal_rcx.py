#!/usr/bin/env python3
"""
Sweep explicit distillation protocol / k / m settings over entanglement rates
and identify RCX regime changes.

This is a light-weight experiment driver: it uses `bucket_simulator -dump-circuit`
to force circuit/noise-model construction, then parses the remote-noise summary
from stderr instead of running Monte Carlo shots.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = EXPERIMENT_DIR / "config.txt"
DEFAULT_SIMULATOR = PROJECT_ROOT / "build" / "bucket_simulator"

MODES = ("xx_merge", "zz_merge", "xx_split", "zz_split")
PROTOCOLS = ("none", "2to1_pumping", "3to1_pumping", "2to1_recurrence", "3to1_recurrence")

REMOTE_PATTERNS = {
    "requested_protocol": re.compile(r"Requested distillation protocol:\s*(\S+)"),
    "effective_protocol": re.compile(r"Effective distillation protocol:\s*(\S+)"),
    "distillation_feasible": re.compile(r"Distillation feasible:\s*(YES|NO)"),
    "distillation_rounds": re.compile(r"Distillation rounds:\s*(\d+)"),
    "backup_batches": re.compile(r"Backup batches \(m\):\s*(\d+)"),
    "remote_cnot_error": re.compile(r"Remote CNOT error:\s*([0-9eE+\-.]+)"),
    "effective_channel_error": re.compile(r"Effective channel error:\s*([0-9eE+\-.]+)"),
    "distillation_qubit_overhead_feasible": re.compile(r"Distillation qubit overhead feasible:\s*(YES|NO)"),
    "distillation_epr_slots_required": re.compile(r"(?:Concurrent )?distillation EPR slots required:\s*(\d+)", re.IGNORECASE),
    "distillation_qubits_required": re.compile(r"Distillation qubits required:\s*(\d+)"),
    "monolithic_equivalent_qubits": re.compile(r"Monolithic equivalent qubits:\s*(\d+)"),
    "remote_cnots_per_cycle": re.compile(r"Remote CNOTs per cycle:\s*(\d+)"),
    "raw_pairs_per_distilled": re.compile(r"Raw pairs per distilled:\s*(\d+)"),
    "distillation_time_ns": re.compile(r"Distillation time:\s*([0-9eE+\-.]+)\s*ns"),
    "measurement_delay_ns": re.compile(r"Measurement delay:\s*([0-9eE+\-.]+)\s*ns"),
    "timing_constraint_satisfied": re.compile(r"Timing constraint satisfied:\s*(YES|NO)"),
}


@dataclass(frozen=True)
class ModeSpec:
    merge_type: str
    experiment_phase: str


MODE_MAP = {
    "xx_merge": ModeSpec("distributed_xx", "merge_only"),
    "zz_merge": ModeSpec("distributed_zz", "merge_only"),
    "xx_split": ModeSpec("distributed_xx", "split_only"),
    "zz_split": ModeSpec("distributed_zz", "split_only"),
}


def parse_bool_token(token: str) -> bool:
    return token.upper() == "YES"


def rate_grid(rate_min_hz: int, rate_max_hz: int, points_per_decade: int) -> list[int]:
    if rate_min_hz <= 0 or rate_max_hz < rate_min_hz:
        raise ValueError("invalid rate range")
    if rate_min_hz == rate_max_hz:
        return [rate_min_hz]
    log_min = math.log10(rate_min_hz)
    log_max = math.log10(rate_max_hz)
    num_points = max(2, int(round((log_max - log_min) * points_per_decade)) + 1)
    values = sorted({int(round(10 ** x)) for x in [log_min + (log_max - log_min) * i / (num_points - 1) for i in range(num_points)]})
    values[0] = rate_min_hz
    values[-1] = rate_max_hz
    return values


def protocol_candidates(max_rounds: int, max_batches: int) -> Iterable[tuple[str, int, int]]:
    for protocol in PROTOCOLS:
        if protocol == "none":
            yield protocol, 0, 1
            continue
        for rounds in range(1, max_rounds + 1):
            for batches in range(1, max_batches + 1):
                yield protocol, rounds, batches


def parse_summary(stderr_text: str) -> dict[str, object]:
    parsed: dict[str, object] = {}
    for key, pattern in REMOTE_PATTERNS.items():
        match = pattern.search(stderr_text)
        if not match:
            continue
        value = match.group(1)
        if key in {
            "distillation_rounds",
            "backup_batches",
            "distillation_epr_slots_required",
            "distillation_qubits_required",
            "monolithic_equivalent_qubits",
            "remote_cnots_per_cycle",
            "raw_pairs_per_distilled",
        }:
            parsed[key] = int(value)
        elif key in {
            "remote_cnot_error",
            "effective_channel_error",
            "distillation_time_ns",
            "measurement_delay_ns",
        }:
            parsed[key] = float(value)
        elif key in {
            "distillation_feasible",
            "distillation_qubit_overhead_feasible",
            "timing_constraint_satisfied",
        }:
            parsed[key] = parse_bool_token(value)
        else:
            parsed[key] = value
    return parsed


def evaluate_candidate(
    simulator: Path,
    config: Path,
    mode: str,
    enr_hz: int,
    protocol: str,
    rounds: int,
    batches: int,
) -> dict[str, object]:
    mode_spec = MODE_MAP[mode]
    with tempfile.TemporaryDirectory(prefix="rcx_regime_") as tmpdir:
        cmd = [
            str(simulator),
            "-config", str(config),
            "-dump-circuit",
            "-output", tmpdir,
            "-merge_type", mode_spec.merge_type,
            "-experiment_phase", mode_spec.experiment_phase,
            "-entanglement_rate", str(enr_hz),
            "-distillation_protocol", protocol,
            "-distillation_rounds", str(rounds),
            "-distillation_backup_batches", str(batches),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"simulator failed for {mode=} {enr_hz=} {protocol=} {rounds=} {batches=}\n{proc.stderr}")

    parsed = parse_summary(proc.stderr)
    parsed.update({
        "mode": mode,
        "entanglement_rate_hz": enr_hz,
        "requested_protocol": protocol,
        "requested_rounds": rounds,
        "requested_backup_batches": batches,
    })
    parsed["fully_feasible"] = bool(
        parsed.get("timing_constraint_satisfied", False) and
        parsed.get("distillation_qubit_overhead_feasible", False)
    )
    return parsed


def choose_best(rows: list[dict[str, object]]) -> dict[str, object]:
    def sort_key(row: dict[str, object]) -> tuple:
        feasible = row.get("fully_feasible", False)
        remote_cnot_error = float(row.get("remote_cnot_error", float("inf")))
        protocol = str(row.get("effective_protocol", row.get("requested_protocol", "none")))
        rounds = int(row.get("requested_rounds", 0))
        batches = int(row.get("requested_backup_batches", 1))
        return (
            0 if feasible else 1,
            remote_cnot_error,
            rounds,
            batches,
            protocol,
        )

    return min(rows, key=sort_key)


def regime_signature(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row.get("effective_protocol", row.get("requested_protocol")),
        row.get("requested_rounds"),
        row.get("requested_backup_batches"),
        row.get("fully_feasible"),
    )


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser(description="Sweep protocol/k/m over ENR and find RCX regimes.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG, help="Base config file")
    parser.add_argument("--simulator", type=Path, default=DEFAULT_SIMULATOR, help="Path to bucket_simulator")
    parser.add_argument("--output-dir", type=Path, default=EXPERIMENT_DIR / "runs" / "manual", help="Output directory")
    parser.add_argument("--modes", nargs="+", default=list(MODES), choices=MODES, help="Modes to sweep")
    parser.add_argument("--rate-min-hz", type=int, default=1_000, help="Minimum ENR in Hz")
    parser.add_argument("--rate-max-hz", type=int, default=100_000, help="Maximum ENR in Hz")
    parser.add_argument("--points-per-decade", type=int, default=12, help="Log-spaced grid density")
    parser.add_argument("--max-rounds", type=int, default=5, help="Maximum k to test")
    parser.add_argument("--max-backup-batches", type=int, default=8, help="Maximum m to test")
    args = parser.parse_args()

    if not args.config.exists():
        raise SystemExit(f"Config not found: {args.config}")
    if not args.simulator.exists():
        raise SystemExit(f"Simulator not found: {args.simulator}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rates = rate_grid(args.rate_min_hz, args.rate_max_hz, args.points_per_decade)

    all_rows: list[dict[str, object]] = []
    best_rows: list[dict[str, object]] = []
    regime_rows: list[dict[str, object]] = []

    for mode in args.modes:
        previous_signature: tuple[object, ...] | None = None
        for enr_hz in rates:
            candidates = [
                evaluate_candidate(
                    args.simulator,
                    args.config,
                    mode,
                    enr_hz,
                    protocol,
                    rounds,
                    batches,
                )
                for protocol, rounds, batches in protocol_candidates(args.max_rounds, args.max_backup_batches)
            ]
            all_rows.extend(candidates)
            best = choose_best(candidates)
            best_rows.append(best)
            signature = regime_signature(best)
            if signature != previous_signature:
                regime_rows.append({
                    "mode": mode,
                    "entanglement_rate_hz": enr_hz,
                    "effective_protocol": best.get("effective_protocol", best.get("requested_protocol")),
                    "requested_rounds": best.get("requested_rounds"),
                    "requested_backup_batches": best.get("requested_backup_batches"),
                    "fully_feasible": best.get("fully_feasible"),
                    "remote_cnot_error": best.get("remote_cnot_error"),
                    "timing_constraint_satisfied": best.get("timing_constraint_satisfied"),
                    "distillation_qubit_overhead_feasible": best.get("distillation_qubit_overhead_feasible"),
                })
                previous_signature = signature

    all_fields = [
        "mode",
        "entanglement_rate_hz",
        "requested_protocol",
        "effective_protocol",
        "requested_rounds",
        "requested_backup_batches",
        "remote_cnot_error",
        "effective_channel_error",
        "distillation_feasible",
        "timing_constraint_satisfied",
        "distillation_qubit_overhead_feasible",
        "fully_feasible",
        "distillation_epr_slots_required",
        "distillation_qubits_required",
        "monolithic_equivalent_qubits",
        "remote_cnots_per_cycle",
        "raw_pairs_per_distilled",
        "distillation_time_ns",
        "measurement_delay_ns",
    ]
    best_fields = [
        "mode",
        "entanglement_rate_hz",
        "effective_protocol",
        "requested_rounds",
        "requested_backup_batches",
        "remote_cnot_error",
        "fully_feasible",
        "timing_constraint_satisfied",
        "distillation_qubit_overhead_feasible",
        "distillation_qubits_required",
        "monolithic_equivalent_qubits",
    ]

    write_csv(args.output_dir / "all_candidates.csv", all_rows, all_fields)
    write_csv(args.output_dir / "best_by_rate.csv", best_rows, best_fields)
    write_csv(args.output_dir / "regime_boundaries.csv", regime_rows, best_fields)

    summary_path = args.output_dir / "summary.txt"
    with summary_path.open("w") as f:
        f.write("RCX Regime Sweep Summary\n")
        f.write("=======================\n")
        f.write(f"config: {args.config}\n")
        f.write(f"simulator: {args.simulator}\n")
        f.write(f"modes: {' '.join(args.modes)}\n")
        f.write(f"rates_hz: {rates}\n")
        f.write(f"max_rounds: {args.max_rounds}\n")
        f.write(f"max_backup_batches: {args.max_backup_batches}\n\n")
        f.write("Regime boundaries:\n")
        for row in regime_rows:
            f.write(
                f"  mode={row['mode']} enr={row['entanglement_rate_hz']}Hz "
                f"protocol={row['effective_protocol']} k={row['requested_rounds']} "
                f"m={row['requested_backup_batches']} rcx={row['remote_cnot_error']} "
                f"feasible={row['fully_feasible']}\n"
            )

    print(f"Wrote all candidates to: {args.output_dir / 'all_candidates.csv'}")
    print(f"Wrote best-by-rate table to: {args.output_dir / 'best_by_rate.csv'}")
    print(f"Wrote regime boundaries to: {args.output_dir / 'regime_boundaries.csv'}")
    print(f"Wrote summary to: {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
