#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

import stim

from distributed_noise import DISTILLATION_PROTOCOLS
from paths import GENERATED_ROOT, INJECTED_ROOT, LER_ROOT, PROJECT_ROOT
from schedule_config import RoundSchedule, circuit_distance_from_k, default_round_schedule_for_k, round_schedule_to_dict


def parse_csv_ints(text: str) -> list[int]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        values.append(int(part))
    if not values:
        raise ValueError("Expected at least one integer value.")
    return values


def parse_csv_floats(text: str) -> list[float]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        values.append(float(part))
    if not values:
        raise ValueError("Expected at least one float value.")
    return values


def parse_csv_strings(text: str) -> list[str]:
    values = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        values.append(part)
    if not values:
        raise ValueError("Expected at least one string value.")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep entanglement rate, backup batches m, and distillation protocol. "
            "For each tuple, maximize circuit_k subject to the timing constraint "
            "(distillation time <= measurement delay), then run the noisy LER experiment."
        )
    )
    parser.add_argument(
        "--source",
        choices=("simple_merge_only", "simple_split_only", "simple_merge_split"),
        default="simple_merge_only",
    )
    parser.add_argument("--basis", choices=("X", "Z"), default="X")
    parser.add_argument(
        "--circuit-k-values",
        type=parse_csv_ints,
        default=[1, 2, 3, 4, 5],
        help="Candidate TQEC circuit_k values to test in ascending order.",
    )
    parser.add_argument(
        "--enr-values",
        type=parse_csv_floats,
        default=[25e6, 50e6, 100e6],
        help="Entanglement-rate values (Hz) to sweep.",
    )
    parser.add_argument(
        "--m-values",
        type=parse_csv_ints,
        default=[1, 2, 4],
        help="Backup batch counts m to sweep.",
    )
    parser.add_argument(
        "--protocols",
        type=parse_csv_strings,
        default=["none", "pumping_2to1", "pumping_3to1", "recurrence_2to1", "recurrence_3to1"],
        help="Comma-separated distillation protocols to sweep.",
    )
    parser.add_argument(
        "--distillation-rounds-values",
        type=parse_csv_ints,
        default=[1, 2, 3],
        help="Comma-separated distillation depth values to sweep (distinct from circuit_k).",
    )
    parser.add_argument(
        "--rounds-per-phase",
        type=int,
        default=None,
        help="Override the default pre/merge schedule with the same round count in each non-post phase.",
    )
    parser.add_argument("--shots", type=int, default=100_000)
    parser.add_argument("--batch-shots", type=int, default=250_000)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--physical-error", type=float, default=0.001)
    parser.add_argument("--measurement-delay", type=float, default=660e-9)
    parser.add_argument("--name", type=str, default="enr_protocol_sweep")
    return parser.parse_args()


def build_env(root: Path) -> dict[str, str]:
    for p in [root / ".tqec", root / ".mplconfig", root / ".cache", root / ".fontconfig"]:
        p.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["XDG_DATA_HOME"] = str(root / ".tqec")
    env["MPLCONFIGDIR"] = str(root / ".mplconfig")
    env["XDG_CACHE_HOME"] = str(root / ".cache")
    env["FONTCONFIG_PATH"] = "/opt/homebrew/etc/fonts"
    env["FONTCONFIG_FILE"] = "/opt/homebrew/etc/fonts/fonts.conf"
    return env


def ensure_clean_circuit(
    py: Path,
    env: dict[str, str],
    root: Path,
    source: str,
    basis: str,
    circuit_k: int,
    schedule: RoundSchedule,
) -> Path:
    base_clean = GENERATED_ROOT / f"{source}_{basis.lower()}_k{circuit_k}_r{schedule.merge_rounds}.stim"
    if not base_clean.exists():
        subprocess.run(
            [
                str(py),
                "python/tqec_import/export_stim.py",
                "--source",
                source,
                "--basis",
                basis,
                "--k",
                str(circuit_k),
                "--syndrome-rounds",
                str(schedule.merge_rounds),
                "--out",
                str(base_clean),
            ],
            cwd=root,
            env=env,
            check=True,
            stdout=subprocess.DEVNULL,
        )
    retimed = GENERATED_ROOT / f"{source}_{basis.lower()}_k{circuit_k}_pre{schedule.pre_rounds}_merge{schedule.merge_rounds}_post{schedule.post_rounds}.stim"
    if not retimed.exists():
        subprocess.run(
            [
                str(py),
                "python/tqec_import/retime_simple_merge.py",
                "--stim",
                str(base_clean),
                "--pre-rounds",
                str(schedule.pre_rounds),
                "--merge-rounds",
                str(schedule.merge_rounds),
                "--post-rounds",
                str(schedule.post_rounds),
                "--out",
                str(retimed),
            ],
            cwd=root,
            env=env,
            check=True,
            stdout=subprocess.DEVNULL,
        )
    return retimed


def validate_deterministic_circuit(path: Path) -> None:
    circuit = stim.Circuit.from_file(path)
    circuit.detector_error_model(allow_gauge_detectors=False)


def schedule_for_circuit_k(circuit_k: int, rounds_per_phase: int | None) -> RoundSchedule:
    if rounds_per_phase is not None:
        return RoundSchedule(
            pre_rounds=rounds_per_phase,
            merge_rounds=rounds_per_phase,
            post_rounds=0,
        )
    return default_round_schedule_for_k(circuit_k)


def main() -> None:
    args = parse_args()
    invalid = [p for p in args.protocols if p not in DISTILLATION_PROTOCOLS or p == "radar"]
    if invalid:
        raise SystemExit(
            "Unsupported protocols in --protocols: "
            + ", ".join(invalid)
            + ". Use values from DISTILLATION_PROTOCOLS excluding radar."
        )

    root = PROJECT_ROOT
    py = root / ".venv-tqec" / "bin" / "python"
    env = build_env(root)

    circuit_k_values = sorted(args.circuit_k_values)
    distillation_rounds_values = sorted(args.distillation_rounds_values)
    results = []

    for protocol in args.protocols:
        for m in args.m_values:
            for distillation_rounds in distillation_rounds_values:
                for enr in args.enr_values:
                    best_result = None
                    best_circuit_k = None
                    best_injected_stim = None
                    best_injected_json = None
                    best_timing_summary = None

                    for circuit_k in circuit_k_values:
                        schedule = schedule_for_circuit_k(circuit_k, args.rounds_per_phase)
                        clean = ensure_clean_circuit(py, env, root, args.source, args.basis, circuit_k, schedule)
                        validate_deterministic_circuit(clean)
                        stem = (
                            f"{args.source}_{args.basis.lower()}_k{circuit_k}_pre{schedule.pre_rounds}"
                            f"_merge{schedule.merge_rounds}_post{schedule.post_rounds}"
                            f".enr{int(enr)}.m{m}.{protocol}.dr{distillation_rounds}"
                        )
                        injected_stim = INJECTED_ROOT / f"{stem}.network_noisy.stim"
                        injected_json = INJECTED_ROOT / f"{stem}.network_noisy.json"

                        subprocess.run(
                            [
                                str(py),
                                "python/tqec_import/inject_noise.py",
                                "--stim",
                                str(clean),
                                "--physical-error",
                                str(args.physical_error),
                                "--accurate-rcx",
                                "--distillation-protocol",
                                protocol,
                                "--distillation-rounds",
                                str(distillation_rounds),
                                "--distillation-backup-batches",
                                str(m),
                                "--entanglement-rate",
                                str(enr),
                                "--measurement-delay",
                                str(args.measurement_delay),
                                "--out-stim",
                                str(injected_stim),
                                "--out-json",
                                str(injected_json),
                            ],
                            cwd=root,
                            env=env,
                            check=True,
                            stdout=subprocess.DEVNULL,
                        )

                        timing_summary = json.loads(injected_json.read_text())
                        timing_ok = all(section["timing_constraint_satisfied"] for section in timing_summary["sections"])
                        if timing_ok:
                            best_circuit_k = circuit_k
                            best_injected_stim = injected_stim
                            best_injected_json = injected_json
                            best_timing_summary = timing_summary
                        else:
                            break

                    combo_summary = {
                        "source": args.source,
                        "basis": args.basis,
                        "protocol": protocol,
                        "m": m,
                        "distillation_rounds": distillation_rounds,
                        "entanglement_rate": enr,
                        "measurement_delay": args.measurement_delay,
                        "candidate_circuit_k_values": circuit_k_values,
                    }

                    if best_circuit_k is None:
                        combo_summary["status"] = "no_feasible_circuit_k"
                        results.append(combo_summary)
                        continue

                    ler_json = LER_ROOT / f"{best_injected_stim.stem}.ler_{args.shots}_mpi.json"
                    proc = subprocess.run(
                        [
                            "mpirun",
                            "-np",
                            str(args.mpi_ranks),
                            str(py),
                            "python/tqec_import/run_ler.py",
                            "--stim",
                            str(best_injected_stim),
                            "--shots",
                            str(args.shots),
                            "--batch-shots",
                            str(args.batch_shots),
                            "--out",
                            str(ler_json),
                        ],
                        cwd=root,
                        env=env,
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                    ler_result = json.loads(proc.stdout)
                    schedule = schedule_for_circuit_k(best_circuit_k, args.rounds_per_phase)
                    combo_summary.update(
                        {
                            "status": "ok",
                            "circuit_k": best_circuit_k,
                            "distance_k": best_circuit_k,
                            "circuit_distance": circuit_distance_from_k(best_circuit_k),
                            "round_schedule": round_schedule_to_dict(schedule),
                            "timing_summary_path": str(best_injected_json),
                            "stim_path": str(best_injected_stim),
                            "timing_summary": {
                                "top_level_remote_pairs": best_timing_summary["top_level_remote_pairs"],
                                "max_idling_time_us": max(section["idling_time_us"] for section in best_timing_summary["sections"]),
                                "max_remote_cnot_error": max(section["remote_cnot_error"] for section in best_timing_summary["sections"]),
                                "timing_constraint_satisfied": all(section["timing_constraint_satisfied"] for section in best_timing_summary["sections"]),
                                "effective_distillation_protocols": [
                                    section["effective_distillation_protocol"] for section in best_timing_summary["sections"]
                                ],
                            },
                            "ler_result": ler_result,
                        }
                    )
                    results.append(combo_summary)

    summary = {
        "source": args.source,
        "basis": args.basis,
        "physical_error": args.physical_error,
        "measurement_delay": args.measurement_delay,
        "mpi_ranks": args.mpi_ranks,
        "batch_shots": args.batch_shots,
        "shots_per_point": args.shots,
        "naming": {
            "circuit_k_meaning": "TQEC spatial scale parameter with distance d = 2*circuit_k + 1",
            "distillation_rounds_meaning": "Distillation depth swept as an independent variable",
        },
        "sweep": {
            "entanglement_rate_values": args.enr_values,
            "m_values": args.m_values,
            "distillation_rounds_values": distillation_rounds_values,
            "protocols": args.protocols,
            "candidate_circuit_k_values": circuit_k_values,
            "rounds_per_phase": args.rounds_per_phase,
        },
        "results": results,
    }
    out = LER_ROOT / f"{args.name}.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
