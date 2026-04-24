#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from paths import PROJECT_ROOT, create_experiment_run_paths
from sweep_enr_protocol import build_env, schedule_for_circuit_k, validate_deterministic_circuit
from sweep_enr_deformation_variants import DeformationVariant, export_clean_variant, parse_variant


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


def parse_concurrency(text: str) -> int:
    value = text.strip().lower()
    if value == "auto":
        return 0
    parsed = int(value)
    if parsed <= 0:
        raise ValueError("--max-concurrent-jobs must be positive or 'auto'.")
    return parsed


def normalize_deformation_mode(mode: str) -> str:
    text = mode.strip().lower()
    if text in {"none", "nosuper"}:
        return "none"
    return text


@dataclass(frozen=True)
class SweepJob:
    variant: DeformationVariant
    clean_stim: Path
    entanglement_rate: float
    injected_stim: Path
    timing_json: Path
    result_json: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep entanglement rate for raw remote pairs only. "
            "The sweep applies optional border deformations, keeps no distillation, "
            "and records the remaining remote-CNOT count for each deformation."
        )
    )
    parser.add_argument("--source", choices=("simple_merge_only", "simple_merge_split"), default="simple_merge_only")
    parser.add_argument("--basis", choices=("X", "Z"), default="X")
    parser.add_argument("--circuit-k", type=int, default=2)
    parser.add_argument(
        "--rounds-per-phase",
        type=int,
        default=6,
        help="Set the pre-merge and merge repeat counts. Post rounds are always removed.",
    )
    parser.add_argument(
        "--deformation-modes",
        type=parse_csv_strings,
        default=["none", "lside", "rside", "twoside"],
        help="Border deformation settings to sweep. 'nosuper' is accepted as an alias for 'none'.",
    )
    parser.add_argument(
        "--enr-values",
        type=parse_csv_floats,
        default=[25e6, 50e6, 75e6, 100e6, 150e6, 200e6, 300e6, 500e6, 750e6, 1e9],
        help="Entanglement-rate values (Hz) to sweep.",
    )
    parser.add_argument("--physical-error", type=float, default=0.001)
    parser.add_argument("--measurement-delay", type=float, default=660e-9)
    parser.add_argument("--shots", type=int, default=100_000)
    parser.add_argument("--batch-shots", type=int, default=250_000)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument(
        "--max-concurrent-jobs",
        type=parse_concurrency,
        default=1,
        help="How many sweep points to run at once. Use 'auto' to derive it from local CPU count and --mpi-ranks.",
    )
    parser.add_argument("--name", type=str, default="enr_raw_pairs_sweep")
    return parser.parse_args()


def compute_required_epr_pairs_per_round(timing_summary: dict, rounds_per_phase: int) -> int:
    sections = timing_summary.get("sections", [])
    if sections:
        return max(section["required_epr_pairs_per_round"] for section in sections)
    total_remote_pairs = int(timing_summary["top_level_remote_pairs"])
    round_cycle_count = rounds_per_phase + 2
    if round_cycle_count <= 0:
        raise ValueError("rounds_per_phase must produce a positive cycle count.")
    if total_remote_pairs % round_cycle_count != 0:
        raise ValueError(
            f"Remote-pair count {total_remote_pairs} is not divisible by the round cycle count {round_cycle_count}."
        )
    return total_remote_pairs // round_cycle_count


def resolve_max_concurrent_jobs(requested: int, mpi_ranks: int) -> int:
    if requested > 0:
        return requested
    cpu_count = os.cpu_count() or 1
    return max(1, cpu_count // max(1, mpi_ranks))


def write_run_metadata(
    *,
    run_dir: Path,
    args: argparse.Namespace,
    normalized_modes: list[str],
    variants: list[DeformationVariant],
    schedule: dict[str, int],
    max_concurrent_jobs: int,
) -> None:
    params_dir = run_dir / "params"
    args_payload = vars(args).copy()
    args_payload["deformation_modes"] = normalized_modes
    args_payload["max_concurrent_jobs_resolved"] = max_concurrent_jobs
    (params_dir / "args.json").write_text(json.dumps(args_payload, indent=2) + "\n", encoding="utf-8")
    manifest = {
        "variants": [{"label": variant.label, "export_args": list(variant.export_args)} for variant in variants],
        "schedule": schedule,
        "naming": {
            "circuit_k_meaning": "TQEC spatial scale parameter with distance d = 2*circuit_k + 1",
            "deformation_baseline_aliases": ["none", "nosuper"],
        },
    }
    (params_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run_job(
    *,
    job: SweepJob,
    py: Path,
    root: Path,
    env: dict[str, str],
    physical_error: float,
    measurement_delay: float,
    shots: int,
    batch_shots: int,
    mpi_ranks: int,
    rounds_per_phase: int,
) -> dict:
    subprocess.run(
        [
            str(py),
            "python/tqec_import/inject_noise.py",
            "--stim",
            str(job.clean_stim),
            "--physical-error",
            str(physical_error),
            "--accurate-rcx",
            "--distillation-protocol",
            "none",
            "--distillation-rounds",
            "0",
            "--distillation-backup-batches",
            "1",
            "--entanglement-rate",
            str(job.entanglement_rate),
            "--measurement-delay",
            str(measurement_delay),
            "--out-stim",
            str(job.injected_stim),
            "--out-json",
            str(job.timing_json),
        ],
        cwd=root,
        env=env,
        check=True,
        stdout=subprocess.DEVNULL,
    )

    timing_summary = json.loads(job.timing_json.read_text())
    required_epr_pairs_per_round = compute_required_epr_pairs_per_round(timing_summary, rounds_per_phase)

    proc = subprocess.run(
        [
            "mpirun",
            "-np",
            str(mpi_ranks),
            str(py),
            "python/tqec_import/run_ler.py",
            "--stim",
            str(job.injected_stim),
            "--shots",
            str(shots),
            "--batch-shots",
            str(batch_shots),
            "--out",
            str(job.result_json),
        ],
        cwd=root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    ler_result = json.loads(proc.stdout)

    return {
        "deformation_mode": job.variant.label,
        "entanglement_rate": job.entanglement_rate,
        "timing_summary_path": str(job.timing_json),
        "stim_path": str(job.injected_stim),
        "ler_json_path": str(job.result_json),
        "required_epr_pairs_per_round": required_epr_pairs_per_round,
        "remaining_remote_cnots_per_section": [
            {
                "section_index": section["section_index"],
                "repeat_count": section["repeat_count"],
                "remote_cnots_per_body": section["remote_cnots_per_body"],
            }
            for section in timing_summary["sections"]
        ],
        "timing_constraint_satisfied": all(section["timing_constraint_satisfied"] for section in timing_summary["sections"])
        if timing_summary["sections"]
        else True,
        "ler_result": ler_result,
    }


def main() -> None:
    args = parse_args()
    root = PROJECT_ROOT
    py = root / ".venv-tqec" / "bin" / "python"
    env = build_env(root)
    max_concurrent_jobs = resolve_max_concurrent_jobs(args.max_concurrent_jobs, args.mpi_ranks)

    normalized_modes = [normalize_deformation_mode(mode) for mode in args.deformation_modes]
    if len(set(normalized_modes)) != len(normalized_modes):
        raise SystemExit("Duplicate deformation modes after normalization. Remove repeated values such as both 'none' and 'nosuper'.")

    schedule = schedule_for_circuit_k(args.circuit_k, args.rounds_per_phase)
    if schedule.post_rounds != 0 or schedule.pre_rounds != args.rounds_per_phase or schedule.merge_rounds != args.rounds_per_phase:
        raise SystemExit("This sweep currently expects a fixed pre=merge=rounds_per_phase and post=0 schedule.")

    run_paths = create_experiment_run_paths("raw_pairs_sweeps", args.name)
    write_run_metadata(
        run_dir=run_paths.root,
        args=args,
        normalized_modes=normalized_modes,
        variants=[parse_variant(mode) for mode in normalized_modes],
        schedule={
            "pre_rounds": schedule.pre_rounds,
            "merge_rounds": schedule.merge_rounds,
            "post_rounds": schedule.post_rounds,
        },
        max_concurrent_jobs=max_concurrent_jobs,
    )

    variants = [parse_variant(mode) for mode in normalized_modes]
    results_by_mode: dict[str, dict] = {}
    jobs: list[SweepJob] = []

    for variant in variants:
        clean = export_clean_variant(
            py=py,
            root=root,
            env=env,
            source=args.source,
            basis=args.basis,
            circuit_k=args.circuit_k,
            rounds_per_phase=args.rounds_per_phase,
            variant=variant,
        )
        validate_deterministic_circuit(clean)

        run_clean = run_paths.clean_circuits_dir / clean.name
        shutil.copyfile(clean, run_clean)

        results_by_mode[variant.label] = {
            "deformation_mode": variant.label,
            "export_args": list(variant.export_args),
            "clean_stim_path": str(run_clean),
            "entanglement_rate_values": args.enr_values,
            "results": [],
        }

        for enr in args.enr_values:
            stem = f"{run_clean.stem}.enr{int(enr)}.rawpairs"
            jobs.append(
                SweepJob(
                    variant=variant,
                    clean_stim=run_clean,
                    entanglement_rate=enr,
                    injected_stim=run_paths.noisy_circuits_dir / f"{stem}.network_noisy.stim",
                    timing_json=run_paths.timing_dir / f"{stem}.network_noisy.json",
                    result_json=run_paths.results_dir / f"{stem}.ler_{args.shots}_mpi.json",
                )
            )

    with ThreadPoolExecutor(max_workers=max_concurrent_jobs) as executor:
        completed = executor.map(
            lambda job: run_job(
                job=job,
                py=py,
                root=root,
                env=env,
                physical_error=args.physical_error,
                measurement_delay=args.measurement_delay,
                shots=args.shots,
                batch_shots=args.batch_shots,
                mpi_ranks=args.mpi_ranks,
                rounds_per_phase=args.rounds_per_phase,
            ),
            jobs,
        )
        for entry in completed:
            results_by_mode[entry["deformation_mode"]]["results"].append(entry)

    results = [results_by_mode[variant.label] for variant in variants]
    summary = {
        "run_dir": str(run_paths.root),
        "source": args.source,
        "basis": args.basis,
        "physical_error": args.physical_error,
        "measurement_delay": args.measurement_delay,
        "mpi_ranks": args.mpi_ranks,
        "max_concurrent_jobs": max_concurrent_jobs,
        "batch_shots": args.batch_shots,
        "shots_per_point": args.shots,
        "circuit_k": args.circuit_k,
        "distance": 2 * args.circuit_k + 1,
        "rounds_per_phase": args.rounds_per_phase,
        "deformation_modes": normalized_modes,
        "results": results,
    }
    out = run_paths.root / "summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
