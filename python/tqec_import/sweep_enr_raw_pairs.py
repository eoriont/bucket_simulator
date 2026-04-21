#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from paths import INJECTED_ROOT, LER_ROOT, PROJECT_ROOT
from sweep_enr_protocol import build_env, ensure_clean_circuit, schedule_for_circuit_k, validate_deterministic_circuit


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
        help="Border deformation settings to sweep.",
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
    parser.add_argument("--name", type=str, default="enr_raw_pairs_sweep")
    return parser.parse_args()


def export_clean_circuit(
    py: Path,
    root: Path,
    env: dict[str, str],
    source: str,
    basis: str,
    circuit_k: int,
    rounds_per_phase: int,
    deformation_mode: str,
) -> Path:
    base_clean = Path("output") / "tqec_import" / "generated" / f"{source}_{basis.lower()}_k{circuit_k}_r{rounds_per_phase}.stim"
    if deformation_mode != "none":
        base_clean = base_clean.with_name(base_clean.stem + f".{deformation_mode}" + base_clean.suffix)
    base_clean = root / base_clean
    if not base_clean.exists():
        cmd = [
            str(py),
            "python/tqec_import/export_stim.py",
            "--source",
            source,
            "--basis",
            basis,
            "--k",
            str(circuit_k),
            "--syndrome-rounds",
            str(rounds_per_phase),
            "--out",
            str(base_clean),
        ]
        if deformation_mode != "none":
            cmd += ["--deformation-mode", deformation_mode]
        subprocess.run(cmd, cwd=root, env=env, check=True, stdout=subprocess.DEVNULL)
    if deformation_mode != "none":
        return base_clean
    retimed = base_clean.with_name(base_clean.stem.replace(f"_r{rounds_per_phase}", f"_pre{rounds_per_phase}_merge{rounds_per_phase}_post0") + base_clean.suffix)
    if not retimed.exists():
        subprocess.run(
            [
                str(py),
                "python/tqec_import/retime_simple_merge.py",
                "--stim",
                str(base_clean),
                "--pre-rounds",
                str(rounds_per_phase),
                "--merge-rounds",
                str(rounds_per_phase),
                "--post-rounds",
                "0",
                "--out",
                str(retimed),
            ],
            cwd=root,
            env=env,
            check=True,
            stdout=subprocess.DEVNULL,
        )
    return retimed


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


def main() -> None:
    args = parse_args()
    root = PROJECT_ROOT
    py = root / ".venv-tqec" / "bin" / "python"
    env = build_env(root)

    results = []
    for deformation_mode in args.deformation_modes:
        clean = export_clean_circuit(
            py=py,
            root=root,
            env=env,
            source=args.source,
            basis=args.basis,
            circuit_k=args.circuit_k,
            rounds_per_phase=args.rounds_per_phase,
            deformation_mode=deformation_mode,
        )
        validate_deterministic_circuit(clean)

        injected_json_for_mode = None
        mode_result = {
            "deformation_mode": deformation_mode,
            "clean_stim_path": str(clean),
            "entanglement_rate_values": args.enr_values,
            "results": [],
        }

        for enr in args.enr_values:
            stem = f"{clean.stem}.enr{int(enr)}.rawpairs"
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
                    "none",
                    "--distillation-rounds",
                    "0",
                    "--distillation-backup-batches",
                    "1",
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
            required_epr_pairs_per_round = compute_required_epr_pairs_per_round(timing_summary, args.rounds_per_phase)

            ler_json = LER_ROOT / f"{stem}.ler_{args.shots}_mpi.json"
            proc = subprocess.run(
                [
                    "mpirun",
                    "-np",
                    str(args.mpi_ranks),
                    str(py),
                    "python/tqec_import/run_ler.py",
                    "--stim",
                    str(injected_stim),
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
            entry = {
                "entanglement_rate": enr,
                "timing_summary_path": str(injected_json),
                "stim_path": str(injected_stim),
                "required_epr_pairs_per_round": required_epr_pairs_per_round,
                "remaining_remote_cnots_per_section": [
                    {
                        "section_index": section["section_index"],
                        "repeat_count": section["repeat_count"],
                        "remote_cnots_per_body": section["remote_cnots_per_body"],
                    }
                    for section in timing_summary["sections"]
                ],
                "timing_constraint_satisfied": all(
                    section["timing_constraint_satisfied"] for section in timing_summary["sections"]
                )
                if timing_summary["sections"]
                else True,
                "ler_result": ler_result,
            }
            mode_result["results"].append(entry)
            injected_json_for_mode = injected_json

        mode_result["timing_summary_path"] = str(injected_json_for_mode) if injected_json_for_mode else None
        results.append(mode_result)

    summary = {
        "source": args.source,
        "basis": args.basis,
        "physical_error": args.physical_error,
        "measurement_delay": args.measurement_delay,
        "mpi_ranks": args.mpi_ranks,
        "batch_shots": args.batch_shots,
        "shots_per_point": args.shots,
        "circuit_k": args.circuit_k,
        "rounds_per_phase": args.rounds_per_phase,
        "deformation_modes": args.deformation_modes,
        "results": results,
    }
    out = LER_ROOT / f"{args.name}.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
