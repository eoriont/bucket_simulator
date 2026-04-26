#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from paths import GENERATED_ROOT, INJECTED_ROOT, LER_ROOT, PROJECT_ROOT
from sweep_enr_protocol import build_env, ensure_clean_circuit, schedule_for_circuit_k, validate_ler_ready_circuit


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


def parse_semicolon_strings(text: str) -> list[str]:
    values = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        values.append(part)
    if not values:
        raise ValueError("Expected at least one string value.")
    return values


@dataclass(frozen=True)
class DeformationVariant:
    label: str
    export_args: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sweep entanglement rate for a mixed set of deformation variants. "
            "Variants may be baseline ('none'), inferred border modes "
            "('lside', 'rside', 'twoside'), or explicit interior superstabilizer "
            "removals such as 'ss@5,11' or 'middle_ss=5,11'."
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
        "--variants",
        type=parse_semicolon_strings,
        default=["none", "lside", "rside", "twoside", "ss_low=5,9", "ss_mid=5,11", "ss_high=5,13"],
        help="Semicolon-separated deformation variants.",
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
    parser.add_argument("--name", type=str, default="enr_deformation_variants_sweep")
    return parser.parse_args()


def parse_variant(spec: str) -> DeformationVariant:
    text = spec.strip()
    if text == "none":
        return DeformationVariant(label="none", export_args=())
    if text in {"lside", "rside", "twoside"}:
        return DeformationVariant(label=text, export_args=("--mode", text))
    if text == "ss_mid":
        return DeformationVariant(label="ss_mid", export_args=("--remove-data", "(5,11)"))

    match = re.fullmatch(r"(?:(?P<label>[A-Za-z0-9_-]+)=)?ss@?(?P<x>-?\d+(?:\.\d+)?),(?P<y>-?\d+(?:\.\d+)?)", text)
    if match:
        label = match.group("label") or f"ss_{match.group('x').replace('.', 'p')}_{match.group('y').replace('.', 'p')}"
        coord = f"({match.group('x')},{match.group('y')})"
        return DeformationVariant(label=label, export_args=("--remove-data", coord))

    match = re.fullmatch(r"(?P<label>[A-Za-z0-9_-]+)=(?P<x>-?\d+(?:\.\d+)?),(?P<y>-?\d+(?:\.\d+)?)", text)
    if match:
        coord = f"({match.group('x')},{match.group('y')})"
        return DeformationVariant(label=match.group("label"), export_args=("--remove-data", coord))

    raise ValueError(f"Unrecognized deformation variant: {spec}")


def export_clean_variant(
    py: Path,
    root: Path,
    env: dict[str, str],
    source: str,
    basis: str,
    circuit_k: int,
    rounds_per_phase: int,
    variant: DeformationVariant,
) -> Path:
    base_name = f"{source}_{basis.lower()}_k{circuit_k}_pre{rounds_per_phase}_merge{rounds_per_phase}_post0"
    out_path = GENERATED_ROOT / f"{base_name}.{variant.label}.stim"
    if out_path.exists():
        return out_path

    schedule = schedule_for_circuit_k(circuit_k, rounds_per_phase)
    base_clean = ensure_clean_circuit(py, env, root, source, basis, circuit_k, schedule)
    if variant.label == "none":
        out_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(base_clean, out_path)
        return out_path

    cmd = [str(py), "python/tqec_import/deform_stim.py", "--stim", str(base_clean), "--out", str(out_path), *variant.export_args]
    subprocess.run(cmd, cwd=root, env=env, check=True, stdout=subprocess.DEVNULL)
    return out_path


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

    schedule = schedule_for_circuit_k(args.circuit_k, args.rounds_per_phase)
    if schedule.post_rounds != 0 or schedule.pre_rounds != args.rounds_per_phase or schedule.merge_rounds != args.rounds_per_phase:
        raise SystemExit("This sweep currently expects a fixed pre=merge=rounds_per_phase and post=0 schedule.")

    variants = [parse_variant(spec) for spec in args.variants]
    results = []

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
        validate_ler_ready_circuit(clean)

        variant_result = {
            "variant": variant.label,
            "export_args": list(variant.export_args),
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
            variant_result["results"].append(
                {
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
            )

        results.append(variant_result)

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
        "variants": args.variants,
        "results": results,
    }
    out = LER_ROOT / f"{args.name}.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
