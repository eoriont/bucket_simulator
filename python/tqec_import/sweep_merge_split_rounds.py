#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

from paths import GENERATED_ROOT, INJECTED_ROOT, LER_ROOT, PROJECT_ROOT
from schedule_config import circuit_distance_from_k, default_round_schedule_for_k, round_schedule_to_dict


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run merge/split round sweeps by retiming the pre-merge, merge, and "
            "post-split top-level repeat blocks in a TQEC simple merge/split circuit."
        )
    )
    parser.add_argument(
        "--source",
        choices=("simple_merge_split", "simple_merge_only", "simple_split_only"),
        default="simple_merge_split",
        help="Which exported TQEC experiment to retime and sweep.",
    )
    parser.add_argument("--basis", choices=("X", "Z"), default="X")
    parser.add_argument(
        "--k",
        dest="circuit_k",
        type=int,
        default=2,
        help="TQEC circuit-size parameter. Physical distance is d = 2*circuit_k + 1.",
    )
    parser.add_argument("--pre-rounds", type=parse_csv_ints)
    parser.add_argument("--merge-rounds", type=parse_csv_ints)
    parser.add_argument("--post-rounds", type=parse_csv_ints)
    parser.add_argument("--shots", type=int, default=10_000_000)
    parser.add_argument("--batch-shots", type=int, default=250_000)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--physical-error", type=float, default=0.001)
    parser.add_argument("--monolithic-baseline", action="store_true")
    parser.add_argument("--name", type=str, default="merge_split_round_sweep")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = PROJECT_ROOT
    py = root / ".venv-tqec" / "bin" / "python"

    default_schedule = default_round_schedule_for_k(args.circuit_k)
    pre_rounds = args.pre_rounds if args.pre_rounds is not None else [default_schedule.pre_rounds]
    merge_rounds = args.merge_rounds if args.merge_rounds is not None else [default_schedule.merge_rounds]
    post_rounds = args.post_rounds if args.post_rounds is not None else [default_schedule.post_rounds]

    for p in [root / ".tqec", root / ".mplconfig", root / ".cache", root / ".fontconfig"]:
        p.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["XDG_DATA_HOME"] = str(root / ".tqec")
    env["MPLCONFIGDIR"] = str(root / ".mplconfig")
    env["XDG_CACHE_HOME"] = str(root / ".cache")
    env["FONTCONFIG_PATH"] = "/opt/homebrew/etc/fonts"
    env["FONTCONFIG_FILE"] = "/opt/homebrew/etc/fonts/fonts.conf"

    base_round = max(max(pre_rounds), max(merge_rounds), max(post_rounds))
    source_tag = args.source.removeprefix("simple_")
    base_clean = GENERATED_ROOT / f"{args.source}_{args.basis.lower()}_k{args.circuit_k}_r{base_round}.stim"
    subprocess.run(
        [
            str(py),
            "python/tqec_import/export_stim.py",
            "--source",
            args.source,
            "--basis",
            args.basis,
            "--k",
            str(args.circuit_k),
            "--syndrome-rounds",
            str(base_round),
            "--out",
            str(base_clean),
        ],
        cwd=root,
        env=env,
        check=True,
        stdout=subprocess.DEVNULL,
    )

    results = []
    for pre in pre_rounds:
        for merge in merge_rounds:
            for post in post_rounds:
                stem = f"{args.source}_{args.basis.lower()}_k{args.circuit_k}_pre{pre}_merge{merge}_post{post}"
                retimed = GENERATED_ROOT / f"{stem}.stim"
                injected = INJECTED_ROOT / f"{stem}.monolithic.stim"
                injected_json = INJECTED_ROOT / f"{stem}.monolithic.json"
                ler_json = LER_ROOT / f"{stem}.monolithic_{args.shots}_mpi.json"

                subprocess.run(
                    [
                        str(py),
                        "python/tqec_import/retime_simple_merge.py",
                        "--stim",
                        str(base_clean),
                        "--pre-rounds",
                        str(pre),
                        "--merge-rounds",
                        str(merge),
                        "--post-rounds",
                        str(post),
                        "--out",
                        str(retimed),
                    ],
                    cwd=root,
                    env=env,
                    check=True,
                    stdout=subprocess.DEVNULL,
                )

                inject_cmd = [
                    str(py),
                    "python/tqec_import/inject_noise.py",
                    "--stim",
                    str(retimed),
                    "--physical-error",
                    str(args.physical_error),
                    "--out-stim",
                    str(injected),
                    "--out-json",
                    str(injected_json),
                ]
                if args.monolithic_baseline:
                    inject_cmd.append("--monolithic-baseline")
                subprocess.run(inject_cmd, cwd=root, env=env, check=True, stdout=subprocess.DEVNULL)

                proc = subprocess.run(
                    [
                        "mpirun",
                        "-np",
                        str(args.mpi_ranks),
                        str(py),
                        "python/tqec_import/run_ler.py",
                        "--stim",
                        str(injected),
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
                data = json.loads(proc.stdout)
                data["pre_rounds"] = pre
                data["merge_rounds"] = merge
                data["post_rounds"] = post
                results.append(data)

    summary = {
        "source": args.source,
        "circuit_k": args.circuit_k,
        "distance_k": args.circuit_k,
        "circuit_distance": circuit_distance_from_k(args.circuit_k),
        "approx_distance": circuit_distance_from_k(args.circuit_k),
        "basis": args.basis,
        "physical_error": args.physical_error,
        "mode": "monolithic_baseline_with_measurement_delay" if args.monolithic_baseline else "network_or_local_noise",
        "mpi_ranks": args.mpi_ranks,
        "batch_shots": args.batch_shots,
        "shots_per_point": args.shots,
        "default_round_schedule": round_schedule_to_dict(default_schedule),
        "sweep": {
            "pre_rounds": pre_rounds,
            "merge_rounds": merge_rounds,
            "post_rounds": post_rounds,
        },
        "naming": {
            "circuit_k_meaning": "TQEC spatial scale parameter with distance d = 2*circuit_k + 1",
            "distillation_k_meaning": "Distillation depth; in code/config this is stored as distillation_rounds",
        },
        "results": results,
    }
    out = LER_ROOT / f"{args.name}.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
