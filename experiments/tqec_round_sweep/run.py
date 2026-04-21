#!/usr/bin/env python3
"""
TQEC XX merge round sweep matching the 500 MHz bucket_simulator experiment.

simple_merge_only (X basis), k=4 (d=9) and k=5 (d=11),
pre=1, merge=1..13, post=1, ENR=500 MHz, p=0.001, raw_epr_fidelity=0.99, 300K shots.
"""
from __future__ import annotations

import json
import os
import subprocess
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PY = PROJECT_ROOT / ".venv-tqec" / "bin" / "python"

ENR = 500_000_000
RAW_EPR_FIDELITY = 0.99
PHYSICAL_ERROR = 0.001
SHOTS = 300_000
MPI_RANKS = 2
BATCH_SHOTS = 150_000

K_VALUES = [4, 5]   # d = 2k+1 → k=4 → d=9, k=5 → d=11
MERGE_ROUNDS = list(range(1, 14))
PRE_ROUNDS = 1
POST_ROUNDS = 1
BASIS = "X"
SOURCE = "simple_merge_only"


def make_env() -> dict:
    env = os.environ.copy()
    env["XDG_DATA_HOME"] = str(PROJECT_ROOT / ".tqec")
    env["MPLCONFIGDIR"] = str(PROJECT_ROOT / ".mplconfig")
    env["XDG_CACHE_HOME"] = str(PROJECT_ROOT / ".cache")
    env["FONTCONFIG_PATH"] = "/opt/homebrew/etc/fonts"
    env["FONTCONFIG_FILE"] = "/opt/homebrew/etc/fonts/fonts.conf"
    return env


def run(cmd, **kwargs):
    return subprocess.run(cmd, check=True, **kwargs)


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = Path(__file__).parent / "runs" / timestamp
    circuits_dir = run_dir / "circuits"
    injected_dir = run_dir / "injected"
    results_dir = run_dir / "results"
    for d in [circuits_dir, injected_dir, results_dir]:
        d.mkdir(parents=True, exist_ok=True)

    env = make_env()
    all_results: list[dict] = []

    for k in K_VALUES:
        dist = 2 * k + 1
        base_max = max(MERGE_ROUNDS)
        base_stim = circuits_dir / f"base_{SOURCE}_x_k{k}_r{base_max}.stim"

        print(f"\n[k={k} d={dist}] Exporting base circuit (max_rounds={base_max})...", flush=True)
        run(
            [str(PY), "python/tqec_import/export_stim.py",
             "--source", SOURCE, "--basis", BASIS,
             "--k", str(k), "--syndrome-rounds", str(base_max),
             "--out", str(base_stim)],
            cwd=PROJECT_ROOT, env=env,
        )

        for merge_r in MERGE_ROUNDS:
            label = f"{SOURCE}_x_k{k}_pre{PRE_ROUNDS}_merge{merge_r}_post{POST_ROUNDS}"
            retimed_stim = circuits_dir / f"{label}.stim"
            injected_stim = injected_dir / f"{label}.stim"
            injected_json = injected_dir / f"{label}.json"
            result_json = results_dir / f"{label}_ler.json"

            if result_json.exists():
                print(f"  [merge_r={merge_r:2d}] SKIP", flush=True)
                data = json.loads(result_json.read_text())
                data.update({"k": k, "distance": dist, "merge_rounds": merge_r})
                all_results.append(data)
                continue

            run(
                [str(PY), "python/tqec_import/retime_simple_merge.py",
                 "--stim", str(base_stim),
                 "--pre-rounds", str(PRE_ROUNDS),
                 "--merge-rounds", str(merge_r),
                 "--post-rounds", str(POST_ROUNDS),
                 "--out", str(retimed_stim)],
                cwd=PROJECT_ROOT, env=env,
                stdout=subprocess.DEVNULL,
            )

            run(
                [str(PY), "python/tqec_import/inject_noise.py",
                 "--stim", str(retimed_stim),
                 "--physical-error", str(PHYSICAL_ERROR),
                 "--entanglement-rate", str(ENR),
                 "--raw-epr-fidelity", str(RAW_EPR_FIDELITY),
                 "--distillation-protocol", "none",
                 "--out-stim", str(injected_stim),
                 "--out-json", str(injected_json)],
                cwd=PROJECT_ROOT, env=env,
                stdout=subprocess.DEVNULL,
            )

            print(f"  [merge_r={merge_r:2d}] Running LER...", flush=True, end="")
            proc = run(
                ["mpirun", "-np", str(MPI_RANKS),
                 str(PY), "python/tqec_import/run_ler.py",
                 "--stim", str(injected_stim),
                 "--shots", str(SHOTS),
                 "--batch-shots", str(BATCH_SHOTS),
                 "--out", str(result_json)],
                cwd=PROJECT_ROOT, env=env,
                capture_output=True, text=True,
            )
            data = json.loads(proc.stdout)
            data.update({"k": k, "distance": dist, "merge_rounds": merge_r})
            all_results.append(data)
            print(f" LER={data['logical_error_rate']:.4f}", flush=True)

    summary = {
        "timestamp": timestamp,
        "source": SOURCE,
        "basis": BASIS,
        "physical_error": PHYSICAL_ERROR,
        "entanglement_rate_hz": ENR,
        "raw_epr_fidelity": RAW_EPR_FIDELITY,
        "pre_rounds": PRE_ROUNDS,
        "post_rounds": POST_ROUNDS,
        "shots_per_point": SHOTS,
        "mpi_ranks": MPI_RANKS,
        "results": all_results,
    }
    summary_path = run_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\nDone. Summary: {summary_path}")
    print(f"To plot: python experiments/tqec_round_sweep/plot.py {run_dir}")


if __name__ == "__main__":
    main()
