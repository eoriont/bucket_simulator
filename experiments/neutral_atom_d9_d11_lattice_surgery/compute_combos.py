#!/usr/bin/env python3
"""
Print all fully-feasible (protocol:rounds:m_star) combos for a given ENR.

Imports the RADar math from rcx_regime_sweep/sweep_optimal_rcx.py so the
formulas stay in one place. Output is space-separated tokens suitable for
direct use in run.sh:

    combos_str=$(python3 compute_combos.py --config config_distributed.txt --enr 300000)
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path


def load_sweep_module():
    sweep_path = (Path(__file__).resolve().parent.parent
                  / "rcx_regime_sweep" / "sweep_optimal_rcx.py")
    if not sweep_path.exists():
        print(f"error: sweep module not found at {sweep_path}", file=sys.stderr)
        sys.exit(1)
    spec = importlib.util.spec_from_file_location("sweep_optimal_rcx", sweep_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["sweep_optimal_rcx"] = mod  # needed for dataclass type resolution
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    parser = argparse.ArgumentParser(
        description="List feasible (protocol:rounds:batches) combos for a given ENR.")
    parser.add_argument("--config", type=Path, required=True,
                        help="Distributed simulator config file")
    parser.add_argument("--enr", type=int, required=True,
                        help="Entanglement rate in Hz")
    parser.add_argument("--max-rounds", type=int, default=5)
    args = parser.parse_args()

    sweep = load_sweep_module()

    file_cfg = sweep.parse_config_file(args.config)
    d = sweep.get_int(file_cfg, "code_distance", 9)
    cfg = sweep.Config(
        code_distance=d,
        raw_epr_fidelity=sweep.get_float(file_cfg, "raw_epr_fidelity", 0.95),
        physical_error=sweep.get_float(file_cfg, "physical_error", 0.005),
        measurement_delay=sweep.get_float(file_cfg, "measurement_delay", 500e-6),
        T1=sweep.get_float(file_cfg, "T1", 119.0),
        T2=sweep.get_float(file_cfg, "T2", 12.6),
        num_remote_cnots=2 * d - 1,
    )

    combos = []
    for protocol, rounds in sweep.protocol_candidates(args.max_rounds):
        result = sweep.evaluate(protocol, rounds, args.enr, cfg)
        if result.fully_feasible:
            # Use rounds=1/batches=1 for "none" so the simulator gets valid args
            r = max(rounds, 1)
            combos.append(f"{protocol}:{r}:{result.m_star}")

    print(" ".join(combos))
    return 0


if __name__ == "__main__":
    sys.exit(main())
