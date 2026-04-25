#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import stim

from distributed_noise import DISTILLATION_PROTOCOLS, NoiseConfig, inject_network_noise
from paths import INJECTED_ROOT


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Inject distributed-network RCX and idling noise into a TQEC-generated Stim circuit."
    )
    parser.add_argument("--stim", type=Path, required=True, help="Input clean .stim circuit.")
    parser.add_argument("--out-stim", type=Path, help="Output noisy .stim path.")
    parser.add_argument("--out-json", type=Path, help="Output JSON summary path.")
    parser.add_argument("--physical-error", type=float, default=0.001)
    parser.add_argument("--interconnect-error", type=float, default=0.0)
    parser.add_argument("--accurate-rcx", action="store_true")
    parser.add_argument("--channel-depolarization-error", type=float)
    parser.add_argument("--raw-epr-fidelity", type=float, default=0.99)
    parser.add_argument("--distillation-protocol", choices=DISTILLATION_PROTOCOLS, default="none")
    parser.add_argument("--distillation-rounds", type=int, default=1)
    parser.add_argument("--distillation-backup-batches", type=int, default=1)
    parser.add_argument("--entanglement-rate", type=float, default=100e6)
    parser.add_argument("--measurement-delay", type=float, default=660e-9)
    parser.add_argument("--t1", type=float, default=125e-6)
    parser.add_argument("--t2", type=float, default=200e-6)
    parser.add_argument("--monolithic-baseline", action="store_true")
    parser.add_argument("--local-baseline", action="store_true")
    parser.add_argument("--split-axis", choices=("auto", "x", "y"), default="auto")
    parser.add_argument("--split-x", type=float)
    return parser.parse_args()


def default_output_paths(stim_path: Path) -> tuple[Path, Path]:
    stem = stim_path.stem
    out_stim = INJECTED_ROOT / f"{stem}.network_noisy.stim"
    out_json = INJECTED_ROOT / f"{stem}.network_noisy.json"
    return out_stim, out_json


def main() -> None:
    args = parse_args()
    circuit = stim.Circuit.from_file(args.stim)
    config = NoiseConfig(
        physical_error=args.physical_error,
        interconnect_error=args.interconnect_error,
        accurate_rcx=args.accurate_rcx,
        channel_depolarization_error=args.channel_depolarization_error,
        raw_epr_fidelity=args.raw_epr_fidelity,
        distillation_protocol=args.distillation_protocol,
        distillation_rounds=args.distillation_rounds,
        distillation_backup_batches=args.distillation_backup_batches,
        entanglement_rate=args.entanglement_rate,
        measurement_delay=args.measurement_delay,
        t1_coherence_time=args.t1,
        t2_coherence_time=args.t2,
        monolithic_baseline=args.monolithic_baseline,
        local_baseline=args.local_baseline,
        split_axis=args.split_axis,
        split_x=args.split_x,
    )

    noisy_circuit, summary = inject_network_noise(circuit, config)
    out_stim, out_json = default_output_paths(args.stim)
    if args.out_stim is not None:
        out_stim = args.out_stim
    if args.out_json is not None:
        out_json = args.out_json

    out_stim.parent.mkdir(parents=True, exist_ok=True)
    out_stim.write_text(str(noisy_circuit), encoding="utf-8")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"stim_path": str(out_stim), "json_path": str(out_json), **summary}, indent=2))


if __name__ == "__main__":
    main()
