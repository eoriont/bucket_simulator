#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import pymatching
import stim

from paths import LER_ROOT

try:
    from mpi4py import MPI
except ImportError:
    MPI = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a batched logical-error-rate experiment on a Stim circuit."
    )
    parser.add_argument("--stim", type=Path, required=True, help="Path to a .stim circuit.")
    parser.add_argument(
        "--shots",
        type=int,
        default=100_000,
        help="Number of Monte Carlo shots.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Optional JSON output path for the summary.",
    )
    parser.add_argument(
        "--batch-shots",
        type=int,
        default=1_000_000,
        help="Maximum number of shots per sampling batch on each rank.",
    )
    return parser.parse_args()


def get_mpi_context() -> tuple[object | None, int, int]:
    launcher_size = (
        os.environ.get("OMPI_COMM_WORLD_SIZE")
        or os.environ.get("PMI_SIZE")
        or os.environ.get("MPI_LOCALNRANKS")
    )
    if MPI is None:
        if launcher_size is not None and int(launcher_size) > 1:
            raise SystemExit(
                "MPI launcher detected but mpi4py is not installed in this environment. "
                "Install mpi4py or run without mpirun."
            )
        return None, 0, 1
    comm = MPI.COMM_WORLD
    return comm, comm.Get_rank(), comm.Get_size()


def split_work(total: int, rank: int, size: int) -> int:
    base = total // size
    remainder = total % size
    return base + (1 if rank < remainder else 0)


def run_local_batches(
    circuit: stim.Circuit,
    matcher: pymatching.Matching,
    local_shots: int,
    batch_shots: int,
) -> tuple[int, int]:
    sampler = circuit.compile_detector_sampler()
    failures = 0
    batches = 0
    remaining = local_shots
    while remaining > 0:
        current_batch = min(batch_shots, remaining)
        detection_events, observable_flips = sampler.sample(
            shots=current_batch,
            separate_observables=True,
        )
        predicted_flips = matcher.decode_batch(detection_events)
        logical_failures = (predicted_flips != observable_flips).any(axis=1)
        failures += int(logical_failures.sum())
        remaining -= current_batch
        batches += 1
    return failures, batches


def main() -> None:
    args = parse_args()
    if args.shots <= 0:
        raise SystemExit("--shots must be positive.")
    if args.batch_shots <= 0:
        raise SystemExit("--batch-shots must be positive.")

    comm, rank, size = get_mpi_context()
    circuit = stim.Circuit.from_file(args.stim)

    t0 = time.time()
    dem = circuit.detector_error_model(
        decompose_errors=True,
        allow_gauge_detectors=False,
    )
    matcher = pymatching.Matching.from_detector_error_model(dem)
    local_shots = split_work(args.shots, rank, size)
    num_failures, local_batches = run_local_batches(
        circuit=circuit,
        matcher=matcher,
        local_shots=local_shots,
        batch_shots=args.batch_shots,
    )
    local_runtime_s = time.time() - t0

    if comm is not None:
        total_failures = comm.reduce(num_failures, op=MPI.SUM, root=0)
        total_batches = comm.reduce(local_batches, op=MPI.SUM, root=0)
        max_runtime_s = comm.reduce(local_runtime_s, op=MPI.MAX, root=0)
    else:
        total_failures = num_failures
        total_batches = local_batches
        max_runtime_s = local_runtime_s

    if rank != 0:
        return

    result = {
        "stim_path": str(args.stim),
        "shots": args.shots,
        "batch_shots": args.batch_shots,
        "qubits": circuit.num_qubits,
        "detectors": circuit.num_detectors,
        "observables": circuit.num_observables,
        "logical_failures": int(total_failures),
        "logical_error_rate": total_failures / args.shots,
        "runtime_seconds": max_runtime_s,
        "detector_validation": "passed",
        "mpi": {
            "enabled": size > 1,
            "ranks": size,
            "mpi4py_available": MPI is not None,
            "total_batches": int(total_batches),
        },
    }

    print(json.dumps(result, indent=2))

    out_path = args.out
    if out_path is None:
        out_path = LER_ROOT / f"{args.stim.stem}_ler.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
