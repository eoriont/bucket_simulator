#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import stim

from stim_annotations import annotate_with_polygons
from stim_layout import deformation_coordinates_for_mode, seam_border_coordinate
from stim_deformation import apply_code_deformation, parse_coordinate_list


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply simple code-deformation edits to a TQEC-exported Stim circuit by "
            "removing selected data qubits and suppressing adjacent repeated ancilla checks."
        )
    )
    parser.add_argument("--stim", type=Path, required=True, help="Input .stim circuit.")
    parser.add_argument("--out", type=Path, required=True, help="Output .stim circuit.")
    parser.add_argument(
        "--remove-data",
        type=str,
        default=None,
        help='Whitespace-separated coordinate list such as "(3,9) (5,9)".',
    )
    parser.add_argument(
        "--remove-seam-border",
        choices=("min", "max"),
        default=None,
        help="Infer the seam orientation and remove one data qubit from the corresponding seam endpoint.",
    )
    parser.add_argument(
        "--mode",
        choices=("lside", "rside", "twoside"),
        default=None,
        help="Infer the seam orientation and remove the left, right, or both seam-edge data qubits.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    circuit = stim.Circuit.from_file(args.stim)
    specified = sum(
        option is not None
        for option in (args.remove_data, args.remove_seam_border, args.mode)
    )
    if specified != 1:
        raise SystemExit("Specify exactly one of --remove-data, --remove-seam-border, or --mode.")
    if args.remove_data is not None:
        removed = parse_coordinate_list(args.remove_data)
    elif args.remove_seam_border is not None:
        removed = [seam_border_coordinate(circuit, args.remove_seam_border)]
    else:
        removed = deformation_coordinates_for_mode(circuit, args.mode)
    result = apply_code_deformation(circuit, removed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(annotate_with_polygons(result.circuit), encoding="utf-8")
    print(
        f"Wrote {args.out} "
        f"(removed_data={len(result.removed_data_qubits)}, suppressed_ancillas={len(result.suppressed_ancillas)}, "
        f"qubits={result.circuit.num_qubits}, detectors={result.circuit.num_detectors})"
    )


if __name__ == "__main__":
    main()
