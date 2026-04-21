#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import stim

from stim_annotations import annotate_with_polygons


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add Crumble/Stim polygon pragmas for stabilizers and seam-crossing interactions."
    )
    parser.add_argument("--stim", type=Path, required=True, help="Input .stim circuit.")
    parser.add_argument("--out", type=Path, required=True, help="Output annotated .stim circuit.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    circuit = stim.Circuit.from_file(args.stim)
    annotated = annotate_with_polygons(circuit)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(annotated, encoding="utf-8")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
