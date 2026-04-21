#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import stim


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Retimes the three top-level repeat sections of a TQEC-exported simple-merge "
            "Stim circuit. This is intended for circuits whose top-level structure is "
            "pre-merge repeat, merge repeat, post-merge repeat."
        )
    )
    parser.add_argument("--stim", type=Path, required=True, help="Input .stim circuit.")
    parser.add_argument("--out", type=Path, required=True, help="Output .stim circuit.")
    parser.add_argument("--pre-rounds", type=int, required=True, help="Repeat count for the first top-level REPEAT block.")
    parser.add_argument("--merge-rounds", type=int, required=True, help="Repeat count for the second top-level REPEAT block.")
    parser.add_argument("--post-rounds", type=int, required=True, help="Repeat count for the third top-level REPEAT block.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for name, value in (
        ("pre-rounds", args.pre_rounds),
        ("merge-rounds", args.merge_rounds),
        ("post-rounds", args.post_rounds),
    ):
        if value < 0:
            raise SystemExit(f"--{name} must be non-negative.")

    circuit = stim.Circuit.from_file(args.stim)
    replacement_counts = [args.pre_rounds, args.merge_rounds, args.post_rounds]
    seen_repeats = 0
    retimed = stim.Circuit()

    for op in circuit:
        if op.name != "REPEAT":
            retimed.append(op)
            continue

        if seen_repeats >= 3:
            raise SystemExit(
                "Expected exactly three top-level REPEAT blocks in the simple-merge circuit; "
                "found more than three."
            )

        repeat_count = replacement_counts[seen_repeats]
        if repeat_count > 0:
            retimed.append(stim.CircuitRepeatBlock(repeat_count, op.body_copy()))
        seen_repeats += 1

    if seen_repeats != 3:
        raise SystemExit(
            "Expected exactly three top-level REPEAT blocks in the simple-merge circuit; "
            f"found {seen_repeats}."
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(str(retimed), encoding="utf-8")


if __name__ == "__main__":
    main()
