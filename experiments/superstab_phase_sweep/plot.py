#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt


MODE_LABELS = {
    "xx_merge": "XX Merge",
    "zz_merge": "ZZ Merge",
    "xx_split": "XX Split",
    "zz_split": "ZZ Split",
}

SUPER_LABELS = {
    "border": "Border",
    "border2": "Border2",
    "middle": "Middle",
    "nosuper": "No superstabilizer",
    "twoside": "Two-sided",
}

SUPER_ORDER = ["nosuper", "border", "border2", "middle", "twoside"]
MODE_ORDER = ["xx_merge", "zz_merge", "xx_split", "zz_split"]


def parse_result_file(path: Path) -> dict | None:
    text = path.read_text()
    rate_match = re.search(r"Entanglement Rate:\s+([\d.]+)\s+MHz", text)
    ler_match = re.search(r"Logical Error Rate:\s+([\d.eE+\-]+)", text)
    name_match = re.match(r"(?P<mode>[a-z]+_[a-z]+)_(?P<super>[a-z0-9]+)_(?P<rate>\d+)MHz_d9_results", path.stem)
    if not rate_match or not ler_match or not name_match:
        print(f"  Warning: could not parse {path.name}")
        return None
    return {
        "mode": name_match.group("mode"),
        "super": name_match.group("super"),
        "rate_mhz": float(rate_match.group(1)),
        "ler": float(ler_match.group(1)),
    }


def plot(run_dir: Path):
    results_dir = run_dir / "results"
    if not results_dir.is_dir():
        sys.exit(f"Error: results directory not found: {results_dir}")

    rows: dict[str, dict[str, list[tuple[float, float]]]] = {}
    for path in sorted(results_dir.glob("*_results.txt")):
        rec = parse_result_file(path)
        if rec is None:
            continue
        rows.setdefault(rec["mode"], {}).setdefault(rec["super"], []).append((rec["rate_mhz"], rec["ler"]))

    if not rows:
        sys.exit("Error: no parseable result files")

    cmap = plt.get_cmap("tab10")
    saved_paths: list[Path] = []

    for mode in MODE_ORDER:
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        mode_data = rows.get(mode, {})
        for idx, super_type in enumerate(SUPER_ORDER):
            pts = sorted(mode_data.get(super_type, []))
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            ax.plot(xs, ys, marker="o", linewidth=1.8, markersize=4, color=cmap(idx),
                    label=SUPER_LABELS.get(super_type, super_type))

        ax.set_title(f"{MODE_LABELS.get(mode, mode)}: ENR vs LER")
        ax.set_yscale("log")
        ax.grid(True, alpha=0.35)
        ax.set_xlabel("ENR (MHz)")
        ax.set_ylabel("Logical Error Rate")
        ax.legend(framealpha=0.9)
        fig.tight_layout()

        out_path = run_dir / f"{mode}_plot.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        saved_paths.append(out_path)

    for path in saved_paths:
        print(f"Saved: {path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python plot.py <run_dir>")
    plot(Path(sys.argv[1]).expanduser().resolve())
