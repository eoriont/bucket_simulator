#!/usr/bin/env python3
"""
Plot tqec round sweep results, optionally overlaid with C++ bucket_simulator results.

Usage:
  python plot.py <tqec_run_dir> [--cpp <cpp_run_dir> ...]
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt


LINESTYLES = ["-", "--", ":", "-."]
COLORS = {
    "tqec_d9":  "#1f77b4",
    "tqec_d11": "#2ca02c",
    "cpp_d9":   "#ff7f0e",
    "cpp_d11":  "#d62728",
}
MARKERS = {"tqec": "o", "cpp": "s"}


def load_tqec(run_dir: Path) -> list[dict]:
    summary = run_dir / "summary.json"
    if not summary.exists():
        sys.exit(f"Error: no summary.json in {run_dir}")
    data = json.loads(summary.read_text())
    rows = []
    for r in data["results"]:
        ler = r["logical_error_rate"]
        shots = r["shots"]
        rows.append({
            "source": "tqec",
            "distance": r["distance"],
            "merge_rounds": r["merge_rounds"],
            "ler": ler,
            "plot_ler": ler if ler > 0 else 0.5 / shots,
        })
    return rows


def load_cpp(run_dir: Path) -> list[dict]:
    results_dir = run_dir / "results"
    if not results_dir.is_dir():
        print(f"  Warning: no results/ in {run_dir}")
        return []
    rows = []
    for path in sorted(results_dir.glob("*_results.txt")):
        text = path.read_text()
        shots_m = re.search(r"Total Shots:\s+(\d+)", text)
        ler_m = re.search(r"Logical Error Rate:\s+([\d.eE+\-]+)", text)
        name_m = re.match(
            r"(?P<mode>[a-z]+_[a-z]+)_(?P<super>[a-z0-9]+)_r(?P<rounds>\d+)_(?P<rate>[0-9A-Za-z]+)_d(?P<dist>\d+)_results",
            path.stem,
        )
        if not shots_m or not ler_m or not name_m:
            continue
        if name_m.group("super") != "nosuper":
            continue
        shots = int(shots_m.group(1))
        ler = float(ler_m.group(1))
        rows.append({
            "source": "cpp",
            "distance": int(name_m.group("dist")),
            "merge_rounds": int(name_m.group("rounds")),
            "ler": ler,
            "plot_ler": ler if ler > 0 else 0.5 / shots,
        })
    return rows


def plot(tqec_dir: Path, cpp_dirs: list[Path]) -> None:
    tqec_rows = load_tqec(tqec_dir)
    cpp_rows: list[dict] = []
    for d in cpp_dirs:
        cpp_rows.extend(load_cpp(d))

    all_rows = tqec_rows + cpp_rows
    if not all_rows:
        sys.exit("No data to plot.")

    distances = sorted(set(r["distance"] for r in all_rows))
    fig, ax = plt.subplots(figsize=(9, 5.8))

    for di, dist in enumerate(distances):
        ls = LINESTYLES[di % len(LINESTYLES)]

        for src in ("tqec", "cpp"):
            pts = sorted(
                (r["merge_rounds"], r["plot_ler"])
                for r in all_rows
                if r["source"] == src and r["distance"] == dist
            )
            if not pts:
                continue
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            key = f"{src}_d{dist}"
            color = COLORS.get(key, "#999999")
            label = f"{'TQEC' if src == 'tqec' else 'C++'} (d={dist})"
            ax.plot(xs, ys,
                    marker=MARKERS[src], linestyle=ls, linewidth=1.8,
                    markersize=5, color=color, label=label)

    ax.set_title("XX merge: merge rounds vs LER\n(TQEC vs C++ bucket_simulator, 500 MHz ENR, no superstabilizer)")
    ax.set_xlabel("Merge rounds")
    ax.set_ylabel("Logical Error Rate")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.35)
    ax.legend(framealpha=0.9, fontsize=8)
    ax.text(0.02, 0.02,
            "TQEC: pre=1, post=1. C++: pre=1, post=0.",
            transform=ax.transAxes, fontsize=7.5, alpha=0.8)
    fig.tight_layout()

    plot_path = tqec_dir / "round_sweep_plot.png"
    fig.savefig(plot_path, dpi=160)
    plt.close(fig)
    print(f"Saved: {plot_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("tqec_run_dir", type=Path)
    parser.add_argument("--cpp", type=Path, nargs="*", default=[])
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    plot(args.tqec_run_dir.expanduser().resolve(),
         [p.expanduser().resolve() for p in args.cpp])
