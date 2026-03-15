#!/usr/bin/env python3
from __future__ import annotations
"""
Plot superstabilizer sweep results.

Usage:
    python plot.py <run_dir>
    python plot.py experiments/superstab_sweep/runs/20260315_120000

Reads all *_results.txt files in <run_dir>/results/, parses LER and
entanglement rate, and plots LER vs rate for each superstabilizer type.
Output is saved as <run_dir>/plot.png (and shown interactively if a
display is available).
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np


# ── Label/color config ────────────────────────────────────────────────────────

SUPER_LABELS = {
    "border":  "Border (y=0.5)",
    "border2": "Border2 (y=4.5)",
    "middle":  "Middle (y=2.5)",
    "nosuper": "No superstabilizer",
    "twoside": "Two-sided",
}

SUPER_ORDER = ["nosuper", "border", "border2", "middle", "twoside"]


# ── Parsing ───────────────────────────────────────────────────────────────────

def parse_result_file(path: Path) -> dict | None:
    """Return a dict with 'super', 'rate_mhz', 'ler', 'ler_num', 'ler_den'
    parsed from a result file, or None if the file can't be parsed."""
    text = path.read_text()

    # Entanglement rate — stored in MHz in the output
    rate_match = re.search(r"Entanglement Rate:\s+([\d.]+)\s+MHz", text)
    # Logical Error Rate line: "  Logical Error Rate: 0.009000 (9/1000)"
    ler_match = re.search(r"Logical Error Rate:\s+([\d.eE+\-]+)(?:\s+\((\d+)/(\d+)\))?", text)

    if not rate_match or not ler_match:
        print(f"  Warning: could not parse {path.name}")
        return None

    rate_mhz = float(rate_match.group(1))
    ler = float(ler_match.group(1))
    ler_num = int(ler_match.group(2)) if ler_match.group(2) else None
    ler_den = int(ler_match.group(3)) if ler_match.group(3) else None

    # Infer superstabilizer type from filename: <type>_<rate>MHz_results.txt
    name = path.stem  # e.g. "border_5MHz_results"
    # strip trailing _results
    name = re.sub(r"_results$", "", name)
    # strip trailing _<N>MHz
    super_type = re.sub(r"_\d+MHz$", "", name)

    return {
        "super": super_type,
        "rate_mhz": rate_mhz,
        "ler": ler,
        "ler_num": ler_num,
        "ler_den": ler_den,
    }


def wilson_ci(successes, total, z=1.96):
    """Wilson score confidence interval for a proportion."""
    if total == 0:
        return 0.0, 0.0
    p = successes / total
    denom = 1 + z**2 / total
    centre = (p + z**2 / (2 * total)) / denom
    half = z * np.sqrt(p * (1 - p) / total + z**2 / (4 * total**2)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot(run_dir: Path):
    results_dir = run_dir / "results"
    if not results_dir.is_dir():
        sys.exit(f"Error: results directory not found: {results_dir}")

    result_files = sorted(results_dir.glob("*_results.txt"))
    if not result_files:
        sys.exit(f"Error: no *_results.txt files in {results_dir}")

    print(f"Parsing {len(result_files)} result files from {results_dir}...")

    # Collect data keyed by (super_type, rate_mhz)
    data: dict[str, dict[float, dict]] = {}
    for path in result_files:
        rec = parse_result_file(path)
        if rec is None:
            continue
        st = rec["super"]
        if st not in data:
            data[st] = {}
        data[st][rec["rate_mhz"]] = rec

    if not data:
        sys.exit("Error: no data could be parsed")

    # Order the super types
    ordered = [s for s in SUPER_ORDER if s in data]
    ordered += [s for s in data if s not in ordered]

    fig, ax = plt.subplots(figsize=(9, 6))

    cmap = plt.get_cmap("tab10")
    for i, super_type in enumerate(ordered):
        rates_data = data[super_type]
        rates = sorted(rates_data.keys())
        lers = [rates_data[r]["ler"] for r in rates]

        # Error bars from Wilson CI if shot counts available
        yerr_lo, yerr_hi = [], []
        for r in rates:
            rec = rates_data[r]
            if rec["ler_num"] is not None and rec["ler_den"] is not None:
                lo, hi = wilson_ci(rec["ler_num"], rec["ler_den"])
                yerr_lo.append(rec["ler"] - lo)
                yerr_hi.append(hi - rec["ler"])
            else:
                yerr_lo.append(0)
                yerr_hi.append(0)

        label = SUPER_LABELS.get(super_type, super_type)
        color = cmap(i)

        has_errors = any(lo > 0 or hi > 0 for lo, hi in zip(yerr_lo, yerr_hi))
        if has_errors:
            ax.errorbar(
                rates, lers,
                yerr=[yerr_lo, yerr_hi],
                label=label, color=color,
                marker="o", markersize=5,
                capsize=3, capthick=1.2,
                linewidth=1.8,
            )
        else:
            ax.plot(rates, lers, label=label, color=color,
                    marker="o", markersize=5, linewidth=1.8)

    ax.set_xlabel("Entanglement Rate (MHz)", fontsize=12)
    ax.set_ylabel("Logical Error Rate", fontsize=12)
    ax.set_title("Superstabilizer Comparison — LER vs Entanglement Rate", fontsize=13)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.35)
    ax.set_yscale("log")

    # x-axis in MHz units
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.0f}"))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5e6))

    # Annotate run ID
    run_id = run_dir.name
    ax.text(0.99, 0.01, f"run: {run_id}", transform=ax.transAxes,
            fontsize=7, color="gray", ha="right", va="bottom")

    fig.tight_layout()

    out_path = run_dir / "plot.png"
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")

    try:
        plt.show()
    except Exception:
        pass  # headless environment


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    plot(Path(sys.argv[1]).expanduser().resolve())
