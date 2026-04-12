#!/usr/bin/env python3
"""
Plot LER results from neutral_atom_d13_lattice_surgery runs.

Two figures:
  1. ler_vs_enr_by_protocol.png  — LER vs ENR, one line per (protocol, k, m),
                                   faceted by mode (xx_merge / zz_merge / etc.)
  2. distributed_vs_monolithic.png — best distributed LER vs ENR compared
                                     against the flat monolithic baseline,
                                     one panel per mode.

Usage:
    python3 plot_results.py runs/<timestamp>
    python3 plot_results.py runs/<timestamp> --modes xx_merge zz_merge
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

LER_RE = re.compile(r"Logical Error Rate:\s*([\d.eE+\-]+)")

PROTOCOL_COLORS = {
    "none":             "#888888",
    "2to1_pumping":     "#4C72B0",
    "3to1_pumping":     "#55A868",
    "2to1_recurrence":  "#C44E52",
    "3to1_recurrence":  "#DD8452",
}

PROTOCOL_LABELS = {
    "none":             "none",
    "2to1_pumping":     "2→1 pump",
    "3to1_pumping":     "3→1 pump",
    "2to1_recurrence":  "2→1 rec",
    "3to1_recurrence":  "3→1 rec",
}

MODES = ["xx_merge", "zz_merge", "xx_split", "zz_split"]


def parse_ler(path: Path) -> float | None:
    try:
        text = path.read_text()
    except OSError:
        return None
    m = LER_RE.search(text)
    if not m:
        return None
    v = float(m.group(1))
    return v if v > 0 else None


def khz(hz: int) -> float:
    return hz / 1e3


def parse_filename(name: str) -> dict | None:
    """
    Parse label from result filename stem, e.g.:
      monolithic_xx_merge_d13_results
      distributed_xx_merge_none_50kHz_d13_results
      distributed_xx_merge_2to1_pumping_k1_m6_650kHz_d13_results
    """
    stem = name.replace("_results.txt", "").replace("_results", "")

    # monolithic
    m = re.fullmatch(r"monolithic_(\w+_(?:merge|split))_d(\d+)", stem)
    if m:
        return dict(arch="monolithic", mode=m.group(1), distance=int(m.group(2)),
                    protocol="monolithic", rounds=0, batches=0, enr_hz=0)

    # distributed none
    m = re.fullmatch(r"distributed_(\w+_(?:merge|split))_none_(\d+)kHz_d(\d+)", stem)
    if m:
        return dict(arch="distributed", mode=m.group(1), distance=int(m.group(3)),
                    protocol="none", rounds=1, batches=1,
                    enr_hz=int(m.group(2)) * 1000)

    # distributed with distillation
    m = re.fullmatch(
        r"distributed_(\w+_(?:merge|split))_(\w+)_k(\d+)_m(\d+)_(\d+)kHz_d(\d+)",
        stem,
    )
    if m:
        return dict(arch="distributed", mode=m.group(1), distance=int(m.group(6)),
                    protocol=m.group(2), rounds=int(m.group(3)), batches=int(m.group(4)),
                    enr_hz=int(m.group(5)) * 1000)

    return None


def load_results(run_dir: Path) -> list[dict]:
    records = []
    results_dir = run_dir / "results"
    if not results_dir.exists():
        print(f"No results/ directory in {run_dir}", file=sys.stderr)
        return records
    for path in sorted(results_dir.glob("*.txt")):
        info = parse_filename(path.name)
        if info is None:
            continue
        ler = parse_ler(path)
        if ler is None:
            continue
        info["ler"] = ler
        info["path"] = path
        records.append(info)
    return records


# ---------------------------------------------------------------------------
# Plot 1: LER vs ENR by protocol, faceted by mode
# ---------------------------------------------------------------------------

def plot_ler_by_protocol(records: list[dict], out_path: Path, modes: list[str]) -> None:
    dist = [r for r in records if r["arch"] == "distributed"]
    if not dist:
        print("No distributed results to plot.", file=sys.stderr)
        return

    ncols = min(2, len(modes))
    nrows = (len(modes) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.5 * nrows),
                             constrained_layout=True, squeeze=False)

    for ax_idx, mode in enumerate(modes):
        ax = axes[ax_idx // ncols][ax_idx % ncols]
        mode_data = [r for r in dist if r["mode"] == mode]

        # Group by (protocol, rounds, batches)
        groups: dict[tuple, list[dict]] = {}
        for r in mode_data:
            key = (r["protocol"], r["rounds"], r["batches"])
            groups.setdefault(key, []).append(r)

        plotted_labels: set[str] = set()
        for (protocol, rounds, batches), pts in sorted(groups.items()):
            pts.sort(key=lambda r: r["enr_hz"])
            xs = [khz(r["enr_hz"]) for r in pts]
            ys = [r["ler"] for r in pts]
            color = PROTOCOL_COLORS.get(protocol, "#333333")

            if protocol == "none":
                label = "none (raw EPR)"
                lw, ls, alpha = 1.8, "--", 0.9
            else:
                label = f"{PROTOCOL_LABELS.get(protocol, protocol)} k={rounds} m={batches}"
                lw, ls, alpha = 1.2, "-", 0.75

            show_label = label not in plotted_labels
            ax.semilogy(xs, ys, color=color, lw=lw, ls=ls, alpha=alpha,
                        marker="o", ms=3, label=label if show_label else "_nolegend_")
            plotted_labels.add(label)

        ax.set_title(mode.replace("_", " "), fontsize=11)
        ax.set_xlabel("ENR (kHz)")
        ax.set_ylabel("LER")
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}"))
        ax.grid(True, which="both", ls=":", alpha=0.4)

        # Legend: deduplicate by protocol color (one entry per protocol family)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, fontsize=7, ncol=1, loc="best",
                  framealpha=0.7, handlelength=1.5)

    # Hide unused axes
    for i in range(len(modes), nrows * ncols):
        axes[i // ncols][i % ncols].set_visible(False)

    fig.suptitle("LER vs ENR by distillation protocol  (d=13)", fontsize=13)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


# ---------------------------------------------------------------------------
# Plot 2: Best distributed vs monolithic
# ---------------------------------------------------------------------------

def plot_distributed_vs_monolithic(records: list[dict], out_path: Path,
                                   modes: list[str]) -> None:
    mono = {r["mode"]: r["ler"] for r in records if r["arch"] == "monolithic"}
    dist = [r for r in records if r["arch"] == "distributed"]
    if not dist:
        print("No distributed results to plot.", file=sys.stderr)
        return

    ncols = min(2, len(modes))
    nrows = (len(modes) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.5 * nrows),
                             constrained_layout=True, squeeze=False)

    for ax_idx, mode in enumerate(modes):
        ax = axes[ax_idx // ncols][ax_idx % ncols]
        mode_dist = [r for r in dist if r["mode"] == mode]

        # All distributed lines, faint
        groups: dict[tuple, list[dict]] = {}
        for r in mode_dist:
            key = (r["protocol"], r["rounds"], r["batches"])
            groups.setdefault(key, []).append(r)

        for (protocol, rounds, batches), pts in groups.items():
            pts.sort(key=lambda r: r["enr_hz"])
            xs = [khz(r["enr_hz"]) for r in pts]
            ys = [r["ler"] for r in pts]
            color = PROTOCOL_COLORS.get(protocol, "#333333")
            ax.semilogy(xs, ys, color=color, lw=0.8, alpha=0.35,
                        marker=".", ms=3, ls="-")

        # Best distributed per ENR (min LER)
        enr_best: dict[int, float] = {}
        for r in mode_dist:
            enr_best[r["enr_hz"]] = min(enr_best.get(r["enr_hz"], 1.0), r["ler"])

        if enr_best:
            xs_best = sorted(enr_best)
            ys_best = [enr_best[e] for e in xs_best]
            ax.semilogy([khz(e) for e in xs_best], ys_best,
                        color="#1a1a1a", lw=2.2, marker="D", ms=5, zorder=5,
                        label="best distributed")

        # Monolithic baseline
        if mode in mono:
            ax.axhline(mono[mode], color="firebrick", lw=1.8, ls="--",
                       label=f"monolithic  (LER={mono[mode]:.2e})")

        ax.set_title(mode.replace("_", " "), fontsize=11)
        ax.set_xlabel("ENR (kHz)")
        ax.set_ylabel("LER")
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}"))
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(fontsize=8, loc="best", framealpha=0.7)

    for i in range(len(modes), nrows * ncols):
        axes[i // ncols][i % ncols].set_visible(False)

    fig.suptitle("Distributed (best per ENR) vs monolithic baseline  (d=13)", fontsize=13)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Plot d=13 lattice surgery LER results.")
    parser.add_argument("run_dir", type=Path, help="Path to a runs/<timestamp> directory")
    parser.add_argument("--modes", nargs="+", default=MODES,
                        help="Modes to plot (default: all four)")
    args = parser.parse_args()

    records = load_results(args.run_dir)
    if not records:
        print("No results found.", file=sys.stderr)
        return 1

    print(f"Loaded {len(records)} results from {args.run_dir}")
    mono_count = sum(1 for r in records if r["arch"] == "monolithic")
    dist_count = sum(1 for r in records if r["arch"] == "distributed")
    print(f"  {mono_count} monolithic,  {dist_count} distributed")

    modes = [m for m in args.modes if m in MODES]

    plot_ler_by_protocol(records, args.run_dir / "ler_vs_enr_by_protocol.png", modes)
    plot_distributed_vs_monolithic(records, args.run_dir / "distributed_vs_monolithic.png", modes)

    return 0


if __name__ == "__main__":
    sys.exit(main())
