#!/usr/bin/env python3
"""
Plot LER results from neutral_atom_d7_lattice_surgery runs.

One PNG per mode (xx_merge, zz_merge, xx_split, zz_split).
Each plot shows:
  - All distributed lines, colored by distillation protocol, faint
  - Best distributed line per ENR highlighted in bold with the protocol color
  - Monolithic baseline as a dashed red line
  - Clear protocol labels in the legend

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
import matplotlib.lines as mlines
import numpy as np

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

LER_RE   = re.compile(r"Logical Error Rate:\s*([\d.eE+\-]+)")
SHOTS_RE = re.compile(r"Total Shots:\s*(\d+)")

PROTOCOL_COLORS = {
    "none":             "#888888",
    "2to1_pumping":     "#4C72B0",
    "3to1_pumping":     "#55A868",
    "2to1_recurrence":  "#C44E52",
    "3to1_recurrence":  "#DD8452",
}

PROTOCOL_LABELS = {
    "none":             "None (raw EPR)",
    "2to1_pumping":     "2→1 pumping",
    "3to1_pumping":     "3→1 pumping",
    "2to1_recurrence":  "2→1 recurrence",
    "3to1_recurrence":  "3→1 recurrence",
}

MODES = ["xx_merge", "zz_merge", "xx_split", "zz_split"]


def parse_ler(path: Path) -> tuple[float | None, bool, int]:
    """Returns (ler, is_upper_bound, total_shots).
    When LER=0, returns 1/shots as an upper bound with is_upper_bound=True."""
    try:
        text = path.read_text()
    except OSError:
        return None, False, 0
    m = LER_RE.search(text)
    if not m:
        return None, False, 0
    v = float(m.group(1))
    ms = SHOTS_RE.search(text)
    shots = int(ms.group(1)) if ms else 0
    if v == 0.0:
        ub = 1.0 / shots if shots > 0 else None
        return ub, True, shots
    return v, False, shots


def khz(hz: int) -> float:
    return hz / 1e3


def parse_filename(name: str) -> dict | None:
    stem = name.replace("_results.txt", "").replace("_results", "")

    m = re.fullmatch(r"monolithic_(\w+_(?:merge|split))_d(\d+)", stem)
    if m:
        return dict(arch="monolithic", mode=m.group(1), distance=int(m.group(2)),
                    protocol="monolithic", rounds=0, batches=0, enr_hz=0)

    m = re.fullmatch(r"distributed_(\w+_(?:merge|split))_none_(\d+)kHz_d(\d+)", stem)
    if m:
        return dict(arch="distributed", mode=m.group(1), distance=int(m.group(3)),
                    protocol="none", rounds=1, batches=1,
                    enr_hz=int(m.group(2)) * 1000)

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
        ler, is_ub, shots = parse_ler(path)
        if ler is None:
            continue
        info["ler"] = ler
        info["is_upper_bound"] = is_ub
        info["shots"] = shots
        info["path"] = path
        records.append(info)
    return records


# ---------------------------------------------------------------------------
# Per-mode plot
# ---------------------------------------------------------------------------

def plot_mode(records: list[dict], mode: str, run_dir: Path) -> None:
    dist = [r for r in records if r["arch"] == "distributed" and r["mode"] == mode]
    mono_records = [r for r in records if r["arch"] == "monolithic" and r["mode"] == mode]
    mono = (mono_records[0]["ler"], mono_records[0]["is_upper_bound"]) if mono_records else None

    d = records[0]["distance"] if records else "?"
    title = f"{mode.replace('_', ' ').title()}  (d={d})"

    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)

    # --- Group by protocol ---
    groups: dict[str, list[dict]] = {}
    for r in dist:
        groups.setdefault(r["protocol"], []).append(r)

    # --- Find best distributed per ENR (min LER) ---
    enr_best: dict[int, tuple[float, str]] = {}  # enr_hz -> (ler, protocol)
    for r in dist:
        enr = r["enr_hz"]
        if enr not in enr_best or r["ler"] < enr_best[enr][0]:
            enr_best[enr] = (r["ler"], r["protocol"])

    # --- Plot all lines faint, grouped by protocol ---
    # Track which protocols actually appear so we can build legend
    seen_protocols: set[str] = set()
    for protocol, pts in sorted(groups.items()):
        color = PROTOCOL_COLORS.get(protocol, "#333333")
        # Group further by (rounds, batches) for individual lines
        sub: dict[tuple, list[dict]] = {}
        for r in pts:
            sub.setdefault((r["rounds"], r["batches"]), []).append(r)
        for (rounds, batches), rpts in sorted(sub.items()):
            rpts.sort(key=lambda r: r["enr_hz"])
            xs = [khz(r["enr_hz"]) for r in rpts]
            ys = [r["ler"] for r in rpts]
            lbl = f"{PROTOCOL_LABELS.get(protocol, protocol)} k={rounds} m={batches}" \
                  if protocol != "none" else PROTOCOL_LABELS["none"]
            ax.semilogy(xs, ys, color=color, lw=0.9, alpha=0.30,
                        marker=".", ms=4, ls="-", zorder=2)
            seen_protocols.add(protocol)

    # --- Highlight best per ENR with thick colored line ---
    # Build segments by protocol so the highlight is colored correctly
    best_by_protocol: dict[str, list[tuple[float, float]]] = {}
    for enr, (ler, prot) in sorted(enr_best.items()):
        best_by_protocol.setdefault(prot, []).append((khz(enr), ler))

    for prot, pts in best_by_protocol.items():
        pts.sort()
        color = PROTOCOL_COLORS.get(prot, "#333333")
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.semilogy(xs, ys, color=color, lw=2.8, alpha=1.0,
                    marker="D", ms=6, ls="-", zorder=5,
                    label=f"★ best: {PROTOCOL_LABELS.get(prot, prot)}")

    # Connect best segments across protocol transitions with a thin black spine
    all_best = sorted(enr_best.items())
    spine_xs = [khz(e) for e, _ in all_best]
    spine_ys = [ler for _, (ler, _) in all_best]
    ax.semilogy(spine_xs, spine_ys, color="#1a1a1a", lw=1.0, alpha=0.5,
                ls="--", zorder=4)

    # --- Monolithic baseline ---
    if mono:
        ler, is_ub = mono
        if is_ub:
            ax.axhline(ler, color="firebrick", lw=1.5, ls=":",
                       label=f"Monolithic  (<{ler:.1e}, 0 errors)")
            xlim = ax.get_xlim()
            xmid = (xlim[0] + xlim[1]) / 2 if xlim[1] > xlim[0] else 500
            ax.annotate("", xy=(xmid, ler * 0.4), xytext=(xmid, ler),
                        arrowprops=dict(arrowstyle="-|>", color="firebrick", lw=1.5))
        else:
            ax.axhline(ler, color="firebrick", lw=2.0, ls="--",
                       label=f"Monolithic  (LER = {ler:.2e})")

    # --- Legend: protocol color patches + best-line entries ---
    legend_handles = []

    # Color patches for each distillation protocol
    for protocol in sorted(seen_protocols, key=lambda p: list(PROTOCOL_COLORS).index(p)
                           if p in PROTOCOL_COLORS else 99):
        color = PROTOCOL_COLORS.get(protocol, "#333333")
        patch = mlines.Line2D([], [], color=color, lw=1.5, alpha=0.6,
                              label=PROTOCOL_LABELS.get(protocol, protocol))
        legend_handles.append(patch)

    # Best-line entries (auto from ax)
    handles, labels = ax.get_legend_handles_labels()
    for h, l in zip(handles, labels):
        legend_handles.append(h)

    ax.legend(handles=legend_handles, fontsize=8, loc="best",
              framealpha=0.8, handlelength=2.0)

    ax.set_title(title, fontsize=13)
    ax.set_xlabel("ENR (kHz)", fontsize=11)
    ax.set_ylabel("Logical Error Rate", fontsize=11)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}"))
    ax.grid(True, which="both", ls=":", alpha=0.4)

    out_path = run_dir / f"{mode}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Plot lattice surgery LER results.")
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
    dist_count  = sum(1 for r in records if r["arch"] == "distributed")
    print(f"  {mono_count} monolithic,  {dist_count} distributed")

    modes = [m for m in args.modes if m in MODES]
    for mode in modes:
        plot_mode(records, mode, args.run_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
