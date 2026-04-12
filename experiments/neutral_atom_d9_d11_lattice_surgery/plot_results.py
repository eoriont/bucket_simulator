#!/usr/bin/env python3
"""
Plot LER results from neutral_atom_d9_d11_lattice_surgery runs.

One PNG per mode (xx_merge, zz_merge), showing d=9 and d=11 on the same axes.
  - All distributed lines faint, colored by distillation protocol
  - Best distributed per ENR highlighted in bold (per distance, distinguished by marker)
  - Monolithic baseline per distance as dashed lines

Usage:
    python3 plot_results.py runs/<timestamp>
    python3 plot_results.py runs/<timestamp> --modes xx_merge
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.lines as mlines

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

# Marker and linestyle per distance
DISTANCE_STYLE = {
    9:  {"marker": "o", "ls": "-",  "mono_ls": "--"},
    11: {"marker": "s", "ls": "--", "mono_ls": ":"},
}

MODES = ["xx_merge", "zz_merge"]

LER_RE   = re.compile(r"Logical Error Rate:\s*([\d.eE+\-]+)")
SHOTS_RE = re.compile(r"Total Shots:\s*(\d+)")


def parse_ler(path: Path) -> tuple[float | None, bool, int]:
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
        r"distributed_(\w+_(?:merge|split))_(\w+)_k(\d+)_m(\d+)_(\d+)kHz_d(\d+)", stem)
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
        records.append(info)
    return records


def plot_mode(records: list[dict], mode: str, run_dir: Path) -> None:
    distances = sorted(set(r["distance"] for r in records))
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)

    seen_protocols: set[str] = set()
    best_legend_entries: list = []

    for d in distances:
        style = DISTANCE_STYLE.get(d, {"marker": "^", "ls": "-", "mono_ls": "--"})

        # --- Monolithic baseline ---
        mono_rec = [r for r in records if r["arch"] == "monolithic"
                    and r["mode"] == mode and r["distance"] == d]
        if mono_rec:
            ler, is_ub = mono_rec[0]["ler"], mono_rec[0]["is_upper_bound"]
            mono_label = f"Monolithic d={d}  (LER={'<' if is_ub else ''}{ler:.2e})"
            color = "firebrick" if d == distances[0] else "darkred"
            ax.axhline(ler, color=color, lw=1.8, ls=style["mono_ls"], label=mono_label, zorder=3)

        # --- Distributed lines ---
        dist = [r for r in records if r["arch"] == "distributed"
                and r["mode"] == mode and r["distance"] == d]

        # Group by (protocol, rounds, batches) → individual lines
        sub: dict[tuple, list[dict]] = defaultdict(list)
        for r in dist:
            sub[(r["protocol"], r["rounds"], r["batches"])].append(r)

        for (protocol, rounds, batches), pts in sorted(sub.items()):
            pts.sort(key=lambda r: r["enr_hz"])
            xs = [khz(r["enr_hz"]) for r in pts]
            ys = [r["ler"] for r in pts]
            color = PROTOCOL_COLORS.get(protocol, "#333333")
            seen_protocols.add(protocol)
            ax.semilogy(xs, ys, color=color, lw=0.8, alpha=0.25,
                        marker=".", ms=3, ls=style["ls"], zorder=2)

        # --- Best distributed per ENR for this distance ---
        enr_best: dict[int, tuple[float, str]] = {}
        for r in dist:
            e = r["enr_hz"]
            if e not in enr_best or r["ler"] < enr_best[e][0]:
                enr_best[e] = (r["ler"], r["protocol"])

        # Draw best line, colored by protocol per segment
        best_by_prot: dict[str, list[tuple]] = defaultdict(list)
        for e, (ler, prot) in sorted(enr_best.items()):
            best_by_prot[prot].append((khz(e), ler))

        for prot, pts in best_by_prot.items():
            pts.sort()
            color = PROTOCOL_COLORS.get(prot, "#333333")
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            line, = ax.semilogy(xs, ys, color=color, lw=2.5, alpha=1.0,
                                marker=style["marker"], ms=7, ls="-", zorder=5)
            best_legend_entries.append(
                mlines.Line2D([], [], color=color, lw=2.0,
                              marker=style["marker"], ms=6,
                              label=f"★ best d={d}: {PROTOCOL_LABELS.get(prot, prot)}"))

        # Thin black spine connecting all best points
        if enr_best:
            xs_all = sorted(enr_best)
            ys_all = [enr_best[e][0] for e in xs_all]
            ax.semilogy([khz(e) for e in xs_all], ys_all,
                        color="#333333", lw=0.8, alpha=0.4, ls="--", zorder=4,
                        marker=style["marker"], ms=4)

    # --- Legend ---
    handles = []
    # Protocol color guide (faint background lines)
    handles.append(mlines.Line2D([], [], color="none", label="— distillation protocols —"))
    for prot in sorted(seen_protocols, key=lambda p: list(PROTOCOL_COLORS).index(p)
                       if p in PROTOCOL_COLORS else 99):
        color = PROTOCOL_COLORS.get(prot, "#333333")
        handles.append(mlines.Line2D([], [], color=color, lw=1.5, alpha=0.7,
                                     label=PROTOCOL_LABELS.get(prot, prot)))

    # Best-line entries per distance
    handles.append(mlines.Line2D([], [], color="none", label="— best per ENR —"))
    handles.extend(best_legend_entries)

    # Monolithic entries (from ax legend)
    ax_handles, ax_labels = ax.get_legend_handles_labels()
    for h, l in zip(ax_handles, ax_labels):
        if "Monolithic" in l:
            handles.append(h)

    # Distance marker guide
    handles.append(mlines.Line2D([], [], color="none", label="— distance marker —"))
    for d in distances:
        style = DISTANCE_STYLE.get(d, {})
        handles.append(mlines.Line2D([], [], color="#555555", lw=1.2,
                                     marker=style.get("marker", "o"), ms=5,
                                     ls=style.get("ls", "-"), label=f"d={d}"))

    ax.legend(handles=handles, fontsize=7.5, loc="best",
              framealpha=0.85, handlelength=2.2, ncol=1)

    title = f"{mode.replace('_', ' ').title()}  —  d={{{', '.join(str(d) for d in distances)}}}"
    ax.set_title(title, fontsize=13)
    ax.set_xlabel("ENR (kHz)", fontsize=11)
    ax.set_ylabel("Logical Error Rate", fontsize=11)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}"))
    ax.grid(True, which="both", ls=":", alpha=0.4)

    out_path = run_dir / f"{mode}.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved {out_path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--modes", nargs="+", default=MODES)
    args = parser.parse_args()

    records = load_results(args.run_dir)
    if not records:
        print("No results found.", file=sys.stderr)
        return 1

    distances = sorted(set(r["distance"] for r in records))
    print(f"Loaded {len(records)} results  —  distances: {distances}")
    for d in distances:
        mc = sum(1 for r in records if r["arch"] == "monolithic" and r["distance"] == d)
        dc = sum(1 for r in records if r["arch"] == "distributed" and r["distance"] == d)
        print(f"  d={d}: {mc} monolithic, {dc} distributed")

    for mode in [m for m in args.modes if m in MODES]:
        plot_mode(records, mode, args.run_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
