#!/usr/bin/env python3
"""
Plot LER results from neutral_atom_d9_d11_lattice_surgery runs.

One PNG per (mode, distance), e.g. xx_merge_d9.png and xx_merge_d11.png.
Each plot shows:
  - All distributed lines faint, colored by distillation protocol
  - Best distributed per ENR highlighted in bold
  - Monolithic baseline as a dashed line

Usage:
    python3 plot_results.py runs/<timestamp>
    python3 plot_results.py runs/<timestamp> --modes xx_merge
"""
from __future__ import annotations

import argparse
import importlib.util
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

MODES = ["xx_merge", "zz_merge"]

LER_RE   = re.compile(r"Logical Error Rate:\s*([\d.eE+\-]+)")
SHOTS_RE = re.compile(r"Total Shots:\s*(\d+)")
SCRIPT_DIR = Path(__file__).resolve().parent
MAX_ROUNDS = 5
REGIME_POINTS_PER_DECADE = 12
REGIME_SHADE_ALPHA = 0.10
REGIME_LABEL_ALPHA = 0.92


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


def load_sweep_module():
    sweep_path = SCRIPT_DIR.parent / "rcx_regime_sweep" / "sweep_optimal_rcx.py"
    if not sweep_path.exists():
        return None
    spec = importlib.util.spec_from_file_location("sweep_optimal_rcx", sweep_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("sweep_optimal_rcx", mod)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


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


def compute_regimes(distance: int, rate_min_hz: int, rate_max_hz: int) -> list[dict]:
    sweep = load_sweep_module()
    if sweep is None:
        return []

    cfg_path = SCRIPT_DIR / f"config_distributed_d{distance}.txt"
    if not cfg_path.exists():
        return []

    file_cfg = sweep.parse_config_file(cfg_path)
    cfg = sweep.Config(
        code_distance=distance,
        raw_epr_fidelity=sweep.get_float(file_cfg, "raw_epr_fidelity", 0.95),
        physical_error=sweep.get_float(file_cfg, "physical_error", 0.005),
        measurement_delay=sweep.get_float(file_cfg, "measurement_delay", 500e-6),
        T1=sweep.get_float(file_cfg, "T1", 119.0),
        T2=sweep.get_float(file_cfg, "T2", 12.6),
        num_remote_cnots=2 * distance - 1,
    )

    rates = sweep.rate_grid(rate_min_hz, rate_max_hz, REGIME_POINTS_PER_DECADE)
    regimes: list[dict] = []
    previous_signature = None

    for enr in rates:
        candidates = [
            sweep.evaluate(protocol, rounds, enr, cfg)
            for protocol, rounds in sweep.protocol_candidates(MAX_ROUNDS)
        ]
        best = sweep.choose_best(candidates)
        signature = sweep.regime_signature(best)
        if signature != previous_signature:
            regimes.append({
                "start_hz": enr,
                "protocol": best.protocol,
                "rounds": best.rounds,
                "m_star": best.m_star,
            })
            previous_signature = signature

    return regimes


def shade_regimes(ax, regimes: list[dict], x_min_khz: float, x_max_khz: float) -> list[mlines.Line2D]:
    if not regimes:
        return []

    handles: list[mlines.Line2D] = []
    seen: set[tuple[str, int]] = set()

    starts_khz = [khz(r["start_hz"]) for r in regimes]
    boundaries = starts_khz + [x_max_khz]

    for i, regime in enumerate(regimes):
        start = max(starts_khz[i], x_min_khz)
        end = min(boundaries[i + 1], x_max_khz)
        if end <= start:
            continue
        color = PROTOCOL_COLORS.get(regime["protocol"], "#333333")
        ax.axvspan(start, end, color=color, alpha=REGIME_SHADE_ALPHA, linewidth=0, zorder=0)
        if i > 0:
            ax.axvline(start, color="#555555", lw=0.9, ls="--", alpha=0.45, zorder=1)

        mid = (start * end) ** 0.5
        label = PROTOCOL_LABELS.get(regime["protocol"], regime["protocol"])
        if regime["protocol"] != "none":
            label = f"{label}  k={regime['rounds']}"
        ax.text(mid, 0.98, label, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8, color="#222222",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                          edgecolor="none", alpha=REGIME_LABEL_ALPHA),
                zorder=6)

        key = (regime["protocol"], regime["rounds"])
        if key not in seen:
            legend_label = PROTOCOL_LABELS.get(regime["protocol"], regime["protocol"])
            if regime["protocol"] != "none":
                legend_label = f"{legend_label}  k={regime['rounds']}"
            handles.append(
                mlines.Line2D([], [], color=color, lw=6, alpha=REGIME_SHADE_ALPHA * 3,
                              label=legend_label)
            )
            seen.add(key)

    return handles


def plot_mode_distance(records: list[dict], mode: str, distance: int, run_dir: Path) -> None:
    records = [r for r in records if r["mode"] == mode and r["distance"] == distance]
    if not records:
        return

    fig, ax = plt.subplots(figsize=(9.4, 6.2), constrained_layout=True)

    seen_protocols: set[str] = set()
    best_legend_entries: list[mlines.Line2D] = []
    style = {"marker": "o", "ls": "-", "mono_ls": "--"} if distance == 9 else {
        "marker": "s", "ls": "-", "mono_ls": "--"
    }
    dist = [r for r in records if r["arch"] == "distributed"]
    if not dist:
        return
    rate_min_khz = min(khz(r["enr_hz"]) for r in dist)
    rate_max_khz = max(khz(r["enr_hz"]) for r in dist)
    regimes = compute_regimes(distance, min(r["enr_hz"] for r in dist), max(r["enr_hz"] for r in dist))
    regime_handles = shade_regimes(ax, regimes, rate_min_khz, rate_max_khz * 1.03)

    # --- Monolithic baseline ---
    mono_rec = [r for r in records if r["arch"] == "monolithic"]
    if mono_rec:
        ler, is_ub = mono_rec[0]["ler"], mono_rec[0]["is_upper_bound"]
        mono_label = f"Monolithic  (LER={'<' if is_ub else ''}{ler:.2e})"
        ax.axhline(ler, color="firebrick", lw=1.8, ls=style["mono_ls"], label=mono_label, zorder=3)

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
        ax.semilogy(xs, ys, color=color, lw=1.6, alpha=0.72,
                    marker="o", ms=4.2, mec="white", mew=0.3,
                    ls=style["ls"], zorder=2)

    # --- Best distributed per ENR ---
    enr_best: dict[int, tuple[float, str]] = {}
    for r in dist:
        e = r["enr_hz"]
        if e not in enr_best or r["ler"] < enr_best[e][0]:
            enr_best[e] = (r["ler"], r["protocol"])

    best_by_prot: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for e, (ler, prot) in sorted(enr_best.items()):
        best_by_prot[prot].append((khz(e), ler))

    for prot, pts in best_by_prot.items():
        pts.sort()
        color = PROTOCOL_COLORS.get(prot, "#333333")
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.semilogy(xs, ys, color="white", lw=4.8, alpha=0.95,
                    marker=style["marker"], ms=8.5, ls="-", zorder=4)
        ax.semilogy(xs, ys, color=color, lw=3.0, alpha=1.0,
                    marker=style["marker"], ms=7.0, mec="white", mew=0.6,
                    ls="-", zorder=5)
        if not any(entry.get_label() == f"Best: {PROTOCOL_LABELS.get(prot, prot)}"
                   for entry in best_legend_entries):
            best_legend_entries.append(
                mlines.Line2D([], [], color=color, lw=3.0,
                              marker=style["marker"], ms=6,
                              label=f"Best: {PROTOCOL_LABELS.get(prot, prot)}"))

    if enr_best:
        xs_all = sorted(enr_best)
        ys_all = [enr_best[e][0] for e in xs_all]
        ax.semilogy([khz(e) for e in xs_all], ys_all,
                    color="#222222", lw=1.0, alpha=0.45, ls="--", zorder=3,
                    marker=style["marker"], ms=3.8)

    # --- Legend ---
    handles = []
    # Protocol color guide
    handles.append(mlines.Line2D([], [], color="none", label="— distillation protocols —"))
    for prot in sorted(seen_protocols, key=lambda p: list(PROTOCOL_COLORS).index(p)
                       if p in PROTOCOL_COLORS else 99):
        color = PROTOCOL_COLORS.get(prot, "#333333")
        handles.append(mlines.Line2D([], [], color=color, lw=2.4, alpha=0.9,
                                     label=PROTOCOL_LABELS.get(prot, prot)))

    if regime_handles:
        handles.append(mlines.Line2D([], [], color="none", label="— RCX best-regime bands —"))
        handles.extend(regime_handles)

    handles.append(mlines.Line2D([], [], color="none", label="— best per ENR —"))
    handles.extend(best_legend_entries)

    # Monolithic entries (from ax legend)
    ax_handles, ax_labels = ax.get_legend_handles_labels()
    for h, l in zip(ax_handles, ax_labels):
        if "Monolithic" in l:
            handles.append(h)

    ax.legend(handles=handles, fontsize=7.5, loc="best",
              framealpha=0.85, handlelength=2.2, ncol=1)

    title = f"{mode.replace('_', ' ').title()}  —  d={distance}"
    ax.set_title(title, fontsize=13)
    ax.set_xlabel("ENR (kHz)", fontsize=11)
    ax.set_ylabel("Logical Error Rate", fontsize=11)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.0f}"))
    ax.set_xlim(rate_min_khz * 0.95, rate_max_khz * 1.05)
    ax.grid(True, which="both", ls=":", alpha=0.4)

    out_path = run_dir / f"{mode}_d{distance}.png"
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
        for distance in distances:
            plot_mode_distance(records, mode, distance, args.run_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
