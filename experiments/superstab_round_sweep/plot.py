#!/usr/bin/env python3
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt


SUPER_LABELS = {
    "nosuper": "No superstabilizer",
    "middle": "Middle superstabilizer",
}

SUPER_ORDER = ["nosuper", "middle"]

SUPER_MARKERS = {"nosuper": "o", "middle": "s"}
SUPER_COLORS   = {"nosuper": "#1f77b4", "middle": "#d62728"}
LINESTYLES     = ["-", "--", ":", "-."]


def parse_result_file(path: Path) -> dict | None:
    text = path.read_text()
    shots_match = re.search(r"Total Shots:\s+(\d+)", text)
    rate_match  = re.search(r"Entanglement Rate:\s+([\d.]+)\s+MHz", text)
    ler_match   = re.search(r"Logical Error Rate:\s+([\d.eE+\-]+)", text)
    name_match  = re.match(
        r"(?P<mode>[a-z]+_[a-z]+)_(?P<super>[a-z0-9]+)_r(?P<rounds>\d+)_(?P<rate>[0-9A-Za-z]+)_d(?P<dist>\d+)_results",
        path.stem,
    )
    if not shots_match or not rate_match or not ler_match or not name_match:
        print(f"  Warning: could not parse {path.name}")
        return None

    shots       = int(shots_match.group(1))
    raw_rounds  = int(name_match.group("rounds"))
    super_type  = name_match.group("super")
    dist        = int(name_match.group("dist"))
    ler         = float(ler_match.group(1))
    eff_rounds  = raw_rounds if super_type == "nosuper" else 2 * raw_rounds
    plot_ler    = ler if ler > 0 else 0.5 / shots
    return {
        "mode":             name_match.group("mode"),
        "super":            super_type,
        "distance":         dist,
        "raw_rounds":       raw_rounds,
        "effective_rounds": eff_rounds,
        "rate_mhz":         float(rate_match.group(1)),
        "ler":              ler,
        "plot_ler":         plot_ler,
        "shots":            shots,
    }


def load_run(run_dir: Path) -> list[dict]:
    results_dir = run_dir / "results"
    if not results_dir.is_dir():
        print(f"  Warning: no results dir in {run_dir}")
        return []
    rows = []
    for path in sorted(results_dir.glob("*_results.txt")):
        rec = parse_result_file(path)
        if rec is not None:
            rows.append(rec)
    return rows


def plot(run_dirs: list[Path]):
    all_rows: list[dict] = []
    for run_dir in run_dirs:
        all_rows.extend(load_run(run_dir))

    if not all_rows:
        sys.exit("Error: no parseable result files")

    distances = sorted(set(r["distance"] for r in all_rows))
    mode      = all_rows[0]["mode"]
    rate_mhz  = all_rows[0]["rate_mhz"]

    fig, ax = plt.subplots(figsize=(8.5, 5.8))

    for super_type in SUPER_ORDER:
        for di, dist in enumerate(distances):
            pts = sorted(
                (r["effective_rounds"], r["plot_ler"], r["raw_rounds"])
                for r in all_rows
                if r["super"] == super_type and r["distance"] == dist
            )
            if not pts:
                continue
            xs     = [p[0] for p in pts]
            ys     = [p[1] for p in pts]
            raws   = [str(p[2]) for p in pts]
            ls     = LINESTYLES[di % len(LINESTYLES)]
            label  = f"{SUPER_LABELS.get(super_type, super_type)} (d={dist})"
            ax.plot(
                xs, ys,
                marker=SUPER_MARKERS[super_type],
                linestyle=ls,
                linewidth=1.8,
                markersize=4,
                color=SUPER_COLORS[super_type],
                label=label,
            )
            if len(distances) == 1:
                for x, y, raw in zip(xs, ys, raws):
                    ax.annotate(raw, (x, y), textcoords="offset points",
                                xytext=(0, 5), ha="center", fontsize=7, alpha=0.75)

    dist_str = "/".join(f"d={d}" for d in distances)
    ax.set_title(f"{mode} at {rate_mhz:g} MHz: effective rounds vs LER\n({dist_str})")
    ax.set_xlabel("Effective merge rounds")
    ax.set_ylabel("Logical Error Rate")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.35)
    ax.legend(framealpha=0.9, fontsize=8)
    ax.text(
        0.02, 0.02,
        "Middle uses 2× effective rounds per raw round. Solid=d smaller, dashed=d larger.",
        transform=ax.transAxes, fontsize=7.5, alpha=0.8,
    )
    fig.tight_layout()

    out_dir   = run_dirs[-1]
    plot_path = out_dir / "round_sweep_plot.png"
    fig.savefig(plot_path, dpi=160)
    plt.close(fig)

    csv_path = out_dir / "round_sweep.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "mode", "super", "distance", "raw_rounds",
            "effective_rounds", "rate_mhz", "shots", "ler",
        ])
        writer.writeheader()
        for row in sorted(all_rows, key=lambda r: (r["distance"], r["super"], r["raw_rounds"])):
            writer.writerow({k: row[k] for k in writer.fieldnames})

    print(f"Saved: {plot_path}")
    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python plot.py <run_dir> [run_dir2 ...]")
    plot([Path(a).expanduser().resolve() for a in sys.argv[1:]])
