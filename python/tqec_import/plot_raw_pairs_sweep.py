#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


MODE_ORDER = ["none", "lside", "rside", "twoside", "ss_mid"]
MODE_COLORS = {
    "none": "#4C4C4C",
    "lside": "#1F77B4",
    "rside": "#D62728",
    "twoside": "#2CA02C",
    "ss_mid": "#9467BD",
}
MODE_LABELS = {
    "none": "nosuper",
    "lside": "lside",
    "rside": "rside",
    "twoside": "twoside",
    "ss_mid": "ss_mid",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a TQEC raw-pairs sweep run directory."
    )
    parser.add_argument("run_dir", type=Path, help="Run directory containing summary.json.")
    return parser.parse_args()


def load_summary(run_dir: Path) -> dict:
    summary_path = run_dir / "summary.json"
    if not summary_path.exists():
        raise SystemExit(f"Missing summary.json in {run_dir}")
    return json.loads(summary_path.read_text())


def mode_sort_key(mode: str) -> tuple[int, str]:
    try:
        return (MODE_ORDER.index(mode), mode)
    except ValueError:
        return (len(MODE_ORDER), mode)


def iter_mode_rows(summary: dict) -> list[dict]:
    rows = []
    for mode_result in summary["results"]:
        mode = mode_result["deformation_mode"]
        for point in mode_result["results"]:
            ler = point["ler_result"]["logical_error_rate"]
            shots = point["ler_result"]["shots"]
            rows.append(
                {
                    "mode": mode,
                    "enr_mhz": point["entanglement_rate"] / 1e6,
                    "ler": ler,
                    "plot_ler": ler if ler > 0 else 0.5 / shots,
                    "is_upper_bound": ler == 0,
                    "required_epr_pairs_per_round": point["required_epr_pairs_per_round"],
                }
            )
    return rows


def add_shared_title(fig: plt.Figure, summary: dict) -> None:
    fig.suptitle(
        (
            f"Raw Pairs Sweep  d={summary['distance']}  "
            f"k={summary['circuit_k']}  rounds={summary['rounds_per_phase']}  "
            f"shots={summary['shots_per_point']}"
        ),
        fontsize=13,
    )


def plot_ler_vs_enr(rows: list[dict], run_dir: Path, summary: dict) -> Path:
    fig, ax = plt.subplots(figsize=(9, 5.8), constrained_layout=True)

    for mode in sorted({row["mode"] for row in rows}, key=mode_sort_key):
        pts = sorted((row for row in rows if row["mode"] == mode), key=lambda row: row["enr_mhz"])
        xs = [row["enr_mhz"] for row in pts]
        ys = [row["plot_ler"] for row in pts]
        color = MODE_COLORS.get(mode, "#333333")
        label = MODE_LABELS.get(mode, mode)

        ax.plot(
            xs,
            ys,
            marker="o",
            markersize=5,
            linewidth=1.9,
            color=color,
            label=label,
        )

        for row in pts:
            if row["is_upper_bound"]:
                ax.annotate(
                    "",
                    xy=(row["enr_mhz"], row["plot_ler"] * 0.55),
                    xytext=(row["enr_mhz"], row["plot_ler"]),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.1),
                )

    ax.set_xlabel("Entanglement Rate (MHz)")
    ax.set_ylabel("Logical Error Rate")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.legend(framealpha=0.9, fontsize=9)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    add_shared_title(fig, summary)
    out_path = run_dir / "ler_vs_enr.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    fig.savefig(run_dir / "ler_vs_enr.svg", bbox_inches="tight")
    fig.savefig(run_dir / "ler_vs_enr.pdf", bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_epr_pairs_vs_enr(rows: list[dict], run_dir: Path, summary: dict) -> Path:
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)

    for mode in sorted({row["mode"] for row in rows}, key=mode_sort_key):
        pts = sorted((row for row in rows if row["mode"] == mode), key=lambda row: row["enr_mhz"])
        xs = [row["enr_mhz"] for row in pts]
        ys = [row["required_epr_pairs_per_round"] for row in pts]
        color = MODE_COLORS.get(mode, "#333333")
        label = MODE_LABELS.get(mode, mode)
        ax.plot(
            xs,
            ys,
            marker="o",
            markersize=5,
            linewidth=1.9,
            color=color,
            label=label,
        )

    ax.set_xlabel("Entanglement Rate (MHz)")
    ax.set_ylabel("Required EPR Pairs / Round")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(framealpha=0.9, fontsize=9)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))

    add_shared_title(fig, summary)
    out_path = run_dir / "required_epr_pairs_vs_enr.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    fig.savefig(run_dir / "required_epr_pairs_vs_enr.svg", bbox_inches="tight")
    fig.savefig(run_dir / "required_epr_pairs_vs_enr.pdf", bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.expanduser().resolve()
    summary = load_summary(run_dir)
    rows = iter_mode_rows(summary)
    if not rows:
        raise SystemExit(f"No raw-pairs points found in {run_dir / 'summary.json'}")

    ler_path = plot_ler_vs_enr(rows, run_dir, summary)
    epr_path = plot_epr_pairs_vs_enr(rows, run_dir, summary)
    print(f"Saved: {ler_path}")
    print(f"Saved: {epr_path}")


if __name__ == "__main__":
    main()
