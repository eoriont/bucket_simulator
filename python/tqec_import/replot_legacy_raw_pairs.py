#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


MODE_ORDER = ["none", "lside", "rside", "twoside"]
MODE_COLORS = {
    "none": "#4C4C4C",
    "lside": "#1F77B4",
    "rside": "#D62728",
    "twoside": "#2CA02C",
}
MODE_LABELS = {
    "none": "No deformation",
    "lside": "Left-side deformation",
    "rside": "Right-side deformation",
    "twoside": "Two-sided deformation",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replot legacy raw-pairs sweep JSON files with thesis-style titles."
    )
    parser.add_argument(
        "json_paths",
        type=Path,
        nargs="*",
        help="Specific legacy JSON files. Defaults to all output/tqec_import/ler/enr_raw_pairs_*.json files.",
    )
    return parser.parse_args()


def mode_sort_key(mode: str) -> tuple[int, str]:
    try:
        return (MODE_ORDER.index(mode), mode)
    except ValueError:
        return (len(MODE_ORDER), mode)


def infer_distance(summary: dict) -> int:
    circuit_k = summary.get("circuit_k")
    if circuit_k is None:
        raise ValueError("Missing circuit_k in summary.")
    return 2 * int(circuit_k) + 1


def load_rows(summary: dict) -> list[dict]:
    rows = []
    for mode_result in summary["results"]:
        mode = mode_result["deformation_mode"]
        for point in mode_result["results"]:
            ler_result = point["ler_result"]
            ler = ler_result["logical_error_rate"]
            shots = ler_result["shots"]
            rows.append(
                {
                    "mode": mode,
                    "enr_mhz": point["entanglement_rate"] / 1e6,
                    "ler": ler,
                    "plot_ler": ler if ler > 0 else 0.5 / shots,
                    "is_upper_bound": ler == 0,
                }
            )
    return rows


def plot_summary(summary: dict, json_path: Path) -> Path:
    rows = load_rows(summary)
    distance = infer_distance(summary)
    basis = summary.get("basis", "X")
    rounds_per_phase = summary.get("rounds_per_phase")
    physical_error = summary.get("physical_error", 0.001)
    shots = summary.get("shots_per_point", 0)

    fig, ax = plt.subplots(figsize=(8.8, 5.8))

    for mode in sorted({row["mode"] for row in rows}, key=mode_sort_key):
        pts = sorted((row for row in rows if row["mode"] == mode), key=lambda row: row["enr_mhz"])
        xs = [row["enr_mhz"] for row in pts]
        ys = [row["plot_ler"] for row in pts]
        color = MODE_COLORS.get(mode, "#333333")

        ax.plot(
            xs,
            ys,
            marker="o",
            markersize=5,
            linewidth=2.0,
            color=color,
            label=MODE_LABELS.get(mode, mode),
        )

        for row in pts:
            if row["is_upper_bound"]:
                ax.annotate(
                    "",
                    xy=(row["enr_mhz"], row["plot_ler"] * 0.55),
                    xytext=(row["enr_mhz"], row["plot_ler"]),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.1),
                )

    ax.set_xlabel("Entanglement Generation Rate (MHz)", fontsize=11)
    ax.set_ylabel("Logical Error Rate", fontsize=11)
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", alpha=0.35)
    ax.legend(framealpha=0.92, fontsize=9)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))

    ax.set_title(
        f"Logical Error Rate vs Entanglement Generation Rate for Distance-{distance} TQEC Merge",
        fontsize=13,
        pad=10,
    )
    fig.subplots_adjust(bottom=0.2, top=0.88)
    fig.text(
        0.5,
        0.055,
        f"X-basis raw-pair sweep, pre/merge rounds = {rounds_per_phase}, "
        f"physical error rate = {physical_error:.0e}, shots per point = {shots:,}",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.text(
        0.5,
        0.025,
        "Zero-failure points shown as 0.5 / shots",
        ha="center",
        va="bottom",
        fontsize=7.5,
        alpha=0.8,
    )

    out_path = json_path.with_name(json_path.stem + "_thesis.png")
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    args = parse_args()
    if args.json_paths:
        json_paths = [path.expanduser().resolve() for path in args.json_paths]
    else:
        json_paths = sorted(Path("output/tqec_import/ler").glob("enr_raw_pairs_*.json"))

    if not json_paths:
        raise SystemExit("No legacy raw-pairs JSON files found.")

    for json_path in json_paths:
        summary = json.loads(json_path.read_text())
        out_path = plot_summary(summary, json_path)
        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
