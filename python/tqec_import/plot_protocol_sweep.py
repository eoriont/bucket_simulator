#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


PROTOCOL_COLORS = {
    "local_baseline": "#222222",
    "none": "#888888",
    "pumping_2to1": "#4C72B0",
    "pumping_3to1": "#55A868",
    "recurrence_2to1": "#C44E52",
    "recurrence_3to1": "#DD8452",
}

PROTOCOL_LABELS = {
    "local_baseline": "Local baseline",
    "none": "No distillation",
    "pumping_2to1": "2→1 pumping",
    "pumping_3to1": "3→1 pumping",
    "recurrence_2to1": "2→1 recurrence",
    "recurrence_3to1": "3→1 recurrence",
}

PROTOCOL_ORDER = ["local_baseline", "none", "pumping_2to1", "pumping_3to1", "recurrence_2to1", "recurrence_3to1"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a TQEC ENR/protocol sweep JSON.")
    parser.add_argument("json_path", type=Path, help="Path to sweep JSON output.")
    return parser.parse_args()


def protocol_sort_key(protocol: str) -> tuple[int, str]:
    try:
        return (PROTOCOL_ORDER.index(protocol), protocol)
    except ValueError:
        return (len(PROTOCOL_ORDER), protocol)


def load_rows(path: Path) -> tuple[dict, list[dict]]:
    raw = json.loads(path.read_text())
    rows = []
    for r in raw["results"]:
        if r.get("status") != "ok":
            continue
        lr = r["ler_result"]
        shots = lr["shots"]
        ler = lr["logical_error_rate"]
        rows.append(
            {
                "protocol": r["protocol"],
                "m": r["m"],
                "dr": r["distillation_rounds"],
                "enr_hz": r["entanglement_rate"],
                "enr_mhz": r["entanglement_rate"] / 1e6,
                "ler": ler,
                "plot_ler": ler if ler > 0 else 0.5 / shots,
                "is_upper_bound": ler == 0.0,
                "distance": r["circuit_distance"],
            }
        )
    return raw, rows


def choose_best_rows(rows: list[dict]) -> list[dict]:
    best_by_point: dict[tuple[str, float], dict] = {}
    for row in rows:
        key = (row["protocol"], row["enr_hz"])
        chosen = best_by_point.get(key)
        if chosen is None:
            best_by_point[key] = row
            continue
        if row["plot_ler"] < chosen["plot_ler"] - 1e-18:
            best_by_point[key] = row
            continue
        if abs(row["plot_ler"] - chosen["plot_ler"]) <= 1e-18:
            if row["m"] < chosen["m"]:
                best_by_point[key] = row
            elif row["m"] == chosen["m"] and row["dr"] < chosen["dr"]:
                best_by_point[key] = row
    return list(best_by_point.values())


def rows_with_distillation_controls(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r["protocol"] != "local_baseline"]


def apply_common_axis_style(ax: plt.Axes) -> None:
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    ax.grid(True, which="both", ls=":", alpha=0.4)


def save_best_ler_plot(summary: dict, rows: list[dict], out_dir: Path, stem: str) -> Path:
    best_rows = choose_best_rows(rows)
    distance = best_rows[0]["distance"]
    shots = summary["shots_per_point"]

    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    for protocol in sorted({r["protocol"] for r in best_rows}, key=protocol_sort_key):
        pts = sorted((r for r in best_rows if r["protocol"] == protocol), key=lambda r: r["enr_mhz"])
        xs = [p["enr_mhz"] for p in pts]
        ys = [p["plot_ler"] for p in pts]
        color = PROTOCOL_COLORS.get(protocol, "#333333")
        ax.semilogy(
            xs,
            ys,
            color=color,
            lw=2.0,
            marker="o",
            ms=5.5,
            mec="white",
            mew=0.4,
            label=PROTOCOL_LABELS.get(protocol, protocol),
        )
        for p in pts:
            if p["is_upper_bound"]:
                ax.annotate(
                    "",
                    xy=(p["enr_mhz"], p["plot_ler"] * 0.55),
                    xytext=(p["enr_mhz"], p["plot_ler"]),
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.1),
                )

    fig.subplots_adjust(bottom=0.2, top=0.88)
    ax.set_title(
        f"Best Logical Error Rate vs Entanglement Generation Rate for Distance-{distance} TQEC Merge",
        fontsize=13,
    )
    ax.set_xlabel("Entanglement Generation Rate (MHz)", fontsize=11)
    ax.set_ylabel("Logical Error Rate", fontsize=11)
    apply_common_axis_style(ax)
    ax.legend(framealpha=0.92, fontsize=9)
    fig.text(
        0.5,
        0.055,
        "Each point shows the best logical error rate across feasible redundancy values m",
        ha="center",
        va="bottom",
        fontsize=8.8,
    )
    fig.text(
        0.5,
        0.028,
        f"Distillation depth is selected automatically; zero-failure points plotted at 0.5 / shots, shots per point = {shots:,}",
        ha="center",
        va="bottom",
        fontsize=8.1,
    )
    out_path = out_dir / f"{stem}_best_ler_thesis.png"
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out_path


def save_selected_m_plot(rows: list[dict], out_dir: Path, stem: str) -> Path:
    best_rows = choose_best_rows(rows_with_distillation_controls(rows))
    distance = best_rows[0]["distance"]
    enr_values = sorted({r["enr_mhz"] for r in best_rows})
    protocols = sorted({r["protocol"] for r in best_rows}, key=protocol_sort_key)
    x_positions = list(range(len(enr_values)))
    width = 0.8 / max(1, len(protocols))

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    for idx, protocol in enumerate(protocols):
        pts = sorted((r for r in best_rows if r["protocol"] == protocol), key=lambda r: r["enr_mhz"])
        color = PROTOCOL_COLORS.get(protocol, "#333333")
        enr_to_m = {p["enr_mhz"]: p["m"] for p in pts}
        xs = []
        ys = []
        for base_x, enr in zip(x_positions, enr_values):
            if enr not in enr_to_m:
                continue
            xs.append(base_x + (idx - (len(protocols) - 1) / 2) * width)
            ys.append(enr_to_m[enr])
        ax.bar(xs, ys, width=width * 0.92, color=color, alpha=0.88, label=PROTOCOL_LABELS.get(protocol, protocol))

    fig.subplots_adjust(bottom=0.18, top=0.88)
    ax.set_title(
        f"Selected Redundancy vs Entanglement Generation Rate for Distance-{distance} TQEC Merge",
        fontsize=13,
    )
    ax.set_xlabel("Entanglement Generation Rate (MHz)", fontsize=11)
    ax.set_ylabel("Selected redundancy m", fontsize=11)
    ax.set_yticks(sorted({r["m"] for r in best_rows}))
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{x:g}" for x in enr_values])
    apply_common_axis_style(ax)
    ax.grid(True, axis="y", ls=":", alpha=0.4)
    ax.legend(framealpha=0.92, fontsize=9)
    fig.text(
        0.5,
        0.03,
        "Local baseline is omitted because it does not use RCX distillation or redundancy selection",
        ha="center",
        va="bottom",
        fontsize=8.2,
    )
    out_path = out_dir / f"{stem}_selected_m_thesis.png"
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out_path


def save_selected_dr_plot(rows: list[dict], out_dir: Path, stem: str) -> Path:
    best_rows = choose_best_rows(rows_with_distillation_controls(rows))
    distance = best_rows[0]["distance"]
    enr_values = sorted({r["enr_mhz"] for r in best_rows})
    protocols = sorted({r["protocol"] for r in best_rows}, key=protocol_sort_key)
    x_positions = list(range(len(enr_values)))
    width = 0.8 / max(1, len(protocols))

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    for idx, protocol in enumerate(protocols):
        pts = sorted((r for r in best_rows if r["protocol"] == protocol), key=lambda r: r["enr_mhz"])
        color = PROTOCOL_COLORS.get(protocol, "#333333")
        enr_to_dr = {p["enr_mhz"]: p["dr"] for p in pts}
        xs = []
        ys = []
        for base_x, enr in zip(x_positions, enr_values):
            if enr not in enr_to_dr:
                continue
            xs.append(base_x + (idx - (len(protocols) - 1) / 2) * width)
            ys.append(enr_to_dr[enr])
        ax.bar(xs, ys, width=width * 0.92, color=color, alpha=0.88, label=PROTOCOL_LABELS.get(protocol, protocol))

    fig.subplots_adjust(bottom=0.18, top=0.88)
    ax.set_title(
        f"Selected Distillation Depth vs Entanglement Generation Rate for Distance-{distance} TQEC Merge",
        fontsize=13,
    )
    ax.set_xlabel("Entanglement Generation Rate (MHz)", fontsize=11)
    ax.set_ylabel("Selected distillation depth (DR)", fontsize=11)
    ax.set_yticks(sorted({r["dr"] for r in best_rows}))
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f"{x:g}" for x in enr_values])
    apply_common_axis_style(ax)
    ax.grid(True, axis="y", ls=":", alpha=0.4)
    ax.legend(framealpha=0.92, fontsize=9)
    fig.text(
        0.5,
        0.03,
        "Local baseline is omitted because it does not use RCX distillation depth selection",
        ha="center",
        va="bottom",
        fontsize=8.2,
    )
    out_path = out_dir / f"{stem}_selected_dr_thesis.png"
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot(json_path: Path) -> list[Path]:
    summary, rows = load_rows(json_path)
    if not rows:
        raise SystemExit(f"No successful results in {json_path}")
    out_dir = json_path.parent
    stem = json_path.stem
    return [
        save_best_ler_plot(summary, rows, out_dir, stem),
        save_selected_m_plot(rows, out_dir, stem),
        save_selected_dr_plot(rows, out_dir, stem),
    ]


def main() -> None:
    args = parse_args()
    out_paths = plot(args.json_path.expanduser().resolve())
    for out_path in out_paths:
        print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
