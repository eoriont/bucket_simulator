#!/usr/bin/env python3
"""
Plot ENR vs LER from simple_merge_only.enr_m_protocol_sweep_100k.json.

Groups lines by protocol; different line styles for m values.
Zero-LER points shown as upper-bound arrows (1/shots).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

SCRIPT_DIR = Path(__file__).resolve().parent
JSON_PATH = SCRIPT_DIR / "simple_merge_only.enr_dr_protocol_sweep_100k.json"
OUT_PATH = SCRIPT_DIR / "simple_merge_only.enr_dr_protocol_sweep.png"

PROTOCOL_COLORS = {
    "none":             "#888888",
    "pumping_2to1":     "#4C72B0",
    "pumping_3to1":     "#55A868",
    "recurrence_2to1":  "#C44E52",
    "recurrence_3to1":  "#DD8452",
}

PROTOCOL_LABELS = {
    "none":             "None (raw EPR)",
    "pumping_2to1":     "2→1 pumping",
    "pumping_3to1":     "3→1 pumping",
    "recurrence_2to1":  "2→1 recurrence",
    "recurrence_3to1":  "3→1 recurrence",
}

DR_STYLES  = {1: "-",  2: "--",  3: ":"}
DR_MARKERS = {1: "o", 2: "s",   3: "^"}


def load_data(path: Path) -> list[dict]:
    with open(path) as f:
        raw = json.load(f)
    rows = []
    for r in raw["results"]:
        if r.get("status") != "ok":
            continue
        lr = r["ler_result"]
        shots = lr["shots"]
        ler = lr["logical_error_rate"]
        is_ub = ler == 0.0
        if is_ub:
            ler = 1.0 / shots
        rows.append({
            "protocol": r["protocol"],
            "dr": r["distillation_rounds"],
            "enr_mhz": r["entanglement_rate"] / 1e6,
            "ler": ler,
            "is_ub": is_ub,
            "distance": r["circuit_distance"],
        })
    return rows


def main() -> None:
    rows = load_data(JSON_PATH)

    protocols = list(dict.fromkeys(r["protocol"] for r in rows))
    dr_values = sorted({r["dr"] for r in rows})

    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)

    for protocol in protocols:
        color = PROTOCOL_COLORS.get(protocol, "#333333")
        for dr in dr_values:
            pts = sorted(
                [r for r in rows if r["protocol"] == protocol and r["dr"] == dr],
                key=lambda r: r["enr_mhz"],
            )
            if not pts:
                continue
            xs = [p["enr_mhz"] for p in pts]
            ys = [p["ler"] for p in pts]
            ls = DR_STYLES.get(dr, ":")
            mk = DR_MARKERS.get(dr, "D")

            ax.semilogy(
                xs, ys,
                color=color,
                lw=1.8,
                ls=ls,
                marker=mk,
                ms=6,
                mec="white",
                mew=0.4,
                alpha=0.85,
                zorder=2,
            )
            for p in pts:
                ax.annotate(
                    f"d={p['distance']}",
                    xy=(p["enr_mhz"], p["ler"]),
                    xytext=(0, 7 if not p["is_ub"] else -12),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=6.5,
                    color=color,
                    alpha=0.75,
                )
                if p["is_ub"]:
                    ax.annotate(
                        "",
                        xy=(p["enr_mhz"], p["ler"] * 0.55),
                        xytext=(p["enr_mhz"], p["ler"]),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.2),
                    )

    # Legend: protocol lines
    proto_handles = [mlines.Line2D([], [], color="none", label="— protocol —")]
    for protocol in protocols:
        color = PROTOCOL_COLORS.get(protocol, "#333333")
        proto_handles.append(
            mlines.Line2D([], [], color=color, lw=2, label=PROTOCOL_LABELS.get(protocol, protocol))
        )

    # Legend: distillation_rounds styles
    dr_handles = [mlines.Line2D([], [], color="none", label="— distillation rounds —")]
    for dr in dr_values:
        dr_handles.append(
            mlines.Line2D([], [], color="#444444", lw=1.8, ls=DR_STYLES.get(dr, ":"),
                          marker=DR_MARKERS.get(dr, "D"), ms=5, label=f"DR={dr}")
        )

    ub_handle = mlines.Line2D([], [], color="#888888", lw=0, marker=r"$\downarrow$",
                               ms=8, label="↓ upper bound (0 failures)")

    ax.legend(
        handles=proto_handles + dr_handles + [ub_handle],
        fontsize=7.5, loc="upper right", framealpha=0.88, handlelength=2.5,
    )

    ax.set_title("Simple Merge Only — ENR vs LER  (p=1e-3, 100k shots, m=1)", fontsize=13)
    ax.set_xlabel("Entanglement Rate (MHz)", fontsize=11)
    ax.set_ylabel("Logical Error Rate", fontsize=11)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    ax.grid(True, which="both", ls=":", alpha=0.4)

    fig.savefig(OUT_PATH, dpi=150)
    plt.close(fig)
    print(f"Saved {OUT_PATH}")


if __name__ == "__main__":
    main()
