#!/usr/bin/env python3
"""
Plot RCX regime-sweep outputs.

Generates:
  - regime_overview.png  — 3-strip shared-axis overview
  - racing_lines.png     — all (protocol, k) candidates with optimal path highlighted
  - qubit_budget.png     — heatmap of qubit usage vs (k, m) with feasibility boundary
"""

from __future__ import annotations

import argparse
from pathlib import Path

import re

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


PROTOCOL_ORDER = ["none", "2to1_pumping", "3to1_pumping", "2to1_recurrence", "3to1_recurrence"]
PROTOCOL_LABELS = {
    "none":            "Raw EPR",
    "2to1_pumping":    "2→1 Pump",
    "3to1_pumping":    "3→1 Pump",
    "2to1_recurrence": "2→1 Recur",
    "3to1_recurrence": "3→1 Recur",
}
PROTOCOL_COLORS = {
    "none":            "#888888",
    "2to1_pumping":    "#4C72B0",
    "3to1_pumping":    "#55A868",
    "2to1_recurrence": "#C44E52",
    "3to1_recurrence": "#DD8452",
}


def khz_formatter(v, _):
    if v >= 1e9:
        return f"{v/1e9:g}THz"
    if v >= 1e6:
        return f"{v/1e6:g}GHz"
    if v >= 1e3:
        return f"{v/1e3:g}MHz"
    return f"{v:g}kHz"


def parse_summary(run_dir: Path) -> dict:
    """Extract code_distance and num_remote_cnots from summary.txt."""
    text = (run_dir / "summary.txt").read_text()
    d = int(re.search(r"code_distance:\s*(\d+)", text).group(1))
    N = int(re.search(r"num_remote_cnots:\s*(\d+)", text).group(1))
    return {"d": d, "N": N}


def load_csvs(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    best = pd.read_csv(run_dir / "best_by_rate.csv")
    regimes = pd.read_csv(run_dir / "regime_boundaries.csv")
    all_cands = pd.read_csv(run_dir / "all_candidates.csv")
    for df in (best, regimes, all_cands):
        df["entanglement_rate_khz"] = df["entanglement_rate_hz"] / 1e3
    best["protocol_color"] = best["protocol"].map(PROTOCOL_COLORS)
    return best, regimes, all_cands


def shade_regimes(ax, best: pd.DataFrame, x_max: float) -> None:
    df = best.sort_values("entanglement_rate_hz")
    boundaries = list(df["entanglement_rate_khz"]) + [x_max]
    for i, (_, row) in enumerate(df.iterrows()):
        ax.axvspan(boundaries[i], boundaries[i + 1],
                   color=PROTOCOL_COLORS[row["protocol"]], alpha=0.10, linewidth=0)


def draw_regime_vlines(ax, regimes: pd.DataFrame) -> None:
    for _, row in regimes.iloc[1:].iterrows():
        ax.axvline(row["entanglement_rate_khz"], color="#444444",
                   linestyle="--", linewidth=0.8, alpha=0.55)


def plot_stepped_line(ax, df: pd.DataFrame, y_col: str, x_max: float) -> None:
    """Draw a protocol-colored stepped line across all data points."""
    df = df.sort_values("entanglement_rate_hz")
    xs = list(df["entanglement_rate_khz"]) + [x_max]
    for i, (_, row) in enumerate(df.iterrows()):
        color = PROTOCOL_COLORS[row["protocol"]]
        x0, x1 = xs[i], xs[i + 1]
        y = row[y_col]
        ax.plot([x0, x1], [y, y], color=color, linewidth=2.5, solid_capstyle="butt")
        if i + 1 < len(df):
            y_next = df.iloc[i + 1][y_col]
            ax.plot([x1, x1], [y, y_next],
                    color=PROTOCOL_COLORS[df.iloc[i + 1]["protocol"]],
                    linewidth=2.5, solid_capstyle="butt")
    ax.scatter(df["entanglement_rate_khz"], df[y_col],
               c=df["protocol_color"], s=28, zorder=5)


# ---------------------------------------------------------------------------
# Plot 1: regime overview (3-strip)
# ---------------------------------------------------------------------------

def save_overview(best: pd.DataFrame, regimes: pd.DataFrame, out_dir: Path) -> None:
    df = best.sort_values("entanglement_rate_hz")
    x_min = df["entanglement_rate_khz"].iloc[0] * 0.7
    x_max = df["entanglement_rate_khz"].iloc[-1] * 1.4

    fig, axes = plt.subplots(
        3, 1, figsize=(10, 7),
        gridspec_kw={"height_ratios": [3, 1, 1]},
        constrained_layout=True,
    )
    fig.get_layout_engine().set(hspace=0.06)

    # ── Strip 1: effective channel error ──────────────────────────────────
    ax = axes[0]
    shade_regimes(ax, df, x_max)
    draw_regime_vlines(ax, regimes)
    plot_stepped_line(ax, df, "effective_channel_error", x_max)
    ax.set_yscale("log")
    ax.set_ylabel("Effective Channel Error", fontsize=11)
    ax.grid(True, which="both", alpha=0.2)
    ax.set_xlim(x_min, x_max)

    legend_patches = [
        mpatches.Patch(color=PROTOCOL_COLORS[p], label=PROTOCOL_LABELS[p])
        for p in PROTOCOL_ORDER if p in df["protocol"].values
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc="upper right",
              framealpha=0.9, ncol=len(legend_patches))
    ax.set_xticklabels([])

    # ── Strip 2: optimal k ───────────────────────────────────────────────
    ax = axes[1]
    shade_regimes(ax, df, x_max)
    draw_regime_vlines(ax, regimes)
    plot_stepped_line(ax, df, "rounds", x_max)
    k_vals = sorted(df["rounds"].unique())
    ax.set_yticks(k_vals)
    ax.set_yticklabels(["raw" if v == 0 else str(v) for v in k_vals])
    ax.set_ylabel("k", fontsize=11)
    ax.grid(True, axis="x", alpha=0.2)
    ax.set_xlim(x_min, x_max)
    ax.set_xticklabels([])

    # ── Strip 3: optimal m* ──────────────────────────────────────────────
    ax = axes[2]
    shade_regimes(ax, df, x_max)
    draw_regime_vlines(ax, regimes)
    plot_stepped_line(ax, df, "m_star", x_max)
    m_vals = sorted(df["m_star"].unique())
    ax.set_yticks(m_vals)
    ax.set_yticklabels([str(v) for v in m_vals])
    ax.set_ylabel("m*", fontsize=11)
    ax.set_xlabel("Entanglement Rate", fontsize=11)
    ax.grid(True, axis="x", alpha=0.2)
    ax.set_xlim(x_min, x_max)

    for ax in axes:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(khz_formatter))

    fig.suptitle("Optimal Distillation Protocol vs Entanglement Rate", fontsize=13)
    fig.savefig(out_dir / "regime_overview.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  regime_overview.png")


# ---------------------------------------------------------------------------
# Plot 2: racing lines
# ---------------------------------------------------------------------------

K_LINESTYLES = {1: "solid", 2: "dashed", 3: "dashdot", 4: (0, (3,1,1,1)), 5: "dotted"}

def save_racing_lines(best: pd.DataFrame, all_cands: pd.DataFrame,
                      regimes: pd.DataFrame, out_dir: Path) -> None:
    best_df = best.sort_values("entanglement_rate_hz")
    # Build a set of optimal (enr, protocol, rounds) for bold path lookup
    optimal_set = set(
        zip(best_df["entanglement_rate_hz"], best_df["protocol"], best_df["rounds"])
    )

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)

    # Plot one line per (protocol, rounds) combination
    seen_labels: set[str] = set()
    protocols_in_data = all_cands["protocol"].unique()

    for protocol in PROTOCOL_ORDER:
        if protocol not in protocols_in_data:
            continue
        color = PROTOCOL_COLORS[protocol]
        proto_df = all_cands[all_cands["protocol"] == protocol]
        k_values = sorted(proto_df["rounds"].unique())

        for k in k_values:
            df_k = proto_df[proto_df["rounds"] == k].sort_values("entanglement_rate_hz")
            if df_k.empty:
                continue

            # Split into feasible / infeasible segments and draw separately
            xs = df_k["entanglement_rate_khz"].values
            ys = df_k["effective_channel_error"].values
            feasible = df_k["fully_feasible"].astype(bool).values

            label = f"{PROTOCOL_LABELS[protocol]}" if protocol not in seen_labels else None
            if label:
                seen_labels.add(protocol)

            ls = K_LINESTYLES.get(k, "solid")

            # Faint background line (full range)
            ax.plot(xs, ys, color=color, linewidth=0.8, linestyle=ls,
                    alpha=0.25, zorder=2)

            # Solid where feasible
            feas_mask = feasible
            if feas_mask.any():
                ax.plot(xs[feas_mask], ys[feas_mask], color=color,
                        linewidth=1.4, linestyle=ls, alpha=0.55, zorder=3,
                        label=label)
                label = None  # only label once per protocol

            # Dashed where infeasible
            infeas_mask = ~feasible
            if infeas_mask.any():
                ax.plot(xs[infeas_mask], ys[infeas_mask], color=color,
                        linewidth=1.0, linestyle="dotted", alpha=0.3, zorder=2)

    # Bold optimal path on top
    x_max = best_df["entanglement_rate_khz"].iloc[-1] * 1.4
    plot_stepped_line_bold(ax, best_df, x_max)

    # Regime vlines
    draw_regime_vlines(ax, regimes)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Entanglement Rate", fontsize=11)
    ax.set_ylabel("Effective Channel Error", fontsize=11)
    ax.set_title("Distillation Protocol Landscape vs Entanglement Rate", fontsize=13)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(khz_formatter))
    ax.grid(True, which="both", alpha=0.2)

    # Protocol color legend
    proto_patches = [
        mpatches.Patch(color=PROTOCOL_COLORS[p], label=PROTOCOL_LABELS[p])
        for p in PROTOCOL_ORDER if p in all_cands["protocol"].values
    ]
    # k linestyle legend
    k_all = sorted(all_cands[all_cands["protocol"] != "none"]["rounds"].unique())
    k_lines = [
        plt.Line2D([0], [0], color="black", linestyle=K_LINESTYLES.get(k, "solid"),
                   linewidth=1.5, label=f"k={k}")
        for k in k_all
    ]
    bold_line = plt.Line2D([0], [0], color="black", linewidth=3, label="Optimal path")

    legend1 = ax.legend(handles=proto_patches, fontsize=9, loc="upper right",
                        framealpha=0.9, title="Protocol")
    ax.add_artist(legend1)
    ax.legend(handles=k_lines + [bold_line], fontsize=9, loc="lower left",
              framealpha=0.9, title="Rounds / path")

    fig.savefig(out_dir / "racing_lines.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  racing_lines.png")


def plot_stepped_line_bold(ax, df: pd.DataFrame, x_max: float) -> None:
    """Draw bold black-outlined optimal path."""
    df = df.sort_values("entanglement_rate_hz")
    xs = list(df["entanglement_rate_khz"]) + [x_max]
    for i, (_, row) in enumerate(df.iterrows()):
        color = PROTOCOL_COLORS[row["protocol"]]
        x0, x1 = xs[i], xs[i + 1]
        y = row["effective_channel_error"]
        # White outline
        ax.plot([x0, x1], [y, y], color="white", linewidth=5.5,
                solid_capstyle="butt", zorder=7)
        ax.plot([x0, x1], [y, y], color=color, linewidth=3.0,
                solid_capstyle="butt", zorder=8)
        if i + 1 < len(df):
            y_next = df.iloc[i + 1]["effective_channel_error"]
            ax.plot([x1, x1], [y, y_next], color="white", linewidth=5.5,
                    solid_capstyle="butt", zorder=7)
            ax.plot([x1, x1], [y, y_next],
                    color=PROTOCOL_COLORS[df.iloc[i + 1]["protocol"]],
                    linewidth=3.0, solid_capstyle="butt", zorder=8)
    ax.scatter(df["entanglement_rate_khz"], df["effective_channel_error"],
               c=df["protocol_color"], s=40, zorder=9,
               edgecolors="white", linewidths=0.8)


# ---------------------------------------------------------------------------
# Plot 3: qubit budget heatmap
# ---------------------------------------------------------------------------

DISTILL_PROTOCOLS = ["2to1_pumping", "3to1_pumping", "2to1_recurrence", "3to1_recurrence"]
DISTILL_TITLES = {
    "2to1_pumping":    "2→1 Pumping\nR(k) = N·2k",
    "3to1_pumping":    "3→1 Pumping\nR(k) = N·3k",
    "2to1_recurrence": "2→1 Recurrence\nR(k) = N·2ᵏ",
    "3to1_recurrence": "3→1 Recurrence\nR(k) = N·3ᵏ",
}


def qubit_slots(protocol: str, k: int, m: int, N: int) -> int:
    """Total EPR slots = N * R(k) * m."""
    n = 3 if "3to1" in protocol else 2
    R_k = n * k if "pumping" in protocol else n ** k
    return N * R_k * m


def save_qubit_budget(run_dir: Path, d: int, N: int,
                      max_k: int = 5, max_m: int = 8) -> None:
    # Budget = monolithic equivalent qubit count = 4d²+d+2 (from lattice surgery circuit).
    # Distillation needs 2 qubits per EPR slot, so max EPR slots = budget / 2.
    qubit_budget = 4 * d * d + d + 2
    epr_budget = qubit_budget // 2   # max EPR slots (each costs 2 qubits)
    ks = list(range(1, max_k + 1))
    ms = list(range(1, max_m + 1))

    # Colormap: green→yellow for feasible (0→1), flat deep red for infeasible (>1)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "budget",
        [(0.0, "#2ecc71"), (0.6, "#f9e04b"), (1.0, "#8b0000")],
    )
    cmap.set_over("#8b0000")
    norm = mcolors.Normalize(vmin=0, vmax=1.0)

    fig, axes = plt.subplots(
        2, 2, figsize=(11, 7.5), constrained_layout=True,
        sharex=True, sharey=True,
    )
    fig.suptitle(
        f"Distillation Qubit Budget  (d={d},  budget={qubit_budget} qubits = 4d²+d+2,  N={N} remote CNOTs)",
        fontsize=12,
    )

    for ax, protocol in zip(axes.flat, DISTILL_PROTOCOLS):
        # Build fraction matrix: rows = m (high→low), cols = k
        data = np.zeros((len(ms), len(ks)))
        for j, k in enumerate(ks):
            for i, m in enumerate(ms):
                slots = qubit_slots(protocol, k, m, N)
                data[i, j] = slots / epr_budget

        # Flip so m=1 is at bottom
        data_plot = data[::-1]
        ms_plot = ms[::-1]

        im = ax.imshow(data_plot, cmap=cmap, norm=norm, aspect="auto",
                       extent=[0.5, max_k + 0.5, 0.5, max_m + 0.5])

        # Annotate each cell with qubit count
        for j, k in enumerate(ks):
            for i, m in enumerate(ms_plot):
                slots = qubit_slots(protocol, k, m, N)
                fraction = slots / epr_budget
                text_color = "white" if fraction > 1.4 or fraction < 0.35 else "black"
                ax.text(k, m, str(slots), ha="center", va="center",
                        fontsize=7.5, color=text_color, fontweight="bold")

        ax.set_title(DISTILL_TITLES[protocol], fontsize=10)
        ax.set_xticks(ks)
        ax.set_xticklabels([str(k) for k in ks])
        ax.set_yticks(ms)
        ax.set_yticklabels([str(m) for m in ms])

    for ax in axes[1]:
        ax.set_xlabel("k  (distillation rounds)", fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel("m  (backup batches)", fontsize=10)

    # Shared colorbar
    cb = fig.colorbar(im, ax=axes, fraction=0.03, pad=0.02, extend="max")
    cb.set_label("Qubit usage", fontsize=10)
    tick_fractions = [0, 0.25, 0.5, 0.75, 1.0]
    cb.set_ticks(tick_fractions)
    cb.set_ticklabels([f"{int(f * epr_budget)}" for f in tick_fractions[:-1]] + [f"≥{epr_budget} (infeasible)"])

    fig.savefig(run_dir / "qubit_budget.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  qubit_budget.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Plot RCX regime sweep outputs.")
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()

    best, regimes, all_cands = load_csvs(args.run_dir)
    params = parse_summary(args.run_dir)
    print(f"Saving plots to: {args.run_dir}")
    save_overview(best, regimes, args.run_dir)
    save_racing_lines(best, all_cands, regimes, args.run_dir)
    save_qubit_budget(args.run_dir, d=params["d"], N=params["N"])
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
