#!/usr/bin/env python3
"""
Plot RCX regime-sweep outputs.

Consumes:
  - best_by_rate.csv
  - regime_boundaries.csv

Generates:
  - rcx_error_vs_enr.png
  - optimal_k_vs_enr.png
  - optimal_m_vs_enr.png
  - protocol_regimes.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


MODE_ORDER = ["xx_merge", "zz_merge", "xx_split", "zz_split"]
MODE_LABELS = {
    "xx_merge": "XX Merge",
    "zz_merge": "ZZ Merge",
    "xx_split": "XX Split",
    "zz_split": "ZZ Split",
}
PROTOCOL_ORDER = ["none", "2to1_pumping", "3to1_pumping", "2to1_recurrence", "3to1_recurrence"]
PROTOCOL_LABELS = {
    "none": "raw",
    "2to1_pumping": "2to1 pump",
    "3to1_pumping": "3to1 pump",
    "2to1_recurrence": "2to1 recur",
    "3to1_recurrence": "3to1 recur",
}
PROTOCOL_TO_INT = {name: idx for idx, name in enumerate(PROTOCOL_ORDER)}


def load_csvs(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    best = pd.read_csv(run_dir / "best_by_rate.csv")
    regimes = pd.read_csv(run_dir / "regime_boundaries.csv")
    best["entanglement_rate_khz"] = best["entanglement_rate_hz"] / 1e3
    regimes["entanglement_rate_khz"] = regimes["entanglement_rate_hz"] / 1e3
    best["protocol_index"] = best["effective_protocol"].map(PROTOCOL_TO_INT)
    return best, regimes


def save_rcx_plot(best: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for ax, mode in zip(axes.flat, MODE_ORDER):
        df = best[best["mode"] == mode].sort_values("entanglement_rate_hz")
        ax.plot(df["entanglement_rate_khz"], df["remote_cnot_error"], marker="o", linewidth=2)
        infeasible = df[~df["fully_feasible"]]
        if not infeasible.empty:
            ax.scatter(infeasible["entanglement_rate_khz"], infeasible["remote_cnot_error"],
                       marker="x", s=60, label="not fully feasible")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(MODE_LABELS[mode])
        ax.set_xlabel("ENR (kHz)")
        ax.set_ylabel("Best Remote CNOT Error")
        ax.grid(True, alpha=0.3)
        if not infeasible.empty:
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "rcx_error_vs_enr.png", dpi=160)
    plt.close(fig)


def save_param_plot(best: pd.DataFrame, out_dir: Path, column: str, ylabel: str, filename: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for ax, mode in zip(axes.flat, MODE_ORDER):
        df = best[best["mode"] == mode].sort_values("entanglement_rate_hz")
        ax.step(df["entanglement_rate_khz"], df[column], where="post", linewidth=2)
        ax.scatter(df["entanglement_rate_khz"], df[column], s=25)
        ax.set_xscale("log")
        ax.set_title(MODE_LABELS[mode])
        ax.set_xlabel("ENR (kHz)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / filename, dpi=160)
    plt.close(fig)


def save_protocol_plot(best: pd.DataFrame, regimes: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for ax, mode in zip(axes.flat, MODE_ORDER):
        df = best[best["mode"] == mode].sort_values("entanglement_rate_hz")
        ax.step(df["entanglement_rate_khz"], df["protocol_index"], where="post", linewidth=2)
        ax.scatter(df["entanglement_rate_khz"], df["protocol_index"], s=25)
        mode_regimes = regimes[regimes["mode"] == mode]
        for _, row in mode_regimes.iterrows():
            ax.axvline(row["entanglement_rate_khz"], color="gray", linestyle=":", alpha=0.4)
        ax.set_xscale("log")
        ax.set_title(MODE_LABELS[mode])
        ax.set_xlabel("ENR (kHz)")
        ax.set_ylabel("Protocol")
        ax.set_yticks(list(PROTOCOL_TO_INT.values()))
        ax.set_yticklabels([PROTOCOL_LABELS[p] for p in PROTOCOL_ORDER])
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "protocol_regimes.png", dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot RCX regime sweep outputs.")
    parser.add_argument("run_dir", type=Path, help="Run directory produced by sweep_optimal_rcx.py")
    args = parser.parse_args()

    best, regimes = load_csvs(args.run_dir)
    save_rcx_plot(best, args.run_dir)
    save_param_plot(best, args.run_dir, "requested_rounds", "Optimal k", "optimal_k_vs_enr.png")
    save_param_plot(best, args.run_dir, "requested_backup_batches", "Optimal m", "optimal_m_vs_enr.png")
    save_protocol_plot(best, regimes, args.run_dir)

    print(f"Wrote plots to: {args.run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
