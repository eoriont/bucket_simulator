#!/usr/bin/env python3
"""
Analytically sweep distillation protocol / k settings over entanglement rates
to identify RCX regime changes, following the RADar algorithm from the paper.

For each (protocol, k), computes the optimal m*(k) as the maximum m satisfying
both the rate constraint and qubit budget constraint (Appendix B, Eq. 19-25).
Then picks k* minimizing effective channel error (proxy for LER).

All calculations are pure Python — no simulator invocations needed.
Formulas mirror simulator.cpp and Appendix B exactly.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


EXPERIMENT_DIR = Path(__file__).resolve().parent

PROTOCOLS = ("none", "2to1_pumping", "3to1_pumping", "2to1_recurrence", "3to1_recurrence")


# ---------------------------------------------------------------------------
# Distillation math (mirrors simulator.cpp)
# ---------------------------------------------------------------------------

def compute_2to1_success_prob(F1: float, F2: float) -> float:
    term1 = (1.0 - 2.0 * (1.0 - F1) / 3.0) * (1.0 - 2.0 * (1.0 - F2) / 3.0)
    term2 = 4.0 * (1.0 - F1) * (1.0 - F2) / 9.0
    return max(0.01, term1 + term2)


def compute_2to1_fidelity(F1: float, F2: float) -> float:
    p_s = compute_2to1_success_prob(F1, F2)
    numerator = (F1 * F2) + ((1.0 - F1) * (1.0 - F2) / 9.0)
    return min(1.0, max(0.0, numerator / p_s))


def compute_3to1_success_prob(F1: float, F2: float) -> float:
    F1_sq = F1 * F1
    one_minus_F1 = 1.0 - F1
    p_s = (F1_sq * F2 + F1_sq * (1.0 - F2) / 9.0
           + 2.0 * F1 * one_minus_F1 * (F2 / 3.0 + (1.0 - F2) / 9.0)
           + one_minus_F1 * one_minus_F1 * (F2 / 9.0 + (1.0 - F2) / 81.0))
    return max(0.01, p_s)


def compute_3to1_fidelity(F1: float, F2: float) -> float:
    p_s = compute_3to1_success_prob(F1, F2)
    F1_sq = F1 * F1
    one_minus_F1 = 1.0 - F1
    numerator = (F1_sq * F2 + F1_sq * (1.0 - F2) / 27.0
                 + 2.0 * F1 * one_minus_F1 * F2 / 9.0
                 + one_minus_F1 * one_minus_F1 * F2 / 81.0)
    return min(1.0, max(0.0, numerator / p_s))


def branching_factor(protocol: str) -> int:
    return 3 if "3to1" in protocol else 2


@dataclass
class DistillationResult:
    """Per-round distillation output, independent of m."""
    output_fidelity: float       # F_k: fidelity after k rounds
    success_probability: float   # P_success(k): per-batch success prob
    raw_pairs_consumed: int      # C(k): raw pairs consumed per distilled pair (timing)
    qubit_slots: int             # R(k): EPR slots active simultaneously (qubit budget)


def compute_distillation(protocol: str, rounds: int, F_raw: float) -> DistillationResult:
    """
    Compute distillation output for k rounds.

    C(k) (raw_pairs_consumed) is used for the rate/timing constraint (Eq. 19).
    R(k) (qubit_slots) is the simultaneous EPR slot count for the qubit budget (Eq. 21):
      - Pumping:    R(k) = n * k   (all k stages active simultaneously, n wide each)
      - Recurrence: R(k) = n^k    (exponential tree)
    These are per-chain values; multiply by num_remote_cnots for total.
    """
    if protocol == "none" or rounds == 0:
        return DistillationResult(F_raw, 1.0, 1, 0)

    n = branching_factor(protocol)
    is_recurrence = "recurrence" in protocol

    F_current = F_raw
    cumulative_success = 1.0

    for k in range(rounds):
        if "2to1" in protocol:
            F_aux = F_raw if not is_recurrence else F_current
            F_out = compute_2to1_fidelity(F_aux, F_current)
            p_s = compute_2to1_success_prob(F_aux, F_current)
        else:
            F_aux = F_raw if not is_recurrence else F_current
            F_out = compute_3to1_fidelity(F_aux, F_current)
            p_s = compute_3to1_success_prob(F_aux, F_current)

        F_current = max(0.0, min(1.0, F_out))
        cumulative_success *= p_s

    # C(k): raw pairs consumed per distilled pair
    if is_recurrence:
        C_k = n ** rounds          # exponential: n^k
    else:
        C_k = n * rounds + 1      # pumping: n*k rounds of auxiliaries + 1 target
        # Note: Eq. 14 gives C_pump = k(n-1)+1; for n=2: k+1, for n=3: 2k+1
        C_k = rounds * (n - 1) + 1

    # R(k): simultaneous EPR slots per chain (qubit budget, Eq. 13/15/16)
    if is_recurrence:
        R_k = n ** rounds          # exponential tree: n^k
    else:
        R_k = n * rounds           # pumping pipeline: n * k stages simultaneously

    return DistillationResult(F_current, cumulative_success, max(1, C_k), R_k)


# ---------------------------------------------------------------------------
# Noise model
# ---------------------------------------------------------------------------

@dataclass
class Config:
    code_distance: int = 9
    raw_epr_fidelity: float = 0.95
    physical_error: float = 0.005
    measurement_delay: float = 500e-6   # seconds (T_QEC - tau_meas window)
    T1: float = 119.0                   # seconds
    T2: float = 12.6                    # seconds
    accurate_rcx: bool = False
    num_remote_cnots: int = 17          # always set to 2d-1 at construction


def raw_channel_error(cfg: Config) -> float:
    p = cfg.physical_error
    return max(0.0, 1.0 - cfg.raw_epr_fidelity * (1.0 - p) * (1.0 - p))


def rcx_error_from_channel(effective_channel_error: float, cfg: Config) -> float:
    if not cfg.accurate_rcx:
        return max(0.0, effective_channel_error)
    p = cfg.physical_error
    return max(0.0, 1.0 - (1.0 - effective_channel_error) * (1.0 - p) * (1.0 - p))


@dataclass
class CandidateResult:
    protocol: str
    rounds: int
    m_star: int                       # optimal m*(k) per paper Eq. 25
    max_m_rate: int                   # rate constraint upper bound on m
    max_m_qubit: int                  # qubit budget upper bound on m
    distilled_fidelity: float
    distillation_success_prob: float
    raw_pairs_consumed: int           # C(k)
    qubit_slots_per_chain: int        # R(k) per chain
    total_qubit_slots: int            # num_remote_cnots * R(k)
    timing_feasible: bool
    qubit_overhead_feasible: bool
    fully_feasible: bool
    effective_channel_error: float
    remote_cnot_error: float


def evaluate(protocol: str, rounds: int, enr_hz: float, cfg: Config) -> CandidateResult:
    """
    Compute optimal m*(k) for a given (protocol, k) at a given ENR,
    following the joint feasibility selection in Appendix B Eq. 22-25.
    """
    distill = compute_distillation(protocol, rounds, cfg.raw_epr_fidelity)
    raw_err = raw_channel_error(cfg)
    d = cfg.code_distance
    d_sq = d * d

    wants_distillation = protocol != "none" and rounds > 0

    if not wants_distillation:
        return CandidateResult(
            protocol=protocol, rounds=rounds, m_star=1,
            max_m_rate=999999, max_m_qubit=999999,
            distilled_fidelity=distill.output_fidelity,
            distillation_success_prob=1.0,
            raw_pairs_consumed=1, qubit_slots_per_chain=0, total_qubit_slots=0,
            timing_feasible=True, qubit_overhead_feasible=True, fully_feasible=True,
            effective_channel_error=raw_err,
            remote_cnot_error=rcx_error_from_channel(raw_err, cfg),
        )

    N = cfg.num_remote_cnots
    C_k = distill.raw_pairs_consumed
    R_k = distill.qubit_slots              # per chain
    P_s = distill.success_probability

    # Rate constraint (Eq. 19-20): m <= r * P_s * T_meas / (N * C(k))
    if enr_hz > 0 and C_k > 0 and N > 0:
        max_m_rate = int(enr_hz * P_s * cfg.measurement_delay / (N * C_k))
    else:
        max_m_rate = 0

    # Qubit budget: distillation qubits (2 * N * R(k) * m) must fit within the
    # monolithic equivalent qubit count = 4d² + d + 2 (from lattice surgery circuit)
    qubit_budget = 4 * d * d + d + 2
    total_r = N * R_k
    if total_r > 0:
        max_m_qubit = (qubit_budget // 2) // total_r   # factor 2: each EPR slot = 2 qubits
    else:
        max_m_qubit = 999999

    # Optimal m*(k) = max feasible m (Eq. 25)
    m_star = min(max_m_rate, max_m_qubit)

    timing_ok = max_m_rate >= 1
    qubit_ok = max_m_qubit >= 1
    fully_ok = m_star >= 1

    if not fully_ok:
        # Return infeasible result with best-effort m=1 for display
        m_star = 1

    # Effective channel error with m_star backup batches
    distilled_err = 1.0 - distill.output_fidelity
    p_all_fail = (1.0 - P_s) ** m_star
    eff_err = (1.0 - p_all_fail) * distilled_err + p_all_fail * raw_err

    rcx_err = rcx_error_from_channel(eff_err, cfg)

    return CandidateResult(
        protocol=protocol,
        rounds=rounds,
        m_star=m_star,
        max_m_rate=max_m_rate,
        max_m_qubit=max_m_qubit,
        distilled_fidelity=distill.output_fidelity,
        distillation_success_prob=P_s,
        raw_pairs_consumed=C_k,
        qubit_slots_per_chain=R_k,
        total_qubit_slots=total_r * m_star if fully_ok else total_r,
        timing_feasible=timing_ok,
        qubit_overhead_feasible=qubit_ok,
        fully_feasible=fully_ok,
        effective_channel_error=eff_err,
        remote_cnot_error=rcx_err,
    )


def protocol_candidates(max_rounds: int) -> Iterable[tuple[str, int]]:
    yield "none", 0
    for protocol in PROTOCOLS[1:]:
        for rounds in range(1, max_rounds + 1):
            yield protocol, rounds


def choose_best(candidates: list[CandidateResult]) -> CandidateResult:
    """Pick k* = argmin effective channel error among fully feasible candidates (Eq. 26).
    Falls back to infeasible candidates if none are feasible."""
    feasible = [c for c in candidates if c.fully_feasible]
    pool = feasible if feasible else candidates
    return min(pool, key=lambda c: (c.effective_channel_error, c.rounds, c.protocol))


def regime_signature(c: CandidateResult) -> tuple:
    return (c.protocol, c.rounds, c.m_star, c.fully_feasible)


# ---------------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------------

def rate_grid(rate_min: int, rate_max: int, points_per_decade: int) -> list[int]:
    if rate_min == rate_max:
        return [rate_min]
    log_min = math.log10(rate_min)
    log_max = math.log10(rate_max)
    num = max(2, int(round((log_max - log_min) * points_per_decade)) + 1)
    vals = sorted({int(round(10 ** (log_min + (log_max - log_min) * i / (num - 1))))
                   for i in range(num)})
    vals[0] = rate_min
    vals[-1] = rate_max
    return vals


# ---------------------------------------------------------------------------
# CSV / output helpers
# ---------------------------------------------------------------------------

ALL_FIELDS = [
    "entanglement_rate_hz",
    "protocol",
    "rounds",
    "m_star",
    "max_m_rate",
    "max_m_qubit",
    "distilled_fidelity",
    "distillation_success_prob",
    "raw_pairs_consumed",
    "qubit_slots_per_chain",
    "total_qubit_slots",
    "timing_feasible",
    "qubit_overhead_feasible",
    "fully_feasible",
    "effective_channel_error",
    "remote_cnot_error",
]

BEST_FIELDS = [
    "entanglement_rate_hz",
    "protocol",
    "rounds",
    "m_star",
    "max_m_rate",
    "max_m_qubit",
    "remote_cnot_error",
    "effective_channel_error",
    "fully_feasible",
    "timing_feasible",
    "qubit_overhead_feasible",
    "total_qubit_slots",
]


def candidate_to_row(enr: int, c: CandidateResult) -> dict:
    return {
        "entanglement_rate_hz": enr,
        "protocol": c.protocol,
        "rounds": c.rounds,
        "m_star": c.m_star,
        "max_m_rate": c.max_m_rate,
        "max_m_qubit": c.max_m_qubit,
        "distilled_fidelity": c.distilled_fidelity,
        "distillation_success_prob": c.distillation_success_prob,
        "raw_pairs_consumed": c.raw_pairs_consumed,
        "qubit_slots_per_chain": c.qubit_slots_per_chain,
        "total_qubit_slots": c.total_qubit_slots,
        "timing_feasible": c.timing_feasible,
        "qubit_overhead_feasible": c.qubit_overhead_feasible,
        "fully_feasible": c.fully_feasible,
        "effective_channel_error": c.effective_channel_error,
        "remote_cnot_error": c.remote_cnot_error,
    }


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


# ---------------------------------------------------------------------------
# Config file parsing
# ---------------------------------------------------------------------------

DEFAULT_CONFIG = EXPERIMENT_DIR / "config_neutral_atom.txt"


def parse_config_file(path: Path) -> dict:
    """Parse key-value pairs from a bucket_simulator config file."""
    result: dict = {}
    for line in path.read_text().splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            result[parts[0]] = parts[1].strip()
    return result


def get_float(cfg: dict, key: str, default: float) -> float:
    return float(cfg[key]) if key in cfg else default


def get_int(cfg: dict, key: str, default: int) -> int:
    return int(cfg[key]) if key in cfg else default


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Analytically sweep distillation params over ENR (RADar algorithm).")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                        help="Config file supplying physical parameters")
    parser.add_argument("--output-dir", type=Path,
                        default=EXPERIMENT_DIR / "runs" / "manual")
    parser.add_argument("--rate-min-hz", type=int, default=1_000)
    parser.add_argument("--rate-max-hz", type=int, default=100_000_000_000)
    parser.add_argument("--points-per-decade", type=int, default=12)
    parser.add_argument("--max-rounds", type=int, default=5)
    parser.add_argument("--code-distance", type=int, default=None,
                        help="Code distance d (overrides config)")
    parser.add_argument("--physical-error", type=float, default=None,
                        help="Physical gate error (overrides config)")
    parser.add_argument("--accurate-rcx", action="store_true")
    args = parser.parse_args()

    file_cfg = parse_config_file(args.config)

    d = args.code_distance if args.code_distance is not None else get_int(file_cfg, "code_distance", 9)
    phys_err = args.physical_error if args.physical_error is not None else get_float(file_cfg, "physical_error", 0.005)

    cfg = Config(
        code_distance=d,
        raw_epr_fidelity=get_float(file_cfg, "raw_epr_fidelity", 0.95),
        physical_error=phys_err,
        measurement_delay=get_float(file_cfg, "measurement_delay", 500e-6),
        T1=get_float(file_cfg, "T1", 119.0),
        T2=get_float(file_cfg, "T2", 12.6),
        accurate_rcx=args.accurate_rcx,
        num_remote_cnots=2 * d - 1,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rates = rate_grid(args.rate_min_hz, args.rate_max_hz, args.points_per_decade)

    all_rows: list[dict] = []
    best_rows: list[dict] = []
    regime_rows: list[dict] = []
    previous_signature = None

    for enr in rates:
        candidates = [
            evaluate(protocol, rounds, enr, cfg)
            for protocol, rounds in protocol_candidates(args.max_rounds)
        ]
        all_rows.extend(candidate_to_row(enr, c) for c in candidates)

        best = choose_best(candidates)
        best_rows.append(candidate_to_row(enr, best))

        sig = regime_signature(best)
        if sig != previous_signature:
            regime_rows.append(candidate_to_row(enr, best))
            previous_signature = sig

    write_csv(args.output_dir / "all_candidates.csv", all_rows, ALL_FIELDS)
    write_csv(args.output_dir / "best_by_rate.csv", best_rows, BEST_FIELDS)
    write_csv(args.output_dir / "regime_boundaries.csv", regime_rows, BEST_FIELDS)

    summary_path = args.output_dir / "summary.txt"
    with summary_path.open("w") as f:
        f.write("RCX Regime Sweep Summary (RADar algorithm)\n")
        f.write("==========================================\n")
        f.write(f"code_distance:      {cfg.code_distance}  (d²={cfg.code_distance**2} qubit budget)\n")
        f.write(f"raw_epr_fidelity:   {cfg.raw_epr_fidelity}\n")
        f.write(f"physical_error:     {cfg.physical_error}\n")
        f.write(f"measurement_delay:  {cfg.measurement_delay*1e6:.1f} us\n")
        f.write(f"T1:                 {cfg.T1} s\n")
        f.write(f"T2:                 {cfg.T2} s\n")
        f.write(f"num_remote_cnots:   {cfg.num_remote_cnots}\n")
        f.write(f"rates_hz:           {rates}\n")
        f.write(f"max_rounds:         {args.max_rounds}\n\n")
        f.write("Regime boundaries:\n")
        for row in regime_rows:
            f.write(
                f"  enr={row['entanglement_rate_hz']}Hz  protocol={row['protocol']}  "
                f"k={row['rounds']}  m*={row['m_star']}  "
                f"(max_m_rate={row['max_m_rate']} max_m_qubit={row['max_m_qubit']})  "
                f"eff_err={row['effective_channel_error']:.6f}  "
                f"feasible={row['fully_feasible']}\n"
            )

    print(f"Wrote {len(all_rows)} candidates, {len(best_rows)} best-by-rate, "
          f"{len(regime_rows)} regime boundaries to: {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
