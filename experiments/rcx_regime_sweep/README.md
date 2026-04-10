# RCX Regime Sweep

Analytically finds optimal distillation protocol (type, rounds k, redundancy m\*) as a function of entanglement rate (ENR), following the RADar algorithm. No simulator required.

## Usage

```bash
# Sweep and plot (reads physical params from config_neutral_atom.txt)
python3 sweep_optimal_rcx.py --output-dir runs/my_run
python3 plot_regimes.py runs/my_run

# Different code distance
python3 sweep_optimal_rcx.py --code-distance 13 --output-dir runs/d13

# Custom config
python3 sweep_optimal_rcx.py --config my_config.txt --output-dir runs/my_run
```

## Parameters

Physical parameters (`raw_epr_fidelity`, `measurement_delay`, `T1`, `T2`, `physical_error`, `code_distance`) are read from the config file. `num_remote_cnots` is computed automatically as `2d - 1`.

| CLI Flag | Default | Description |
|---|---|---|
| `--config` | `config_neutral_atom.txt` | Physical parameter source |
| `--code-distance` | from config | Overrides config code distance |
| `--rate-min-hz` | 1,000 | Sweep start ENR |
| `--rate-max-hz` | 100,000,000,000 | Sweep end ENR |
| `--max-rounds` | 5 | Max distillation rounds k |

## Outputs

Each run directory contains:

- `summary.txt` — regime boundaries in plain text
- `best_by_rate.csv` — optimal (protocol, k, m\*) at each ENR
- `all_candidates.csv` — all evaluated (protocol, k) combinations
- `regime_overview.png` — 3-strip plot: effective channel error, k, m\* vs ENR
- `racing_lines.png` — full candidate landscape with optimal path highlighted
- `qubit_budget.png` — heatmap of qubit usage vs (k, m) per protocol

## Qubit budget

Feasibility threshold = `4d² + d + 2` (monolithic equivalent qubit count). Rate constraint: `m ≤ ENR · P_success · T_meas / (N · C(k))`. Optimal m\* = max satisfying both.
