# RADar Paper Reproduction Experiment

This experiment reproduces the key results from:

> **RADar: Resource-Aware Entanglement-Distillation for Distributed Surface Code**
> ISCA 2026 Submission #24

## Overview

The experiment consists of four sub-experiments:

| Sub-experiment | Description | Configs | Paper Reference |
|----------------|-------------|---------|-----------------|
| `monolithic_baseline` | Monolithic surface codes (no distribution) | 3 | Figure 4 baseline |
| `table1_low_rate` | Low entanglement rate regime (80-100 MHz) | 18 | Table I |
| `experiment` | High entanglement rate regime (150-500 MHz) | 105 | Table II |
| `rate_sweep` | Rate sweep showing timing threshold | 9 | Figure 6 |

**Total: 135 configurations**

## Parameters (from Paper Section III.B)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Physical error rate (p) | 10⁻³ | Gate/measurement error |
| T1 coherence time | 125 μs | Relaxation time |
| T2 coherence time | 200 μs | Dephasing time |
| Measurement delay (t_meas) | 660 ns | QEC cycle timing bound |
| Raw EPR fidelity | 0.99 | Initial entanglement quality |

## Code Distances

- **d = 13**: 25 remote CNOTs per cycle
- **d = 15**: 29 remote CNOTs per cycle
- **d = 17**: 33 remote CNOTs per cycle

## Distillation Protocols

1. **None** - No distillation (baseline distributed)
2. **2→1 Pumping** - Linear resource cost O(2k)
3. **3→1 Pumping** - Linear resource cost O(3k)
4. **2→1 Recurrence** - Exponential resource cost O(2^k)
5. **3→1 Recurrence** - Exponential resource cost O(3^k)

## Usage

### Generate configs only
```bash
./run.sh --generate
```

### Run full experiment
```bash
./run.sh -n 8           # 8 MPI processes
```

### Run specific sub-experiment
```bash
./run.sh --sub monolithic_baseline -n 4
./run.sh --sub table1_low_rate -n 4
```

### List available sub-experiments
```bash
./run.sh --list
```

## Expected Results

### Table I (Low Rate Regime)
At 80-100 MHz, only shallow distillation (1 round of 2→1) is timing-safe.
Expected LER gap reduction: ~30-60%

### Table II (High Rate Regime)
At 150-500 MHz, deeper distillation becomes feasible.
- 2→1 Recurrence provides best results at medium rates
- 3→1 Recurrence optimal at high rates (≥400 MHz)
- Expected LER gap reduction: up to 98%

## Output Structure

```
runs/
└── YYYYMMDD_HHMMSS/
    ├── experiment_info.txt     # Experiment metadata
    ├── results/
    │   ├── monolithic_baseline/
    │   ├── table1_low_rate/
    │   ├── experiment/
    │   └── rate_sweep/
    ├── logs/
    │   └── (same structure)
    └── results_summary.csv
```

## Analysis

After running, analyze results with:
```bash
python3 ../scripts/analyze_results.py runs/<timestamp> --csv --plot
```

## Estimated Runtime

With 10M shots per config and 8 MPI processes:
- Monolithic baseline: ~5-10 minutes
- Table I (18 configs): ~1-2 hours
- Table II (105 configs): ~6-12 hours
- Rate sweep (9 configs): ~30-60 minutes

**Total: ~8-15 hours** (varies significantly with hardware)

For faster iteration, reduce `total_shots` in the JSON specs.
