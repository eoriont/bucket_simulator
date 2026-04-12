#!/bin/bash
#
# Neutral-atom d=7 lattice-surgery LER sweep.
# Covers:
#   - monolithic baseline (ENR-independent, once per mode)
#   - distributed: all feasible (protocol, k, m) combinations at each ENR point
#
# Feasible combos are computed at runtime via compute_combos.py using the RADar
# math from rcx_regime_sweep/sweep_optimal_rcx.py — parameters come from the
# distributed config file, so changing the config automatically changes the sweep.
#
# Usage:
#   ./run.sh [options]
#
# Options:
#   -n, --num-procs NUM    MPI processes per sim (default: 4)
#   -j, --parallel NUM     Simulations in parallel (default: 1)
#   --dry-run              Print commands without executing
#   --circuits-only        Dump circuits instead of running simulations
#   --resume RUN_DIR       Resume an existing run directory
#   --shots SHOTS          Override total_shots, e.g. 10K

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SIMULATOR="$PROJECT_ROOT/build/bucket_simulator"
MONOLITHIC_CONFIG="$SCRIPT_DIR/config_monolithic.txt"
DISTRIBUTED_CONFIG="$SCRIPT_DIR/config_distributed.txt"

NUM_PROCS=4
PARALLEL_JOBS=1
DRY_RUN=false
CIRCUITS_ONLY=false
SHOTS=""
RESUME_RUN_DIR=""
DISTANCE=7
MAX_ROUNDS=5

# ENR sweep points (Hz), bracketing all regime transitions for d=7
# at F=0.99, p=0.001:
#   ~56kHz  (none → 2to1_pumping k=1,       eff 1.20% → 0.68%)
#   ~83kHz  (→ 3to1_pumping k=1,            eff 0.68% → 0.54%)
#  ~121kHz  (→ 2to1_recurrence k=2,         eff 0.54% → 0.47%)
#
# Use extra density around the three transition windows and only a few
# plateau points afterward.
RATES_HZ=(20000 50000 70000 80000 90000 110000 130000 180000 300000 500000)
MODES=(xx_merge zz_merge)

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)   NUM_PROCS="$2"; shift 2 ;;
        -j|--parallel)    PARALLEL_JOBS="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=true; shift ;;
        --circuits-only)  CIRCUITS_ONLY=true; shift ;;
        --resume)         RESUME_RUN_DIR="$2"; shift 2 ;;
        --shots)          SHOTS="$2"; shift 2 ;;
        -h|--help)        head -27 "$0" | tail -25; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------- helpers --------------------------------------------------------

enr_combos() {
    python3 "$SCRIPT_DIR/compute_combos.py" \
        --config "$DISTRIBUTED_CONFIG" \
        --enr "$1" \
        --max-rounds "$MAX_ROUNDS"
}

mode_merge_type() {
    local arch="$1"
    local mode="$2"
    case "$arch:$mode" in
        monolithic:xx_merge|monolithic:xx_split) echo "xx" ;;
        monolithic:zz_merge|monolithic:zz_split) echo "zz" ;;
        distributed:xx_merge|distributed:xx_split) echo "distributed_xx" ;;
        distributed:zz_merge|distributed:zz_split) echo "distributed_zz" ;;
    esac
}

mode_phase() {
    case "$1" in
        xx_merge|zz_merge) echo "merge_only" ;;
        xx_split|zz_split) echo "split_only" ;;
    esac
}

rate_label() {
    local rate_hz="$1"
    if [[ "$rate_hz" -ge 1000 ]]; then
        echo "$((rate_hz / 1000))kHz"
    else
        echo "${rate_hz}Hz"
    fi
}

job_label() {
    local arch="$1"
    local mode="$2"
    local rate_hz="$3"
    local protocol="$4"
    local rounds="$5"
    local batches="$6"
    if [[ "$arch" == "monolithic" ]]; then
        echo "${arch}_${mode}_d${DISTANCE}"
    elif [[ "$protocol" == "none" ]]; then
        echo "${arch}_${mode}_none_$(rate_label "$rate_hz")_d${DISTANCE}"
    else
        echo "${arch}_${mode}_${protocol}_k${rounds}_m${batches}_$(rate_label "$rate_hz")_d${DISTANCE}"
    fi
}

job_is_done() {
    local label="$1"
    if [[ "$CIRCUITS_ONLY" == "true" ]]; then
        [[ -f "$CIRCUITS_DIR/${label}.stim" ]]
    else
        [[ -f "$RESULTS_DIR/${label}_results.txt" ]]
    fi
}

# ---------- pre-flight checks ----------------------------------------------

if [[ ! -x "$SIMULATOR" ]]; then
    echo "Error: simulator not found at $SIMULATOR"
    echo "Build first: cd $PROJECT_ROOT/build && make"
    exit 1
fi

if [[ ! -f "$MONOLITHIC_CONFIG" || ! -f "$DISTRIBUTED_CONFIG" ]]; then
    echo "Error: expected base configs not found in $SCRIPT_DIR"
    exit 1
fi

# ---------- run directory --------------------------------------------------

if [[ -n "$RESUME_RUN_DIR" ]]; then
    RUN_DIR="$RESUME_RUN_DIR"
    TIMESTAMP=$(basename "$RUN_DIR")
else
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    RUN_DIR="$SCRIPT_DIR/runs/$TIMESTAMP"
fi

RESULTS_DIR="$RUN_DIR/results"
CIRCUITS_DIR="$RUN_DIR/circuits"
LOGS_DIR="$RUN_DIR/logs"
mkdir -p "$RESULTS_DIR" "$CIRCUITS_DIR" "$LOGS_DIR"

TOTAL_MONO=${#MODES[@]}
TOTAL_DIST=0
for rate_hz in "${RATES_HZ[@]}"; do
    combos_str=$(enr_combos "$rate_hz")
    n_combos=$(echo "$combos_str" | wc -w)
    TOTAL_DIST=$(( TOTAL_DIST + n_combos * ${#MODES[@]} ))
done
TOTAL_CONFIGS=$(( TOTAL_MONO + TOTAL_DIST ))

FAIL_FILE="$RUN_DIR/failed.txt"
: > "$FAIL_FILE"

echo "=============================================="
echo "Neutral-Atom d=7 Lattice Surgery LER Sweep"
echo "=============================================="
echo "Run ID:    $TIMESTAMP"
echo "Output:    $RUN_DIR"
echo "Distance:  $DISTANCE"
echo "Modes:     ${MODES[*]}"
echo "Rates:     ${RATES_HZ[*]}"
echo "Shots:     ${SHOTS:-from config}"
echo "MaxRounds: $MAX_ROUNDS"
echo "Circuits:  $CIRCUITS_ONLY"
echo "MPI:       $NUM_PROCS processes"
echo "Parallel:  $PARALLEL_JOBS"
echo "Jobs:      $TOTAL_CONFIGS  (${TOTAL_MONO} monolithic + ${TOTAL_DIST} distributed)"
[[ -n "$RESUME_RUN_DIR" ]] && echo "Resume:    yes"
echo "=============================================="

cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: neutral_atom_d7_lattice_surgery
Run ID: $TIMESTAMP
Started: $(date)
Monolithic config: $MONOLITHIC_CONFIG
Distributed config: $DISTRIBUTED_CONFIG
Distance: $DISTANCE
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Shots override: ${SHOTS:-none}
Max distillation rounds: $MAX_ROUNDS
Circuits only: $CIRCUITS_ONLY
Modes: ${MODES[*]}
Distributed entanglement rates (Hz): ${RATES_HZ[*]}
Total jobs: $TOTAL_CONFIGS
Host: $(hostname)
Simulator: $SIMULATOR
EOF

# ---------- run_job --------------------------------------------------------

run_job() {
    local idx="$1"
    local total="$2"
    local arch="$3"
    local mode="$4"
    local rate_hz="$5"
    local protocol="${6:-}"
    local rounds="${7:-}"
    local batches="${8:-}"
    local config_file
    local label
    local merge_type
    local phase
    local result_file
    local log_file
    local circuit_file
    local dump_dir
    local extra_args=()

    if [[ "$arch" == "monolithic" ]]; then
        config_file="$MONOLITHIC_CONFIG"
    else
        config_file="$DISTRIBUTED_CONFIG"
    fi

    label=$(job_label "$arch" "$mode" "$rate_hz" "$protocol" "$rounds" "$batches")
    merge_type=$(mode_merge_type "$arch" "$mode")
    phase=$(mode_phase "$mode")
    result_file="$RESULTS_DIR/${label}_results.txt"
    log_file="$LOGS_DIR/${label}.log"
    circuit_file="$CIRCUITS_DIR/${label}.stim"
    dump_dir="$CIRCUITS_DIR/.tmp_${label}"

    [[ -n "$SHOTS" ]] && extra_args+=(-total_shots "$SHOTS")
    if [[ "$arch" == "distributed" ]]; then
        extra_args+=(-entanglement_rate "$rate_hz")
        extra_args+=(-distillation_protocol "$protocol")
        extra_args+=(-distillation_rounds "$rounds")
        extra_args+=(-distillation_backup_batches "$batches")
    fi

    if job_is_done "$label"; then
        echo "[SKIP $idx/$total] $label"
        return 0
    fi

    if [[ "$arch" == "monolithic" ]]; then
        echo "[RUN $idx/$total] $label  (merge_type=$merge_type phase=$phase)"
    else
        echo "[RUN $idx/$total] $label  (merge_type=$merge_type phase=$phase enr=$(rate_label "$rate_hz"))"
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  mpirun -n $NUM_PROCS $SIMULATOR \\"
        echo "    -config $config_file \\"
        echo "    -merge_type $merge_type \\"
        echo "    -experiment_phase $phase \\"
        if [[ "$arch" == "distributed" ]]; then
            echo "    -entanglement_rate $rate_hz \\"
            echo "    -distillation_protocol $protocol \\"
            echo "    -distillation_rounds $rounds \\"
            echo "    -distillation_backup_batches $batches \\"
        fi
        [[ -n "$SHOTS" ]] && echo "    -total_shots $SHOTS \\"
        if [[ "$CIRCUITS_ONLY" == "true" ]]; then
            echo "    -dump-circuit \\"
            echo "    -output $CIRCUITS_DIR"
        else
            echo "    -output $RESULTS_DIR"
        fi
        return 0
    fi

    rm -rf "$dump_dir"
    mkdir -p "$dump_dir"

    mpirun -n 1 "$SIMULATOR" \
        -config "$config_file" \
        -merge_type "$merge_type" \
        -experiment_phase "$phase" \
        "${extra_args[@]}" \
        -dump-circuit \
        -output "$dump_dir" >> "$log_file" 2>&1

    local dumped="$dump_dir/$(basename "$config_file" .txt).stim"
    if [[ -f "$dumped" ]]; then
        mv "$dumped" "$circuit_file"
    fi
    rmdir "$dump_dir" 2>/dev/null || true

    if [[ "$CIRCUITS_ONLY" == "true" ]]; then
        return 0
    fi

    # Use a per-job temp output dir so the simulator's timestamped filename
    # is unambiguous even when jobs run in parallel.
    local tmp_out="$RUN_DIR/.tmp_out_${label}"
    rm -rf "$tmp_out"
    mkdir -p "$tmp_out"

    mpirun -n "$NUM_PROCS" "$SIMULATOR" \
        -config "$config_file" \
        -merge_type "$merge_type" \
        -experiment_phase "$phase" \
        "${extra_args[@]}" \
        -output "$tmp_out" >> "$log_file" 2>&1

    local outfile
    outfile=$(ls "$tmp_out"/*.txt 2>/dev/null | head -1)
    if [[ -n "$outfile" ]]; then
        mv "$outfile" "$result_file"
    else
        echo "$label" >> "$FAIL_FILE"
        echo "  Warning: no result file found for $label" >> "$log_file"
    fi
    rmdir "$tmp_out" 2>/dev/null || true
}

# ---------- main loop ------------------------------------------------------

TASK_INDEX=0

for mode in "${MODES[@]}"; do
    TASK_INDEX=$((TASK_INDEX + 1))
    run_job "$TASK_INDEX" "$TOTAL_CONFIGS" monolithic "$mode" 0
done

for rate_hz in "${RATES_HZ[@]}"; do
    combos_str=$(enr_combos "$rate_hz")
    for combo in $combos_str; do
        IFS=':' read -r protocol rounds batches <<< "$combo"
        for mode in "${MODES[@]}"; do
            TASK_INDEX=$((TASK_INDEX + 1))
            run_job "$TASK_INDEX" "$TOTAL_CONFIGS" distributed "$mode" "$rate_hz" "$protocol" "$rounds" "$batches"
        done
    done
done

echo ""
echo "Done. Results in $RUN_DIR"
