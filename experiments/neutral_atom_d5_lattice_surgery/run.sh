#!/bin/bash
#
# Neutral-atom d=5 lattice-surgery sweep.
# Covers:
#   - monolithic baseline
#   - distributed lattice surgery
# Across:
#   - xx_merge
#   - zz_merge
#   - xx_split
#   - zz_split
#
# The distributed sweep spans entanglement generation rates from 1-100 kHz.
# The monolithic baseline is ENR-independent, so it is run once per mode.
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
#   --protocol NAME        Override distillation protocol for distributed runs
#   --dist-rounds K        Override distillation rounds for distributed runs
#   --backup-batches M     Override distillation backup batches for distributed runs

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
DISTILLATION_PROTOCOL=""
DISTILLATION_ROUNDS=""
DISTILLATION_BACKUP_BATCHES=""
DISTANCE=5

# 1-100 kHz sweep, spaced to cover the low-rate regime without exploding run count.
RATES_HZ=(1000 2000 5000 10000 20000 50000 100000)
MODES=(xx_merge zz_merge xx_split zz_split)
ARCHS=(monolithic distributed)

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)   NUM_PROCS="$2"; shift 2 ;;
        -j|--parallel)    PARALLEL_JOBS="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=true; shift ;;
        --circuits-only)  CIRCUITS_ONLY=true; shift ;;
        --resume)         RESUME_RUN_DIR="$2"; shift 2 ;;
        --shots)          SHOTS="$2"; shift 2 ;;
        --protocol)       DISTILLATION_PROTOCOL="$2"; shift 2 ;;
        --dist-rounds)    DISTILLATION_ROUNDS="$2"; shift 2 ;;
        --backup-batches) DISTILLATION_BACKUP_BATCHES="$2"; shift 2 ;;
        -h|--help)        head -33 "$0" | tail -31; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

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
    if [[ "$arch" == "monolithic" ]]; then
        echo "${arch}_${mode}_d${DISTANCE}"
    else
        echo "${arch}_${mode}_$(rate_label "$rate_hz")_d${DISTANCE}"
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

if [[ ! -x "$SIMULATOR" ]]; then
    echo "Error: simulator not found at $SIMULATOR"
    echo "Build first: cd $PROJECT_ROOT/build && make"
    exit 1
fi

if [[ ! -f "$MONOLITHIC_CONFIG" || ! -f "$DISTRIBUTED_CONFIG" ]]; then
    echo "Error: expected base configs not found in $SCRIPT_DIR"
    exit 1
fi

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

TOTAL_CONFIGS=$(( ${#MODES[@]} + (${#MODES[@]} * ${#RATES_HZ[@]}) ))
FAIL_FILE="$RUN_DIR/failed.txt"
: > "$FAIL_FILE"

echo "=============================================="
echo "Neutral-Atom d=5 Lattice Surgery Sweep"
echo "=============================================="
echo "Run ID:    $TIMESTAMP"
echo "Output:    $RUN_DIR"
echo "Distance:  $DISTANCE"
echo "Modes:     ${MODES[*]}"
echo "Rates:     ${RATES_HZ[*]}"
echo "Shots:     ${SHOTS:-from config}"
echo "Protocol:  ${DISTILLATION_PROTOCOL:-from distributed config}"
echo "Rounds k:  ${DISTILLATION_ROUNDS:-from distributed config}"
echo "Batches m: ${DISTILLATION_BACKUP_BATCHES:-from distributed config}"
echo "Circuits:  $CIRCUITS_ONLY"
echo "MPI:       $NUM_PROCS processes"
echo "Parallel:  $PARALLEL_JOBS"
echo "Configs:   $TOTAL_CONFIGS"
[[ -n "$RESUME_RUN_DIR" ]] && echo "Resume:    yes"
echo "=============================================="

cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: neutral_atom_d5_lattice_surgery
Run ID: $TIMESTAMP
Started: $(date)
Monolithic config: $MONOLITHIC_CONFIG
Distributed config: $DISTRIBUTED_CONFIG
Distance: $DISTANCE
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Shots override: ${SHOTS:-none}
Distillation protocol override: ${DISTILLATION_PROTOCOL:-none}
Distillation rounds override: ${DISTILLATION_ROUNDS:-none}
Distillation backup batches override: ${DISTILLATION_BACKUP_BATCHES:-none}
Circuits only: $CIRCUITS_ONLY
Modes: ${MODES[*]}
Distributed entanglement rates (Hz): ${RATES_HZ[*]}
Host: $(hostname)
Simulator: $SIMULATOR
EOF

run_job() {
    local idx="$1"
    local total="$2"
    local arch="$3"
    local mode="$4"
    local rate_hz="$5"
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

    label=$(job_label "$arch" "$mode" "$rate_hz")
    merge_type=$(mode_merge_type "$arch" "$mode")
    phase=$(mode_phase "$mode")
    result_file="$RESULTS_DIR/${label}_results.txt"
    log_file="$LOGS_DIR/${label}.log"
    circuit_file="$CIRCUITS_DIR/${label}.stim"
    dump_dir="$CIRCUITS_DIR/.tmp_${label}"

    [[ -n "$SHOTS" ]] && extra_args+=(-total_shots "$SHOTS")
    if [[ "$arch" == "distributed" ]]; then
        extra_args+=(-entanglement_rate "$rate_hz")
        [[ -n "$DISTILLATION_PROTOCOL" ]] && extra_args+=(-distillation_protocol "$DISTILLATION_PROTOCOL")
        [[ -n "$DISTILLATION_ROUNDS" ]] && extra_args+=(-distillation_rounds "$DISTILLATION_ROUNDS")
        [[ -n "$DISTILLATION_BACKUP_BATCHES" ]] && extra_args+=(-distillation_backup_batches "$DISTILLATION_BACKUP_BATCHES")
    fi

    if job_is_done "$label"; then
        echo "[SKIP $idx/$total] $label (already complete)"
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
            [[ -n "$DISTILLATION_PROTOCOL" ]] && echo "    -distillation_protocol $DISTILLATION_PROTOCOL \\"
            [[ -n "$DISTILLATION_ROUNDS" ]] && echo "    -distillation_rounds $DISTILLATION_ROUNDS \\"
            [[ -n "$DISTILLATION_BACKUP_BATCHES" ]] && echo "    -distillation_backup_batches $DISTILLATION_BACKUP_BATCHES \\"
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

    mpirun -n "$NUM_PROCS" "$SIMULATOR" \
        -config "$config_file" \
        -merge_type "$merge_type" \
        -experiment_phase "$phase" \
        "${extra_args[@]}" \
        -output "$RESULTS_DIR" >> "$log_file" 2>&1

    local latest
    latest=$(ls -t "$RESULTS_DIR"/results_*.txt 2>/dev/null | head -1)
    if [[ -n "$latest" ]]; then
        mv "$latest" "$result_file"
    else
        echo "$label" >> "$FAIL_FILE"
        echo "  Warning: no result file found for $label" >> "$log_file"
    fi
}

TASK_INDEX=0
TOTAL_JOBS=$TOTAL_CONFIGS

for mode in "${MODES[@]}"; do
    TASK_INDEX=$((TASK_INDEX + 1))
    run_job "$TASK_INDEX" "$TOTAL_JOBS" monolithic "$mode" 0
done

for mode in "${MODES[@]}"; do
    for rate_hz in "${RATES_HZ[@]}"; do
        TASK_INDEX=$((TASK_INDEX + 1))
        run_job "$TASK_INDEX" "$TOTAL_JOBS" distributed "$mode" "$rate_hz"
    done
done

echo "Prepared experiment in $RUN_DIR"
