#!/bin/bash
#
# Neutral-atom d=9 and d=11 lattice-surgery LER sweep.
# Covers:
#   - monolithic baseline (ENR-independent, once per mode and distance)
#   - distributed: all feasible (protocol, k, m) combinations at each ENR point
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

NUM_PROCS=4
PARALLEL_JOBS=1
DRY_RUN=false
CIRCUITS_ONLY=false
SHOTS=""
RESUME_RUN_DIR=""
MAX_ROUNDS=5

DISTANCES=(9 11)
MODES=(xx_merge zz_merge)

# ENR sweep points (Hz), bracketing all regime transitions for d=9 and d=11
# at F=0.99, p=0.001.
RATES_HZ=(30000 80000 120000 160000 250000 400000 600000 1000000)

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)   NUM_PROCS="$2"; shift 2 ;;
        -j|--parallel)    PARALLEL_JOBS="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=true; shift ;;
        --circuits-only)  CIRCUITS_ONLY=true; shift ;;
        --resume)         RESUME_RUN_DIR="$2"; shift 2 ;;
        --shots)          SHOTS="$2"; shift 2 ;;
        -h|--help)        head -25 "$0" | tail -23; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

config_path() {
    local arch="$1"
    local distance="$2"
    echo "$SCRIPT_DIR/config_${arch}_d${distance}.txt"
}

enr_combos() {
    local distance="$1"
    local rate_hz="$2"
    python3 "$SCRIPT_DIR/compute_combos.py" \
        --config "$(config_path distributed "$distance")" \
        --enr "$rate_hz" \
        --max-rounds "$MAX_ROUNDS"
}

mode_merge_type() {
    local arch="$1"
    local mode="$2"
    case "$arch:$mode" in
        monolithic:xx_merge) echo "xx" ;;
        monolithic:zz_merge) echo "zz" ;;
        distributed:xx_merge) echo "distributed_xx" ;;
        distributed:zz_merge) echo "distributed_zz" ;;
    esac
}

mode_phase() {
    case "$1" in
        xx_merge|zz_merge) echo "merge_only" ;;
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
    local distance="$3"
    local rate_hz="$4"
    local protocol="$5"
    local rounds="$6"
    local batches="$7"
    if [[ "$arch" == "monolithic" ]]; then
        echo "${arch}_${mode}_d${distance}"
    elif [[ "$protocol" == "none" ]]; then
        echo "${arch}_${mode}_none_$(rate_label "$rate_hz")_d${distance}"
    else
        echo "${arch}_${mode}_${protocol}_k${rounds}_m${batches}_$(rate_label "$rate_hz")_d${distance}"
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

for distance in "${DISTANCES[@]}"; do
    if [[ ! -f "$(config_path monolithic "$distance")" || ! -f "$(config_path distributed "$distance")" ]]; then
        echo "Error: expected base configs not found for d=$distance in $SCRIPT_DIR"
        exit 1
    fi
done

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

TOTAL_MONO=$(( ${#DISTANCES[@]} * ${#MODES[@]} ))
TOTAL_DIST=0
for distance in "${DISTANCES[@]}"; do
    for rate_hz in "${RATES_HZ[@]}"; do
        combos_str=$(enr_combos "$distance" "$rate_hz")
        n_combos=$(echo "$combos_str" | wc -w)
        TOTAL_DIST=$(( TOTAL_DIST + n_combos * ${#MODES[@]} ))
    done
done
TOTAL_CONFIGS=$(( TOTAL_MONO + TOTAL_DIST ))

FAIL_FILE="$RUN_DIR/failed.txt"
: > "$FAIL_FILE"

echo "=================================================="
echo "Neutral-Atom d=9/d=11 Lattice Surgery LER Sweep"
echo "=================================================="
echo "Run ID:    $TIMESTAMP"
echo "Output:    $RUN_DIR"
echo "Distances: ${DISTANCES[*]}"
echo "Modes:     ${MODES[*]}"
echo "Rates:     ${RATES_HZ[*]}"
echo "Shots:     ${SHOTS:-from config}"
echo "MaxRounds: $MAX_ROUNDS"
echo "Circuits:  $CIRCUITS_ONLY"
echo "MPI:       $NUM_PROCS processes"
echo "Parallel:  $PARALLEL_JOBS"
echo "Jobs:      $TOTAL_CONFIGS  (${TOTAL_MONO} monolithic + ${TOTAL_DIST} distributed)"
[[ -n "$RESUME_RUN_DIR" ]] && echo "Resume:    yes"
echo "=================================================="

cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: neutral_atom_d9_d11_lattice_surgery
Run ID: $TIMESTAMP
Started: $(date)
Distances: ${DISTANCES[*]}
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

run_job() {
    local idx="$1"
    local total="$2"
    local distance="$3"
    local arch="$4"
    local mode="$5"
    local rate_hz="$6"
    local protocol="${7:-}"
    local rounds="${8:-}"
    local batches="${9:-}"
    local config_file
    local label
    local merge_type
    local phase
    local result_file
    local log_file
    local circuit_file
    local dump_dir
    local extra_args=()

    config_file=$(config_path "$arch" "$distance")
    label=$(job_label "$arch" "$mode" "$distance" "$rate_hz" "$protocol" "$rounds" "$batches")
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
    outfile=$(find "$tmp_out" -maxdepth 1 -name 'results_*.txt' | head -1)
    if [[ -n "$outfile" ]]; then
        mv "$outfile" "$result_file"
    else
        echo "$label" >> "$FAIL_FILE"
        echo "  Warning: no result file found for $label" >> "$log_file"
    fi
    rmdir "$tmp_out" 2>/dev/null || true
}

wait_for_slot() {
    while (( $(jobs -rp | wc -l) >= PARALLEL_JOBS )); do
        wait -n 2>/dev/null || true
    done
}

TASK_INDEX=0

for distance in "${DISTANCES[@]}"; do
    echo ""
    echo "========================================"
    echo "Distance d=$distance"
    echo "========================================"

    for mode in "${MODES[@]}"; do
        TASK_INDEX=$((TASK_INDEX + 1))
        wait_for_slot
        run_job "$TASK_INDEX" "$TOTAL_CONFIGS" "$distance" monolithic "$mode" 0 &
    done

    for rate_hz in "${RATES_HZ[@]}"; do
        combos_str=$(enr_combos "$distance" "$rate_hz")
        for combo in $combos_str; do
            IFS=':' read -r protocol rounds batches <<< "$combo"
            for mode in "${MODES[@]}"; do
                TASK_INDEX=$((TASK_INDEX + 1))
                wait_for_slot
                run_job "$TASK_INDEX" "$TOTAL_CONFIGS" "$distance" distributed "$mode" "$rate_hz" "$protocol" "$rounds" "$batches" &
            done
        done
    done
done

wait

result_count=$(find "$RESULTS_DIR" -maxdepth 1 -name '*.txt' | wc -l)
fail_count=$(wc -l < "$FAIL_FILE" 2>/dev/null || echo 0)

echo ""
echo "========================================"
echo "Run complete: $RUN_DIR"
echo "  Results: $result_count / $TOTAL_CONFIGS"
echo "  Failed:  $fail_count"
[[ "$fail_count" -gt 0 ]] && echo "--- failures ---" && cat "$FAIL_FILE"
echo "========================================"
