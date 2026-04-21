#!/bin/bash
#
# Compare round-count trends for a no-superstabilizer circuit and a middle
# superstabilizer circuit. The middle circuit is plotted/run with the same raw
# merge_rounds input, but each merge round is implemented as two gauge half-rounds.
#
# Usage:
#   bash run.sh [options]
#
# Options:
#   -n, --num-procs NUM    MPI processes per sim (default: 2)
#   -j, --parallel NUM     Simulations in parallel (default: 5)
#   --dry-run              Print commands without executing
#   --circuits-only        Dump circuits for all configs without running simulations
#   --resume RUN_DIR       Resume an existing run directory
#   --shots SHOTS          Override total_shots, e.g. 300K
#   --rate HZ              Override entanglement rate in Hz (default: 10000000)
#   --mode MODE            xx_merge or zz_merge (default: xx_merge)
#   --rounds A,B,...       Comma-separated raw merge_rounds values (default: 1..13)
#   --distance D           Code distance, must be odd >= 3 (default: 9)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BASE_CONFIG="$SCRIPT_DIR/config.txt"
SIMULATOR="$PROJECT_ROOT/build/bucket_simulator"

NUM_PROCS=2
PARALLEL_JOBS=5
DRY_RUN=false
CIRCUITS_ONLY=false
SHOTS=""
RATE_HZ=10000000
MODE="xx_merge"
DISTANCE=9
RESUME_RUN_DIR=""
ROUND_VALUES=(1 2 3 4 5 6 7 8 9 10 11 12 13)
SUPERS=(nosuper middle)

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs) NUM_PROCS="$2"; shift 2 ;;
        -j|--parallel) PARALLEL_JOBS="$2"; shift 2 ;;
        --dry-run) DRY_RUN=true; shift ;;
        --circuits-only) CIRCUITS_ONLY=true; shift ;;
        --resume) RESUME_RUN_DIR="$2"; shift 2 ;;
        --shots) SHOTS="$2"; shift 2 ;;
        --rate) RATE_HZ="$2"; shift 2 ;;
        --mode)
            MODE="$2"
            case "$MODE" in
                xx_merge|zz_merge) ;;
                *) echo "Unknown mode: $MODE"; exit 1 ;;
            esac
            shift 2
            ;;
        --rounds)
            IFS=',' read -ra ROUND_VALUES <<< "$2"
            shift 2
            ;;
        --distance) DISTANCE="$2"; shift 2 ;;
        -h|--help)
            head -22 "$0" | tail -20
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

get_super_coords() {
    local type="$1"
    local d="$2"
    local x
    x=$(echo "$d + 0.5" | bc)
    local y_mid
    y_mid=$(echo "scale=1; ($d - 1) / 2 + 0.5" | bc)
    case "$type" in
        nosuper) echo "none" ;;
        middle) echo "($x,$y_mid)" ;;
    esac
}

mode_merge_type() {
    case "$1" in
        xx_merge) echo "distributed_xx" ;;
        zz_merge) echo "distributed_zz" ;;
    esac
}

mode_phase() {
    echo "merge_only"
}

rate_label() {
    local rate_hz="$1"
    if [[ "$rate_hz" -ge 1000000 ]]; then
        echo "$((rate_hz / 1000000))MHz"
    elif [[ "$rate_hz" -ge 1000 ]]; then
        echo "$((rate_hz / 1000))kHz"
    else
        echo "${rate_hz}Hz"
    fi
}

job_label() {
    local super_type="$1"
    local rounds="$2"
    echo "${MODE}_${super_type}_r${rounds}_$(rate_label "$RATE_HZ")_d${DISTANCE}"
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
if [[ ! -f "$BASE_CONFIG" ]]; then
    echo "Error: base config not found at $BASE_CONFIG"
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

TOTAL_CONFIGS=$(( ${#SUPERS[@]} * ${#ROUND_VALUES[@]} ))
TASKS_FILE="$RUN_DIR/tasks.txt"
FAIL_FILE="$RUN_DIR/failed.txt"
: > "$FAIL_FILE"

echo "=============================================="
echo "Superstabilizer Round Sweep"
echo "=============================================="
echo "Run ID:    $TIMESTAMP"
echo "Output:    $RUN_DIR"
echo "Config:    $BASE_CONFIG"
echo "Distance:  $DISTANCE"
echo "Mode:      $MODE"
echo "Supers:    ${SUPERS[*]}"
echo "Rounds:    ${ROUND_VALUES[*]}"
echo "ENR:       $RATE_HZ Hz"
echo "Shots:     ${SHOTS:-from config}"
echo "Circuits:  $CIRCUITS_ONLY"
echo "MPI:       $NUM_PROCS processes"
echo "Parallel:  $PARALLEL_JOBS"
echo "Configs:   $TOTAL_CONFIGS"
[[ -n "$RESUME_RUN_DIR" ]] && echo "Resume:    yes"
echo "=============================================="

cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: superstab_round_sweep
Run ID: $TIMESTAMP
Started: $(date)
Base config: $BASE_CONFIG
Distance: $DISTANCE
Mode: $MODE
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Shots override: ${SHOTS:-none}
Entanglement rate (Hz): $RATE_HZ
Superstabilizers: ${SUPERS[*]}
Raw merge_rounds values: ${ROUND_VALUES[*]}
Middle effective rounds: 2 * merge_rounds
Host: $(hostname)
Simulator: $SIMULATOR
EOF

run_job() {
    local idx="$1"
    local total="$2"
    local super_type="$3"
    local rounds="$4"
    local label
    label=$(job_label "$super_type" "$rounds")
    local result_file="$RESULTS_DIR/${label}_results.txt"
    local log_file="$LOGS_DIR/${label}.log"
    local circuit_file="$CIRCUITS_DIR/${label}.stim"
    local dump_dir="$CIRCUITS_DIR/.tmp_${label}"
    local result_dir="$RESULTS_DIR/.tmp_${label}"
    local coords
    coords=$(get_super_coords "$super_type" "$DISTANCE")
    local merge_type
    merge_type=$(mode_merge_type "$MODE")
    local phase
    phase=$(mode_phase "$MODE")
    local extra_args=()

    [[ -n "$SHOTS" ]] && extra_args+=(-total_shots "$SHOTS")
    extra_args+=(-code_distance "$DISTANCE")

    if job_is_done "$label"; then
        echo "[SKIP $idx/$total] $label (already complete)"
        return 0
    fi

    local effective_rounds="$rounds"
    if [[ "$super_type" == "middle" ]]; then
        effective_rounds=$(( 2 * rounds ))
    fi

    echo "[RUN $idx/$total] $label  (raw_rounds=$rounds effective_rounds=$effective_rounds superstabilizers=\"$coords\")"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  mpirun -n $NUM_PROCS $SIMULATOR \\"
        echo "    -config $BASE_CONFIG \\"
        echo "    -merge_type $merge_type \\"
        echo "    -experiment_phase $phase \\"
        echo "    -superstabilizers \"$coords\" \\"
        echo "    -merge_rounds $rounds \\"
        echo "    -entanglement_rate $RATE_HZ \\"
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
        -config "$BASE_CONFIG" \
        -merge_type "$merge_type" \
        -experiment_phase "$phase" \
        -superstabilizers "$coords" \
        -merge_rounds "$rounds" \
        -entanglement_rate "$RATE_HZ" \
        "${extra_args[@]}" \
        -dump-circuit \
        -output "$dump_dir" >> "$log_file" 2>&1

    local dumped="$dump_dir/config.stim"
    if [[ -f "$dumped" ]]; then
        mv "$dumped" "$circuit_file"
    fi
    rmdir "$dump_dir" 2>/dev/null || true

    if [[ "$CIRCUITS_ONLY" == "true" ]]; then
        echo "[DONE $idx/$total] $label (circuit only)"
        return 0
    fi

    rm -rf "$result_dir"
    mkdir -p "$result_dir"
    local start_time
    start_time=$(date +%s)
    if mpirun -n "$NUM_PROCS" "$SIMULATOR" \
        -config "$BASE_CONFIG" \
        -merge_type "$merge_type" \
        -experiment_phase "$phase" \
        -superstabilizers "$coords" \
        -merge_rounds "$rounds" \
        -entanglement_rate "$RATE_HZ" \
        "${extra_args[@]}" \
        -output "$result_dir" >> "$log_file" 2>&1; then
        local end_time
        end_time=$(date +%s)
        echo "[DONE $idx/$total] $label ($((end_time - start_time))s)"
        local produced
        produced=$(find "$result_dir" -maxdepth 1 -name 'results_*.txt' | head -1)
        if [[ -n "$produced" ]]; then
            mv "$produced" "$result_file"
        fi
        rmdir "$result_dir" 2>/dev/null || true
    else
        rmdir "$result_dir" 2>/dev/null || true
        echo "$label" >> "$FAIL_FILE"
        echo "[FAIL $idx/$total] $label — see $log_file"
        return 1
    fi
}

run_job_record() {
    local record="$1"
    local idx total super_type rounds
    IFS='|' read -r idx total super_type rounds <<< "$record"
    run_job "$idx" "$total" "$super_type" "$rounds"
}

export -f run_job run_job_record get_super_coords mode_merge_type mode_phase rate_label job_label job_is_done
export RESULTS_DIR CIRCUITS_DIR LOGS_DIR SIMULATOR NUM_PROCS DRY_RUN CIRCUITS_ONLY SHOTS BASE_CONFIG DISTANCE MODE RATE_HZ FAIL_FILE

echo ""
echo "Running simulations..."
START=$(date +%s)
COMPLETED=0

task_idx=0
> "$TASKS_FILE"
for super_type in "${SUPERS[@]}"; do
    for rounds in "${ROUND_VALUES[@]}"; do
        label=$(job_label "$super_type" "$rounds")
        if job_is_done "$label"; then
            COMPLETED=$((COMPLETED + 1))
            continue
        fi
        task_idx=$((task_idx + 1))
        printf '%s|%s|%s|%s\n' "$task_idx" "$TOTAL_CONFIGS" "$super_type" "$rounds" >> "$TASKS_FILE"
    done
done

PENDING=$task_idx
echo "Progress:  $COMPLETED/$TOTAL_CONFIGS complete"
if [[ "$PENDING" -eq 0 ]]; then
    echo "Nothing to do."
fi

if [[ "$PENDING" -gt 0 && ( "$DRY_RUN" == "true" || "$PARALLEL_JOBS" -le 1 ) ]]; then
    while IFS= read -r record; do
        run_job_record "$record" || true
    done < "$TASKS_FILE"
elif [[ "$PENDING" -gt 0 ]]; then
    echo "Launching up to $PARALLEL_JOBS MPI jobs in parallel..."
    xargs -P "$PARALLEL_JOBS" -I {} bash -lc 'run_job_record "$1"' _ {} < "$TASKS_FILE" || true
fi

FAILED=0
if [[ -f "$FAIL_FILE" ]]; then
    FAILED=$(wc -l < "$FAIL_FILE" | tr -d ' ')
fi

END=$(date +%s)
TOTAL=$((END - START))

echo ""
echo "=============================================="
echo "Sweep Complete"
echo "=============================================="
echo "Total time: ${TOTAL}s"
echo "Failed:     $FAILED"
echo "Results in: $RESULTS_DIR"

SUMMARY="$RUN_DIR/summary.txt"
{
    echo "=============================================="
    echo "Superstabilizer Round Sweep Summary"
    echo "=============================================="
    echo "Run ID:    $TIMESTAMP"
    echo "Completed: $(date)"
    echo "Runtime:   ${TOTAL}s"
    echo "Failed:    $FAILED"
    echo ""
    printf "%-40s  %-6s  %-9s  %s\n" "Config" "RawR" "EffR" "LER"
    echo "------------------------------------------------------------------"
    for r in "$RESULTS_DIR"/*_results.txt; do
        [[ -f "$r" ]] || continue
        cfg=$(basename "$r" _results.txt)
        raw=$(echo "$cfg" | sed -E 's/.*_r([0-9]+)_.*/\1/')
        eff="$raw"
        if [[ "$cfg" == *"_middle_"* ]]; then
            eff=$(( 2 * raw ))
        fi
        ler=$(grep "Logical Error Rate:" "$r" | head -1 | awk '{print $4}')
        printf "%-40s  %-6s  %-9s  %s\n" "$cfg" "$raw" "$eff" "$ler"
    done | sort
} | tee "$SUMMARY"

echo ""
echo "To plot:  python $SCRIPT_DIR/plot.py $RUN_DIR"
