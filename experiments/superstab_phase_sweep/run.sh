#!/bin/bash
#
# Superstabilizer phase/basis sweep for d=9.
#
# Sweeps all superstabilizer placements across:
#   - xx_merge
#   - zz_merge
#   - xx_split
#   - zz_split
#
# Usage:
#   ./run.sh [options]
#
# Options:
#   -n, --num-procs NUM    MPI processes per sim (default: 4)
#   -j, --parallel NUM     Simulations in parallel (default: 1)
#   --dry-run              Print commands without executing
#   --circuits-only        Dump circuits for all configs without running simulations
#   --shots SHOTS          Override total_shots, e.g. 10K (default: use config value)
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BASE_CONFIG="$SCRIPT_DIR/config.txt"
SIMULATOR="$PROJECT_ROOT/build/bucket_simulator"

NUM_PROCS=4
PARALLEL_JOBS=1
DRY_RUN=false
CIRCUITS_ONLY=false
SHOTS=""
DISTANCE=9

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)   NUM_PROCS="$2"; shift 2 ;;
        -j|--parallel)    PARALLEL_JOBS="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=true; shift ;;
        --circuits-only)  CIRCUITS_ONLY=true; shift ;;
        --shots)          SHOTS="$2"; shift 2 ;;
        -h|--help)        head -22 "$0" | tail -20; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

get_super_coords() {
    local type="$1"
    local d="$2"
    local x
    x=$(echo "$d + 0.5" | bc)
    local y_top="0.5"
    local y_bot
    y_bot=$(echo "$d - 0.5" | bc)
    local y_mid
    y_mid=$(echo "scale=1; ($d - 1) / 2 + 0.5" | bc)
    case "$type" in
        border)  echo "($x,$y_top)" ;;
        border2) echo "($x,$y_bot)" ;;
        middle)  echo "($x,$y_mid)" ;;
        nosuper) echo "none" ;;
        twoside) echo "($x,$y_top) ($x,$y_bot)" ;;
    esac
}

get_merge_rounds() {
    case "$1" in
        middle) echo "1" ;;
        *)      echo "$2" ;;
    esac
}

mode_merge_type() {
    case "$1" in
        xx_merge|xx_split) echo "distributed_xx" ;;
        zz_merge|zz_split) echo "distributed_zz" ;;
    esac
}

mode_phase() {
    case "$1" in
        xx_merge|zz_merge) echo "merge_only" ;;
        xx_split|zz_split) echo "split_only" ;;
    esac
}

SUPERS=(border border2 middle nosuper twoside)
MODES=(xx_merge zz_merge xx_split zz_split)
RATES=(5000000 10000000 15000000 20000000 25000000 30000000 35000000 40000000 45000000 50000000)

if [[ ! -x "$SIMULATOR" ]]; then
    echo "Error: simulator not found at $SIMULATOR"
    echo "Build first: cd $PROJECT_ROOT/build && make"
    exit 1
fi
if [[ ! -f "$BASE_CONFIG" ]]; then
    echo "Error: base config not found at $BASE_CONFIG"
    exit 1
fi

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="$SCRIPT_DIR/runs/$TIMESTAMP"
RESULTS_DIR="$RUN_DIR/results"
CIRCUITS_DIR="$RUN_DIR/circuits"
LOGS_DIR="$RUN_DIR/logs"
mkdir -p "$RESULTS_DIR" "$CIRCUITS_DIR" "$LOGS_DIR"

TOTAL_CONFIGS=$(( ${#SUPERS[@]} * ${#MODES[@]} * ${#RATES[@]} ))
TASKS_FILE="$RUN_DIR/tasks.txt"
FAIL_FILE="$RUN_DIR/failed.txt"
touch "$FAIL_FILE"

echo "=============================================="
echo "Superstabilizer Phase/Basis Sweep"
echo "=============================================="
echo "Run ID:    $TIMESTAMP"
echo "Output:    $RUN_DIR"
echo "Config:    $BASE_CONFIG"
echo "Distance:  $DISTANCE"
echo "Modes:     ${MODES[*]}"
echo "Supers:    ${SUPERS[*]}"
echo "Rates:     ${RATES[*]}"
echo "Shots:     ${SHOTS:-from config}"
echo "Circuits:  $CIRCUITS_ONLY"
echo "MPI:       $NUM_PROCS processes"
echo "Parallel:  $PARALLEL_JOBS"
echo "Configs:   $TOTAL_CONFIGS"
echo "=============================================="

cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: superstab_phase_sweep
Run ID: $TIMESTAMP
Started: $(date)
Base config: $BASE_CONFIG
Distance: $DISTANCE
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Shots override: ${SHOTS:-none}
Circuits only: $CIRCUITS_ONLY
Modes: ${MODES[*]}
Superstabilizer types: ${SUPERS[*]}
Entanglement rates (Hz): ${RATES[*]}
Host: $(hostname)
Simulator: $SIMULATOR
EOF

run_job() {
    local idx="$1"
    local total="$2"
    local mode="$3"
    local name="$4"
    local rate="$5"
    local rate_mhz
    rate_mhz=$(echo "scale=0; $rate / 1000000" | bc)
    local label="${mode}_${name}_${rate_mhz}MHz_d${DISTANCE}"
    local result_file="$RESULTS_DIR/${label}_results.txt"
    local log_file="$LOGS_DIR/${label}.log"
    local circuit_file="$CIRCUITS_DIR/${label}.stim"
    local dump_dir="$CIRCUITS_DIR/.tmp_${label}"
    local coords
    coords=$(get_super_coords "$name" "$DISTANCE")
    local mrounds
    mrounds=$(get_merge_rounds "$name" "$DISTANCE")
    local merge_type
    merge_type=$(mode_merge_type "$mode")
    local phase
    phase=$(mode_phase "$mode")

    local extra_args=()
    [[ -n "$SHOTS" ]] && extra_args+=(-total_shots "$SHOTS")

    echo "[RUN $idx/$total] $label  (merge_type=$merge_type phase=$phase superstabilizers=\"$coords\" merge_rounds=$mrounds)"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  mpirun -n $NUM_PROCS $SIMULATOR \\"
        echo "    -config $BASE_CONFIG \\"
        echo "    -merge_type $merge_type \\"
        echo "    -experiment_phase $phase \\"
        echo "    -superstabilizers \"$coords\" \\"
        echo "    -merge_rounds $mrounds \\"
        echo "    -entanglement_rate $rate \\"
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
        -merge_rounds "$mrounds" \
        -entanglement_rate "$rate" \
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

    local start_time
    start_time=$(date +%s)
    if mpirun -n "$NUM_PROCS" "$SIMULATOR" \
        -config "$BASE_CONFIG" \
        -merge_type "$merge_type" \
        -experiment_phase "$phase" \
        -superstabilizers "$coords" \
        -merge_rounds "$mrounds" \
        -entanglement_rate "$rate" \
        "${extra_args[@]}" \
        -output "$RESULTS_DIR" >> "$log_file" 2>&1; then
        local end_time
        end_time=$(date +%s)
        echo "[DONE $idx/$total] $label ($((end_time - start_time))s)"
        local latest
        latest=$(ls -t "$RESULTS_DIR"/results_*.txt 2>/dev/null | head -1)
        if [[ -n "$latest" && "$latest" != "$result_file" ]]; then
            mv "$latest" "$result_file"
        fi
    else
        echo "$label" >> "$FAIL_FILE"
        echo "[FAIL $idx/$total] $label — see $log_file"
        return 1
    fi
}

run_job_record() {
    local record="$1"
    local idx total mode name rate
    IFS='|' read -r idx total mode name rate <<< "$record"
    run_job "$idx" "$total" "$mode" "$name" "$rate"
}

export -f run_job run_job_record get_super_coords get_merge_rounds mode_merge_type mode_phase
export RESULTS_DIR CIRCUITS_DIR LOGS_DIR SIMULATOR NUM_PROCS DRY_RUN CIRCUITS_ONLY SHOTS BASE_CONFIG DISTANCE FAIL_FILE

echo ""
echo "Running simulations..."
START=$(date +%s)
FAILED=0
echo "Progress:  0/$TOTAL_CONFIGS complete"

task_idx=0
> "$TASKS_FILE"
for mode in "${MODES[@]}"; do
    for name in "${SUPERS[@]}"; do
        for rate in "${RATES[@]}"; do
            task_idx=$((task_idx + 1))
            printf '%s|%s|%s|%s|%s\n' "$task_idx" "$TOTAL_CONFIGS" "$mode" "$name" "$rate" >> "$TASKS_FILE"
        done
    done
done

if [[ "$DRY_RUN" == "true" || "$PARALLEL_JOBS" -le 1 ]]; then
    while IFS= read -r record; do
        run_job_record "$record" || true
    done < "$TASKS_FILE"
else
    echo "Launching up to $PARALLEL_JOBS MPI jobs in parallel..."
    xargs -P "$PARALLEL_JOBS" -I {} bash -lc 'run_job_record "$1"' _ {} < "$TASKS_FILE" || true
fi

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
    echo "Superstabilizer Phase/Basis Sweep Summary"
    echo "=============================================="
    echo "Run ID:    $TIMESTAMP"
    echo "Completed: $(date)"
    echo "Runtime:   ${TOTAL}s"
    echo "Failed:    $FAILED"
    echo ""
    printf "%-48s  %s\n" "Config" "LER"
    echo "------------------------------------------------------------------"
    for r in "$RESULTS_DIR"/*_results.txt; do
        [[ -f "$r" ]] || continue
        cfg=$(basename "$r" _results.txt)
        ler=$(grep "Logical Error Rate:" "$r" | head -1 | awk '{print $4}')
        printf "%-48s  %s\n" "$cfg" "$ler"
    done
} | tee "$SUMMARY"

echo ""
echo "To plot:  python $SCRIPT_DIR/plot.py $RUN_DIR"
