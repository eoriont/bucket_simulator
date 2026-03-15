#!/bin/bash
#
# Superstabilizer sweep experiment
#
# Sweeps all 5 superstabilizer types across entanglement rates 5–50 MHz.
# Uses a single base config; superstabilizers, entanglement_rate, and
# merge_rounds are passed as CLI overrides.
#
# Usage:
#   ./run.sh [options]
#
# Options:
#   -n, --num-procs NUM    MPI processes per sim (default: 4)
#   -j, --parallel NUM     Simulations in parallel (default: 1)
#   -d, --distance NUM     Code distance (default: use config value)
#   --dry-run              Print commands without executing
#   --shots SHOTS          Override total_shots, e.g. 10K (default: use config value)
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BASE_CONFIG="$SCRIPT_DIR/config.txt"
SIMULATOR="$PROJECT_ROOT/build/bucket_simulator"

# Defaults
NUM_PROCS=4
PARALLEL_JOBS=1
DRY_RUN=false
SHOTS=""
DISTANCE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)  NUM_PROCS="$2"; shift 2 ;;
        -j|--parallel)   PARALLEL_JOBS="$2"; shift 2 ;;
        -d|--distance)   DISTANCE="$2"; shift 2 ;;
        --dry-run)       DRY_RUN=true; shift ;;
        --shots)         SHOTS="$2"; shift 2 ;;
        -h|--help)       head -20 "$0" | tail -18; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Superstabilizer type definitions ──────────────────────────────────────────
get_super_coords() {
    case "$1" in
        border)  echo "(5.5,0.5)" ;;
        border2) echo "(5.5,4.5)" ;;
        middle)  echo "(5.5,2.5)" ;;
        nosuper) echo "none" ;;
        twoside) echo "(5.5,0.5) (5.5,4.5)" ;;
    esac
}

get_merge_rounds() {
    case "$1" in
        middle) echo "1" ;;
        *)      echo "5" ;;
    esac
}

SUPERS=(border border2 middle nosuper twoside)

# Entanglement rates in Hz (5 MHz to 50 MHz)
RATES=(5000000 10000000 15000000 20000000 25000000 30000000 35000000 40000000 45000000 50000000)

# ── Validation ────────────────────────────────────────────────────────────────
if [[ ! -x "$SIMULATOR" ]]; then
    echo "Error: simulator not found at $SIMULATOR"
    echo "Build first: cd $PROJECT_ROOT/build && make"
    exit 1
fi
if [[ ! -f "$BASE_CONFIG" ]]; then
    echo "Error: base config not found at $BASE_CONFIG"
    exit 1
fi

# ── Setup run directory ───────────────────────────────────────────────────────
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="$SCRIPT_DIR/runs/$TIMESTAMP"
RESULTS_DIR="$RUN_DIR/results"
CIRCUITS_DIR="$RUN_DIR/circuits"
LOGS_DIR="$RUN_DIR/logs"
mkdir -p "$RESULTS_DIR" "$CIRCUITS_DIR" "$LOGS_DIR"

TOTAL_CONFIGS=$(( ${#SUPERS[@]} * ${#RATES[@]} ))

echo "=============================================="
echo "Superstabilizer Sweep"
echo "=============================================="
echo "Run ID:   $TIMESTAMP"
echo "Output:   $RUN_DIR"
echo "Config:   $BASE_CONFIG"
echo "Supers:   ${SUPERS[*]}"
echo "Rates:    ${RATES[*]}"
echo "Distance: ${DISTANCE:-from config}"
echo "Shots:    ${SHOTS:-from config}"
echo "MPI:      $NUM_PROCS processes"
echo "Parallel: $PARALLEL_JOBS"
echo "Configs:  $TOTAL_CONFIGS"
echo "=============================================="

cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: superstab_sweep
Run ID: $TIMESTAMP
Started: $(date)
Base config: $BASE_CONFIG
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Distance override: ${DISTANCE:-none}
Shots override: ${SHOTS:-none}
Superstabilizer types: ${SUPERS[*]}
Entanglement rates (Hz): ${RATES[*]}
Host: $(hostname)
Simulator: $SIMULATOR
EOF

# ── Runner function ───────────────────────────────────────────────────────────
run_job() {
    local name="$1"
    local rate="$2"
    local rate_mhz
    rate_mhz=$(echo "scale=0; $rate / 1000000" | bc)
    local label="${name}_${rate_mhz}MHz${DISTANCE:+_d${DISTANCE}}"
    local result_file="$RESULTS_DIR/${label}_results.txt"
    local log_file="$LOGS_DIR/${label}.log"

    local coords
    coords=$(get_super_coords "$name")
    local mrounds
    mrounds=$(get_merge_rounds "$name")

    local extra_args=()
    [[ -n "$DISTANCE" ]] && extra_args+=(-code_distance "$DISTANCE")
    [[ -n "$SHOTS" ]] && extra_args+=(-total_shots "$SHOTS")

    echo "[RUN] $label  (superstabilizers=\"$coords\" merge_rounds=$mrounds)"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  mpirun -n $NUM_PROCS $SIMULATOR \\"
        echo "    -config $BASE_CONFIG \\"
        echo "    -superstabilizers \"$coords\" \\"
        echo "    -merge_rounds $mrounds \\"
        echo "    -entanglement_rate $rate \\"
        [[ -n "$DISTANCE" ]] && echo "    -code_distance $DISTANCE \\"
        [[ -n "$SHOTS" ]] && echo "    -total_shots $SHOTS \\"
        echo "    -output $RESULTS_DIR"
        return 0
    fi

    # Dump Stim circuit (single rank, no simulation)
    mpirun -n 1 "$SIMULATOR" \
            -config "$BASE_CONFIG" \
            -superstabilizers "$coords" \
            -merge_rounds "$mrounds" \
            -entanglement_rate "$rate" \
            "${extra_args[@]}" \
            -dump-circuit \
            -output "$CIRCUITS_DIR" >> "$log_file" 2>&1
    # Rename the generic "config.stim" to the labelled name
    local dumped="$CIRCUITS_DIR/config.stim"
    if [[ -f "$dumped" ]]; then
        mv "$dumped" "$CIRCUITS_DIR/${label}.stim"
    fi

    local start_time
    start_time=$(date +%s)
    if mpirun -n "$NUM_PROCS" "$SIMULATOR" \
            -config "$BASE_CONFIG" \
            -superstabilizers "$coords" \
            -merge_rounds "$mrounds" \
            -entanglement_rate "$rate" \
            "${extra_args[@]}" \
            -output "$RESULTS_DIR" >> "$log_file" 2>&1; then
        local end_time
        end_time=$(date +%s)
        echo "[DONE] $label ($((end_time - start_time))s)"
        local latest
        latest=$(ls -t "$RESULTS_DIR"/results_*.txt 2>/dev/null | head -1)
        if [[ -n "$latest" && "$latest" != "$result_file" ]]; then
            mv "$latest" "$result_file"
        fi
    else
        echo "[FAIL] $label — see $log_file"
        return 1
    fi
}
export -f run_job get_super_coords get_merge_rounds
export RESULTS_DIR CIRCUITS_DIR LOGS_DIR SIMULATOR NUM_PROCS DRY_RUN SHOTS DISTANCE BASE_CONFIG

echo ""
echo "Running simulations..."
START=$(date +%s)
FAILED=0

for name in "${SUPERS[@]}"; do
    for rate in "${RATES[@]}"; do
        run_job "$name" "$rate" || ((FAILED++)) || true
    done
done

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
    echo "Superstabilizer Sweep Summary"
    echo "=============================================="
    echo "Run ID:    $TIMESTAMP"
    echo "Completed: $(date)"
    echo "Runtime:   ${TOTAL}s"
    echo "Failed:    $FAILED"
    echo ""
    printf "%-40s  %s\n" "Config" "LER"
    echo "--------------------------------------------------------------"
    for r in "$RESULTS_DIR"/*_results.txt; do
        [[ -f "$r" ]] || continue
        cfg=$(basename "$r" _results.txt)
        ler=$(grep "Logical Error Rate:" "$r" | head -1 | awk '{print $4}')
        printf "%-40s  %s\n" "$cfg" "$ler"
    done
} | tee "$SUMMARY"

echo ""
echo "To plot:  python $SCRIPT_DIR/plot.py $RUN_DIR"

[[ $FAILED -eq 0 ]] || exit 1
