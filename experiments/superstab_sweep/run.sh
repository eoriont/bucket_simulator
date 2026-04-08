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
#   -d, --distances D,...  Comma-separated distances to sweep, e.g. 5,7,9 (default: use config value)
#   -p, --phase PHASE      Experiment phase: merge_and_split, merge_only, split_only (default: merge_and_split)
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
DISTANCES=()
PHASE="merge_and_split"

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)   NUM_PROCS="$2"; shift 2 ;;
        -j|--parallel)    PARALLEL_JOBS="$2"; shift 2 ;;
        -d|--distances)   IFS=',' read -ra DISTANCES <<< "$2"; shift 2 ;;
        -p|--phase)
            PHASE="$2"
            case "$PHASE" in
                merge_and_split|merge_only|split_only|full) ;;
                *) echo "Unknown phase: $PHASE"; exit 1 ;;
            esac
            shift 2
            ;;
        --dry-run)        DRY_RUN=true; shift ;;
        --shots)          SHOTS="$2"; shift 2 ;;
        -h|--help)        head -20 "$0" | tail -18; exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Superstabilizer type definitions ──────────────────────────────────────────
# Args: <type> <distance>
# Seam data qubits are at x=d+0.5, y=0.5,1.5,...,d-0.5
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

# Args: <type> <distance>
get_merge_rounds() {
    case "$1" in
        middle) echo "1" ;;
        *)      echo "$2" ;;
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

NUM_DISTANCES=$(( ${#DISTANCES[@]} > 0 ? ${#DISTANCES[@]} : 1 ))
TOTAL_CONFIGS=$(( ${#SUPERS[@]} * ${#RATES[@]} * NUM_DISTANCES ))

echo "=============================================="
echo "Superstabilizer Sweep"
echo "=============================================="
echo "Run ID:    $TIMESTAMP"
echo "Output:    $RUN_DIR"
echo "Config:    $BASE_CONFIG"
echo "Supers:    ${SUPERS[*]}"
echo "Rates:     ${RATES[*]}"
echo "Distances: ${DISTANCES[*]:-from config}"
echo "Shots:     ${SHOTS:-from config}"
echo "Phase:     $PHASE"
echo "MPI:       $NUM_PROCS processes"
echo "Parallel:  $PARALLEL_JOBS"
echo "Configs:   $TOTAL_CONFIGS"
echo "=============================================="

cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: superstab_sweep
Run ID: $TIMESTAMP
Started: $(date)
Base config: $BASE_CONFIG
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Distances: ${DISTANCES[*]:-from config}
Shots override: ${SHOTS:-none}
Experiment phase: $PHASE
Superstabilizer types: ${SUPERS[*]}
Entanglement rates (Hz): ${RATES[*]}
Host: $(hostname)
Simulator: $SIMULATOR
EOF

# ── Runner function ───────────────────────────────────────────────────────────
# Args: <type> <rate> <distance>  (distance="" means use config value)
run_job() {
    local name="$1"
    local rate="$2"
    local dist="$3"
    local rate_mhz
    rate_mhz=$(echo "scale=0; $rate / 1000000" | bc)
    local label="${name}_${rate_mhz}MHz${dist:+_d${dist}}"
    local result_file="$RESULTS_DIR/${label}_results.txt"
    local log_file="$LOGS_DIR/${label}.log"

    local coords
    coords=$(get_super_coords "$name" "${dist:-5}")
    local mrounds
    mrounds=$(get_merge_rounds "$name" "${dist:-5}")

    local extra_args=()
    [[ -n "$dist" ]] && extra_args+=(-code_distance "$dist")
    [[ -n "$SHOTS" ]] && extra_args+=(-total_shots "$SHOTS")

    echo "[RUN] $label  (superstabilizers=\"$coords\" merge_rounds=$mrounds)"
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  mpirun -n $NUM_PROCS $SIMULATOR \\"
        echo "    -config $BASE_CONFIG \\"
        echo "    -experiment_phase $PHASE \\"
        echo "    -superstabilizers \"$coords\" \\"
        echo "    -merge_rounds $mrounds \\"
        echo "    -entanglement_rate $rate \\"
        [[ -n "$dist" ]] && echo "    -code_distance $dist \\"
        [[ -n "$SHOTS" ]] && echo "    -total_shots $SHOTS \\"
        echo "    -output $RESULTS_DIR"
        return 0
    fi

    # Dump Stim circuit (single rank, no simulation)
    mpirun -n 1 "$SIMULATOR" \
            -config "$BASE_CONFIG" \
            -experiment_phase "$PHASE" \
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
            -experiment_phase "$PHASE" \
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
export RESULTS_DIR CIRCUITS_DIR LOGS_DIR SIMULATOR NUM_PROCS DRY_RUN SHOTS BASE_CONFIG PHASE

echo ""
echo "Running simulations..."
START=$(date +%s)
FAILED=0

if [[ ${#DISTANCES[@]} -gt 0 ]]; then
    for dist in "${DISTANCES[@]}"; do
        for name in "${SUPERS[@]}"; do
            for rate in "${RATES[@]}"; do
                run_job "$name" "$rate" "$dist" || ((FAILED++)) || true
            done
        done
    done
else
    for name in "${SUPERS[@]}"; do
        for rate in "${RATES[@]}"; do
            run_job "$name" "$rate" "" || ((FAILED++)) || true
        done
    done
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
