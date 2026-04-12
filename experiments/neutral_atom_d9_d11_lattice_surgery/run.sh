#!/bin/bash
#
# Neutral-atom d=9 and d=11 lattice-surgery LER sweep.
# Runs:
#   - Monolithic baseline for d=9 and d=11 (2 modes × 2 distances = 4 jobs)
#   - Distributed: all feasible (protocol, k, m) combinations at each ENR
#     for both distances (combos computed at runtime via compute_combos.py)
#
# Note: split modes are excluded — the split observable requires separate
# investigation. Only xx_merge and zz_merge are run here.
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
#   --shots SHOTS          Override total_shots, e.g. 500K or 1M

set -euo pipefail

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

# Distances to sweep
DISTANCES=(9 11)

# Modes (split modes excluded — observable not validated)
MODES=(xx_merge zz_merge)

# ENR sweep points (Hz), bracketing all regime transitions for both d=9 and d=11
# at F=0.99, p=0.001:
#   ~75kHz  d=9 first transition (none → 2to1_pump)
#  ~100kHz  d=11 first transition
#  ~133kHz  both: → 3to1_pump
#  ~178kHz  both: → 2to1_rec k=2
#  ~316kHz  d=9:  → 2to1_rec k=3
#  ~422kHz  both: → 3to1_rec k=2  (best achievable: eff=0.30%)
RATES_HZ=(30000 80000 120000 160000 250000 400000 600000 1000000)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)  NUM_PROCS="$2";    shift 2 ;;
        -j|--parallel)   PARALLEL_JOBS="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=true;       shift   ;;
        --circuits-only)  CIRCUITS_ONLY=true; shift   ;;
        --shots)          SHOTS="$2";         shift 2 ;;
        --resume)         RESUME_RUN_DIR="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Run directory setup
# ---------------------------------------------------------------------------
if [[ -n "$RESUME_RUN_DIR" ]]; then
    RUN_DIR="$RESUME_RUN_DIR"
    echo "Resuming run: $RUN_DIR"
else
    RUN_ID="$(date +%Y%m%d_%H%M%S)"
    RUN_DIR="$SCRIPT_DIR/runs/$RUN_ID"
    mkdir -p "$RUN_DIR"/{results,logs,circuits}
fi

RESULTS_DIR="$RUN_DIR/results"
LOG_DIR="$RUN_DIR/logs"
CIRCUIT_DIR="$RUN_DIR/circuits"
FAIL_FILE="$RUN_DIR/failed.txt"

# ---------------------------------------------------------------------------
# Compute total jobs for metadata
# ---------------------------------------------------------------------------
total_jobs=0
for DISTANCE in "${DISTANCES[@]}"; do
    DIST_CONFIG="$SCRIPT_DIR/config_distributed_d${DISTANCE}.txt"
    for mode in "${MODES[@]}"; do
        # monolithic
        total_jobs=$((total_jobs + 1))
        # distributed
        for enr in "${RATES_HZ[@]}"; do
            combos=$(python3 "$SCRIPT_DIR/compute_combos.py" \
                --config "$DIST_CONFIG" --enr "$enr" --max-rounds "$MAX_ROUNDS" 2>/dev/null)
            for _ in $combos; do
                total_jobs=$((total_jobs + 1))
            done
        done
    done
done

# ---------------------------------------------------------------------------
# Write metadata
# ---------------------------------------------------------------------------
if [[ -z "$RESUME_RUN_DIR" ]]; then
    cat > "$RUN_DIR/metadata.txt" <<EOF
Experiment: neutral_atom_d9_d11_lattice_surgery
Run ID: $RUN_ID
Started: $(date)
Distances: ${DISTANCES[*]}
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Shots override: ${SHOTS:-none}
Max distillation rounds: $MAX_ROUNDS
Circuits only: $CIRCUITS_ONLY
Modes: ${MODES[*]}
Distributed entanglement rates (Hz): ${RATES_HZ[*]}
Total jobs: $total_jobs
Host: $(hostname)
Simulator: $SIMULATOR
EOF
    echo "Run directory: $RUN_DIR"
    echo "Total jobs: $total_jobs"
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
ACTIVE_JOBS=0

# Build base simulator args
sim_args() {
    local config="$1" label="$2"
    local args=(-config "$config" -output "$RUN_DIR/.tmp_out_${label}")
    [[ -n "$SHOTS" ]]         && args+=(-shots "$SHOTS")
    $CIRCUITS_ONLY             && args+=(--circuits-only -circuit-output "$CIRCUIT_DIR/${label}.stim")
    ! $CIRCUITS_ONLY           && args+=(-circuit-output "$CIRCUIT_DIR/${label}.stim")
    echo "${args[@]}"
}

run_job() {
    local label="$1" config="$2"
    local result_file="$RESULTS_DIR/${label}_results.txt"
    local log_file="$LOG_DIR/${label}.log"

    # Skip if already done (resume mode)
    if [[ -f "$result_file" ]]; then
        return 0
    fi

    local tmp_out="$RUN_DIR/.tmp_out_${label}"
    rm -rf "$tmp_out"
    mkdir -p "$tmp_out"

    local args
    args=$(sim_args "$config" "$label")

    local cmd="mpirun -n $NUM_PROCS $SIMULATOR $args"
    echo "[$(date +%H:%M:%S)] $label"

    if $DRY_RUN; then
        echo "  DRY-RUN: $cmd"
        rmdir "$tmp_out" 2>/dev/null || true
        return 0
    fi

    $cmd >> "$log_file" 2>&1 &
    local pid=$!

    # Wait and collect
    wait "$pid"
    local exit_code=$?

    if [[ $exit_code -ne 0 ]]; then
        echo "$label (exit $exit_code)" >> "$FAIL_FILE"
        rmdir "$tmp_out" 2>/dev/null || true
        return 0
    fi

    if $CIRCUITS_ONLY; then
        rmdir "$tmp_out" 2>/dev/null || true
        return 0
    fi

    local outfile
    outfile=$(ls "$tmp_out"/*.txt 2>/dev/null | head -1)
    if [[ -n "$outfile" ]]; then
        mv "$outfile" "$result_file"
    else
        echo "$label (no output)" >> "$FAIL_FILE"
    fi
    rmdir "$tmp_out" 2>/dev/null || true
}

wait_for_slot() {
    while [[ $(jobs -rp | wc -l) -ge $PARALLEL_JOBS ]]; do
        wait -n 2>/dev/null || true
    done
}

# ---------------------------------------------------------------------------
# Write a temporary config with specific overrides applied
# ---------------------------------------------------------------------------
make_config() {
    local base_config="$1"
    local out_config="$2"
    shift 2
    # Copy base, then append overrides
    cp "$base_config" "$out_config"
    while [[ $# -ge 2 ]]; do
        local key="$1" val="$2"
        shift 2
        # Remove existing line for this key (if any), then append
        grep -v "^${key}\b" "$out_config" > "${out_config}.tmp" && mv "${out_config}.tmp" "$out_config"
        echo "$key $val" >> "$out_config"
    done
}

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
job_num=0
for DISTANCE in "${DISTANCES[@]}"; do
    MONO_CONFIG="$SCRIPT_DIR/config_monolithic_d${DISTANCE}.txt"
    DIST_CONFIG="$SCRIPT_DIR/config_distributed_d${DISTANCE}.txt"

    echo ""
    echo "========================================"
    echo "  Distance d=$DISTANCE"
    echo "========================================"

    # --- Monolithic jobs ---
    for mode in "${MODES[@]}"; do
        merge_type="${mode//_merge/}"  # xx or zz
        label="monolithic_${mode}_d${DISTANCE}"
        result_file="$RESULTS_DIR/${label}_results.txt"

        if [[ -f "$result_file" ]]; then
            echo "[skip] $label (already done)"
            continue
        fi

        job_num=$((job_num + 1))
        echo "[JOB $job_num/$total_jobs] $label"

        # Write a per-job config with the correct merge_type
        tmp_cfg="$RUN_DIR/.cfg_${label}.txt"
        make_config "$MONO_CONFIG" "$tmp_cfg" \
            merge_type "${merge_type}_merge"

        wait_for_slot
        run_job "$label" "$tmp_cfg" &
    done

    # --- Distributed jobs ---
    for enr in "${RATES_HZ[@]}"; do
        enr_khz=$((enr / 1000))

        # Get feasible combos for this distance+ENR
        combos=$(python3 "$SCRIPT_DIR/compute_combos.py" \
            --config "$DIST_CONFIG" --enr "$enr" --max-rounds "$MAX_ROUNDS" 2>/dev/null)

        for combo in $combos; do
            IFS=':' read -r protocol rounds batches <<< "$combo"

            for mode in "${MODES[@]}"; do
                merge_type="${mode//_merge/}"

                if [[ "$protocol" == "none" ]]; then
                    label="distributed_${mode}_none_${enr_khz}kHz_d${DISTANCE}"
                else
                    label="distributed_${mode}_${protocol}_k${rounds}_m${batches}_${enr_khz}kHz_d${DISTANCE}"
                fi

                result_file="$RESULTS_DIR/${label}_results.txt"
                if [[ -f "$result_file" ]]; then
                    echo "[skip] $label (already done)"
                    continue
                fi

                job_num=$((job_num + 1))
                echo "[JOB $job_num/$total_jobs] $label"

                tmp_cfg="$RUN_DIR/.cfg_${label}.txt"
                make_config "$DIST_CONFIG" "$tmp_cfg" \
                    merge_type "distributed_${merge_type}_merge" \
                    entanglement_rate "$enr" \
                    distillation_protocol "$protocol" \
                    distillation_rounds "$rounds" \
                    distillation_backup_batches "$batches"

                wait_for_slot
                run_job "$label" "$tmp_cfg" &
            done
        done
    done
done

# Wait for all remaining jobs
wait

# Clean up temp config files
rm -f "$RUN_DIR"/.cfg_*.txt

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
result_count=$(ls "$RESULTS_DIR"/*.txt 2>/dev/null | wc -l)
fail_count=$(wc -l < "$FAIL_FILE" 2>/dev/null || echo 0)

echo ""
echo "========================================"
echo "Run complete: $RUN_DIR"
echo "  Results: $result_count"
echo "  Failed:  $fail_count"
[[ $fail_count -gt 0 ]] && cat "$FAIL_FILE"
echo "========================================"
