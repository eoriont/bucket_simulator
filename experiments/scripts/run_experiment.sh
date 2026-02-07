#!/bin/bash
#
# Experiment Runner for bucket_simulator
#
# Usage:
#   ./run_experiment.sh <experiment_dir> [options]
#
# Options:
#   -n, --num-procs NUM    Number of MPI processes (default: 4)
#   -j, --parallel NUM     Run NUM simulations in parallel (default: 1)
#   --dry-run              Print commands without executing
#   --resume               Skip already completed configs
#
# Example:
#   ./run_experiment.sh radar_sweep -n 8 -j 2

set -e

# Default values
NUM_PROCS=4
PARALLEL_JOBS=1
DRY_RUN=false
RESUME=false

# Parse arguments
EXPERIMENT_DIR=""
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--num-procs)
            NUM_PROCS="$2"
            shift 2
            ;;
        -j|--parallel)
            PARALLEL_JOBS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --resume)
            RESUME=true
            shift
            ;;
        -h|--help)
            head -20 "$0" | tail -18
            exit 0
            ;;
        *)
            if [[ -z "$EXPERIMENT_DIR" ]]; then
                EXPERIMENT_DIR="$1"
            else
                echo "Unknown option: $1"
                exit 1
            fi
            shift
            ;;
    esac
done

if [[ -z "$EXPERIMENT_DIR" ]]; then
    echo "Error: experiment directory required"
    echo "Usage: $0 <experiment_dir> [options]"
    exit 1
fi

# Resolve paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPERIMENTS_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$EXPERIMENTS_ROOT")"
SIMULATOR="$PROJECT_ROOT/build/bucket_simulator"

# Check if experiment dir is absolute or relative to experiments/
if [[ "$EXPERIMENT_DIR" = /* ]]; then
    EXPERIMENT_PATH="$EXPERIMENT_DIR"
else
    EXPERIMENT_PATH="$EXPERIMENTS_ROOT/$EXPERIMENT_DIR"
fi

if [[ ! -d "$EXPERIMENT_PATH" ]]; then
    echo "Error: Experiment directory not found: $EXPERIMENT_PATH"
    exit 1
fi

# Check simulator exists
if [[ ! -x "$SIMULATOR" ]]; then
    echo "Error: Simulator not found at $SIMULATOR"
    echo "Please build the project first: cd build && cmake .. && make"
    exit 1
fi

# Create timestamped run directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="$EXPERIMENT_PATH/runs/$TIMESTAMP"
RESULTS_DIR="$RUN_DIR/results"
LOGS_DIR="$RUN_DIR/logs"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

echo "=============================================="
echo "Bucket Simulator Experiment Runner"
echo "=============================================="
echo "Experiment: $EXPERIMENT_DIR"
echo "Run ID: $TIMESTAMP"
echo "MPI processes: $NUM_PROCS"
echo "Parallel jobs: $PARALLEL_JOBS"
echo "Output: $RUN_DIR"
echo "=============================================="

# Find config files
CONFIG_DIR="$EXPERIMENT_PATH/configs"
if [[ ! -d "$CONFIG_DIR" ]]; then
    echo "Error: No configs directory found at $CONFIG_DIR"
    echo "Run the config generator first or create configs manually."
    exit 1
fi

CONFIG_FILES=($(find "$CONFIG_DIR" -name "*.txt" -type f | sort))
TOTAL_CONFIGS=${#CONFIG_FILES[@]}

if [[ $TOTAL_CONFIGS -eq 0 ]]; then
    echo "Error: No config files found in $CONFIG_DIR"
    exit 1
fi

echo "Found $TOTAL_CONFIGS configuration(s)"
echo ""

# Save experiment metadata
cat > "$RUN_DIR/metadata.txt" << EOF
Experiment: $EXPERIMENT_DIR
Run ID: $TIMESTAMP
Started: $(date)
MPI Processes: $NUM_PROCS
Parallel Jobs: $PARALLEL_JOBS
Total Configs: $TOTAL_CONFIGS
Host: $(hostname)
User: $(whoami)
Simulator: $SIMULATOR
EOF

# Function to run a single config
run_config() {
    local config_file="$1"
    local config_name=$(basename "$config_file" .txt)
    local result_file="$RESULTS_DIR/${config_name}_results.txt"
    local log_file="$LOGS_DIR/${config_name}.log"

    # Skip if resuming and result already exists
    if [[ "$RESUME" == "true" && -f "$result_file" ]]; then
        echo "[SKIP] $config_name (already completed)"
        return 0
    fi

    echo "[RUN] $config_name"

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  Would run: mpirun -n $NUM_PROCS $SIMULATOR -config $config_file -output $RESULTS_DIR"
        return 0
    fi

    # Run simulation
    local start_time=$(date +%s)
    if mpirun -n "$NUM_PROCS" "$SIMULATOR" -config "$config_file" -output "$RESULTS_DIR" > "$log_file" 2>&1; then
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        echo "[DONE] $config_name (${duration}s)"

        # Rename result file to include config name
        local latest_result=$(ls -t "$RESULTS_DIR"/results_*.txt 2>/dev/null | head -1)
        if [[ -n "$latest_result" && "$latest_result" != "$result_file" ]]; then
            mv "$latest_result" "$result_file"
        fi
    else
        echo "[FAIL] $config_name - check $log_file"
        return 1
    fi
}

export -f run_config
export RESULTS_DIR LOGS_DIR SIMULATOR NUM_PROCS DRY_RUN RESUME

# Run configs
START_TIME=$(date +%s)
FAILED=0

if [[ $PARALLEL_JOBS -gt 1 ]] && command -v parallel &> /dev/null; then
    # Use GNU parallel if available
    echo "Running $PARALLEL_JOBS simulations in parallel..."
    printf '%s\n' "${CONFIG_FILES[@]}" | parallel -j "$PARALLEL_JOBS" run_config {} || FAILED=$?
else
    # Sequential execution
    for config_file in "${CONFIG_FILES[@]}"; do
        run_config "$config_file" || ((FAILED++))
    done
fi

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

# Generate summary
echo ""
echo "=============================================="
echo "Experiment Complete"
echo "=============================================="
echo "Total time: ${TOTAL_TIME}s"
echo "Configs run: $TOTAL_CONFIGS"
echo "Failed: $FAILED"
echo "Results in: $RESULTS_DIR"

# Create summary file
SUMMARY_FILE="$RUN_DIR/summary.txt"
cat > "$SUMMARY_FILE" << EOF
============================================
Experiment Summary
============================================
Experiment: $EXPERIMENT_DIR
Run ID: $TIMESTAMP
Completed: $(date)
Total Runtime: ${TOTAL_TIME}s
Configs Run: $TOTAL_CONFIGS
Failed: $FAILED

Results:
EOF

# Append LER from each result
for result_file in "$RESULTS_DIR"/*_results.txt; do
    if [[ -f "$result_file" ]]; then
        config_name=$(basename "$result_file" _results.txt)
        ler=$(grep "Logical Error Rate:" "$result_file" | head -1 | awk '{print $4}')
        echo "  $config_name: LER = $ler" >> "$SUMMARY_FILE"
    fi
done

echo ""
cat "$SUMMARY_FILE"

# Exit with error if any configs failed
if [[ $FAILED -gt 0 ]]; then
    exit 1
fi
