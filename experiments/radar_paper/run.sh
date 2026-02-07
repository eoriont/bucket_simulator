#!/bin/bash
#
# RADar Paper Reproduction Experiment
# ====================================
# Reproduces key results from "RADar: Resource-Aware Entanglement-Distillation
# for Distributed Surface Code" (ISCA 2026)
#
# Sub-experiments:
#   1. monolithic_baseline  - Monolithic baselines (Figure 4 reference)
#   2. table1_low_rate      - Table I: Low rate regime (80-100 MHz)
#   3. experiment           - Table II: High rate regime (150-500 MHz)
#   4. rate_sweep           - Figure 6 style rate sweep
#
# Usage:
#   ./run.sh                    # Run all sub-experiments
#   ./run.sh --generate         # Only generate configs
#   ./run.sh --sub <name>       # Run specific sub-experiment
#   ./run.sh -n 8               # Use 8 MPI processes
#
# Estimated configs: 3 + 18 + 105 + 9 = 135 total
# Estimated runtime: Several hours with 10M shots each

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Sub-experiments to run
SUB_EXPERIMENTS=("monolithic_baseline" "table1_low_rate" "experiment" "rate_sweep")

GENERATE_ONLY=false
RUN_ONLY=false
SPECIFIC_SUB=""
EXTRA_ARGS=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --generate)
            GENERATE_ONLY=true
            shift
            ;;
        --run)
            RUN_ONLY=true
            shift
            ;;
        --sub)
            SPECIFIC_SUB="$2"
            shift 2
            ;;
        --list)
            echo "Available sub-experiments:"
            for sub in "${SUB_EXPERIMENTS[@]}"; do
                echo "  - $sub"
            done
            exit 0
            ;;
        -h|--help)
            head -25 "$0" | tail -23
            exit 0
            ;;
        *)
            EXTRA_ARGS="$EXTRA_ARGS $1"
            shift
            ;;
    esac
done

# Filter to specific sub-experiment if requested
if [[ -n "$SPECIFIC_SUB" ]]; then
    SUB_EXPERIMENTS=("$SPECIFIC_SUB")
fi

echo "=============================================="
echo "RADar Paper Reproduction Experiment"
echo "=============================================="
echo ""

# Generate configs for each sub-experiment
if [[ "$RUN_ONLY" != "true" ]]; then
    echo "Generating configurations..."
    echo ""

    TOTAL_CONFIGS=0
    for sub in "${SUB_EXPERIMENTS[@]}"; do
        spec_file="${sub}.json"
        if [[ ! -f "$spec_file" ]]; then
            echo "Warning: $spec_file not found, skipping"
            continue
        fi

        output_dir="configs/${sub}"
        mkdir -p "$output_dir"

        echo "  [$sub]"
        python3 ../scripts/generate_configs.py "$spec_file" "$output_dir" 2>&1 | sed 's/^/    /'

        # Count configs
        count=$(find "$output_dir" -name "*.txt" | wc -l | tr -d ' ')
        TOTAL_CONFIGS=$((TOTAL_CONFIGS + count))
        echo ""
    done

    echo "Total configs generated: $TOTAL_CONFIGS"
    echo ""
fi

# Run experiments
if [[ "$GENERATE_ONLY" != "true" ]]; then
    echo "Running experiments..."
    echo ""

    # Create unified run directory
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    RUN_DIR="runs/$TIMESTAMP"
    mkdir -p "$RUN_DIR"

    # Save experiment metadata
    cat > "$RUN_DIR/experiment_info.txt" << EOF
RADar Paper Reproduction
========================
Started: $(date)
Sub-experiments: ${SUB_EXPERIMENTS[*]}
MPI args: $EXTRA_ARGS

Paper Reference:
  Title: RADar: Resource-Aware Entanglement-Distillation for Distributed Surface Code
  Venue: ISCA 2026

Key Parameters (from paper Section III.B):
  Physical error rate: p = 10^-3
  T1 coherence time: 125 μs
  T2 coherence time: 200 μs
  Measurement delay: 660 ns
  Raw EPR fidelity: 0.99
EOF

    # Run each sub-experiment
    for sub in "${SUB_EXPERIMENTS[@]}"; do
        config_dir="configs/${sub}"
        if [[ ! -d "$config_dir" ]]; then
            echo "Warning: $config_dir not found, skipping"
            continue
        fi

        echo "=============================================="
        echo "Running sub-experiment: $sub"
        echo "=============================================="

        # Create output directories for this sub-experiment
        sub_results="$RUN_DIR/results/${sub}"
        sub_logs="$RUN_DIR/logs/${sub}"
        mkdir -p "$sub_results" "$sub_logs"

        # Run each config in this sub-experiment
        for config_file in "$config_dir"/*.txt; do
            [[ -f "$config_file" ]] || continue
            [[ "$(basename "$config_file")" == "manifest.json" ]] && continue

            config_name=$(basename "$config_file" .txt)
            log_file="$sub_logs/${config_name}.log"

            echo -n "  [$config_name] "

            start_time=$(date +%s)
            if mpirun $EXTRA_ARGS ../../build/bucket_simulator \
                -config "$config_file" \
                -output "$sub_results" > "$log_file" 2>&1; then

                end_time=$(date +%s)
                duration=$((end_time - start_time))
                echo "done (${duration}s)"

                # Rename result file
                latest=$(ls -t "$sub_results"/results_*.txt 2>/dev/null | head -1)
                if [[ -n "$latest" ]]; then
                    mv "$latest" "$sub_results/${config_name}_results.txt"
                fi
            else
                echo "FAILED (see $log_file)"
            fi
        done
        echo ""
    done

    # Generate summary
    echo "=============================================="
    echo "Generating Summary"
    echo "=============================================="

    python3 ../scripts/analyze_results.py "$RUN_DIR" --csv 2>&1 || true

    echo ""
    echo "=============================================="
    echo "Experiment Complete"
    echo "=============================================="
    echo "Results in: $RUN_DIR"
    echo ""
    echo "To analyze results:"
    echo "  python3 ../scripts/analyze_results.py $RUN_DIR --csv --plot"
fi
