#!/bin/bash
# Quick runner for radar_sweep experiment
#
# Usage:
#   ./run.sh              # Generate configs and run with defaults
#   ./run.sh --generate   # Only generate configs
#   ./run.sh --run        # Only run (configs must exist)
#   ./run.sh -n 8         # Run with 8 MPI processes

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

GENERATE=true
RUN=true
EXTRA_ARGS=""

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --generate)
            RUN=false
            shift
            ;;
        --run)
            GENERATE=false
            shift
            ;;
        *)
            EXTRA_ARGS="$EXTRA_ARGS $1"
            shift
            ;;
    esac
done

# Generate configs
if [[ "$GENERATE" == "true" ]]; then
    echo "Generating configs..."
    python3 ../scripts/generate_configs.py experiment.json configs/
    echo ""
fi

# Run experiment
if [[ "$RUN" == "true" ]]; then
    echo "Running experiment..."
    ../scripts/run_experiment.sh radar_sweep $EXTRA_ARGS
fi
