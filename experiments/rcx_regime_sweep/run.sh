#!/bin/bash
# Run the RCX regime sweep experiment.
#
# Usage:
#   ./run.sh [python args...]
#
# Example:
#   ./run.sh --output-dir runs/$(date +%Y%m%d_%H%M%S)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ ! -x ../../build/bucket_simulator ]]; then
    echo "Error: simulator not found at ../../build/bucket_simulator"
    echo "Build first: cd ../../build && cmake .. && make"
    exit 1
fi

python3 ./sweep_optimal_rcx.py "$@"
