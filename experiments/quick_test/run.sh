#!/bin/bash
# Quick test runner
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Generating configs..."
python3 ../scripts/generate_configs.py experiment.json configs/

echo ""
echo "Running experiment..."
../scripts/run_experiment.sh quick_test "$@"
