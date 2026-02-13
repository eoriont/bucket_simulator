#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${ROOT_DIR}/build"
SIM_BIN="${BUILD_DIR}/bucket_simulator"
INV_TEST_BIN="${BUILD_DIR}/distributed_lattice_surgery_invariants_test"

CONFIG_D3="${ROOT_DIR}/configs/test_distributed_ls_d3.txt"
CONFIG_D5="${ROOT_DIR}/configs/test_distributed_ls_d5.txt"
OUT_DIR="${ROOT_DIR}/output/orchestrated_distributed_ls"

SKIP_RUN=0
if [[ "${1:-}" == "--skip-run" ]]; then
  SKIP_RUN=1
fi

log() {
  printf "[%s] %s\n" "$1" "$2"
}

log "orchestrator" "start distributed lattice-surgery pipeline"

log "build-agent" "configure and build required targets"
cmake -S "${ROOT_DIR}" -B "${BUILD_DIR}"
cmake --build "${BUILD_DIR}" -j8 --target distributed_lattice_surgery_invariants_test bucket_simulator

log "spec-agent" "validate seam-column invariants (d=3,d=5)"
"${INV_TEST_BIN}"
ctest --test-dir "${BUILD_DIR}" -R distributed_lattice_surgery_invariants --output-on-failure

if [[ "${SKIP_RUN}" -eq 1 ]]; then
  log "orchestrator" "tests passed; simulator run skipped by flag"
  exit 0
fi

mkdir -p "${OUT_DIR}"

log "run-agent" "execute distributed d=3 and d=5 smoke runs in parallel"
"${SIM_BIN}" -config "${CONFIG_D3}" -output "${OUT_DIR}" > "${OUT_DIR}/d3.log" 2>&1 &
pid_d3=$!
"${SIM_BIN}" -config "${CONFIG_D5}" -output "${OUT_DIR}" > "${OUT_DIR}/d5.log" 2>&1 &
pid_d5=$!

wait "${pid_d3}"
wait "${pid_d5}"

log "report-agent" "runs complete"
log "report-agent" "logs: ${OUT_DIR}/d3.log and ${OUT_DIR}/d5.log"
