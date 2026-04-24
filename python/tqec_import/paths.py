from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "output" / "tqec_import"
GENERATED_ROOT = OUTPUT_ROOT / "generated"
INJECTED_ROOT = OUTPUT_ROOT / "injected"
LER_ROOT = OUTPUT_ROOT / "ler"


@dataclass(frozen=True)
class ExperimentRunPaths:
    root: Path
    circuits_dir: Path
    clean_circuits_dir: Path
    noisy_circuits_dir: Path
    timing_dir: Path
    results_dir: Path
    params_dir: Path

    def create(self) -> "ExperimentRunPaths":
        for path in (
            self.root,
            self.circuits_dir,
            self.clean_circuits_dir,
            self.noisy_circuits_dir,
            self.timing_dir,
            self.results_dir,
            self.params_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)
        return self


def create_experiment_run_paths(experiment_name: str, run_name: str) -> ExperimentRunPaths:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = OUTPUT_ROOT / experiment_name / f"{run_name}_{timestamp}"
    circuits_dir = root / "circuits"
    return ExperimentRunPaths(
        root=root,
        circuits_dir=circuits_dir,
        clean_circuits_dir=circuits_dir / "clean",
        noisy_circuits_dir=circuits_dir / "noisy",
        timing_dir=root / "timing",
        results_dir=root / "results",
        params_dir=root / "params",
    ).create()
