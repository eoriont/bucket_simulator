from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import stim

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "python" / "tqec_import"))
TQEC_PYTHON = PROJECT_ROOT / ".venv-tqec" / "bin" / "python"

from distributed_noise import NoiseConfig, inject_network_noise
from sweep_enr_raw_pairs import compute_required_epr_pairs_per_round, normalize_deformation_mode, resolve_max_concurrent_jobs


class TqecRawPairsSweepTest(unittest.TestCase):
    def test_baseline_alias_and_auto_concurrency_helpers(self) -> None:
        self.assertEqual(normalize_deformation_mode("nosuper"), "none")
        self.assertEqual(normalize_deformation_mode("none"), "none")
        self.assertEqual(normalize_deformation_mode("lside"), "lside")
        self.assertGreaterEqual(resolve_max_concurrent_jobs(0, 8), 1)
        self.assertEqual(resolve_max_concurrent_jobs(3, 8), 3)

    def test_required_epr_pairs_track_border_deformation(self) -> None:
        export_script = PROJECT_ROOT / "python" / "tqec_import" / "export_stim.py"
        retime_script = PROJECT_ROOT / "python" / "tqec_import" / "retime_simple_merge.py"

        counts: dict[str, int] = {}
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            for mode in ("none", "lside", "rside", "twoside"):
                clean = temp_dir / f"{mode}.stim"
                export_cmd = [
                    str(TQEC_PYTHON),
                    str(export_script),
                    "--source",
                    "simple_merge_only",
                    "--basis",
                    "X",
                    "--k",
                    "2",
                    "--syndrome-rounds",
                    "6",
                    "--out",
                    str(clean),
                ]
                if mode != "none":
                    export_cmd += ["--deformation-mode", mode]
                subprocess.run(export_cmd, cwd=PROJECT_ROOT, check=True)
                if mode == "none":
                    retimed = temp_dir / f"{mode}.retimed.stim"
                    subprocess.run(
                        [
                            str(TQEC_PYTHON),
                            str(retime_script),
                            "--stim",
                            str(clean),
                            "--pre-rounds",
                            "6",
                            "--merge-rounds",
                            "6",
                            "--post-rounds",
                            "0",
                            "--out",
                            str(retimed),
                        ],
                        cwd=PROJECT_ROOT,
                        check=True,
                    )
                else:
                    retimed = clean

                circuit = stim.Circuit.from_file(retimed)
                injected_circuit, summary = inject_network_noise(
                    circuit,
                    NoiseConfig(
                        distillation_protocol="none",
                        distillation_rounds=0,
                        distillation_backup_batches=1,
                        entanglement_rate=1e9,
                        measurement_delay=660e-9,
                    ),
                )
                counts[mode] = compute_required_epr_pairs_per_round(summary, 6)
                self.assertGreater(str(injected_circuit).count("PAULI_CHANNEL_1"), 0)
                self.assertEqual(counts[mode], summary["required_epr_pairs_per_round"])

        self.assertGreater(len(set(counts.values())), 1)


if __name__ == "__main__":
    unittest.main()
