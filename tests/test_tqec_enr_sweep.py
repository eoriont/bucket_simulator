from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "python" / "tqec_import"))
TQEC_PYTHON = PROJECT_ROOT / ".venv-tqec" / "bin" / "python"

from sweep_enr_protocol import validate_deterministic_circuit


class TqecEnrSweepTest(unittest.TestCase):
    def test_default_d5_merge_schedule_drops_post_phase(self) -> None:
        export_script = PROJECT_ROOT / "python" / "tqec_import" / "export_stim.py"
        retime_script = PROJECT_ROOT / "python" / "tqec_import" / "retime_simple_merge.py"

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            source = temp_dir / "base.stim"
            out = temp_dir / "d5_merge_6_6_0.stim"

            subprocess.run(
                [
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
                    str(source),
                ],
                cwd=PROJECT_ROOT,
                check=True,
            )
            subprocess.run(
                [
                    str(TQEC_PYTHON),
                    str(retime_script),
                    "--stim",
                    str(source),
                    "--pre-rounds",
                    "6",
                    "--merge-rounds",
                    "6",
                    "--post-rounds",
                    "0",
                    "--out",
                    str(out),
                ],
                cwd=PROJECT_ROOT,
                check=True,
            )

            validate_deterministic_circuit(out)
            text = out.read_text()
            self.assertNotIn("post", text.lower())

    def test_remote_cnot_noise_changes_with_entanglement_rate(self) -> None:
        stim_path = PROJECT_ROOT / "output" / "tqec_import" / "generated" / "simple_merge_only_x_k2_pre6_merge6_post0.ss_mid.stim"
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
            low_json = temp_dir / "low.json"
            high_json = temp_dir / "high.json"
            low_stim = temp_dir / "low.stim"
            high_stim = temp_dir / "high.stim"

            for enr, out_stim, out_json in [(1_000_000, low_stim, low_json), (50_000_000, high_stim, high_json)]:
                subprocess.run(
                    [
                        str(TQEC_PYTHON),
                        str(PROJECT_ROOT / "python" / "tqec_import" / "inject_noise.py"),
                        "--stim",
                        str(stim_path),
                        "--physical-error",
                        "0.001",
                        "--accurate-rcx",
                        "--distillation-protocol",
                        "none",
                        "--distillation-rounds",
                        "0",
                        "--distillation-backup-batches",
                        "1",
                        "--entanglement-rate",
                        str(enr),
                        "--measurement-delay",
                        "6.6e-07",
                        "--out-stim",
                        str(out_stim),
                        "--out-json",
                        str(out_json),
                    ],
                    cwd=PROJECT_ROOT,
                    check=True,
                    stdout=subprocess.DEVNULL,
                )

            low = json.loads(low_json.read_text())
            high = json.loads(high_json.read_text())
            low_rcx = max(section["remote_cnot_error"] for section in low["sections"])
            high_rcx = max(section["remote_cnot_error"] for section in high["sections"])
            self.assertGreater(low_rcx, high_rcx)


if __name__ == "__main__":
    unittest.main()
