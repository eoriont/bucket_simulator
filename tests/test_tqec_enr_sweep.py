from __future__ import annotations

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


if __name__ == "__main__":
    unittest.main()
