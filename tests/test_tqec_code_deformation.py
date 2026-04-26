from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import stim

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "python" / "tqec_import"))

from stim_deformation import apply_code_deformation


class TqecCodeDeformationTest(unittest.TestCase):
    def test_simple_merge_code_deformation_stays_deterministic(self) -> None:
        circuit = stim.Circuit.from_file(
            PROJECT_ROOT / "output" / "tqec_import" / "generated" / "simple_merge_x_k1_r5.stim"
        )

        result = apply_code_deformation(circuit, [(5.0, 7.0)])
        self.assertEqual(result.removed_data_qubits, frozenset({41}))
        self.assertEqual(result.suppressed_ancillas, frozenset({33, 48}))
        self.assertLess(result.circuit.num_detectors, circuit.num_detectors)
        self.assertEqual(result.dropped_detector_coords, ())

        dem = result.circuit.detector_error_model(allow_gauge_detectors=False)
        self.assertGreater(dem.num_detectors, 0)

    def test_horizontal_seam_edge_deformation_keeps_same_type_weight_three_check(self) -> None:
        circuit = stim.Circuit.from_file(
            PROJECT_ROOT / "output" / "tqec_import" / "generated" / "simple_merge_x_k2_r5.stim"
        )

        # Removing (1,11) collapses the Z-type border check at (0,12) to weight-1
        # (auto-suppressed) and shrinks the two neighboring weight-4 checks to
        # weight-3: the Z-type at (2,10) commutes with its same-type neighbor and
        # is kept, while the X-type at (2,12) anti-commutes and must be suppressed.
        result = apply_code_deformation(circuit, [(1.0, 11.0)])
        self.assertEqual(result.removed_data_qubits, frozenset({17}))
        self.assertEqual(result.suppressed_ancillas, frozenset({6, 29}))
        self.assertEqual(result.dropped_detector_coords, ())

        transformed = result.circuit.flattened()
        coords = {}
        for op in transformed:
            if op.name == "QUBIT_COORDS":
                coords[op.targets_copy()[0].qubit_value] = tuple(op.gate_args_copy()[:2])

        supports = {}
        for op in transformed:
            if op.name not in {"CX", "CZ"}:
                continue
            targets = [t.qubit_value for t in op.targets_copy() if t.is_qubit_target]
            for i in range(0, len(targets), 2):
                a = targets[i]
                b = targets[i + 1]
                if a in coords and coords[a] == (2.0, 10.0) and b in coords:
                    supports.setdefault(coords[a], set()).add(coords[b])
                if b in coords and coords[b] == (2.0, 10.0) and a in coords:
                    supports.setdefault(coords[b], set()).add(coords[a])

        self.assertEqual(len(supports[(2.0, 10.0)]), 3)
        self.assertNotIn((2.0, 12.0), set(coords.values()))
        circuit_dem = result.circuit.detector_error_model(allow_gauge_detectors=False)
        self.assertGreater(circuit_dem.num_detectors, 0)

    def test_horizontal_seam_right_and_two_side_deformation(self) -> None:
        circuit = stim.Circuit.from_file(
            PROJECT_ROOT / "output" / "tqec_import" / "generated" / "simple_merge_x_k2_r5.stim"
        )

        right = apply_code_deformation(circuit, [(9.0, 11.0)])
        self.assertEqual(right.removed_data_qubits, frozenset({109}))
        self.assertEqual(right.suppressed_ancillas, frozenset({97, 120}))
        self.assertEqual(right.dropped_detector_coords, ())
        self.assertGreater(right.circuit.detector_error_model(allow_gauge_detectors=False).num_detectors, 0)

        both = apply_code_deformation(circuit, [(1.0, 11.0), (9.0, 11.0)])
        self.assertEqual(both.removed_data_qubits, frozenset({17, 109}))
        self.assertEqual(both.suppressed_ancillas, frozenset({6, 29, 97, 120}))
        self.assertEqual(both.dropped_detector_coords, ())
        self.assertGreater(both.circuit.detector_error_model(allow_gauge_detectors=False).num_detectors, 0)

    def test_ss_mid_interior_deformation_keeps_one_logical_observable(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test_ss_mid.stim"
            subprocess.run(
                [
                    str(PROJECT_ROOT / ".venv-tqec" / "bin" / "python"),
                    str(PROJECT_ROOT / "python" / "tqec_import" / "export_stim.py"),
                    "--source",
                    "simple_merge_only",
                    "--basis",
                    "X",
                    "--k",
                    "2",
                    "--syndrome-rounds",
                    "6",
                    "--remove-data",
                    "(5,11)",
                    "--out",
                    str(output_path),
                ],
                cwd=PROJECT_ROOT,
                check=True,
            )

            circuit = stim.Circuit.from_file(output_path)
            self.assertEqual(circuit.num_observables, 1)
            self.assertIn("OBSERVABLE_INCLUDE", output_path.read_text())
            dem = circuit.detector_error_model(allow_gauge_detectors=False)
            self.assertGreater(dem.num_detectors, 0)


if __name__ == "__main__":
    unittest.main()
