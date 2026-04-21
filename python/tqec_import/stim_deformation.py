from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import re
from typing import Iterable

import stim

from stim_layout import ancilla_interactions, repeated_measure_qubits


Pair = tuple[float, float]

_PAIR_GATES = {"CX", "CZ"}
_MEASUREMENT_GATES = {"M", "MX", "MY", "MR", "MRX", "MRY"}
_SINGLE_QUBIT_GATES = {
    "R",
    "RX",
    "RY",
    "H",
    "H_XY",
    "H_XZ",
    "H_YZ",
    "S",
    "SQRT_X",
    "SQRT_Y",
    "SQRT_Z",
    "X",
    "Y",
    "Z",
}


@dataclass(frozen=True)
class CodeDeformationResult:
    circuit: stim.Circuit
    removed_data_qubits: frozenset[int]
    suppressed_ancillas: frozenset[int]
    dropped_detector_coords: tuple[tuple[float, ...], ...] = ()


def parse_coordinate_list(text: str) -> list[Pair]:
    coords: list[Pair] = []
    for part in text.split():
        chunk = part.strip()
        if not chunk:
            continue
        if chunk[0] == "(" and chunk[-1] == ")":
            chunk = chunk[1:-1]
        x_text, y_text = [piece.strip() for piece in chunk.split(",", 1)]
        coords.append((float(x_text), float(y_text)))
    return coords


def apply_code_deformation(circuit: stim.Circuit, removed_coords: Iterable[Pair]) -> CodeDeformationResult:
    flat = circuit.flattened()
    coord_to_qubit, qubit_to_coord = _build_coordinate_maps(flat)
    removed_data_qubits = _resolve_removed_data_qubits(removed_coords, coord_to_qubit)
    if not removed_data_qubits:
        return CodeDeformationResult(circuit=flat, removed_data_qubits=frozenset(), suppressed_ancillas=frozenset())

    ancilla_candidates = repeated_measure_qubits(flat)
    suppressed_ancillas = _find_suppressed_ancillas(flat, removed_data_qubits, ancilla_candidates)
    suppressed_coord_prefixes = {qubit_to_coord[a] for a in suppressed_ancillas if a in qubit_to_coord}
    transformed = _rewrite_circuit(flat, removed_data_qubits | suppressed_ancillas, suppressed_coord_prefixes)
    transformed, dropped_coords = _drop_nondeterministic_detectors(transformed)

    return CodeDeformationResult(
        circuit=transformed,
        removed_data_qubits=frozenset(removed_data_qubits),
        suppressed_ancillas=frozenset(suppressed_ancillas),
        dropped_detector_coords=dropped_coords,
    )


def _build_coordinate_maps(circuit: stim.Circuit) -> tuple[dict[Pair, int], dict[int, Pair]]:
    coord_to_qubit: dict[Pair, int] = {}
    qubit_to_coord: dict[int, Pair] = {}
    for op in circuit:
        if op.name != "QUBIT_COORDS":
            continue
        args = op.gate_args_copy()
        if len(args) < 2:
            continue
        target = op.targets_copy()[0]
        if not target.is_qubit_target:
            continue
        pair = (float(args[0]), float(args[1]))
        qubit = target.qubit_value
        coord_to_qubit[pair] = qubit
        qubit_to_coord[qubit] = pair
    return coord_to_qubit, qubit_to_coord


def _resolve_removed_data_qubits(removed_coords: Iterable[Pair], coord_to_qubit: dict[Pair, int]) -> set[int]:
    missing: list[Pair] = []
    removed: set[int] = set()
    for coord in removed_coords:
        qubit = coord_to_qubit.get((float(coord[0]), float(coord[1])))
        if qubit is None:
            missing.append(coord)
            continue
        removed.add(qubit)
    if missing:
        formatted = ", ".join(f"({x},{y})" for x, y in missing)
        raise ValueError(f"No TQEC qubits found at deformation coordinates: {formatted}")
    return removed


def _find_suppressed_ancillas(circuit: stim.Circuit, removed_data_qubits: set[int], ancilla_candidates: set[int]) -> set[int]:
    families, touched = ancilla_interactions(circuit, ancilla_candidates)
    suppressed: set[int] = set()

    affected: list[tuple[int, set[int], set[int], str | None]] = []
    for ancilla in ancilla_candidates:
        support = touched.get(ancilla, set())
        overlap = support & removed_data_qubits
        if not overlap:
            continue
        remaining = support - removed_data_qubits
        affected.append((ancilla, remaining, overlap, families.get(ancilla)))

    # Boundary checks that collapse to weight ≤ 1 are removed. Their family tells
    # us the local border stabilizer type at each removed data qubit.
    border_type_by_data: dict[int, set[str]] = {}
    for ancilla, remaining, overlap, family in affected:
        if len(remaining) <= 1:
            suppressed.add(ancilla)
            if family is not None:
                for data in overlap:
                    border_type_by_data.setdefault(data, set()).add(family)

    # A weight-4 check dropping to weight-3 of the OPPOSITE family from the border
    # now shares a single data qubit with its same-direction neighbor (which also
    # dropped to weight-3), so the pair anti-commutes during measurement. Suppress
    # the opposite-type one — matches the C++ gauge-layering fix for this case.
    for ancilla, remaining, overlap, family in affected:
        if ancilla in suppressed or family is None or len(remaining) != 3:
            continue
        for data in overlap:
            border_types = border_type_by_data.get(data)
            if border_types and family not in border_types:
                suppressed.add(ancilla)
                break

    return suppressed


def _rewrite_circuit(
    circuit: stim.Circuit,
    removed_qubits: set[int],
    suppressed_detector_coord_prefixes: set[Pair] = frozenset(),
) -> stim.Circuit:
    out = stim.Circuit()
    old_to_new: list[int | None] = []
    old_measure_count = 0
    new_measure_count = 0

    for op in circuit:
        name = op.name
        targets = op.targets_copy()
        args = op.gate_args_copy()

        if name == "QUBIT_COORDS":
            kept = [t for t in targets if not (t.is_qubit_target and t.qubit_value in removed_qubits)]
            if kept:
                out.append(name, kept, args)
            continue

        if name in _PAIR_GATES:
            kept: list[stim.GateTarget] = []
            qubit_targets = [t for t in targets if t.is_qubit_target]
            for i in range(0, len(qubit_targets), 2):
                a = qubit_targets[i]
                b = qubit_targets[i + 1]
                if a.qubit_value in removed_qubits or b.qubit_value in removed_qubits:
                    continue
                kept.extend([a, b])
            if kept:
                out.append(name, kept, args)
            continue

        if name in _MEASUREMENT_GATES:
            kept: list[stim.GateTarget] = []
            for target in targets:
                if not target.is_qubit_target:
                    continue
                if target.qubit_value in removed_qubits:
                    old_to_new.append(None)
                else:
                    kept.append(target)
                    old_to_new.append(new_measure_count)
                    new_measure_count += 1
                old_measure_count += 1
            if kept:
                out.append(name, kept, args)
            continue

        if name in _SINGLE_QUBIT_GATES:
            kept = [t for t in targets if t.is_qubit_target and t.qubit_value not in removed_qubits]
            if kept:
                out.append(name, kept, args)
            continue

        if name == "DETECTOR":
            if len(args) >= 2:
                prefix = (float(args[0]), float(args[1]))
                if prefix in suppressed_detector_coord_prefixes:
                    continue
            remapped = _remap_record_targets(targets, old_to_new, old_measure_count, new_measure_count)
            if remapped:
                out.append(name, remapped, args)
            continue

        if name == "OBSERVABLE_INCLUDE":
            remapped = _remap_record_targets(targets, old_to_new, old_measure_count, new_measure_count)
            if remapped:
                out.append(name, remapped, args)
            continue

        out.append(op)

    return out


def _remap_record_targets(
    targets: list[stim.GateTarget],
    old_to_new: list[int | None],
    old_measure_count: int,
    new_measure_count: int,
) -> list[stim.GateTarget]:
    remapped: list[stim.GateTarget] = []
    for target in targets:
        if not target.is_measurement_record_target:
            remapped.append(target)
            continue
        old_index = old_measure_count + target.value
        if old_index < 0 or old_index >= len(old_to_new):
            raise ValueError(f"Measurement record target {target} points outside the seen measurement range.")
        new_index = old_to_new[old_index]
        if new_index is None:
            continue
        remapped.append(stim.target_rec(new_index - new_measure_count))
    return remapped


def _drop_nondeterministic_detectors(
    circuit: stim.Circuit,
    max_passes: int = 32,
) -> tuple[stim.Circuit, tuple[tuple[float, ...], ...]]:
    """Remove detectors flagged as nondeterministic by stim's flow analysis.

    After a code deformation, certain transition-round detectors become gauge
    operators — typically those referencing an ancilla whose stabilizer support
    changed because one of its data neighbors was removed. Stim only surfaces
    one failure at a time (it stops at the first collapse that anti-commutes
    with a detector), so we iterate until the circuit is clean or give up.
    Returns the cleaned circuit and the coordinates of every dropped detector.
    """
    dropped_coords: list[tuple[float, ...]] = []
    current = circuit
    for _ in range(max_passes):
        try:
            current.detector_error_model(allow_gauge_detectors=False)
            return current, tuple(dropped_coords)
        except ValueError as exc:
            detector_ids = sorted({int(value) for value in re.findall(r"D(\d+)", str(exc))})
            if not detector_ids:
                raise
            target_ids = set(detector_ids)
            seen = 0
            for op in current.flattened():
                if op.name == "DETECTOR":
                    if seen in target_ids:
                        dropped_coords.append(tuple(op.gate_args_copy()))
                    seen += 1
            current = _remove_detector_ids(current, detector_ids)
    raise ValueError("Exceeded detector-pruning passes while cleaning a deformed Stim circuit.")


def _remove_detector_ids(circuit: stim.Circuit, detector_ids: list[int]) -> stim.Circuit:
    to_remove = set(detector_ids)
    out = stim.Circuit()
    detector_index = 0
    for op in circuit:
        if op.name == "DETECTOR":
            if detector_index not in to_remove:
                out.append(op)
            detector_index += 1
            continue
        out.append(op)
    return out
