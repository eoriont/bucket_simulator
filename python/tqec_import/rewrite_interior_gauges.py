#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import stim
from stim_deformation import parse_coordinate_list
from stim_layout import ancilla_interactions, qubit_coords, repeated_measure_qubits


@dataclass(frozen=True)
class RewriteSummary:
    removed_data_qubits: tuple[int, ...]
    removed_coords: tuple[tuple[float, float], ...]
    gauge_x_ancillas: tuple[int, ...]
    gauge_z_ancillas: tuple[int, ...]
    gauge_x_clusters: tuple[tuple[int, ...], ...]
    gauge_z_clusters: tuple[tuple[int, ...], ...]
    x_entry_transition_pre_ancillas: tuple[int, ...]
    z_entry_transition_pre_ancillas: tuple[int, ...]
    original_repeat_count: int
    rewritten_repeat_count: int
    original_measurements_per_round: int
    half1_measurements: int
    half2_measurements: int


_PAIR_GATES = {"CX", "CZ"}
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
_MEASUREMENT_GATES = {"M", "MX", "MY", "MR", "MRX", "MRY"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prototype rewrite for interior seam gauge superstabilizers. "
            "This rewrites the merge REPEAT body into two subrounds per logical "
            "round, halves the repeat count, and rebuilds detector history for "
            "the split merge schedule."
        )
    )
    parser.add_argument("--stim", type=Path, required=True, help="Input .stim circuit.")
    parser.add_argument("--out", type=Path, required=True, help="Output .stim circuit.")
    parser.add_argument(
        "--remove-data",
        type=str,
        required=True,
        help='Whitespace-separated coordinate list such as "(5,9) (5,11)".',
    )
    return parser.parse_args()


def _resolve_removed_data_qubits(
    circuit: stim.Circuit,
    removed_coords: Iterable[tuple[float, float]],
) -> tuple[int, ...]:
    flat = circuit.flattened()
    coord_to_qubit = {(float(x), float(y)): q for q, (x, y) in qubit_coords(flat).items()}
    missing: list[tuple[float, float]] = []
    removed: list[int] = []
    for coord in removed_coords:
        qubit = coord_to_qubit.get((float(coord[0]), float(coord[1])))
        if qubit is None:
            missing.append(coord)
            continue
        removed.append(qubit)
    if missing:
        formatted = ", ".join(f"({x},{y})" for x, y in missing)
        raise ValueError(f"No TQEC qubits found at deformation coordinates: {formatted}")
    return tuple(sorted(set(removed)))


def _find_merge_repeat(circuit: stim.Circuit) -> tuple[int, stim.CircuitRepeatBlock]:
    repeats: list[tuple[int, stim.CircuitRepeatBlock]] = []
    for i, op in enumerate(circuit):
        if isinstance(op, stim.CircuitRepeatBlock):
            repeats.append((i, op))
    if not repeats:
        raise ValueError("Expected at least one top-level REPEAT block; found none.")
    if len(repeats) == 1:
        return repeats[0]
    # Retimed simple-merge circuits use top-level REPEAT blocks in
    # pre-merge, merge, post-merge order. The interior gauge rewrite must
    # target the merge body, which is the second top-level repeat.
    return repeats[1]


def _find_preceding_seed_block_start(circuit: stim.Circuit, repeat_index: int) -> int:
    """Returns the index of the explicit seeded round preceding a merge repeat.

    TQEC-exported simple-merge circuits emit one explicit round before the
    repeated body of each phase. That explicit seed begins immediately after the
    previous top-level SHIFT_COORDS.
    """
    shift_indices: list[int] = []
    for i in range(repeat_index - 1, -1, -1):
        op = circuit[i]
        if isinstance(op, stim.CircuitRepeatBlock):
            continue
        if op.name == "SHIFT_COORDS":
            shift_indices.append(i)
            if len(shift_indices) == 2:
                return i + 1
    return 0 if not shift_indices else shift_indices[-1] + 1


def _find_merge_phase_end(circuit: stim.Circuit, repeat_index: int) -> int:
    """Returns the index of the top-level SHIFT_COORDS ending the merge phase."""
    for i in range(repeat_index + 1, len(circuit)):
        op = circuit[i]
        if isinstance(op, stim.CircuitRepeatBlock):
            continue
        if op.name == "SHIFT_COORDS":
            return i
    return repeat_index


def _measurement_targets(body: stim.Circuit) -> list[int]:
    measured: list[int] = []
    for op in body:
        if op.name not in _MEASUREMENT_GATES:
            continue
        measured.extend(t.qubit_value for t in op.targets_copy() if t.is_qubit_target)
    return measured


def _is_interior_data_coordinate(circuit: stim.Circuit, coord: tuple[float, float]) -> bool:
    flat = circuit.flattened()
    coords = qubit_coords(flat)
    ancillas = repeated_measure_qubits(flat)
    data_coords = [xy for qubit, xy in coords.items() if qubit not in ancillas]
    xs = [x for x, _ in data_coords]
    ys = [y for _, y in data_coords]
    x, y = coord
    return min(xs) < x < max(xs) and min(ys) < y < max(ys)


def supports_interior_gauge_rewrite(
    circuit: stim.Circuit,
    removed_coords: Iterable[tuple[float, float]],
) -> bool:
    coords = [tuple(map(float, coord)) for coord in removed_coords]
    return bool(coords) and all(_is_interior_data_coordinate(circuit, coord) for coord in coords)


def _active_measurements(
    measured_qubits: list[int],
    gauge_x: set[int],
    gauge_z: set[int],
    active_family: str,
) -> list[int]:
    out: list[int] = []
    for qubit in measured_qubits:
        if qubit in gauge_x:
            if active_family == "X":
                out.append(qubit)
            continue
        if qubit in gauge_z:
            if active_family == "Z":
                out.append(qubit)
            continue
        out.append(qubit)
    return out


def _connected_removed_clusters(
    coords_by_qubit: dict[int, tuple[float, float]],
    removed_data_qubits: set[int],
) -> list[set[int]]:
    remaining = set(removed_data_qubits)
    clusters: list[set[int]] = []
    while remaining:
        root = remaining.pop()
        cluster = {root}
        stack = [root]
        while stack:
            current = stack.pop()
            cx, cy = coords_by_qubit[current]
            neighbors = [
                other
                for other in list(remaining)
                if abs(coords_by_qubit[other][0] - cx) <= 2.0 and abs(coords_by_qubit[other][1] - cy) <= 2.0
            ]
            for other in neighbors:
                remaining.remove(other)
                cluster.add(other)
                stack.append(other)
        clusters.append(cluster)
    return clusters


def _gauge_clusters_for_family(
    touched: dict[int, set[int]],
    gauge_ancillas: set[int],
    removed_clusters: list[set[int]],
) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for cluster in removed_clusters:
        members = sorted(
            ancilla
            for ancilla in gauge_ancillas
            if touched.get(ancilla, set()) & cluster
        )
        if members:
            out.append(tuple(members))
    return out


def _cluster_center(
    cluster: set[int],
    coords_by_qubit: dict[int, tuple[float, float]],
) -> tuple[float, float]:
    xs = [coords_by_qubit[q][0] for q in cluster]
    ys = [coords_by_qubit[q][1] for q in cluster]
    return (sum(xs) / len(xs), sum(ys) / len(ys))


def _xor_support(acc: set[int], nxt: Iterable[int]) -> set[int]:
    out = set(acc)
    for q in nxt:
        if q in out:
            out.remove(q)
        else:
            out.add(q)
    return out


def _find_entry_transition_subset(
    *,
    family: str,
    ancillas_by_family: dict[str, set[int]],
    touched: dict[int, set[int]],
    removed_data_qubits: set[int],
    gauge_ancillas: set[int],
) -> tuple[int, ...]:
    if not gauge_ancillas:
        return ()

    merge_support_xor: set[int] = set()
    for anc in gauge_ancillas:
        merge_support_xor = _xor_support(merge_support_xor, touched.get(anc, set()) - removed_data_qubits)
    if not merge_support_xor:
        return ()

    family_ancillas = ancillas_by_family.get(family, set())
    pre_candidates: list[int] = []
    for anc in sorted(family_ancillas):
        if anc in gauge_ancillas:
            continue
        support = touched.get(anc, set())
        if not support:
            continue
        if (support & merge_support_xor) or (support & removed_data_qubits):
            pre_candidates.append(anc)

    max_subset = min(4, len(pre_candidates))
    for subset_size in range(1, max_subset + 1):
        picks = list(range(subset_size))
        while True:
            support_xor: set[int] = set()
            ancs: list[int] = []
            for pick in picks:
                anc = pre_candidates[pick]
                ancs.append(anc)
                support_xor = _xor_support(support_xor, touched.get(anc, set()))
            if support_xor == merge_support_xor:
                return tuple(ancs)

            pos = subset_size - 1
            while pos >= 0 and picks[pos] == len(pre_candidates) - subset_size + pos:
                pos -= 1
            if pos < 0:
                break
            picks[pos] += 1
            for k in range(pos + 1, subset_size):
                picks[k] = picks[k - 1] + 1

    return ()


def _append_round_detectors(
    out: stim.Circuit,
    *,
    measured_qubits: list[int],
    previous_normal_measured_qubits: list[int] | None,
    prev_normal_intermediate: int,
    previous_gauge_measured_qubits: list[int] | None,
    prev_gauge_intermediate: int,
    coords_by_qubit: dict[int, tuple[float, float]],
    normal_ancillas: Iterable[int],
    gauge_clusters: Iterable[tuple[int, ...]],
    gauge_cluster_centers: Iterable[tuple[float, float]],
    first_normal: bool = False,
    first_gauge: bool = False,
) -> None:
    measured_index = {qubit: i for i, qubit in enumerate(measured_qubits)}
    previous_normal_index = {qubit: i for i, qubit in enumerate(previous_normal_measured_qubits or [])}
    previous_gauge_index = {qubit: i for i, qubit in enumerate(previous_gauge_measured_qubits or [])}
    n_curr = len(measured_qubits)

    for ancilla in normal_ancillas:
        idx = measured_index.get(ancilla)
        if idx is None:
            continue
        curr = -(n_curr - idx)
        x, y = coords_by_qubit[ancilla]
        if first_normal:
            continue
        prev_idx = previous_normal_index.get(ancilla)
        if prev_idx is None:
            continue
        prev_block_len = len(previous_normal_measured_qubits or [])
        prev = -(n_curr + prev_normal_intermediate + (prev_block_len - prev_idx))
        out.append("DETECTOR", [stim.target_rec(curr), stim.target_rec(prev)], [x, y, 0.0])

    for cluster, center in zip(gauge_clusters, gauge_cluster_centers):
        if first_gauge:
            continue
        tgts: list[stim.GateTarget] = []
        for ancilla in cluster:
            idx = measured_index.get(ancilla)
            if idx is None:
                continue
            curr = -(n_curr - idx)
            tgts.append(stim.target_rec(curr))
            prev_idx = previous_gauge_index.get(ancilla)
            if prev_idx is None:
                continue
            prev_block_len = len(previous_gauge_measured_qubits or [])
            prev = -(n_curr + prev_gauge_intermediate + (prev_block_len - prev_idx))
            tgts.append(stim.target_rec(prev))
        if tgts:
            x, y = center
            out.append("DETECTOR", tgts, [x, y, 0.0])


def _append_half_round_detectors(
    out: stim.Circuit,
    *,
    measured_qubits: list[int],
    previous_normal_measured_qubits: list[int] | None,
    previous_gauge_measured_qubits: list[int] | None,
    coords_by_qubit: dict[int, tuple[float, float]],
    normal_ancillas: Iterable[int],
    gauge_clusters: Iterable[tuple[int, ...]],
    gauge_cluster_centers: Iterable[tuple[float, float]],
    family: str,
    half1_measurements: int,
    half2_measurements: int,
    original_measurements_per_round: int,
    first: bool = False,
) -> None:
    measured_index = {qubit: i for i, qubit in enumerate(measured_qubits)}
    previous_normal_index = {qubit: i for i, qubit in enumerate(previous_normal_measured_qubits or [])}
    previous_gauge_index = {qubit: i for i, qubit in enumerate(previous_gauge_measured_qubits or [])}
    for ancilla in normal_ancillas:
        idx = measured_index.get(ancilla)
        if idx is None:
            continue
        if family == "X":
            curr = -(half1_measurements - idx)
        else:
            curr = -(half2_measurements - idx)
        x, y = coords_by_qubit[ancilla]
        prev_idx = previous_normal_index.get(ancilla)
        if prev_idx is None:
            continue
        prev_block_len = len(previous_normal_measured_qubits or [])
        prev = -(len(measured_qubits) + prev_block_len - prev_idx)
        out.append("DETECTOR", [stim.target_rec(curr), stim.target_rec(prev)], [x, y, 0.0])

    for cluster, center in zip(gauge_clusters, gauge_cluster_centers):
        tgts: list[stim.GateTarget] = []
        for ancilla in cluster:
            idx = measured_index.get(ancilla)
            if idx is None:
                continue
            if family == "X":
                curr = -(half1_measurements - idx)
                prev_span = 0 if first else half2_measurements
            else:
                curr = -(half2_measurements - idx)
                prev_span = half1_measurements
            tgts.append(stim.target_rec(curr))
            if not first:
                prev_idx = previous_gauge_index.get(ancilla)
                if prev_idx is None:
                    continue
                prev = -(len(measured_qubits) + prev_span + len(previous_gauge_measured_qubits or []) - prev_idx)
                tgts.append(stim.target_rec(prev))
        if tgts:
            if first:
                continue
            x, y = center
            out.append("DETECTOR", tgts, [x, y, 0.0])


def _filter_destruction_round(
    op: stim.CircuitInstruction,
    removed_data_qubits: set[int],
) -> stim.CircuitInstruction | None:
    """Strip removed data qubits from a destruction-round op.

    Returns None for DETECTORs (they have stale rec offsets) and for ops that
    become empty. Any other operation has its target list filtered down.
    """
    name = op.name
    if name == "DETECTOR":
        # Keep the original DETECTOR so _rebuild_post_transition_layer can
        # locate the layer (its rec offsets are no longer correct because we
        # filtered the round, but the rebuild will replace them anyway).
        return op
    targets = op.targets_copy()
    args = op.gate_args_copy()
    if name in _PAIR_GATES:
        qubits = [t.qubit_value for t in targets if t.is_qubit_target]
        kept: list[stim.GateTarget] = []
        for i in range(0, len(qubits), 2):
            a = qubits[i]
            b = qubits[i + 1]
            if a in removed_data_qubits or b in removed_data_qubits:
                continue
            kept.extend([stim.GateTarget(a), stim.GateTarget(b)])
        if not kept:
            return None
        return stim.CircuitInstruction(name, kept, args)
    if name in _SINGLE_QUBIT_GATES:
        kept = [t for t in targets if t.is_qubit_target and t.qubit_value not in removed_data_qubits]
        if not kept:
            return None
        return stim.CircuitInstruction(name, kept, args)
    if name in _MEASUREMENT_GATES:
        # Keep removed data qubits in the destructive measurement so that the
        # downstream rec offsets in post-merge detectors stay valid. The qubit
        # is in some unentangled state because the merge half-rounds and the
        # destruction-round gates no longer touch it.
        return op
    return op


def _filter_body_for_family(
    body: stim.Circuit,
    removed_data_qubits: set[int],
    gauge_x: set[int],
    gauge_z: set[int],
    active_family: str,
) -> tuple[stim.Circuit, list[int]]:
    out = stim.Circuit()
    inactive_gauges = gauge_z if active_family == "X" else gauge_x
    active_gauges = gauge_x if active_family == "X" else gauge_z
    measured_qubits: list[int] = []

    for op in body:
        name = op.name
        targets = op.targets_copy()
        args = op.gate_args_copy()

        if name == "SHIFT_COORDS" or name == "DETECTOR":
            continue

        if name in _PAIR_GATES:
            qubits = [t.qubit_value for t in targets if t.is_qubit_target]
            kept: list[stim.GateTarget] = []
            for i in range(0, len(qubits), 2):
                a = qubits[i]
                b = qubits[i + 1]
                if a in removed_data_qubits or b in removed_data_qubits:
                    continue
                if a in inactive_gauges or b in inactive_gauges:
                    continue
                kept.extend([stim.GateTarget(a), stim.GateTarget(b)])
            if kept:
                out.append(name, kept, args)
            continue

        if name in _SINGLE_QUBIT_GATES:
            kept = []
            for target in targets:
                if not target.is_qubit_target:
                    continue
                qubit = target.qubit_value
                if qubit in removed_data_qubits or qubit in inactive_gauges:
                    continue
                kept.append(target)
            if kept:
                out.append(name, kept, args)
            continue

        if name in _MEASUREMENT_GATES:
            kept = []
            for target in targets:
                if not target.is_qubit_target:
                    continue
                qubit = target.qubit_value
                if qubit in removed_data_qubits or qubit in inactive_gauges:
                    continue
                kept.append(target)
                measured_qubits.append(qubit)
            if kept:
                out.append(name, kept, args)
            continue

        out.append(op)

    return out, measured_qubits


def _drop_detectors_at_coords(
    circuit: stim.Circuit,
    blocked_xy: set[tuple[float, float]],
) -> stim.Circuit:
    if not blocked_xy:
        return circuit
    out = stim.Circuit()
    for op in circuit:
        if op.name == "DETECTOR" and len(op.gate_args_copy()) >= 2:
            xy = (float(op.gate_args_copy()[0]), float(op.gate_args_copy()[1]))
            if xy in blocked_xy:
                continue
        out.append(op)
    return out


def _is_detector_at_coords(
    op: stim.CircuitInstruction | stim.CircuitRepeatBlock,
    blocked_xy: set[tuple[float, float]],
) -> bool:
    if not isinstance(op, stim.CircuitInstruction):
        return False
    if op.name != "DETECTOR" or len(op.gate_args_copy()) < 2:
        return False
    xy = (float(op.gate_args_copy()[0]), float(op.gate_args_copy()[1]))
    return xy in blocked_xy


def _drop_detectors_at_coords_from_z(
    circuit: stim.Circuit,
    blocked_xy: set[tuple[float, float]],
    min_z: float,
) -> stim.Circuit:
    out = stim.Circuit()
    for op in circuit:
        if isinstance(op, stim.CircuitRepeatBlock):
            out.append(stim.CircuitRepeatBlock(op.repeat_count, _drop_detectors_at_coords_from_z(op.body_copy(), blocked_xy, min_z)))
            continue
        if op.name == "DETECTOR" and len(op.gate_args_copy()) >= 3:
            x, y, z = map(float, op.gate_args_copy()[:3])
            if (x, y) in blocked_xy and z >= min_z:
                continue
        out.append(op)
    return out


def _measurement_blocks(flat_circuit: stim.Circuit) -> list[tuple[int, list[int]]]:
    blocks: list[tuple[int, list[int]]] = []
    for i, op in enumerate(flat_circuit):
        if op.name in _MEASUREMENT_GATES:
            blocks.append((i, [t.qubit_value for t in op.targets_copy() if t.is_qubit_target]))
    return blocks


def _count_measurement_blocks(circuit: stim.Circuit) -> int:
    count = 0
    for op in circuit:
        if isinstance(op, stim.CircuitRepeatBlock):
            count += op.repeat_count * _count_measurement_blocks(op.body_copy())
        elif op.name in _MEASUREMENT_GATES:
            count += 1
    return count


def _count_pre_merge_blocks(circuit: stim.Circuit, seed_start_index: int) -> int:
    count = 0
    for i in range(seed_start_index):
        op = circuit[i]
        if isinstance(op, stim.CircuitRepeatBlock):
            count += op.repeat_count * _count_measurement_blocks(op.body_copy())
        elif op.name in _MEASUREMENT_GATES:
            count += 1
    return count


def _measurement_qubits_for_observable(
    flat_circuit: stim.Circuit,
) -> tuple[list[int], list[tuple[int, int, list[int], int, int]]]:
    measurements: list[int] = []
    blocks: list[tuple[int, int, list[int], int, int]] = []
    block_idx = 0
    for i, op in enumerate(flat_circuit):
        if op.name not in _MEASUREMENT_GATES:
            continue
        qubits = [t.qubit_value for t in op.targets_copy() if t.is_qubit_target]
        start = len(measurements)
        measurements.extend(qubits)
        blocks.append((block_idx, i, qubits, start, len(measurements) - 1))
        block_idx += 1
    return measurements, blocks


def _coords_by_qubit_from_circuit(flat_circuit: stim.Circuit) -> dict[int, tuple[float, float]]:
    out: dict[int, tuple[float, float]] = {}
    for op in flat_circuit:
        if op.name != "QUBIT_COORDS":
            continue
        t = op.targets_copy()[0]
        out[t.qubit_value] = (float(op.gate_args_copy()[0]), float(op.gate_args_copy()[1]))
    return out


def _original_observable_pieces(
    flat_circuit: stim.Circuit,
) -> dict[int, list[tuple[int, list[int]]]]:
    measurements, blocks = _measurement_qubits_for_observable(flat_circuit)
    abs_to_block_and_qubit: dict[int, tuple[int, int]] = {}
    for block_idx, _op_i, qubits, start, _end in blocks:
        for offset, qubit in enumerate(qubits):
            abs_to_block_and_qubit[start + offset] = (block_idx, qubit)

    pieces: dict[int, list[tuple[int, list[int]]]] = {}
    for op in flat_circuit:
        if op.name != "OBSERVABLE_INCLUDE":
            continue
        obs_index = int(op.gate_args_copy()[0]) if op.gate_args_copy() else 0
        by_block: dict[int, list[int]] = {}
        for target in op.targets_copy():
            if not target.is_measurement_record_target:
                continue
            abs_idx = len(measurements) + target.value
            block_idx, qubit = abs_to_block_and_qubit[abs_idx]
            by_block.setdefault(block_idx, []).append(qubit)
        for block_idx, qubits in by_block.items():
            pieces.setdefault(obs_index, []).append((block_idx, qubits))
    return pieces


def _rebuild_observables(
    original_flat: stim.Circuit,
    rewritten_flat: stim.Circuit,
    *,
    gauge_x_ancillas: tuple[int, ...],
    gauge_x_clusters: tuple[tuple[int, ...], ...] = (),
) -> stim.Circuit:
    original_pieces = _original_observable_pieces(original_flat)
    if not original_pieces:
        print("_rebuild_observables: no original observables", flush=True)
        return rewritten_flat

    measurements, blocks = _measurement_qubits_for_observable(rewritten_flat)
    coords_by_qubit = _coords_by_qubit_from_circuit(rewritten_flat)
    idx_by_block_qubit = {
        (block_idx, qubit): start + offset
        for block_idx, _op_i, qubits, start, _end in blocks
        for offset, qubit in enumerate(qubits)
    }

    base = stim.Circuit()
    for op in rewritten_flat:
        if op.name != "OBSERVABLE_INCLUDE":
            base.append(op)

    gauge_x_set = set(gauge_x_ancillas)
    total_measurements = len(measurements)
    print(
        f"_rebuild_observables: originals={len(original_pieces)} blocks={len(blocks)} "
        f"measurements={total_measurements} gauge_x={len(gauge_x_set)} "
        f"cluster_candidates={len(gauge_x_clusters)}",
        flush=True,
    )

    def describe_targets(targets: list[stim.GateTarget]) -> list[tuple[int, tuple[float, float] | None]]:
        described: list[tuple[int, tuple[float, float] | None]] = []
        for target in targets:
            if not target.is_measurement_record_target:
                continue
            abs_idx = total_measurements + target.value
            if abs_idx < 0 or abs_idx >= len(measurements):
                described.append((abs_idx, None))
                continue
            qubit = measurements[abs_idx]
            described.append((qubit, coords_by_qubit.get(qubit)))
        return described

    def append_observable_with_log(obs_index: int, targets: list[stim.GateTarget], reason: str) -> None:
        print(
            f"_rebuild_observables[{obs_index}]: selected {reason} "
            f"targets={describe_targets(targets)}",
            flush=True,
        )
        base.append("OBSERVABLE_INCLUDE", targets, [float(obs_index)])

    def try_border_observable(obs_index: int) -> list[stim.GateTarget] | None:
        def solve_candidates(candidates: list[list[stim.GateTarget]]) -> list[stim.GateTarget] | None:
            print(
                f"_rebuild_observables[{obs_index}]: trying {len(candidates)} border candidates",
                flush=True,
            )
            for tgts in candidates:
                trial = stim.Circuit(str(base))
                trial.append("OBSERVABLE_INCLUDE", tgts, [float(obs_index)])
                try:
                    trial.detector_error_model(allow_gauge_detectors=True)
                    return tgts
                except ValueError:
                    continue
            return None

        final_block = blocks[-1]
        final_qubits = final_block[2]
        final_xs = sorted(
            {
                coords_by_qubit[q][0]
                for q in final_qubits
                if q in coords_by_qubit and int(round(coords_by_qubit[q][0])) % 2 == 1
            }
        )
        final_ys = sorted(
            {
                coords_by_qubit[q][1]
                for q in final_qubits
                if q in coords_by_qubit and int(round(coords_by_qubit[q][1])) % 2 == 1
            }
        )

        # Export-generated schedules often want a final-readout-only border
        # representative. Try one or two outer data columns/rows first.
        final_candidates: list[list[stim.GateTarget]] = []
        for x in final_xs:
            col = [q for q in final_qubits if coords_by_qubit.get(q, (None, None))[0] == x]
            col = [q for _, q in sorted((coords_by_qubit[q][1], q) for q in col)]
            if col:
                final_candidates.append(
                    [stim.target_rec(idx_by_block_qubit[(final_block[0], q)] - total_measurements) for q in col]
                )
        for y in final_ys:
            row = [q for q in final_qubits if coords_by_qubit.get(q, (None, None))[1] == y]
            row = [q for _, q in sorted((coords_by_qubit[q][0], q) for q in row)]
            if row:
                final_candidates.append(
                    [stim.target_rec(idx_by_block_qubit[(final_block[0], q)] - total_measurements) for q in row]
                )
        edge_x_pairs = [(0, 1), (-2, -1)] if len(final_xs) >= 2 else []
        edge_y_pairs = [(0, 1), (-2, -1)] if len(final_ys) >= 2 else []
        for a, b in edge_x_pairs:
            xs = {final_xs[a], final_xs[b]}
            cols = [q for q in final_qubits if coords_by_qubit.get(q, (None, None))[0] in xs]
            cols = [q for _, _, q in sorted((coords_by_qubit[q][0], coords_by_qubit[q][1], q) for q in cols)]
            if cols:
                final_candidates.append(
                    [stim.target_rec(idx_by_block_qubit[(final_block[0], q)] - total_measurements) for q in cols]
                )
        for a, b in edge_y_pairs:
            ys = {final_ys[a], final_ys[b]}
            rows = [q for q in final_qubits if coords_by_qubit.get(q, (None, None))[1] in ys]
            rows = [q for _, _, q in sorted((coords_by_qubit[q][1], coords_by_qubit[q][0], q) for q in rows)]
            if rows:
                final_candidates.append(
                    [stim.target_rec(idx_by_block_qubit[(final_block[0], q)] - total_measurements) for q in rows]
                )
        solved = solve_candidates(final_candidates)
        if solved is not None:
            print(f"_rebuild_observables[{obs_index}]: solved by final-block border candidate", flush=True)
            return solved

        if len(blocks) < 3:
            return None
        late_block = blocks[-3]
        late_qubits = late_block[2]
        late_xs = sorted(
            {
                coords_by_qubit[q][0]
                for q in late_qubits
                if q in coords_by_qubit and int(round(coords_by_qubit[q][0])) % 2 == 1
            }
        )
        if not late_xs or not final_xs:
            return None

        candidates: list[list[stim.GateTarget]] = []
        for x in (max(late_xs[0], final_xs[0]), min(late_xs[-1], final_xs[-1])):
            late_col = [q for q in late_qubits if coords_by_qubit.get(q, (None, None))[0] == x]
            final_col = [q for q in final_qubits if coords_by_qubit.get(q, (None, None))[0] == x]
            late_col = [q for _, q in sorted((coords_by_qubit[q][1], q) for q in late_col)]
            final_col = [q for _, q in sorted((coords_by_qubit[q][1], q) for q in final_col)]
            if not late_col or not final_col:
                continue
            tgts = [stim.target_rec(idx_by_block_qubit[(late_block[0], q)] - total_measurements) for q in late_col]
            tgts += [stim.target_rec(idx_by_block_qubit[(final_block[0], q)] - total_measurements) for q in final_col]
            candidates.append(tgts)

        return solve_candidates(candidates)

    def try_left_column_split_observable(obs_index: int) -> list[stim.GateTarget] | None:
        if len(blocks) < 3:
            return None
        final_block = blocks[-1]
        late_block = blocks[-3]
        final_qubits = final_block[2]
        late_qubits = late_block[2]
        odd_xs = sorted(
            {
                coords_by_qubit[q][0]
                for q in final_qubits
                if q in coords_by_qubit and int(round(coords_by_qubit[q][0])) % 2 == 1
            }
        )
        if not odd_xs:
            return None
        left_x = odd_xs[0]
        final_col = [
            q for q in final_qubits
            if coords_by_qubit.get(q, (None, None))[0] == left_x
        ]
        late_col = [
            q for q in late_qubits
            if coords_by_qubit.get(q, (None, None))[0] == left_x
        ]
        final_col = [q for _, q in sorted((coords_by_qubit[q][1], q) for q in final_col)]
        late_col = [q for _, q in sorted((coords_by_qubit[q][1], q) for q in late_col)]
        if not final_col or not late_col:
            return None
        tgts = [
            stim.target_rec(idx_by_block_qubit[(late_block[0], q)] - total_measurements)
            for q in late_col
        ]
        tgts += [
            stim.target_rec(idx_by_block_qubit[(final_block[0], q)] - total_measurements)
            for q in final_col
        ]
        trial = stim.Circuit(str(base))
        trial.append("OBSERVABLE_INCLUDE", tgts, [float(obs_index)])
        try:
            trial.detector_error_model(allow_gauge_detectors=True)
            print(
                f"_rebuild_observables[{obs_index}]: solved by left-column split candidate x={left_x}",
                flush=True,
            )
            return tgts
        except ValueError:
            print(
                f"_rebuild_observables[{obs_index}]: left-column split candidate rejected x={left_x}",
                flush=True,
            )
            return None

    _, original_blocks = _measurement_qubits_for_observable(original_flat)
    original_block_count = len(original_blocks)
    rewritten_block_count = len(blocks)

    def remap_block(orig_block_idx: int) -> int:
        if orig_block_idx == original_block_count - 1:
            return rewritten_block_count - 1
        return orig_block_idx

    # Build cluster-pair candidates: per (block_idx, cluster), include the XOR
    # of all gauge X measurements in that cluster within that block. Each is
    # a true X-superstabilizer in that block — it commutes with Z errors and
    # so adding it to an X observable cannot introduce Z anti-commutation.
    cluster_pair_candidates: list[list[stim.GateTarget]] = []
    if gauge_x_clusters:
        for block_idx in range(rewritten_block_count):
            for cluster in gauge_x_clusters:
                tgts: list[stim.GateTarget] = []
                for q in sorted(cluster):
                    if (block_idx, q) in idx_by_block_qubit:
                        tgts.append(
                            stim.target_rec(idx_by_block_qubit[(block_idx, q)] - total_measurements)
                        )
                if len(tgts) == len(cluster):
                    cluster_pair_candidates.append(tgts)

    def try_with_cluster_pairs(seed_targets: list[stim.GateTarget], obs_index: int) -> list[stim.GateTarget] | None:
        max_r = min(len(cluster_pair_candidates), 8)
        print(
            f"_rebuild_observables[{obs_index}]: trying cluster-pair expansions "
            f"seed_len={len(seed_targets)} candidates={len(cluster_pair_candidates)} max_r={max_r}",
            flush=True,
        )
        for r in range(max_r + 1):
            print(f"_rebuild_observables[{obs_index}]: cluster expansion size {r}", flush=True)
            for subset in combinations(cluster_pair_candidates, r):
                tgts = list(seed_targets)
                for pair in subset:
                    tgts.extend(pair)
                trial = stim.Circuit(str(base))
                trial.append("OBSERVABLE_INCLUDE", tgts, [float(obs_index)])
                try:
                    trial.detector_error_model(allow_gauge_detectors=True)
                    return tgts
                except ValueError:
                    continue
        return None

    for obs_index, pieces in sorted(original_pieces.items()):
        print(f"_rebuild_observables[{obs_index}]: processing {len(pieces)} pieces", flush=True)
        fixed_targets: list[stim.GateTarget] = []
        missing_pieces: list[tuple[int, list[int]]] = []
        for block_idx, qubits in pieces:
            mapped = remap_block(block_idx)
            kept = [q for q in qubits if (mapped, q) in idx_by_block_qubit]
            for q in kept:
                fixed_targets.append(stim.target_rec(idx_by_block_qubit[(mapped, q)] - total_measurements))
            if len(kept) != len(qubits):
                missing_pieces.append((mapped, kept))

        split_targets = try_left_column_split_observable(obs_index)
        if split_targets is not None:
            append_observable_with_log(obs_index, split_targets, "left-column split candidate")
            continue

        if fixed_targets:
            solved_fixed = try_with_cluster_pairs(fixed_targets, obs_index)
            if solved_fixed is not None:
                append_observable_with_log(obs_index, solved_fixed, "fixed-target cluster repair")
                continue

        border_targets = try_border_observable(obs_index)
        if border_targets is not None:
            append_observable_with_log(obs_index, border_targets, "border candidate")
            continue

        if fixed_targets:
            print(f"_rebuild_observables[{obs_index}]: falling back to fixed targets", flush=True)
            append_observable_with_log(obs_index, fixed_targets, "fixed targets fallback")
            continue

        raise ValueError("Failed to rebuild a deterministic observable for the deformed circuit.")

    return base


def _rebuild_post_transition_layer(
    flat_circuit: stim.Circuit,
    *,
    coords_by_qubit: dict[int, tuple[float, float]],
    gauge_ancilla_coords: set[tuple[float, float]],
    rewritten_repeat_count: int,
    pre_merge_block_count: int,
    block_offset: int = 0,
) -> stim.Circuit:
    blocks = _measurement_blocks(flat_circuit)
    post_block_index = pre_merge_block_count + 2 * rewritten_repeat_count + block_offset
    if len(blocks) <= post_block_index:
        return flat_circuit

    curr_op_index, curr_qubits = blocks[post_block_index]
    _, prev_qubits = blocks[post_block_index - 1]
    next_measure_index = blocks[post_block_index + 1][0] if post_block_index + 1 < len(blocks) else len(flat_circuit)

    detector_indices: list[int] = []
    layer_z: float | None = None
    for i in range(curr_op_index + 1, next_measure_index):
        op = flat_circuit[i]
        if op.name != "DETECTOR":
            continue
        detector_indices.append(i)
        if layer_z is None and len(op.gate_args_copy()) >= 3:
            layer_z = float(op.gate_args_copy()[2])
    if not detector_indices or layer_z is None:
        return flat_circuit

    curr_index = {q: i for i, q in enumerate(curr_qubits)}
    prev_index = {q: i for i, q in enumerate(prev_qubits)}
    rebuilt: list[stim.CircuitInstruction] = []
    for q in curr_qubits:
        xy = (float(coords_by_qubit[q][0]), float(coords_by_qubit[q][1]))
        if xy in gauge_ancilla_coords:
            continue
        pi = prev_index.get(q)
        ci = curr_index.get(q)
        if pi is None or ci is None:
            continue
        curr = -(len(curr_qubits) - ci)
        prev = -(len(curr_qubits) + len(prev_qubits) - pi)
        inst = stim.CircuitInstruction("DETECTOR", [stim.target_rec(curr), stim.target_rec(prev)], [xy[0], xy[1], layer_z])
        rebuilt.append(inst)

    out = stim.Circuit()
    rebuilt_iter = iter(rebuilt)
    first = detector_indices[0]
    last = detector_indices[-1]
    for i, op in enumerate(flat_circuit):
        if i == first:
            for det in rebuilt_iter:
                out.append(det)
        if first <= i <= last and op.name == "DETECTOR":
            continue
        out.append(op)
    return out


def rewrite_merge_repeat_for_interior_gauges(
    circuit: stim.Circuit,
    removed_coords: Iterable[tuple[float, float]],
) -> tuple[stim.Circuit, RewriteSummary]:
    removed_coords = [tuple(map(float, coord)) for coord in removed_coords]
    print(f"rewrite_merge_repeat_for_interior_gauges: start removed_coords={removed_coords}", flush=True)
    if not removed_coords:
        raise ValueError("At least one interior seam data qubit must be provided.")
    if not supports_interior_gauge_rewrite(circuit, removed_coords):
        raise ValueError("This prototype only supports interior data qubits, not border removals.")

    removed_data_qubits = set(_resolve_removed_data_qubits(circuit, removed_coords))
    print(f"rewrite_merge_repeat_for_interior_gauges: resolved removed_data_qubits={sorted(removed_data_qubits)}", flush=True)
    flat = circuit.flattened()
    coords_by_qubit = qubit_coords(flat)
    ancillas = repeated_measure_qubits(flat)
    families, touched = ancilla_interactions(flat, ancillas)
    ancillas_by_family: dict[str, set[int]] = {"X": set(), "Z": set()}
    for ancilla, family in families.items():
        ancillas_by_family.setdefault(family, set()).add(ancilla)

    gauge_x = {
        ancilla
        for ancilla, support in touched.items()
        if families.get(ancilla) == "X" and support & removed_data_qubits
    }
    gauge_z = {
        ancilla
        for ancilla, support in touched.items()
        if families.get(ancilla) == "Z" and support & removed_data_qubits
    }
    print(
        "rewrite_merge_repeat_for_interior_gauges: gauge_ancillas "
        f"X={len(gauge_x)} Z={len(gauge_z)}",
        flush=True,
    )
    if not gauge_x and not gauge_z:
        raise ValueError("No repeated ancilla checks were turned into interior gauges by the selected data qubits.")

    removed_clusters = _connected_removed_clusters(coords_by_qubit, removed_data_qubits)
    gauge_x_clusters = _gauge_clusters_for_family(touched, gauge_x, removed_clusters)
    gauge_z_clusters = _gauge_clusters_for_family(touched, gauge_z, removed_clusters)
    gauge_cluster_centers = tuple(_cluster_center(cluster, coords_by_qubit) for cluster in removed_clusters)
    x_cluster_centers = gauge_cluster_centers[: len(gauge_x_clusters)]
    z_cluster_centers = gauge_cluster_centers[: len(gauge_z_clusters)]
    gauge_ancilla_coords = {
        (float(coords_by_qubit[anc][0]), float(coords_by_qubit[anc][1]))
        for anc in sorted(gauge_x | gauge_z)
    }
    x_entry_transition_pre_ancillas = _find_entry_transition_subset(
        family="X",
        ancillas_by_family=ancillas_by_family,
        touched=touched,
        removed_data_qubits=removed_data_qubits,
        gauge_ancillas=gauge_x,
    )
    z_entry_transition_pre_ancillas = _find_entry_transition_subset(
        family="Z",
        ancillas_by_family=ancillas_by_family,
        touched=touched,
        removed_data_qubits=removed_data_qubits,
        gauge_ancillas=gauge_z,
    )

    repeat_index, repeat_block = _find_merge_repeat(circuit)
    print(
        "rewrite_merge_repeat_for_interior_gauges: located merge repeat "
        f"index={repeat_index} count={repeat_block.repeat_count}",
        flush=True,
    )
    seed_start_index = _find_preceding_seed_block_start(circuit, repeat_index)
    merge_end_index = _find_merge_phase_end(circuit, repeat_index)
    if repeat_block.repeat_count % 2 != 0:
        raise ValueError(
            f"Expected an even merge repeat count so it can be halved after subround splitting; got {repeat_block.repeat_count}."
        )
    body = repeat_block.body_copy()
    original_measured = _measurement_targets(body)

    body_init_qubits: dict[str, set[int]] = {}
    for op in body:
        if op.name in {"R", "RX", "RY"}:
            body_init_qubits.setdefault(op.name, set()).update(
                t.qubit_value for t in op.targets_copy() if t.is_qubit_target
            )

    seam_data_inits: list[tuple[str, list[int]]] = []
    for j in range(seed_start_index, repeat_index):
        op = circuit[j]
        if isinstance(op, stim.CircuitRepeatBlock):
            continue
        if op.name in {"R", "RX", "RY"}:
            qs = [t.qubit_value for t in op.targets_copy() if t.is_qubit_target]
            extras = [
                q for q in qs
                if q not in body_init_qubits.get(op.name, set()) and q not in removed_data_qubits
            ]
            if extras:
                seam_data_inits.append((op.name, extras))

    pre_merge_tail_measured: list[int] = []
    for j in range(repeat_index - 1, -1, -1):
        prev_op = circuit[j]
        if isinstance(prev_op, stim.CircuitRepeatBlock):
            pre_merge_tail_measured = _measurement_targets(prev_op.body_copy())
            break
    half1_measured = _active_measurements(original_measured, gauge_x, gauge_z, active_family="X")
    half2_measured = _active_measurements(original_measured, gauge_x, gauge_z, active_family="Z")

    half1_seed, half1_order = _filter_body_for_family(body, removed_data_qubits, gauge_x, gauge_z, active_family="X")
    half2_seed, half2_order = _filter_body_for_family(body, removed_data_qubits, gauge_x, gauge_z, active_family="Z")
    print(
        "rewrite_merge_repeat_for_interior_gauges: filtered half rounds "
        f"half1_measured={len(half1_order)} half2_measured={len(half2_order)}",
        flush=True,
    )
    half1_base = stim.Circuit(str(half1_seed))
    half2_base = stim.Circuit(str(half2_seed))

    normal_half1 = sorted(set(half1_order) - gauge_x)
    normal_half2 = sorted(set(half2_order) - gauge_z)

    _append_half_round_detectors(
        half1_seed,
        measured_qubits=half1_order,
        previous_normal_measured_qubits=pre_merge_tail_measured,
        previous_gauge_measured_qubits=pre_merge_tail_measured,
        coords_by_qubit=coords_by_qubit,
        normal_ancillas=normal_half1,
        gauge_clusters=gauge_x_clusters,
        gauge_cluster_centers=x_cluster_centers,
        family="X",
        half1_measurements=len(half1_order),
        half2_measurements=len(half2_order),
        original_measurements_per_round=len(original_measured),
        first=True,
    )
    half1_seed = _drop_detectors_at_coords(half1_seed, gauge_ancilla_coords)
    half1_seed.append("SHIFT_COORDS", [], [0.0, 0.0, 1.0])
    _append_half_round_detectors(
        half2_seed,
        measured_qubits=half2_order,
        previous_normal_measured_qubits=half1_order,
        previous_gauge_measured_qubits=original_measured,
        coords_by_qubit=coords_by_qubit,
        normal_ancillas=normal_half2,
        gauge_clusters=gauge_z_clusters,
        gauge_cluster_centers=z_cluster_centers,
        family="Z",
        half1_measurements=len(half1_order),
        half2_measurements=len(half2_order),
        original_measurements_per_round=len(original_measured),
        first=True,
    )
    half2_seed = _drop_detectors_at_coords(half2_seed, gauge_ancilla_coords)
    half2_seed.append("SHIFT_COORDS", [], [0.0, 0.0, 1.0])

    half1 = stim.Circuit(str(half1_base))
    half2 = stim.Circuit(str(half2_base))
    _append_half_round_detectors(
        half1,
        measured_qubits=half1_order,
        previous_normal_measured_qubits=half2_order,
        previous_gauge_measured_qubits=half1_order,
        coords_by_qubit=coords_by_qubit,
        normal_ancillas=normal_half1,
        gauge_clusters=gauge_x_clusters,
        gauge_cluster_centers=x_cluster_centers,
        family="X",
        half1_measurements=len(half1_order),
        half2_measurements=len(half2_order),
        original_measurements_per_round=len(original_measured),
    )
    half1 = _drop_detectors_at_coords(half1, gauge_ancilla_coords)
    half1.append("SHIFT_COORDS", [], [0.0, 0.0, 1.0])
    _append_half_round_detectors(
        half2,
        measured_qubits=half2_order,
        previous_normal_measured_qubits=half1_order,
        previous_gauge_measured_qubits=half2_order,
        coords_by_qubit=coords_by_qubit,
        normal_ancillas=normal_half2,
        gauge_clusters=gauge_z_clusters,
        gauge_cluster_centers=z_cluster_centers,
        family="Z",
        half1_measurements=len(half1_order),
        half2_measurements=len(half2_order),
        original_measurements_per_round=len(original_measured),
    )
    half2 = _drop_detectors_at_coords(half2, gauge_ancilla_coords)
    half2.append("SHIFT_COORDS", [], [0.0, 0.0, 1.0])

    seeded_body = stim.Circuit()
    for init_name, init_qubits in seam_data_inits:
        seeded_body.append(init_name, init_qubits, [])
    seeded_body += half1_seed
    seeded_body += half2_seed

    rewritten_body = stim.Circuit()
    rewritten_body += half1
    rewritten_body += half2

    out = stim.Circuit()
    inserted_rewrite = False
    print("rewrite_merge_repeat_for_interior_gauges: assembling rewritten circuit", flush=True)
    for i, op in enumerate(circuit):
        if seed_start_index <= i <= repeat_index:
            if not inserted_rewrite and i == seed_start_index:
                logical_rounds = repeat_block.repeat_count // 2
                if logical_rounds > 0:
                    out += seeded_body
                    if logical_rounds > 1:
                        out.append(stim.CircuitRepeatBlock(logical_rounds - 1, rewritten_body))
                inserted_rewrite = True
            continue
        if repeat_index < i <= merge_end_index:
            filtered = _filter_destruction_round(op, removed_data_qubits)
            if filtered is None:
                continue
            out.append(filtered)
            continue
        if i > merge_end_index and _is_detector_at_coords(op, gauge_ancilla_coords):
            continue
        out.append(op)

    flat_out = out.flattened()
    cluster_zs = [
        float(op.gate_args_copy()[2])
        for op in flat_out
        if op.name == "DETECTOR"
        and len(op.gate_args_copy()) >= 3
        and (float(op.gate_args_copy()[0]), float(op.gate_args_copy()[1])) in set(x_cluster_centers) | set(z_cluster_centers)
    ]
    local_gauge_zs = [
        float(op.gate_args_copy()[2])
        for op in flat_out
        if op.name == "DETECTOR"
        and len(op.gate_args_copy()) >= 3
        and (float(op.gate_args_copy()[0]), float(op.gate_args_copy()[1])) in gauge_ancilla_coords
    ]
    if cluster_zs and local_gauge_zs:
        later_local_gauge_zs = [z for z in local_gauge_zs if z > max(cluster_zs)]
        if later_local_gauge_zs:
            cutoff_z = min(later_local_gauge_zs)
            out = _drop_detectors_at_coords_from_z(flat_out, gauge_ancilla_coords, cutoff_z)
        else:
            out = flat_out
    else:
        out = flat_out

    pre_merge_block_count = _count_pre_merge_blocks(circuit, seed_start_index)
    print("rewrite_merge_repeat_for_interior_gauges: rebuilding post-transition layer", flush=True)
    out = _rebuild_post_transition_layer(
        out.flattened(),
        coords_by_qubit=coords_by_qubit,
        gauge_ancilla_coords=gauge_ancilla_coords,
        rewritten_repeat_count=repeat_block.repeat_count // 2,
        pre_merge_block_count=pre_merge_block_count,
    )
    print("rewrite_merge_repeat_for_interior_gauges: rebuilding observables", flush=True)
    try:
        out = _rebuild_observables(
            flat,
            out.flattened(),
            gauge_x_ancillas=tuple(sorted(gauge_x)),
            gauge_x_clusters=tuple(gauge_x_clusters),
        )
    except ValueError:
        # Keep the rewritten circuit exportable for inspection even when the
        # logical observable has not yet been reconstructed successfully.
        out = out.flattened()
    print("rewrite_merge_repeat_for_interior_gauges: done", flush=True)

    summary = RewriteSummary(
        removed_data_qubits=tuple(sorted(removed_data_qubits)),
        removed_coords=tuple(removed_coords),
        gauge_x_ancillas=tuple(sorted(gauge_x)),
        gauge_z_ancillas=tuple(sorted(gauge_z)),
        gauge_x_clusters=tuple(gauge_x_clusters),
        gauge_z_clusters=tuple(gauge_z_clusters),
        x_entry_transition_pre_ancillas=x_entry_transition_pre_ancillas,
        z_entry_transition_pre_ancillas=z_entry_transition_pre_ancillas,
        original_repeat_count=repeat_block.repeat_count,
        rewritten_repeat_count=repeat_block.repeat_count // 2,
        original_measurements_per_round=len(original_measured),
        half1_measurements=len(half1_order),
        half2_measurements=len(half2_order),
    )
    return out, summary


def main() -> None:
    args = parse_args()
    circuit = stim.Circuit.from_file(args.stim)
    rewritten, summary = rewrite_merge_repeat_for_interior_gauges(
        circuit,
        parse_coordinate_list(args.remove_data),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(str(rewritten), encoding="utf-8")
    print(f"Wrote {args.out}")
    print(
        "Summary: "
        f"removed_data={list(summary.removed_data_qubits)} "
        f"gauge_x={list(summary.gauge_x_ancillas)} "
        f"gauge_z={list(summary.gauge_z_ancillas)} "
        f"x_clusters={[list(cluster) for cluster in summary.gauge_x_clusters]} "
        f"z_clusters={[list(cluster) for cluster in summary.gauge_z_clusters]} "
        f"x_entry_pre={list(summary.x_entry_transition_pre_ancillas)} "
        f"z_entry_pre={list(summary.z_entry_transition_pre_ancillas)} "
        f"repeat={summary.original_repeat_count}->{summary.rewritten_repeat_count} "
        f"mx={summary.original_measurements_per_round}->{summary.half1_measurements}/{summary.half2_measurements}"
    )


if __name__ == "__main__":
    main()
