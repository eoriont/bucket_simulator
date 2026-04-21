from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass

import stim


@dataclass(frozen=True)
class InferredSeam:
    orientation: str
    seam_coordinate: float
    orthogonal_min: float
    orthogonal_max: float


def qubit_coords(circuit: stim.Circuit) -> dict[int, tuple[float, float]]:
    coords: dict[int, tuple[float, float]] = {}
    for op in circuit:
        if op.name != "QUBIT_COORDS":
            continue
        target = op.targets_copy()[0]
        if not target.is_qubit_target:
            continue
        x, y = op.gate_args_copy()[:2]
        coords[target.qubit_value] = (float(x), float(y))
    return coords


def repeated_measure_qubits(circuit: stim.Circuit) -> set[int]:
    counts: Counter[int] = Counter()
    for op in circuit:
        if op.name not in {"M", "MX", "MY", "MR", "MRX", "MRY"}:
            continue
        for target in op.targets_copy():
            if target.is_qubit_target:
                counts[target.qubit_value] += 1
    return {q for q, n in counts.items() if n > 1}


def ancilla_interactions(
    circuit: stim.Circuit,
    ancillas: set[int],
) -> tuple[dict[int, str], dict[int, set[int]]]:
    family_votes: dict[int, set[str]] = defaultdict(set)
    touched: dict[int, set[int]] = defaultdict(set)

    for op in circuit:
        if op.name not in {"CX", "CZ"}:
            continue
        qubits = [t.qubit_value for t in op.targets_copy() if t.is_qubit_target]
        for i in range(0, len(qubits), 2):
            a = qubits[i]
            b = qubits[i + 1]
            if a in ancillas and b not in ancillas:
                touched[a].add(b)
                family_votes[a].add("Z" if op.name == "CZ" else "X")
            if b in ancillas and a not in ancillas:
                touched[b].add(a)
                family_votes[b].add("Z" if op.name == "CZ" else "X")

    families: dict[int, str] = {}
    for ancilla, votes in family_votes.items():
        if len(votes) == 1:
            families[ancilla] = next(iter(votes))
    return families, touched


def infer_seam(circuit: stim.Circuit) -> InferredSeam:
    coords = qubit_coords(circuit.flattened())
    ancillas = repeated_measure_qubits(circuit.flattened())
    data = [xy for q, xy in coords.items() if q not in ancillas]
    if not data:
        raise ValueError("Unable to infer seam without data-qubit coordinates.")

    xs = sorted({x for x, _ in data})
    ys = sorted({y for _, y in data})
    x_span = xs[-1] - xs[0]
    y_span = ys[-1] - ys[0]

    if y_span >= x_span:
        seam_coordinate = ys[len(ys) // 2]
        seam_band = [(x, y) for x, y in data if abs(y - seam_coordinate) < 1e-9]
        return InferredSeam(
            orientation="horizontal",
            seam_coordinate=seam_coordinate,
            orthogonal_min=min(x for x, _ in seam_band),
            orthogonal_max=max(x for x, _ in seam_band),
        )

    seam_coordinate = xs[len(xs) // 2]
    seam_band = [(x, y) for x, y in data if abs(x - seam_coordinate) < 1e-9]
    return InferredSeam(
        orientation="vertical",
        seam_coordinate=seam_coordinate,
        orthogonal_min=min(y for _, y in seam_band),
        orthogonal_max=max(y for _, y in seam_band),
    )


def seam_border_coordinate(circuit: stim.Circuit, side: str) -> tuple[float, float]:
    seam = infer_seam(circuit)
    if seam.orientation == "horizontal":
        if side == "min":
            return (seam.orthogonal_min, seam.seam_coordinate)
        if side == "max":
            return (seam.orthogonal_max, seam.seam_coordinate)
    else:
        if side == "min":
            return (seam.seam_coordinate, seam.orthogonal_min)
        if side == "max":
            return (seam.seam_coordinate, seam.orthogonal_max)
    raise ValueError(f"Unknown seam border side: {side}")


def deformation_coordinates_for_mode(circuit: stim.Circuit, mode: str) -> list[tuple[float, float]]:
    if mode == "lside":
        return [seam_border_coordinate(circuit, "min")]
    if mode == "rside":
        return [seam_border_coordinate(circuit, "max")]
    if mode == "twoside":
        return [seam_border_coordinate(circuit, "min"), seam_border_coordinate(circuit, "max")]
    raise ValueError(f"Unknown deformation mode: {mode}")
