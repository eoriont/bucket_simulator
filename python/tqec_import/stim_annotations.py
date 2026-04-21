from __future__ import annotations

from collections import Counter, defaultdict
from math import atan2

import stim

from stim_layout import ancilla_interactions, infer_seam, qubit_coords, repeated_measure_qubits

def annotate_with_polygons(circuit: stim.Circuit) -> str:
    flat = circuit.flattened()
    coords = qubit_coords(flat)
    ancillas = repeated_measure_qubits(flat)
    families, touched = ancilla_interactions(flat, ancillas)
    seam = infer_seam(flat)

    stabilizer_lines: list[str] = []
    remote_lines: list[str] = []

    for ancilla in sorted(ancillas):
        support = sorted(touched.get(ancilla, set()))
        if not support:
            continue
        family = families.get(ancilla)
        if family == "X":
            stabilizer_lines.append(_polygon_line("0,0,1,0.25", ancilla, support, coords))
        elif family == "Z":
            stabilizer_lines.append(_polygon_line("1,0,0,0.25", ancilla, support, coords))

        for data in support:
            if _crosses_seam(coords[ancilla], coords[data], seam):
                remote_lines.append(f"#!pragma POLYGON(1,0,1,0.5) {ancilla} {data}")

    lines = [line for line in stabilizer_lines + remote_lines if line]
    if not lines:
        return str(circuit)

    rendered: list[str] = []
    inserted = False
    for line in str(circuit).splitlines():
        rendered.append(line)
        if not inserted and not line.startswith("QUBIT_COORDS"):
            rendered[-1:] = lines + [line]
            inserted = True
    if not inserted:
        rendered.extend(lines)
    return "\n".join(rendered) + "\n"


def _polygon_line(color: str, ancilla: int, support: list[int], coords: dict[int, tuple[float, float]]) -> str:
    ax, ay = coords[ancilla]
    ordered = sorted(
        support,
        key=lambda q: (
            atan2(coords[q][1] - ay, coords[q][0] - ax),
            coords[q][0],
            coords[q][1],
        ),
    )
    return f"#!pragma POLYGON({color}) " + " ".join(str(q) for q in ordered)


def _crosses_seam(
    ancilla_xy: tuple[float, float],
    data_xy: tuple[float, float],
    seam,
) -> bool:
    if seam.orientation == "horizontal":
        ay = ancilla_xy[1]
        dy = data_xy[1]
        return (ay < seam.seam_coordinate <= dy) or (dy < seam.seam_coordinate <= ay)
    ax = ancilla_xy[0]
    dx = data_xy[0]
    return (ax < seam.seam_coordinate <= dx) or (dx < seam.seam_coordinate <= ax)
