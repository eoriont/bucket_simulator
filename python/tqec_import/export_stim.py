#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from tqec import BlockGraph, NoiseModel, compile_block_graph
from tqec.computation.cube import ZXCube
from tqec.computation.open_graph import FilledGraph
from tqec.gallery.cnot import cnot
from tqec.utils.enums import Basis
from tqec.utils.position import Position3D
from tqec.utils.scale import LinearFunction

from paths import GENERATED_ROOT
from rewrite_interior_gauges import rewrite_merge_repeat_for_interior_gauges, supports_interior_gauge_rewrite
from stim_annotations import annotate_with_polygons
from stim_layout import deformation_coordinates_for_mode, seam_border_coordinate
from stim_deformation import apply_code_deformation, parse_coordinate_list


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a verified Stim circuit from TQEC."
    )
    parser.add_argument(
        "--source",
        choices=("gallery_cnot", "simple_merge", "simple_merge_only", "simple_split_only", "simple_merge_split", "dae"),
        required=True,
        help="Where to load the TQEC computation from.",
    )
    parser.add_argument(
        "--basis",
        choices=("X", "Z"),
        default="X",
        help="Boundary basis used for the gallery logical CNOT.",
    )
    parser.add_argument(
        "--dae",
        type=Path,
        help="Path to a TQEC .dae computation when --source dae is used.",
    )
    parser.add_argument(
        "--k",
        type=int,
        default=1,
        help="TQEC spatial scale factor. Physical distance is typically d = 2k + 1.",
    )
    parser.add_argument(
        "--syndrome-rounds",
        type=int,
        default=None,
        help=(
            "Override the block temporal height with a constant number of stabilizer rounds. "
            "Use this to keep distance fixed while changing the number of syndrome cycles."
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Output .stim path. Open-port DAE imports may emit multiple suffixed files.",
    )
    parser.add_argument(
        "--manhattan-radius",
        type=int,
        default=2,
        help="Detector search radius passed to TQEC.",
    )
    parser.add_argument(
        "--with-uniform-depolarizing-noise",
        type=float,
        default=None,
        help="Optional TQEC noise level. Leave unset for a clean circuit export.",
    )
    parser.add_argument(
        "--remove-data",
        type=str,
        default=None,
        help=(
            "Optional simple code-deformation list applied after TQEC export, "
            'for example "(3,9) (5,9)".'
        ),
    )
    parser.add_argument(
        "--remove-seam-border",
        choices=("min", "max"),
        default=None,
        help="Infer the seam orientation and remove one data qubit from the corresponding seam endpoint.",
    )
    parser.add_argument(
        "--deformation-mode",
        choices=("lside", "rside", "twoside"),
        default=None,
        help="Infer the seam orientation and remove the left, right, or both seam-edge data qubits.",
    )
    return parser.parse_args()


def make_noise_model(noise: float | None) -> NoiseModel | None:
    if noise is None:
        return None
    return NoiseModel.uniform_depolarizing(noise)


def build_gallery_cnot(basis_name: str) -> tuple[BlockGraph, list | None]:
    basis = Basis(basis_name)
    graph = cnot(basis)
    observables = graph.find_correlation_surfaces()
    return graph, observables


def build_simple_merge(basis_name: str) -> tuple[BlockGraph, list | None]:
    basis = Basis(basis_name)
    fill_kind = ZXCube.from_str("ZXX" if basis == Basis.X else "ZXZ")

    graph = BlockGraph(f"Two-patch {basis.value}-basis merge")
    nodes = [
        (Position3D(0, 0, 0), "P", "In_A"),
        (Position3D(0, 1, 0), "P", "In_B"),
        (Position3D(0, 0, 1), "ZXX", ""),
        (Position3D(0, 1, 1), "ZXX", ""),
        (Position3D(0, 0, 2), "P", "Out_M"),
    ]

    for pos, kind, label in nodes:
        graph.add_cube(pos, kind, label)
    for start, end in [(0, 2), (1, 3), (2, 3), (2, 4)]:
        graph.add_pipe(nodes[start][0], nodes[end][0])

    graph.fill_ports({"In_A": fill_kind, "In_B": fill_kind, "Out_M": fill_kind})
    observables = graph.find_correlation_surfaces()
    return graph, observables


def build_simple_split(basis_name: str) -> tuple[BlockGraph, list | None]:
    basis = Basis(basis_name)
    fill_kind = ZXCube.from_str("ZXX" if basis == Basis.X else "ZXZ")

    graph = BlockGraph(f"One-patch {basis.value}-basis split")
    nodes = [
        (Position3D(0, 0, 0), "P", "In_M"),
        (Position3D(0, 0, 1), "ZXX", ""),
        (Position3D(0, 1, 1), "ZXX", ""),
        (Position3D(0, 0, 2), "P", "Out_A"),
        (Position3D(0, 1, 2), "P", "Out_B"),
    ]

    for pos, kind, label in nodes:
        graph.add_cube(pos, kind, label)
    for start, end in [(0, 1), (1, 2), (1, 3), (2, 4)]:
        graph.add_pipe(nodes[start][0], nodes[end][0])

    graph.fill_ports({"In_M": fill_kind, "Out_A": fill_kind, "Out_B": fill_kind})
    observables = graph.find_correlation_surfaces()
    return graph, observables


def build_simple_merge_split(basis_name: str) -> tuple[BlockGraph, list | None]:
    # The existing "simple_merge" graph has three top-level temporal sections in the
    # emitted Stim circuit, so the retiming/sweep tooling uses it as a convenient
    # merge/split schedule proxy even though the logical graph is a merge-only
    # computation.
    return build_simple_merge(basis_name)


def build_dae_graph(dae_path: Path) -> tuple[list[FilledGraph], bool]:
    graph = BlockGraph.from_dae_file(dae_path.resolve(), graph_name=dae_path.stem)
    if graph.num_ports == 0:
        return [FilledGraph(graph=graph, stabilizers=[], observables=graph.find_correlation_surfaces())], False
    return graph.fill_ports_for_minimal_simulation(), True


def compile_and_write(
    graph: BlockGraph,
    observables: list | None,
    out_path: Path,
    k: int,
    syndrome_rounds: int | None,
    manhattan_radius: int,
    noise: float | None,
    remove_data: str | None,
    remove_seam_border: str | None,
    deformation_mode: str | None,
) -> None:
    compile_kwargs = {
        "observables": observables if observables else "auto",
    }
    if syndrome_rounds is not None:
        # TQEC requires cube blocks to remain scalable in time. This affine form
        # evaluates to the requested round count at the chosen spatial scale k.
        compile_kwargs["block_temporal_height"] = LinearFunction(1, syndrome_rounds - k)
    compiled = compile_block_graph(graph, **compile_kwargs)
    circuit = compiled.generate_stim_circuit(
        k=k,
        noise_model=make_noise_model(noise),
        manhattan_radius=manhattan_radius,
    )
    if remove_data:
        removed_coords = parse_coordinate_list(remove_data)
        if supports_interior_gauge_rewrite(circuit, removed_coords):
            circuit, _summary = rewrite_merge_repeat_for_interior_gauges(circuit, removed_coords)
        else:
            circuit = apply_code_deformation(circuit, removed_coords).circuit
    elif remove_seam_border:
        circuit = apply_code_deformation(circuit, [seam_border_coordinate(circuit, remove_seam_border)]).circuit
    elif deformation_mode:
        circuit = apply_code_deformation(circuit, deformation_coordinates_for_mode(circuit, deformation_mode)).circuit
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(annotate_with_polygons(circuit), encoding="utf-8")
    print(
        f"Wrote {out_path} "
        f"(qubits={circuit.num_qubits}, detectors={circuit.num_detectors}, observables={circuit.num_observables})"
    )


def default_output_path(args: argparse.Namespace) -> Path:
    if args.source == "dae":
        stem = args.dae.stem if args.dae is not None else "from_dae"
        return GENERATED_ROOT / f"{stem}_k{args.k}.stim"
    rounds_suffix = f"_r{args.syndrome_rounds}" if args.syndrome_rounds is not None else ""
    return GENERATED_ROOT / f"{args.source}_{args.basis.lower()}_k{args.k}{rounds_suffix}.stim"


def main() -> None:
    args = parse_args()
    specified = sum(
        option is not None
        for option in (args.remove_data, args.remove_seam_border, args.deformation_mode)
    )
    if specified > 1:
        raise SystemExit("Specify at most one of --remove-data, --remove-seam-border, or --deformation-mode.")
    if args.out is None:
        args.out = default_output_path(args)

    if args.source == "gallery_cnot":
        graph, observables = build_gallery_cnot(args.basis)
        compile_and_write(
            graph=graph,
            observables=observables,
            out_path=args.out,
            k=args.k,
            syndrome_rounds=args.syndrome_rounds,
            manhattan_radius=args.manhattan_radius,
            noise=args.with_uniform_depolarizing_noise,
            remove_data=args.remove_data,
            remove_seam_border=args.remove_seam_border,
            deformation_mode=args.deformation_mode,
        )
        return

    if args.source == "simple_merge":
        graph, observables = build_simple_merge(args.basis)
        compile_and_write(
            graph=graph,
            observables=observables,
            out_path=args.out,
            k=args.k,
            syndrome_rounds=args.syndrome_rounds,
            manhattan_radius=args.manhattan_radius,
            noise=args.with_uniform_depolarizing_noise,
            remove_data=args.remove_data,
            remove_seam_border=args.remove_seam_border,
            deformation_mode=args.deformation_mode,
        )
        return

    if args.source == "simple_merge_only":
        graph, observables = build_simple_merge(args.basis)
        compile_and_write(
            graph=graph,
            observables=observables,
            out_path=args.out,
            k=args.k,
            syndrome_rounds=args.syndrome_rounds,
            manhattan_radius=args.manhattan_radius,
            noise=args.with_uniform_depolarizing_noise,
            remove_data=args.remove_data,
            remove_seam_border=args.remove_seam_border,
            deformation_mode=args.deformation_mode,
        )
        return

    if args.source == "simple_split_only":
        graph, observables = build_simple_split(args.basis)
        compile_and_write(
            graph=graph,
            observables=observables,
            out_path=args.out,
            k=args.k,
            syndrome_rounds=args.syndrome_rounds,
            manhattan_radius=args.manhattan_radius,
            noise=args.with_uniform_depolarizing_noise,
            remove_data=args.remove_data,
            remove_seam_border=args.remove_seam_border,
            deformation_mode=args.deformation_mode,
        )
        return

    if args.source == "simple_merge_split":
        graph, observables = build_simple_merge_split(args.basis)
        compile_and_write(
            graph=graph,
            observables=observables,
            out_path=args.out,
            k=args.k,
            syndrome_rounds=args.syndrome_rounds,
            manhattan_radius=args.manhattan_radius,
            noise=args.with_uniform_depolarizing_noise,
            remove_data=args.remove_data,
            remove_seam_border=args.remove_seam_border,
            deformation_mode=args.deformation_mode,
        )
        return

    if args.dae is None:
        raise SystemExit("--dae is required when --source dae is used.")

    filled_graphs, used_auto_fill = build_dae_graph(args.dae)
    if not used_auto_fill:
        compile_and_write(
            graph=filled_graphs[0].graph,
            observables=filled_graphs[0].observables,
            out_path=args.out,
            k=args.k,
            syndrome_rounds=args.syndrome_rounds,
            manhattan_radius=args.manhattan_radius,
            noise=args.with_uniform_depolarizing_noise,
            remove_data=args.remove_data,
            remove_seam_border=args.remove_seam_border,
            deformation_mode=args.deformation_mode,
        )
        return

    stem = args.out.stem
    suffix = args.out.suffix or ".stim"
    for index, filled in enumerate(filled_graphs):
        variant_path = args.out.with_name(f"{stem}.variant_{index}{suffix}")
        compile_and_write(
            graph=filled.graph,
            observables=filled.observables,
            out_path=variant_path,
            k=args.k,
            syndrome_rounds=args.syndrome_rounds,
            manhattan_radius=args.manhattan_radius,
            noise=args.with_uniform_depolarizing_noise,
            remove_data=args.remove_data,
            remove_seam_border=args.remove_seam_border,
            deformation_mode=args.deformation_mode,
        )
        print(f"  variant_{index} stabilizers={filled.stabilizers}")


if __name__ == "__main__":
    main()
