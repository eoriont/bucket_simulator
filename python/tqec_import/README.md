# TQEC Import Path

This directory is a separate Python-based path for generating verified lattice-surgery Stim circuits from TQEC without touching the existing C++ implementation.

Generated files are organized under `output/tqec_import/`:

- `output/tqec_import/generated/`
- `output/tqec_import/injected/`
- `output/tqec_import/ler/`

## What Works Today

- Export a verified logical CNOT lattice-surgery circuit from `tqec.gallery.cnot`.
- Export a minimal custom two-patch merge computation from a hand-built TQEC `BlockGraph`.
- Import a custom TQEC `.dae` computation with `BlockGraph.from_dae_file(...)`.
- Write the resulting Stim circuit to `output/` for inspection in Crumble.

## Important Constraint

TQEC does not expose a standalone "merge" or "split" helper in the public gallery API used here. The supported paths are:

- Use a built-in logical computation such as `tqec.gallery.cnot`, which is implemented via lattice surgery.
- Import a custom computation from `.dae` and compile that into a Stim circuit.

For the built-in logical CNOT, the TQEC graph must have its open ports filled before compilation. In practice that means either:

- call `cnot(Basis.X)` or `cnot(Basis.Z)`, or
- call `fill_ports_for_minimal_simulation()` on an open graph and export one or more filled variants.

## Setup

`k` in the export/sweep commands below always means the TQEC circuit-size parameter
(`circuit_k`), where the corresponding code distance is `d = 2k + 1`.

The distillation-depth parameter is a different quantity and is always referred to
as `distillation_rounds` in the code and command line. If you want a short name for
it in discussions, use `distillation_k`, not `k`.

```bash
uv venv .venv-tqec --python 3.12
uv pip install --python .venv-tqec/bin/python -r python/tqec_import/requirements.txt
```

For multicore LER runs, install an MPI runtime such as OpenMPI or MPICH and use
`mpirun` with `python/tqec_import/run_ler.py`.

## Export The Gallery Logical CNOT

```bash
.venv-tqec/bin/python python/tqec_import/export_stim.py \
  --source gallery_cnot \
  --basis X \
  --k 1
```

## Export A Minimal Two-Patch Merge

```bash
.venv-tqec/bin/python python/tqec_import/export_stim.py \
  --source simple_merge \
  --basis X \
  --k 1
```

Explicit merge/split alias:

```bash
.venv-tqec/bin/python python/tqec_import/export_stim.py \
  --source simple_merge_split \
  --basis X \
  --k 1
```

## Export From A `.dae` File

```bash
.venv-tqec/bin/python python/tqec_import/export_stim.py \
  --source dae \
  --dae path/to/computation.dae \
  --k 1
```

If the imported `.dae` graph still has open ports, the script will use `fill_ports_for_minimal_simulation()` and emit one Stim file per filled variant.

## Inject Distributed-Network Noise

This mirrors the existing C++ distributed-noise model:

- RCX / remote-CNOT noise from channel error, distillation protocol, backup batches `m`, and optional accurate-RCX folding
- entanglement-limited idling noise from remote-CNOT demand, entanglement rate, measurement delay, `T1`, and `T2`

Example:

```bash
.venv-tqec/bin/python python/tqec_import/inject_noise.py \
  --stim output/tqec_import/generated/simple_merge_x_k1_r5.stim \
  --accurate-rcx \
  --distillation-protocol pumping_2to1 \
  --distillation-rounds 2 \
  --distillation-backup-batches 2 \
  --entanglement-rate 100e6
```

This writes:

- `output/tqec_import/injected/<name>.network_noisy.stim`
- `output/tqec_import/injected/<name>.network_noisy.json`

## Run Batched LER

Single-process batched run:

```bash
.venv-tqec/bin/python python/tqec_import/run_ler.py \
  --stim output/tqec_import/injected/simple_merge_x_k1_r5.network_noisy.stim \
  --shots 1000000 \
  --batch-shots 250000
```

MPI batched run across 8 ranks:

```bash
mpirun -np 8 .venv-tqec/bin/python python/tqec_import/run_ler.py \
  --stim output/tqec_import/injected/simple_merge_x_k1_r5.network_noisy.stim \
  --shots 20000000 \
  --batch-shots 250000
```

The runner splits total shots across MPI ranks, samples and decodes in batches on
each rank, then reduces the logical-failure counts back to rank 0 for reporting.

## Sweep Merge/Split Schedule Parameters

This path can retime the top-level repeat sections of the generated merge/split
circuit to independently control pre-merge, merged, and post-split round counts.

Example monolithic sweep:

```bash
.venv-tqec/bin/python python/tqec_import/sweep_merge_split_rounds.py \
  --basis X \
  --k 2 \
  --pre-rounds 0,1,2,3,4,5 \
  --merge-rounds 5 \
  --post-rounds 0,1,2,3,4,5 \
  --shots 10000000 \
  --batch-shots 250000 \
  --mpi-ranks 8 \
  --monolithic-baseline \
  --name example_merge_split_sweep
```

If `--pre-rounds`, `--merge-rounds`, and `--post-rounds` are omitted, the script
uses the default schedule policy:

- `r_pre = d + 1`
- `r_merge = d + 1`
- `r_post = 0`

where `d = 2k + 1`.

To pin a deterministic d=5 merge schedule with six rounds before and during the
merge and then sweep entanglement rate, pass `--circuit-k-values 2 --rounds-per-phase 6`.

## Raw-Pairs ENR Sweep With Border Deformations

This sweep keeps distillation off (`protocol=none`) and compares the border
deformation modes `none`, `lside`, `rside`, and `twoside`. It also records the
remaining remote-CNOT count as the required raw EPR pairs per round.

```bash
.venv-tqec/bin/python python/tqec_import/sweep_enr_raw_pairs.py \
  --circuit-k 2 \
  --rounds-per-phase 6 \
  --deformation-modes none,lside,rside,twoside \
  --enr-values 25000000,50000000,75000000,100000000,150000000,200000000,300000000,500000000,750000000,1000000000
```

The sweep output includes the remaining remote-CNOT count for each deformation
and the derived `required_epr_pairs_per_round` field.

## Apply Simple Code Deformation

This path now includes a Stim-level post-processing pass for the simple case
where code deformation is modeled by removing specific data qubits and
suppressing the adjacent repeated ancilla checks.

Standalone pass on an already exported circuit:

```bash
.venv-tqec/bin/python python/tqec_import/deform_stim.py \
  --stim output/tqec_import/generated/simple_merge_x_k1_r5.stim \
  --out output/tqec_import/generated/simple_merge_x_k1_r5.deformed.stim \
  --remove-data "(5,7)"
```

Or let the tooling infer the seam orientation and remove one endpoint of the seam:

```bash
.venv-tqec/bin/python python/tqec_import/deform_stim.py \
  --stim output/tqec_import/generated/simple_merge_x_k2_r5.stim \
  --out output/tqec_import/generated/simple_merge_x_k2_r5.deformed_edge_example.stim \
  --mode lside
```

Available inferred seam-endpoint modes are:

- `lside`: left endpoint of the inferred seam
- `rside`: right endpoint of the inferred seam
- `twoside`: both seam endpoints

Or apply the same pass directly during export:

```bash
.venv-tqec/bin/python python/tqec_import/export_stim.py \
  --source simple_merge \
  --basis X \
  --k 1 \
  --remove-data "(5,7)"
```

The current deformation pass is intentionally limited to the simpler
"remove data qubits / suppress nearby checks" case. Superstabilizer gauge-round
alternation is not implemented in the TQEC import path yet.

## Add Polygon Annotations

Deformation/export outputs are now written annotated by default, with
`#!pragma POLYGON(...)` lines for inferred X/Z stabilizers plus magenta
seam-crossing interaction edges. The standalone annotator still exists if you
want to re-annotate an already generated circuit:

```bash
.venv-tqec/bin/python python/tqec_import/annotate_stim.py \
  --stim output/tqec_import/generated/simple_merge_x_k2_r5.deformed_example.stim \
  --out output/tqec_import/generated/simple_merge_x_k2_r5.deformed_example.annotated.stim
```

For the TQEC path these magenta edges are inferred seam-crossing interactions,
not literal remote-CNOT instructions as in the distributed C++ simulator.
