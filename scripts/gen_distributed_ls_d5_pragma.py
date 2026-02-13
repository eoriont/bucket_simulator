#!/usr/bin/env python3
"""Generate d=5 distributed lattice surgery Stim circuit with pragma coloring.

Layout: Two distance-5 rotated surface code patches side by side,
with a single column of merge ancillas bridging them. No merge data qubits.

Pragma colors:
  Blue   (0,0,1,0.25)   = Patch X-stabilizer
  Red    (1,0,0,0.25)   = Patch Z-stabilizer
  Green  (0,0.7,0,0.3)  = Merge X-stabilizer (cross-boundary)
  Orange (1,0.5,0,0.3)  = Merge Z-stabilizer (cross-boundary)
"""

import sys
import os

d = 5
num_merge_rounds = d
interconnect_error = 0.01

# ============================================================
# 1. Build qubit layout (matches C++ initialize_layout)
# ============================================================
DATA, X_ANC, Z_ANC = 'DATA', 'X_ANC', 'Z_ANC'
PATCH_A, PATCH_B, MERGE = 'A', 'B', 'M'

# Each qubit: (index, type, patch, x, y)
qubits = []

def add_qubit(qt, patch, x, y):
    idx = len(qubits)
    qubits.append((idx, qt, patch, x, y))
    return idx

# -- Patch A --
add_qubit(X_ANC, PATCH_A, 0.0, float(d - 1))  # corner

for ix in range(d - 1):
    xd = ix + 0.5
    xa = ix + 1
    for j in range(d):
        add_qubit(DATA, PATCH_A, xd, j + 0.5)
    y_min = 1 if (xa % 2 == 0) else 0
    for j in range(d):
        y = y_min + j
        is_x = ((y - y_min) % 2 == 1)
        add_qubit(X_ANC if is_x else Z_ANC, PATCH_A, float(xa), float(y))

for j in range(d):
    add_qubit(DATA, PATCH_A, d - 0.5, j + 0.5)

# -- Merge column at x = d --
mx = d
my_min = 1 if (mx % 2 == 0) else 0
for j in range(d):
    y = my_min + j
    is_x = ((y - my_min) % 2 == 1)
    add_qubit(X_ANC if is_x else Z_ANC, MERGE, float(mx), float(y))

# -- Patch B --
pb_start = d + 0.5
for ix in range(d - 1):
    xd = pb_start + ix
    xa_int = d + 1 + ix
    for j in range(d):
        add_qubit(DATA, PATCH_B, xd, j + 0.5)
    y_min = 1 if (xa_int % 2 == 0) else 0
    for j in range(d):
        y = y_min + j
        is_x = ((y - y_min) % 2 == 1)
        add_qubit(X_ANC if is_x else Z_ANC, PATCH_B, float(xa_int), float(y))

for j in range(d):
    add_qubit(DATA, PATCH_B, pb_start + d - 1, j + 0.5)

add_qubit(X_ANC, PATCH_B, float(2 * d), 1.0)  # corner

# ============================================================
# 2. Build stabilizers
# ============================================================

def find_data_at(x, y):
    for idx, qt, p, qx, qy in qubits:
        if qt == DATA and abs(qx - x) < 0.1 and abs(qy - y) < 0.1:
            return idx
    return -1

# (ancilla_idx, is_x, patch, [data_qubits], crosses_boundary)
all_stabilizers = []

for idx, qt, patch, x, y in qubits:
    if qt == DATA:
        continue
    is_x = (qt == X_ANC)
    dqs = []
    has_a, has_b = False, False
    for dx, dy in [(-0.5, -0.5), (0.5, -0.5), (-0.5, 0.5), (0.5, 0.5)]:
        di = find_data_at(x + dx, y + dy)
        if di >= 0:
            dqs.append(di)
            if qubits[di][2] == PATCH_A: has_a = True
            if qubits[di][2] == PATCH_B: has_b = True
    if dqs:
        all_stabilizers.append((idx, is_x, patch, dqs, has_a and has_b))

patch_stabs = [s for s in all_stabilizers if s[2] != MERGE]
merge_stab_list = [s for s in all_stabilizers if s[2] == MERGE]

# ============================================================
# 3. Classify qubits
# ============================================================

all_data = [i for i, qt, p, x, y in qubits if qt == DATA]
z_anc = [i for i, qt, p, x, y in qubits if qt == Z_ANC and p != MERGE]
x_anc = [i for i, qt, p, x, y in qubits if qt == X_ANC and p != MERGE]
mz_anc = [i for i, qt, p, x, y in qubits if qt == Z_ANC and p == MERGE]
mx_anc = [i for i, qt, p, x, y in qubits if qt == X_ANC and p == MERGE]

# ============================================================
# 4. Build CX layers (NE, SE, NW, SW from ancilla to data)
# ============================================================
offsets = [(0.5, -0.5), (0.5, 0.5), (-0.5, -0.5), (-0.5, 0.5)]

cx_patch = [[] for _ in range(4)]
cx_merge = [[] for _ in range(4)]

for idx, qt, patch, x, y in qubits:
    if qt == DATA:
        continue
    is_x = (qt == X_ANC)
    is_m = (patch == MERGE)

    for li, (dx, dy) in enumerate(offsets):
        di = find_data_at(x + dx, y + dy)
        if di < 0:
            continue
        if is_x:
            ctrl, tgt = idx, di
        else:
            ctrl, tgt = di, idx
        target = cx_merge[li] if is_m else cx_patch[li]
        target.extend([ctrl, tgt])

# ============================================================
# 5. Generate circuit text
# ============================================================

lines = []

class MeasTracker:
    """Track measurement record indices for detector construction."""
    def __init__(self):
        self.total = 0
        self._m = {}

    def record(self, label, q):
        self._m[(label, q)] = self.total
        self.total += 1

    def rec(self, label, q):
        return f"rec[{self._m[(label, q)] - self.total}]"

mt = MeasTracker()

def fv(v):
    """Format coordinate value."""
    if isinstance(v, float) and v == int(v):
        return str(int(v))
    return str(v)

def ql(indices):
    return " ".join(str(i) for i in indices)

def make_polygon(anc_idx, is_x, data_qs, is_merge=False):
    """Generate a POLYGON pragma with data qubits in clockwise order."""
    if is_merge:
        color = "0,0.7,0,0.3" if is_x else "1,0.5,0,0.3"
    else:
        color = "0,0,1,0.25" if is_x else "1,0,0,0.25"

    ax, ay = qubits[anc_idx][3], qubits[anc_idx][4]
    pos_map = {}
    for dq in data_qs:
        dx = round((qubits[dq][3] - ax) * 2) / 2
        dy = round((qubits[dq][4] - ay) * 2) / 2
        pos_map[(dx, dy)] = dq

    # Clockwise: SW, SE, NE, NW
    cw = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
    ordered = [pos_map[k] for k in cw if k in pos_map]
    return f"#!pragma POLYGON({color}) {ql(ordered)}"

def emit_patch_pragmas():
    for anc, is_x, patch, dqs, cross in patch_stabs:
        lines.append(make_polygon(anc, is_x, dqs))

def emit_merge_pragmas():
    for anc, is_x, patch, dqs, cross in merge_stab_list:
        lines.append(make_polygon(anc, is_x, dqs, True))

# -- QUBIT_COORDS --
for idx, qt, p, x, y in qubits:
    lines.append(f"QUBIT_COORDS({fv(x)}, {fv(y)}) {idx}")
lines.append("")

# -- Initial stabilizer coloring (patches only) --
emit_patch_pragmas()
lines.append("")

# -- Reset all data qubits --
lines.append("TICK")
lines.append(f"R {ql(all_data)}")

# -- Round 0: Pre-merge (patches only, no merge ancillas) --
lines.append("TICK")
lines.append(f"R {ql(z_anc)}")
lines.append(f"RX {ql(x_anc)}")
for li in range(4):
    lines.append("TICK")
    if cx_patch[li]:
        lines.append(f"CX {ql(cx_patch[li])}")
lines.append("TICK")
lines.append(f"M {ql(z_anc)}")
for q in z_anc:
    mt.record('r0', q)
lines.append(f"MX {ql(x_anc)}")
for q in x_anc:
    mt.record('r0', q)

# Round 0 detectors: Z-ancillas only (deterministic after |0> init)
for q in z_anc:
    x, y = qubits[q][3], qubits[q][4]
    lines.append(f"DETECTOR({fv(x)}, {fv(y)}, 0) {mt.rec('r0', q)}")
lines.append("TICK")

# -- Merge rounds --
for r in range(num_merge_rounds):
    rl = f"m{r}"
    pl = 'r0' if r == 0 else f"m{r-1}"

    # Pragmas showing stabilizers for this round
    emit_patch_pragmas()
    emit_merge_pragmas()
    lines.append("")

    lines.append("TICK")
    all_z = z_anc + mz_anc
    all_x = x_anc + mx_anc
    lines.append(f"R {ql(all_z)}")
    lines.append(f"RX {ql(all_x)}")

    for li in range(4):
        lines.append("TICK")
        combined = cx_patch[li] + cx_merge[li]
        if combined:
            lines.append(f"CX {ql(combined)}")
        if interconnect_error > 0 and cx_merge[li]:
            lines.append(f"DEPOLARIZE2({interconnect_error}) {ql(cx_merge[li])}")

    lines.append("TICK")
    lines.append(f"M {ql(all_z)}")
    for q in all_z:
        mt.record(rl, q)
    lines.append(f"MX {ql(all_x)}")
    for q in all_x:
        mt.record(rl, q)

    # Patch Z-ancilla detectors (compare with previous round)
    for q in z_anc:
        x, y = qubits[q][3], qubits[q][4]
        lines.append(f"DETECTOR({fv(x)}, {fv(y)}, {r+1}) {mt.rec(rl, q)} {mt.rec(pl, q)}")

    # Merge Z-ancilla detectors (from merge round 1 onward)
    if r >= 1:
        pml = f"m{r-1}"
        for q in mz_anc:
            x, y = qubits[q][3], qubits[q][4]
            lines.append(f"DETECTOR({fv(x)}, {fv(y)}, {r+1}) {mt.rec(rl, q)} {mt.rec(pml, q)}")

    # Patch X-ancilla detectors
    for q in x_anc:
        x, y = qubits[q][3], qubits[q][4]
        lines.append(f"DETECTOR({fv(x)}, {fv(y)}, {r+1}) {mt.rec(rl, q)} {mt.rec(pl, q)}")

    # Merge X-ancilla detectors (from merge round 1 onward)
    if r >= 1:
        pml = f"m{r-1}"
        for q in mx_anc:
            x, y = qubits[q][3], qubits[q][4]
            lines.append(f"DETECTOR({fv(x)}, {fv(y)}, {r+1}) {mt.rec(rl, q)} {mt.rec(pml, q)}")

    lines.append("TICK")

# -- Final round: Post-merge (patches only) --
fl = 'final'
lml = f"m{num_merge_rounds - 1}"

emit_patch_pragmas()
lines.append("")

lines.append("TICK")
lines.append(f"R {ql(z_anc)}")
lines.append(f"RX {ql(x_anc)}")
for li in range(4):
    lines.append("TICK")
    if cx_patch[li]:
        lines.append(f"CX {ql(cx_patch[li])}")
lines.append("TICK")
lines.append(f"M {ql(z_anc)}")
for q in z_anc:
    mt.record(fl, q)
lines.append(f"MX {ql(x_anc)}")
for q in x_anc:
    mt.record(fl, q)

final_round = num_merge_rounds + 1
for q in z_anc:
    x, y = qubits[q][3], qubits[q][4]
    lines.append(f"DETECTOR({fv(x)}, {fv(y)}, {final_round}) {mt.rec(fl, q)} {mt.rec(lml, q)}")
for q in x_anc:
    x, y = qubits[q][3], qubits[q][4]
    lines.append(f"DETECTOR({fv(x)}, {fv(y)}, {final_round}) {mt.rec(fl, q)} {mt.rec(lml, q)}")
lines.append("TICK")

# -- Final data measurement --
dl = 'data'
lines.append(f"M {ql(all_data)}")
for q in all_data:
    mt.record(dl, q)

data_round = final_round + 1
for anc, is_x, patch, dqs, cross in patch_stabs:
    if is_x:
        continue  # Only Z-type stabilizers for Z-basis data measurement
    ax_c, ay_c = qubits[anc][3], qubits[anc][4]
    targets = [mt.rec(dl, dq) for dq in dqs] + [mt.rec(fl, anc)]
    lines.append(f"DETECTOR({fv(ax_c)}, {fv(ay_c)}, {data_round}) {' '.join(targets)}")

# -- Logical observables --
# Patch A Z-logical: leftmost data column (x=0.5)
obs_a = [mt.rec(dl, q) for q in all_data
         if qubits[q][2] == PATCH_A and abs(qubits[q][3] - 0.5) < 0.1]
if obs_a:
    lines.append(f"OBSERVABLE_INCLUDE(0) {' '.join(obs_a)}")

# Patch B Z-logical: rightmost data column
max_bx = max(qubits[q][3] for q in all_data if qubits[q][2] == PATCH_B)
obs_b = [mt.rec(dl, q) for q in all_data
         if qubits[q][2] == PATCH_B and abs(qubits[q][3] - max_bx) < 0.1]
if obs_b:
    lines.append(f"OBSERVABLE_INCLUDE(1) {' '.join(obs_b)}")

# ============================================================
# 6. Write output
# ============================================================
outpath = os.path.join(os.path.dirname(__file__), '..', 'output', 'distributed_ls_d5_pragma.stim')
outpath = os.path.abspath(outpath)
os.makedirs(os.path.dirname(outpath), exist_ok=True)

with open(outpath, 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"Generated {outpath}")
print(f"  Qubits: {len(qubits)} ({len(all_data)} data, {len(z_anc)+len(x_anc)} patch anc, {len(mz_anc)+len(mx_anc)} merge anc)")
print(f"  Patch stabilizers: {len(patch_stabs)}")
print(f"  Merge stabilizers: {len(merge_stab_list)}")
print(f"  Total measurements: {mt.total}")
print(f"  Lines: {len(lines)}")

# Print stabilizer summary
for anc, is_x, patch, dqs, cross in merge_stab_list:
    x, y = qubits[anc][3], qubits[anc][4]
    print(f"  Merge {'X' if is_x else 'Z'} at ({fv(x)},{fv(y)}) weight={len(dqs)} cross={cross}")
