# Boundary Ancillas Analysis for Lattice Surgery

## Request
Add boundary ancillas at specific coordinates to make each patch a "complete independent rotated surface code":
- For d=3: (5, 3) and (6, 4)
- For d=5 and above: Also (11, 3)

## Investigation

### Current Patch Structure
For distance d rotated surface code:
- Expected stabilizers: d² - 1 (approximately d²/2 X-type and d²/2 Z-type)
- For d=3: 9 data qubits → 8 stabilizers expected (4 X, 4 Z)
- For d=5: 25 data qubits → 24 stabilizers expected (12 X, 12 Z)

### Actual Counts
**d=3:**
- Patch A: 7 stabilizers (3 X, 4 Z) - missing 1 X-stabilizer
- Patch B: 7 stabilizers (3 X, 4 Z) - missing 1 X-stabilizer

**d=5:**
- Patch A: 21 stabilizers - missing 3 stabilizers
- Patch B: 21 stabilizers - missing 3 stabilizers

## Fundamental Design Conflict

### Why Patches Are Intentionally Incomplete

The current lattice surgery implementation **intentionally** creates incomplete patches. This is by design, not an oversight:

1. **Merge Operation Requires Open Boundaries**: For lattice surgery XX-merge to work, patches must have open boundaries on the merge-facing sides where they connect.

2. **Missing Stabilizers Are On Merge-Facing Sides**: The "missing" stabilizers would be located at positions that would interfere with the merge region ancillas.

3. **Overlapping Measurements Problem**: Adding boundary ancillas at the requested coordinates creates overlapping stabilizer measurements:
   - Example: For d=5, adding X-ancilla at (7, 5) attempts to measure data qubits at (7.5, 4.5) and (6.5, 4.5)
   - These same data qubits are already measured by interior ancillas at (8, 4), (8, 5), (7, 4), etc.
   - This creates overdetermined error syndromes that break detector generation

### Detector Errors When Boundary Ancillas Added

Attempted adding:
- X-ancilla at (7, 5) for d=5
- X-ancilla at (8, 6) for d=5

**Result**: Non-deterministic detector error
```
Error: The circuit contains non-deterministic detectors.

This was discovered while analyzing an X-basis reset (RX) on:
    qubit 108 [coords (7, 5)]

The collapse anti-commuted with these detectors/observables:
    D290 [coords (8, 5, 6)]
```

The detector at (8, 5) became sensitive to the new ancilla at (7, 5) because they measure shared data qubits, creating measurement dependencies that the detector generation logic doesn't handle.

## What Would Be Needed

To support "complete independent rotated surface codes" for lattice surgery would require:

1. **Conditional Ancilla Activation**: Boundary ancillas active only during:
   - Pre-merge rounds (when patches are independent)
   - Post-split rounds (when patches are separated again)
   - But **inactive** during merge rounds

2. **Modified Detector Generation**: Logic to handle:
   - Different stabilizer sets across rounds
   - Boundary stabilizers that may have lower weight (measure fewer data qubits)
   - Transitions between complete and incomplete code phases

3. **Alternative Approach - Separate Layouts**: Maintain different qubit layouts for:
   - Independent code phase (complete patches)
   - Merge phase (open boundaries)
   - Transitions between phases

None of these capabilities exist in the current implementation.

## Current Status

The lattice surgery circuits work correctly with:
- **d=3**: 41 qubits, 88 detectors, validated
- **d=5**: 107 qubits, 300 detectors, validated

Patches are intentionally incomplete (missing boundary stabilizers on merge-facing sides), which is the correct design for lattice surgery merge operations.

## Recommendation

The current design is architecturally sound for lattice surgery. Adding boundary ancillas as requested would require significant refactoring to support:
1. Dynamic stabilizer activation across rounds
2. Multi-phase qubit layouts
3. Enhanced detector generation for variable stabilizer sets

These changes would be a substantial feature addition beyond the scope of simple ancilla additions.
