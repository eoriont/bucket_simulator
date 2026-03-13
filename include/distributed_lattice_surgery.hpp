#pragma once

#include "config.hpp"
#include <stim.h>
#include <cstdint>
#include <vector>
#include <map>
#include <unordered_set>

namespace bucket_sim {

// Qubit types in the distributed lattice surgery layout
enum class DQubitType {
    DATA,           // Data qubit
    X_ANCILLA,      // X-stabilizer ancilla
    Z_ANCILLA,      // Z-stabilizer ancilla
};

// Which patch a qubit belongs to
enum class DPatchID {
    PATCH_A,    // Left patch (QPU 1)
    PATCH_B,    // Right patch (QPU 2)
    SEAM_A,     // Seam boundary X ancillas owned by Patch A (at x=d)
    SEAM_B,     // Seam boundary X ancillas owned by Patch B (at x=d)
    MERGE,      // Merge boundary ancillas (on interconnect)
};

// Information about a qubit in the layout
struct DQubitInfo {
    uint32_t index;
    DQubitType type;
    DPatchID patch;
    double x, y;     // Coordinates
};

// Stabilizer information
struct DStabilizer {
    bool is_x_type;                     // X or Z stabilizer
    std::vector<uint32_t> data_qubits;  // Data qubits involved
    uint32_t ancilla;                   // Measurement ancilla
    DPatchID patch;                     // Which patch (or MERGE for merge stabilizers)
    bool crosses_boundary;              // True if touches data on both patches
};

// Distributed lattice surgery: two patches side-by-side, Patch B flipped,
// merge uses a single column of ancillas bridging existing data qubits
// via remote CNOTs. No new data qubits are introduced during merge.
class DistributedLatticeSurgeryCircuit {
private:
    Config config_;
    uint32_t distance_;
    uint32_t merge_rounds_;
    MergeType merge_type_;

    // Qubit layout
    std::vector<DQubitInfo> qubits_;

    // Stabilizers
    std::vector<DStabilizer> patch_a_stabilizers_;
    std::vector<DStabilizer> patch_b_stabilizers_;
    std::vector<DStabilizer> seam_a_stabilizers_;  // Seam X stabs owned by Patch A
    std::vector<DStabilizer> seam_b_stabilizers_;  // Seam X stabs owned by Patch B
    std::vector<DStabilizer> merge_stabilizers_;   // Boundary stabilizers during merge

    // Removed data qubits (data qubit removal superstabilizer approach)
    std::unordered_set<uint32_t> removed_data_qubits_;
    // Merge-Z ancillas suppressed during merge (weight-1 after data qubit removal)
    std::unordered_set<uint32_t> merge_suppressed_z_ancillas_;

    // Circuit
    stim::Circuit circuit_;

    // Layout
    void initialize_layout();
    uint32_t add_qubit(DQubitType type, DPatchID patch, double x, double y);

    // Stabilizer construction
    void build_stabilizers();

    // Circuit generation
    void generate_circuit();
    void generate_general_circuit();

public:
    DistributedLatticeSurgeryCircuit(const Config& config);

    // Generate the full circuit
    stim::Circuit generate();

    // Accessors
    size_t num_data_qubits() const;
    size_t num_detectors() const;
    const std::vector<DQubitInfo>& qubits() const { return qubits_; }
    const std::vector<DStabilizer>& patch_a_stabs() const { return patch_a_stabilizers_; }
    const std::vector<DStabilizer>& patch_b_stabs() const { return patch_b_stabilizers_; }
    const std::vector<DStabilizer>& seam_a_stabs() const { return seam_a_stabilizers_; }
    const std::vector<DStabilizer>& seam_b_stabs() const { return seam_b_stabilizers_; }
    const std::vector<DStabilizer>& merge_stabs() const { return merge_stabilizers_; }

    // Generate stim text with #!pragma POLYGON annotations for stabilizer coloring
    std::string annotated_stim_str() const;
};

} // namespace bucket_sim
