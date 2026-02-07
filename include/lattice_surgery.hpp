#pragma once

#include "config.hpp"
#include <stim.h>
#include <cstdint>
#include <vector>
#include <map>

namespace bucket_sim {

// Qubit types in the lattice surgery layout
enum class QubitType {
    DATA,           // Data qubit
    X_ANCILLA,      // X-stabilizer ancilla
    Z_ANCILLA,      // Z-stabilizer ancilla
    MERGE_ANCILLA   // Ancilla for merge boundary stabilizers
};

// Which patch a qubit belongs to
enum class PatchID {
    PATCH_A,    // Left patch
    PATCH_B,    // Right patch
    BOUNDARY    // Merge boundary (shared)
};

// Information about a qubit in the layout
struct QubitInfo {
    uint32_t index;
    QubitType type;
    PatchID patch;
    double x, y;     // Coordinates
    bool is_boundary; // True if on merge boundary
};

// Stabilizer information
struct Stabilizer {
    bool is_x_type;                     // X or Z stabilizer
    std::vector<uint32_t> data_qubits;  // Data qubits involved
    uint32_t ancilla;                   // Measurement ancilla
    PatchID patch;                      // Which patch (or BOUNDARY for merge stabilizers)
};

class LatticeSurgeryCircuit {
private:
    Config config_;
    uint32_t distance_;
    uint32_t merge_rounds_;
    MergeType merge_type_;
    bool distributed_;

    // Qubit layout
    std::vector<QubitInfo> qubits_;
    std::map<uint32_t, QubitInfo> qubit_map_;

    // Stabilizers
    std::vector<Stabilizer> patch_a_stabilizers_;
    std::vector<Stabilizer> patch_b_stabilizers_;
    std::vector<Stabilizer> merge_stabilizers_;  // Boundary stabilizers during merge

    // Circuit
    stim::Circuit circuit_;

    // Layout helpers
    void initialize_layout();
    void initialize_rotated_layout();
    uint32_t add_qubit(QubitType type, PatchID patch, double x, double y, bool boundary = false);
    int32_t find_qubit_at(double x, double y, double tolerance = 0.1) const;
    void build_stabilizers();
    void add_patch_stabilizers(PatchID patch, double x_offset);

    // Circuit generation helpers
    void generate_circuit();
    void generate_d3_circuit();      // Hardcoded d=3 circuit matching reference
    void generate_general_circuit(); // General algorithm for any distance
    void add_reset_layer();
    void add_stabilizer_round(bool is_merge_round, uint32_t round_num);
    void add_stabilizer_round_v2(bool is_merge_round, uint32_t round_num, size_t& measurement_count);
    void add_cnot_layer(const std::vector<Stabilizer>& stabilizers, bool first_half);
    void add_measurement_layer(const std::vector<Stabilizer>& stabilizers);
    void add_detector(uint32_t ancilla, uint32_t round, bool is_x_type);
    void add_final_measurements();
    void add_logical_observable();

    // Error injection
    void add_data_noise();
    void add_measurement_noise();
    void add_reset_noise();
    void add_distributed_noise(const std::vector<Stabilizer>& stabilizers);

    // Coordinate calculation
    std::pair<double, double> data_qubit_coords(uint32_t row, uint32_t col, PatchID patch);
    std::pair<double, double> ancilla_coords(uint32_t row, uint32_t col, PatchID patch, bool is_x);

public:
    LatticeSurgeryCircuit(const Config& config);

    // Generate the full circuit
    stim::Circuit generate();

    // Accessors
    size_t num_qubits() const { return qubits_.size(); }
    size_t num_data_qubits() const;
    size_t num_detectors() const;
    const std::vector<QubitInfo>& get_qubits() const { return qubits_; }
    const std::vector<Stabilizer>& get_merge_stabilizers() const { return merge_stabilizers_; }

    // Check if a stabilizer crosses QPU boundary (for distributed mode)
    bool stabilizer_crosses_boundary(const Stabilizer& stab) const;
};

} // namespace bucket_sim
