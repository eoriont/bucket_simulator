/**
 * Ported to cpp from the Lin et al paper
 * Codesign of quantum error-correcting codes and modular chiplets in the presence of defects
 */
#include <cassert>
#include <cstdint>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <algorithm>
#include <string_view>
using namespace std::literals;
#include <fstream>
#include <iostream>
#include <string_view>
using namespace std::literals;

#include "stim/circuit/circuit.h"
#include "stim/circuit/gate_target.h"
#include "superstabilizer_circuit.hpp"

namespace superstabilizer {
    // ---------- Graph that can count # of shortest paths via BFS (unweighted, undirected) ----------
    struct Graph {
        // node ids are int64_t (you use -1 and -2 as virtual boundary nodes)
        std::unordered_map<int64_t, std::vector<int64_t>> adj;
        void add_node(int64_t v) {
            adj.emplace(v, std::vector<int64_t>{});
        }
        void add_edge(int64_t u, int64_t v) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }
        // returns {num_shortest_paths, length}
        std::pair<uint64_t,int> bfs_shortest_paths(int64_t start, int64_t target) const {
            std::unordered_map<int64_t,int> dist;
            std::unordered_map<int64_t,uint64_t> ways;
            std::deque<int64_t> q;

            dist[start] = 0;
            ways[start] = 1;
            q.push_back(start);

            while (!q.empty()) {
                int64_t u = q.front(); q.pop_front();
                if (u == target) return {ways.at(u), dist.at(u)};
                auto it = adj.find(u);
                if (it == adj.end()) continue;
                for (int64_t v : it->second) {
                    if (!dist.count(v)) {
                        dist[v] = dist[u] + 1;
                        ways[v] = ways[u];
                        q.push_back(v);
                    } else if (dist[v] == dist[u] + 1) {
                        ways[v] += ways[u];
                    }
                }
            }
            return {0, 0};
        }
    };

    // ---------- DSU / Union-Find ----------
    struct DSU {
        std::vector<int> p, r;
        explicit DSU(int n): p(n), r(n,1) { for (int i=0;i<n;i++) p[i]=i; }
        int find(int x){ return p[x]==x?x:p[x]=find(p[x]); }
        bool unite(int a,int b){
            a=find(a); b=find(b);
            if(a==b) return false;
            if(r[a]<r[b]) std::swap(a,b);
            p[b]=a;
            if(r[a]==r[b]) r[a]++;
            return true;
        }
        bool connected(int a,int b){ return find(a)==find(b); }
        // O(n^2) components extraction (fine for moderate n)
        std::vector<std::vector<int>> components() {
            std::vector<std::vector<int>> out;
            std::vector<int> reps(p.size());
            for (int i=0;i<(int)p.size();++i) reps[i]=find(i);
            std::unordered_map<int,int> idx;
            for (int i=0;i<(int)p.size();++i){
                int r=reps[i];
                if(!idx.count(r)){ idx[r]=(int)out.size(); out.push_back({}); }
                out[idx[r]].push_back(i);
            }
            return out;
        }
    };

    // ---------- Basic data types mirroring the Python classes ----------
    struct DataQubit {
        int name;
        std::pair<int,int> coords;  // (x,y)
        std::string repr() const {
            return std::to_string(name) + ", Coords: (" + std::to_string(coords.first) + "," + std::to_string(coords.second) + ")";
        }
    };

    struct MeasureQubit {
        int name;
        std::pair<int,int> coords;
        // up to 4 neighbors; nulls in Python -> use -1 here
        std::array<int,4> data_qubits; // names of data qubits, or -1
        char basis; // 'X' or 'Z'
        std::string repr() const {
            return "|" + std::to_string(name) + ", Coords: (" + std::to_string(coords.first) + "," + std::to_string(coords.second) + "), Basis: " + basis + "|";
        }
    };

    struct Defect {
        // cluster of disabled data qubits
        std::vector<int> names;                     // data-qubit names
        std::vector<std::pair<int,int>> coords;     // data-qubit coords
        std::vector<MeasureQubit> x_gauges;
        std::vector<MeasureQubit> z_gauges;
        int diameter=0, horizontal_len=0, vertical_len=0;

        Defect() = default;
        Defect(std::vector<int> n, std::vector<std::pair<int,int>> c,
            std::vector<MeasureQubit> xg, std::vector<MeasureQubit> zg)
            : names(std::move(n)), coords(std::move(c)), x_gauges(std::move(xg)), z_gauges(std::move(zg)) {
            int max_x=0,max_y=0,min_x=INT32_MAX,min_y=INT32_MAX;
            for (auto &xy: coords){
                max_x=std::max(max_x, xy.first);
                max_y=std::max(max_y, xy.second);
                min_x=std::min(min_x, xy.first);
                min_y=std::min(min_y, xy.second);
            }
            diameter = (int)(std::max(max_x - min_x, max_y - min_y)/2) + 1;
            horizontal_len = (int)((max_x - min_x)/2) + 1;
            vertical_len   = (int)((max_y - min_y)/2) + 1;
        }
        void change_diameter(int d){ diameter = d; }
    };

    // Helper: convert qubit indices to GateTargets.
    static inline std::vector<uint32_t> qubit_targets(const std::vector<int> &qs){
        std::vector<uint32_t> t;
        t.reserve(qs.size());
        for (int q: qs) t.push_back(stim::GateTarget::qubit(q).data);
        return t;
    }
    static inline std::vector<stim::GateTarget> rec_targets(const std::vector<int64_t> &rec_lookbacks){
        std::vector<stim::GateTarget> t;
        t.reserve(rec_lookbacks.size());
        for (int64_t r: rec_lookbacks) t.push_back(stim::GateTarget::rec((int32_t)r));
        return t;
    }

    // ---------- LogicalQubit skeleton ----------
    struct LogicalQubit {
        int d;
        double readout_err, gate1_err, gate2_err;
        bool percolated=false, too_many_qubits_lost=false;
        int vertical_distance=0, horizontal_distance=0;
        int num_inactive_data=0, num_inactive_syn=0, num_data_superstabilizer=0;

        // geometry / bookkeeping
        std::vector<DataQubit> data;               // active data qubits
        std::vector<MeasureQubit> x_ancilla, z_ancilla;
        std::vector<MeasureQubit> x_gauges, z_gauges;
        std::vector<Defect> defect;
        std::vector<int> all_qubits;               // active qubit names (data + syn)
        std::array<double,4> boundary_deformation{}; // L,R,Top,Bottom
        std::vector<int> observable;               // names of data qubits along logical Z (top->bottom) path
        // dynamic boundaries: [left, right, top, bottom], each a list of data-qubit names
        std::array<std::vector<int>,4> dynamic_boundaries{};
        // measurement bookkeeping
        // Each round: map qubit-> lookback index (negative offsets)
        std::vector<std::unordered_map<int,int64_t>> meas_record;
        int64_t total_meas_count = 0;

        LogicalQubit(int d_,
                    double readout_err_, double gate1_err_, double gate2_err_,
                    const std::vector<std::pair<int,int>> &missing_coords = {},
                    double data_loss_tolerance = 1.0,
                    bool verbose=false,
                    bool get_metrics=false)
            : d(d_), readout_err(readout_err_), gate1_err(gate1_err_), gate2_err(gate2_err_) {

            // --- small helper to pack 2D coords into a hashable key ---
            auto pack = [](int x, int y) -> int64_t {
                return ( (int64_t)x << 32 ) ^ ( (int64_t)y & 0xffffffffLL );
            };

            // --- locals that mirror your Python variables ---
            const int DD = d * d;

            // Active(0) / Disabled(1) / Deleted(-1). We'll start at all 0 then mark as we go.
            // Size matches Python: 2*d^2 - 1 (data + syndrome).
            std::vector<int8_t> is_disabled(2 * d * d - 1, 0);

            // Data lookup tables.
            std::vector<DataQubit> data_list(DD);                    // index -> DataQubit
            std::vector<std::vector<int>> data_matching(2 * d, std::vector<int>(2 * d, -1)); // (2x,2y) -> data name

            // Syndrome ID by coord; syn_info and data_info like in Python.
            std::unordered_map<int64_t, int> syndrome_matching;      // (2x+1,2y+1)-> syn name
            std::vector<MeasureQubit> syn_info;                      // all syndrome qubits
            std::vector<std::vector<int>> data_info(DD);             // data_name -> list of neighbor syn names

            // Clear edges we'll build (left, right, top, bottom)
            std::array<std::vector<int>, 4> edges;  // 0:L,1:R,2:T,3:B

            // Fast membership for missing coords
            std::unordered_set<int64_t> missing_set;
            missing_set.reserve(missing_coords.size() * 2 + 1);
            for (auto &c : missing_coords) {
                // coords are *physical* even grid (e.g., (4,4)), as used below
                auto k = pack(c.first, c.second);
                missing_set.insert(k);
            }

            // Build data qubits and initial dynamic boundaries (like Python)
            int disabled_count = 0;
            for (int x = 0; x < d; ++x) {
                for (int y = 0; y < d; ++y) {
                    int name = d * x + y;
                    std::pair<int,int> coords = {2 * x, 2 * y};

                    // Boundaries by *grid indices*, matching Python
                    if (x == 0)        edges[2].push_back(name); // top
                    else if (x == d-1) edges[3].push_back(name); // bottom
                    if (y == 0)        edges[0].push_back(name); // left
                    else if (y == d-1) edges[1].push_back(name); // right

                    // Disable if listed
                    if (missing_set.count(pack(coords.first, coords.second))) {
                        is_disabled[name] = 1;  // needs superstabilizer
                        ++disabled_count;
                    }

                    // Record data qubit
                    data_list[name] = DataQubit{name, coords};
                    data_matching[coords.first][coords.second] = name;
                }
            }


            // Move into the member to avoid self-assign bugs
            this->dynamic_boundaries = std::move(edges);

            // ------------------------
            // Build syndrome (X / Z)
            // ------------------------
            int q = DD; // syndrome indices start at d*d (like Python)
            syn_info.clear();

            for (int x = -1; x < d; ++x) {
                for (int y = -1; y < d; ++y) {
                    // X syndrome: (x+y) odd, and x != -1, x != d-1
                    if ( ((x + y) & 1) && x != -1 && x != d - 1 ) {
                        int sx = 2 * x + 1, sy = 2 * y + 1;
                        syndrome_matching[pack(sx, sy)] = q;

                        MeasureQubit mq;
                        mq.name = q;
                        mq.coords = {sx, sy};
                        mq.basis = 'X';
                        mq.data_qubits = {-1, -1, -1, -1};

                        // neighbors: up to 4 data qubits (use -1 when Python had None)
                        if (y != d - 1) {
                            mq.data_qubits[0] = data_matching[sx + 1][sy + 1];
                            mq.data_qubits[1] = data_matching[sx - 1][sy + 1];
                        }
                        if (y != -1) {
                            mq.data_qubits[2] = data_matching[sx + 1][sy - 1];
                            mq.data_qubits[3] = data_matching[sx - 1][sy - 1];
                        }

                        // Map this syn to each present neighbor data
                        for (int k = 0; k < 4; ++k) {
                            int dn = mq.data_qubits[k];
                            if (dn != -1) data_info[dn].push_back(q);
                        }

                        if (missing_set.count(pack(sx, sy))) {
                            is_disabled[q] = 1; // defective syndrome
                        }

                        syn_info.push_back(mq);
                        ++q;

                    // Z syndrome: (x+y) even, and y != -1, y != d-1
                    } else if ( ((x + y) & 1) == 0 && y != -1 && y != d - 1 ) {
                        int sx = 2 * x + 1, sy = 2 * y + 1;
                        syndrome_matching[pack(sx, sy)] = q;

                        MeasureQubit mq;
                        mq.name = q;
                        mq.coords = {sx, sy};
                        mq.basis = 'Z';
                        mq.data_qubits = {-1, -1, -1, -1};

                        if (x != d - 1) {
                            mq.data_qubits[0] = data_matching[sx + 1][sy + 1];
                            mq.data_qubits[1] = data_matching[sx + 1][sy - 1];
                        }
                        if (x != -1) {
                            mq.data_qubits[2] = data_matching[sx - 1][sy + 1];
                            mq.data_qubits[3] = data_matching[sx - 1][sy - 1];
                        }

                        for (int k = 0; k < 4; ++k) {
                            int dn = mq.data_qubits[k];
                            if (dn != -1) data_info[dn].push_back(q);
                        }

                        if (missing_set.count(pack(sx, sy))) {
                            is_disabled[q] = 1; // defective syndrome
                        }

                        syn_info.push_back(mq);
                        ++q;
                    }
                }
            }

            // (Optional sanity check; Python expects total = 2*d^2 - 1)
            #ifndef NDEBUG
            if (q != 2 * d * d - 1) {
                std::cerr << "[warn] total qubits q=" << q << " != 2*d^2-1=" << (2*d*d - 1) << "\n";
            }
            #endif

            // -------------------------------
            // Corner list (same indices as py)
            // -------------------------------
            std::array<int, 8> corners = {
                0,
                DD + d/2,
                d - 1,
                DD + d/2 - 1,
                DD - d,
                2*DD - 1 - d/2,
                DD - 1,
                2*DD - 2 - d/2
            };

            // -------------------------------
            // Helpers for boundary detection
            // -------------------------------
            auto is_boundary = [&](int qname) -> bool {
                if (qname < DD) {
                    // data
                    auto [x,y] = data_list[qname].coords;
                    return x==0 || x==2*(d-1) || y==0 || y==2*(d-1);
                } else {
                    // syndrome (weight-2 live at -1 / 2d-1 “ghost” coords in one axis)
                    const auto &mq = syn_info[qname - DD];
                    int x = mq.coords.first;
                    int y = mq.coords.second;
                    return x==-1 || x==2*d-1 || y==-1 || y==2*d-1;
                }
            };

            auto near_blue_boundary = [&](int qname) -> bool {
                // blue edges are x==1 or x==2d-3 in syndrome coords
                const auto &mq = syn_info[qname - DD];
                int x = mq.coords.first;
                return x==1 || x==2*d-3;
            };

            auto near_red_boundary = [&](int qname) -> bool {
                // red edges are y==1 or y==2d-3 in syndrome coords
                const auto &mq = syn_info[qname - DD];
                int y = mq.coords.second;
                return y==1 || y==2*d-3;
            };

            // ---------------------------------------
            // Handle a broken DATA on the orig border
            // Delete: one weight-2 syn (boundary) and
            //         its two adjacent data; then also
            //         the one weight-4 syn consistent.
            // ---------------------------------------
            auto handle_boundary_data = [&](int qname) {
                // Given a *data* qubit on original boundary, delete the 4 around it.
                //  - find its 3 syndromes: 2 weight-4 (interior), 1 weight-2 (boundary)
                std::vector<int> weight4_syn;
                int weight2_syn = -1;
                for (int synq : data_info[qname]) {
                    if (!is_boundary(synq)) weight4_syn.push_back(synq);
                    else                     weight2_syn = synq;
                }
                assert((int)weight4_syn.size()==2 && weight2_syn!=-1);

                // delete the boundary weight-2 syn
                is_disabled[weight2_syn] = -1;
                if (verbose) std::cerr << "Delete weight-2 syndrome " << weight2_syn << "\n";

                // delete its two adjacent data qubits
                std::vector<int> data_to_delete;
                for (int k=0;k<4;k++){
                    int dn = syn_info[weight2_syn - DD].data_qubits[k];
                    if (dn != -1) {
                        is_disabled[dn] = -1;
                        data_to_delete.push_back(dn);
                    }
                }
                if (verbose) {
                    std::cerr << "Delete data qubits [";
                    for (size_t i=0;i<data_to_delete.size();++i){
                        if (i) std::cerr << ", ";
                        std::cerr << data_to_delete[i];
                    }
                    std::cerr << "]\n";
                }

                // of the two weight-4 syn, delete the one that touched both deleted data
                auto touches_both = [&](int synq)->bool{
                    std::array<int,4> neigh = syn_info[synq - DD].data_qubits;
                    int count=0;
                    for (int dn : neigh) if (dn!=-1 && (dn==data_to_delete[0] || dn==data_to_delete[1])) count++;
                    return count==2;
                };
                if (touches_both(weight4_syn[0])) {
                    is_disabled[weight4_syn[0]] = -1;
                    if (verbose) std::cerr << "Delete weight-4 syndrome " << weight4_syn[0] << "\n";
                } else {
                    is_disabled[weight4_syn[1]] = -1;
                    if (verbose) std::cerr << "Delete weight-4 syndrome " << weight4_syn[1] << "\n";
                }
            };

            // -----------------------------------------------------------
            // Handle a broken *qubit* (data or syn) on/near orig boundary
            // Returns true if it actually handled something (i.e. deleted)
            // -----------------------------------------------------------
            auto handle_boundary_defects = [&](int qname) -> bool {
                if (qname < DD) {
                    // data qubit
                    if (is_boundary(qname)) {
                        handle_boundary_data(qname);
                        return true;
                    }
                    return false;
                } else {
                    // syndrome qubit
                    const auto &mq = syn_info[qname - DD];
                    int x = mq.coords.first, y = mq.coords.second;
                    char syn_basis = mq.basis; // 'X' or 'Z'

                    if (is_boundary(qname)) {
                        // weight-2 syndrome: delete it and the 4-qubit block via its adjacent data
                        // (pick one of the adjacent data; either triggers same 4-qubit pattern)
                        for (int k=0;k<4;k++){
                            int dn = mq.data_qubits[k];
                            if (dn != -1) { handle_boundary_data(dn); break; }
                        }
                        return true;
                    }

                    // Near boundary (one step inside)
                    auto near_boundary = [&](int q)->bool { return near_blue_boundary(q) || near_red_boundary(q); };
                    if (!near_boundary(qname)) return false;

                    // delete this near-boundary syn (always):
                    is_disabled[qname] = -1;
                    if (verbose) std::cerr << "Delete near-boundary syndrome " << qname << "\n";

                    // Easy cases:
                    // - X basis and near blue boundary (x==1 or 2d-3): delete neighbor weight-2 at (x1,y)
                    // - Z basis and near red boundary  (y==1 or 2d-3): delete neighbor weight-2 at (x,y1)
                    bool easy=false;
                    int neighbor_weight2 = -1;
                    if (syn_basis=='X' && near_blue_boundary(qname)) {
                        int x1 = (x==1) ? -1 : (2*d - 1);
                        auto it = syndrome_matching.find(((int64_t)x1<<32) ^ ( (int64_t)y & 0xffffffffLL ));
                        if (it != syndrome_matching.end()) {
                            neighbor_weight2 = it->second;
                            easy = true;
                        }
                    } else if (syn_basis=='Z' && near_red_boundary(qname)) {
                        int y1 = (y==1) ? -1 : (2*d - 1);
                        auto it = syndrome_matching.find(((int64_t)x<<32) ^ ( (int64_t)y1 & 0xffffffffLL ));
                        if (it != syndrome_matching.end()) {
                            neighbor_weight2 = it->second;
                            easy = true;
                        }
                    }

                    if (easy) {
                        // delete the boundary weight-2 and its two adjacent data
                        is_disabled[neighbor_weight2] = -1;
                        if (verbose) std::cerr << "Delete boundary weight-2 syndrome " << neighbor_weight2 << "\n";
                        std::vector<int> data_to_delete;
                        for (int k=0;k<4;k++){
                            int dn = syn_info[neighbor_weight2 - DD].data_qubits[k];
                            if (dn != -1) {
                                is_disabled[dn] = -1;
                                data_to_delete.push_back(dn);
                            }
                        }
                        if (verbose) {
                            std::cerr << "Delete data qubits [";
                            for (size_t i=0;i<data_to_delete.size();++i){
                                if (i) std::cerr << ", ";
                                std::cerr << data_to_delete[i];
                            }
                            std::cerr << "]\n";
                        }
                        return true;
                    }

                    // Hard cases (same color as boundary, close to corner):
                    // Delete a set of nearby same/different-color syns as in Python logic.
                    // Build the exact neighbor list by position:
                    std::vector<int> syn_to_delete;
                    auto push_if_live = [&](int cx, int cy){
                        auto it = syndrome_matching.find(((int64_t)cx<<32) ^ ( (int64_t)cy & 0xffffffffLL ));
                        if (it != syndrome_matching.end()) {
                            int s = it->second;
                            if (is_disabled[s] != -1) syn_to_delete.push_back(s);
                        }
                    };

                    if (x == 1) { // upper edge
                        push_if_live(x, y-2);
                        push_if_live(x, y+2);
                        push_if_live(x+2, y);
                        push_if_live(x-2, y-2);
                        push_if_live(x-2, y+2);
                    } else if (x == 2*d - 3) { // lower edge
                        push_if_live(x, y-2);
                        push_if_live(x, y+2);
                        push_if_live(x-2, y);
                        push_if_live(x+2, y-2);
                        push_if_live(x+2, y+2);
                    } else if (y == 1) { // left edge
                        push_if_live(x-2, y);
                        push_if_live(x+2, y);
                        push_if_live(x, y+2);
                        push_if_live(x-2, y-2);
                        push_if_live(x+2, y-2);
                    } else {
                        // right edge: y == 2*d-3
                        push_if_live(x-2, y);
                        push_if_live(x+2, y);
                        push_if_live(x, y-2);
                        push_if_live(x-2, y+2);
                        push_if_live(x+2, y+2);
                    }

                    if (verbose) {
                        std::cerr << "Hard near-boundary case " << qname << " delete syns [";
                        for (size_t i=0;i<syn_to_delete.size();++i){
                            if (i) std::cerr << ", ";
                            std::cerr << syn_to_delete[i];
                        }
                        std::cerr << "]\n";
                    }
                    for (int s: syn_to_delete) is_disabled[s] = -1;
                    return true;
                }
            };

            // ----------------------------------------
            // Delete any corner where either member is
            // already disabled (-> delete both members)
            // ----------------------------------------
            for (int i=0;i<4;i++){
                int dq = corners[2*i];
                int sq = corners[2*i+1];
                if (is_disabled[dq] || is_disabled[sq]) {
                    is_disabled[dq] = -1;
                    is_disabled[sq] = -1;
                    if (verbose) std::cerr << "Delete corner data="<<dq<<" syn="<<sq<<"\n";
                }
            }

            // -----------------------------------------------------
            // Count inactive, set flags, and decide if we must stop
            // -----------------------------------------------------
            auto check_remaining_qubits = [&]() -> bool {
                // count in data section [0..DD), syn section [DD..2*DD-1)
                num_inactive_data = 0;
                num_inactive_syn  = 0;
                num_data_superstabilizer = 0;
                for (int i=0;i<DD;i++){
                    if (is_disabled[i] != 0) {
                        num_inactive_data++;
                        if (is_disabled[i] == 1) num_data_superstabilizer++;
                    }
                }
                for (int i=DD;i<2*DD-1;i++){
                    if (is_disabled[i] != 0) num_inactive_syn++;
                }
                too_many_qubits_lost = (num_inactive_data > (int)(DD * data_loss_tolerance));
                return too_many_qubits_lost;
            };

            // -----------------------------------------------------
            // Disconnect cleanup (boundary-only mode)
            // -----------------------------------------------------
            auto handle_disconnected_syn_boundary = [&]() -> bool {
                bool any=false;
                for (size_t si=0; si<syn_info.size(); ++si){
                    int sname = (int)si + DD;
                    // active or (marked disabled==1 but treatable in boundary pass)
                    if (is_disabled[sname]==0 || is_disabled[sname]==1){
                        int active_data=0;
                        bool near_deleted=false;
                        for (int k=0;k<4;k++){
                            int dn = syn_info[si].data_qubits[k];
                            if (dn==-1) continue;
                            if (is_disabled[dn]==0 || is_disabled[dn]==1) active_data++;
                            else if (is_disabled[dn]==-1) near_deleted=true;
                        }
                        if (active_data<=1){
                            if (near_deleted){
                                // boundary pass: delete this syn outright
                                is_disabled[sname] = -1;
                                any=true;
                                if (verbose) std::cerr<<"Boundary cleanup: delete syn "<<sname<<"\n";
                            } else if (is_disabled[sname]==0){
                                // interior-ish but touched by earlier boundary marks: mark disabled (tentative)
                                is_disabled[sname] = 1;
                                any=true;
                                if (verbose) std::cerr<<"Boundary cleanup: disable syn "<<sname<<"\n";
                            }
                        }
                    }
                }
                return any;
            };

            auto handle_disconnected_data_boundary = [&]() -> bool {
                bool any=false;
                for (int dn=0; dn<DD; ++dn){
                    if (is_disabled[dn]==0 || is_disabled[dn]==1){
                        int active_x=0, active_z=0;
                        bool near_deleted=false;
                        for (int sname : data_info[dn]){
                            if (is_disabled[sname]==0 || is_disabled[sname]==1){
                                char b = syn_info[sname - DD].basis;
                                if (b=='X') active_x++; else active_z++;
                            } else if (is_disabled[sname]==-1){
                                near_deleted=true;
                            }
                        }
                        if (active_x==0 || active_z==0){
                            if (near_deleted){
                                is_disabled[dn] = -1;
                                any=true;
                                if (verbose) std::cerr<<"Boundary cleanup: delete data "<<dn<<"\n";
                            } else if (is_disabled[dn]==0){
                                is_disabled[dn] = 1;
                                any=true;
                                if (verbose) std::cerr<<"Boundary cleanup: disable data "<<dn<<"\n";
                            }
                        }
                    }
                }
                return any;
            };

            auto remove_disconnected_component_boundary = [&]() -> bool {
                // Build list of currently "not missing" nodes for boundary mode:
                auto not_missing = [&](int qname)->bool {
                    return is_disabled[qname] != -1;
                };

                std::vector<int> active_data_names;
                active_data_names.reserve(DD);
                for (int i=0;i<DD;i++) if (not_missing(i)) active_data_names.push_back(i);
                int n = (int)active_data_names.size();
                if (n==0) return false;

                DSU dsu(n);
                auto index_of = [&](int data_name)->int {
                    // small helper: linear scan is fine here
                    for (int i=0;i<n;i++) if (active_data_names[i]==data_name) return i;
                    return -1;
                };

                // connect data nodes that share any not-missing syndrome
                for (int i=0;i<n;i++){
                    int di = active_data_names[i];
                    std::unordered_set<int> syn_i;
                    for (int s : data_info[di]) if (not_missing(s)) syn_i.insert(s);
                    for (int j=i+1;j<n;j++){
                        int dj = active_data_names[j];
                        bool share=false;
                        for (int s : data_info[dj]){
                            if (not_missing(s) && syn_i.count(s)){ share=true; break; }
                        }
                        if (share) dsu.unite(i,j);
                    }
                }

                // get components
                auto comps = dsu.components();
                if (comps.size() <= 1) return false;

                // keep largest; delete others
                std::sort(comps.begin(), comps.end(),
                        [](const std::vector<int>&a,const std::vector<int>&b){return a.size()>b.size();});
                if (verbose) std::cerr<<"Boundary cleanup: multiple comps, remove all but largest\n";
                for (size_t k=1;k<comps.size();++k){
                    for (int idx : comps[k]){
                        int dn = active_data_names[idx];
                        is_disabled[dn] = -1;
                        if (verbose) std::cerr<<"\tRemove data "<<dn<<"\n";
                    }
                }
                return true;
            };

            auto handle_disconnected_qubits_boundary = [&]() -> bool {
                bool any=false;
                bool changed=true;
                while (changed){
                    changed=false;
                    if (handle_disconnected_syn_boundary()){ changed=true; any=true; }
                    if (handle_disconnected_data_boundary()){ changed=true; any=true; }
                    if (remove_disconnected_component_boundary()){ changed=true; any=true; }
                }
                return any;
            };

            // -----------------------------------------------------
            // Boundary graph neighbors for rebuilding boundaries
            // -----------------------------------------------------
            auto neighbors_on_boundary = [&](int qname, char boundary_type) -> std::vector<int> {
                // Return data neighbors of qname along a boundary of given color ('X' for red, 'Z' for blue)
                auto [x0,y0] = data_list[qname].coords;

                // syndrome of that color adjacent to qname that are NOT deleted
                std::vector<int> active_syn;
                for (int sname : data_info[qname]){
                    if (is_disabled[sname] != -1){
                        if (syn_info[sname - DD].basis == boundary_type){
                            active_syn.push_back(sname);
                        }
                    }
                }
                std::vector<int> out;
                out.reserve(2);

                for (int sname : active_syn){
                    // split remaining vs deleted data around this syn
                    std::vector<int> remaining, deleted;
                    for (int k=0;k<4;k++){
                        int dn = syn_info[sname - DD].data_qubits[k];
                        if (dn==-1) continue;
                        if (is_disabled[dn]==-1) deleted.push_back(dn);
                        else                     remaining.push_back(dn);
                    }

                    if ((int)remaining.size() == 2){
                        // simple 2-neighbor case
                        if (remaining[0]==qname) out.push_back(remaining[1]);
                        else                     out.push_back(remaining[0]);

                    } else if ((int)remaining.size() == 3){
                        // pick the diagonal neighbor (differs by 2 in both coords)
                        for (int dn : remaining){
                            auto [x1,y1] = data_list[dn].coords;
                            if (std::abs(x1-x0)==2 && std::abs(y1-y0)==2){
                                out.push_back(dn);
                                break;
                            }
                        }

                    } else {
                        // 4 neighbors case; decide using neighbor syn presence like Python
                        auto [sx,sy] = syn_info[sname - DD].coords;
                        std::pair<int,int> neighbor_syn0 = {sx, sy + 2*(y0 - sy)};
                        std::pair<int,int> neighbor_data0 = {x0 + 2*(sx - x0), y0};

                        std::pair<int,int> neighbor_syn1 = {sx + 2*(x0 - sx), sy};
                        std::pair<int,int> neighbor_data1 = {x0, y0 + 2*(sy - y0)};

                        auto key0 = ((int64_t)neighbor_syn0.first<<32) ^ ( (int64_t)neighbor_syn0.second & 0xffffffffLL );
                        auto key1 = ((int64_t)neighbor_syn1.first<<32) ^ ( (int64_t)neighbor_syn1.second & 0xffffffffLL );

                        bool syn0_missing = (syndrome_matching.find(key0)==syndrome_matching.end()) ||
                                            (is_disabled[syndrome_matching[key0]]==-1);
                        bool syn1_missing = (syndrome_matching.find(key1)==syndrome_matching.end()) ||
                                            (is_disabled[syndrome_matching[key1]]==-1);

                        if (syn0_missing) {
                            out.push_back( data_matching[neighbor_data0.first][neighbor_data0.second] );
                        }
                        if (syn1_missing) {
                            out.push_back( data_matching[neighbor_data1.first][neighbor_data1.second] );
                        }
                    }
                }
                // dedup
                std::sort(out.begin(), out.end());
                out.erase(std::unique(out.begin(), out.end()), out.end());
                return out;
            };

            // recursively extend [ .... ] as in Python
            std::function<std::vector<int>(const std::vector<int>&, char)> extend_boundary_list =
            [&](const std::vector<int> &boundary_list, char boundary_type) -> std::vector<int> {
                if (!boundary_list.empty() && boundary_list.front()==-1 && boundary_list.back()==-1){
                    // both ends terminated -> strip sentinels
                    return std::vector<int>(boundary_list.begin()+1, boundary_list.end()-1);
                }
                std::vector<int> ext = boundary_list;

                auto extend_front = [&](){
                    if (ext.front()==-1) return;
                    auto nbrs = neighbors_on_boundary(ext.front(), boundary_type);
                    if (nbrs.size()==1){
                        if (ext.size()==1){
                            ext.insert(ext.begin(), -1);
                            ext.push_back(nbrs[0]);
                        } else {
                            ext.insert(ext.begin(), -1);
                        }
                    } else if (nbrs.size()==2){
                        if (ext.size()==1){
                            ext.insert(ext.begin(), nbrs[0]);
                            ext.push_back(nbrs[1]);
                        } else {
                            // choose the one that's not the current second
                            if (nbrs[0]==ext[1]) ext.insert(ext.begin(), nbrs[1]);
                            else                 ext.insert(ext.begin(), nbrs[0]);
                        }
                    }
                };

                auto extend_back = [&](){
                    if (ext.size()>1 && ext.back()==-1) return;
                    auto nbrs = neighbors_on_boundary(ext.back(), boundary_type);
                    if (nbrs.size()==1){
                        // must equal previous
                        ext.push_back(-1);
                    } else if (nbrs.size()==2){
                        if (nbrs[0]==ext[ext.size()-2]) ext.push_back(nbrs[1]);
                        else                             ext.push_back(nbrs[0]);
                    }
                };

                extend_front();
                if (ext.size()>1) extend_back();

                return extend_boundary_list(ext, boundary_type);
            };

            // recompute all four boundaries; return true if invalid/percolated
            auto update_dynamic_boundary = [&]() -> bool {
                for (int i=0;i<4;i++){
                    // unchanged nodes from current boundary
                    std::vector<int> unchanged_idx;
                    for (int j=0;j<(int)dynamic_boundaries[i].size(); ++j){
                        int dn = dynamic_boundaries[i][j];
                        if (is_disabled[dn] != -1) unchanged_idx.push_back(j);
                    }
                    if ((int)unchanged_idx.size() == (int)dynamic_boundaries[i].size()){
                        continue; // no change for this edge
                    } else if (unchanged_idx.empty()){
                        throw std::runtime_error("Entire boundary shifted away: edge "+std::to_string(i));
                    }

                    // new boundary from a surviving seed
                    int seed = dynamic_boundaries[i][unchanged_idx[0]];
                    char btype = (i<=1) ? 'X' : 'Z'; // red edges 0/1 use 'X', blue 2/3 use 'Z'
                    auto nb = extend_boundary_list({seed}, btype);

                    // orient from smaller to larger axis coordinate like Python
                    auto [xs,ys] = data_list[nb.front()].coords;
                    auto [xe,ye] = data_list[nb.back()].coords;
                    if ( (i<=1 && xs>xe) || (i>1 && ys>ye) ) std::reverse(nb.begin(), nb.end());

                    dynamic_boundaries[i] = std::move(nb);
                }

                if (verbose){
                    std::cerr<<"New boundaries:\n";
                    for (int i=0;i<4;i++){
                        std::cerr<<"  edge "<<i<<": size="<<dynamic_boundaries[i].size()<<"\n";
                    }
                }

                // invalid if red boundaries touch
                {
                    std::unordered_set<int> R(dynamic_boundaries[0].begin(), dynamic_boundaries[0].end());
                    for (int dn : dynamic_boundaries[1]){
                        if (R.count(dn)) {
                            if (verbose) std::cerr<<"Red edges collapse.\n";
                            return true;
                        }
                    }
                }
                // invalid if blue boundaries touch
                {
                    std::unordered_set<int> B(dynamic_boundaries[2].begin(), dynamic_boundaries[2].end());
                    for (int dn : dynamic_boundaries[3]){
                        if (B.count(dn)) {
                            if (verbose) std::cerr<<"Blue edges collapse.\n";
                            return true;
                        }
                    }
                }
                // corner connectivity checks (top-left, bottom-left, top-right, bottom-right)
                if (dynamic_boundaries[0].empty() || dynamic_boundaries[1].empty() ||
                    dynamic_boundaries[2].empty() || dynamic_boundaries[3].empty()){
                    return true; // degenerate
                }
                if ( dynamic_boundaries[0].front() != dynamic_boundaries[2].front() ||
                    dynamic_boundaries[0].back()  != dynamic_boundaries[3].front() ||
                    dynamic_boundaries[2].back()  != dynamic_boundaries[1].front() ||
                    dynamic_boundaries[1].back()  != dynamic_boundaries[3].back() ) {
                    if (verbose) std::cerr<<"Boundaries don’t connect properly.\n";
                    return true;
                }
                return false;
            };


            // -----------------------------------------------------------------------
            // Pass 1 over currently disabled (==1) qubits to handle boundary defects.
            // (Python: loop over is_disabled, skip corners, call handle_boundary_defects)
            // NOTE: In Python, after handling each, they call boundary/disconnect updates.
            // We’ll do the deletion here; you’ll plug in the update calls in the next step.
            // -----------------------------------------------------------------------
            for (int i = 0; i < (int)is_disabled.size(); ++i) {
                if (is_disabled[i] == 1) {
                    // skip the 4 corners (Python asserts "i not in corners")
                    bool is_corner = false;
                    for (int k = 0; k < 4; k++) {
                        if (i == corners[2 * k] || i == corners[2 * k + 1]) {
                            is_corner = true;
                            break;
                        }
                    }
                    if (is_corner) continue;

                    bool handled = handle_boundary_defects(i);
                    if (handled) {
                        // boundary-only disconnect cleanup
                        handle_disconnected_qubits_boundary();//(/*only_handle_boundary=*/true);

                        // too many qubits lost?
                        if (check_remaining_qubits()) {
                            return;  // flags already set (too_many_qubits_lost)
                        }

                        // recompute dynamic boundaries; if invalid -> percolation
                        if (update_dynamic_boundary()) {
                            percolated = true;
                            return;
                        }
                    }
                }
            }


            // ========== Helpers used by the new-boundary & interior passes ==========

            // Count remaining (non-deleted) data neighbors of a syndrome.
            auto count_remaining_dataq = [&](int sname)->int {
                assert(sname >= DD);
                int c=0;
                for (int k=0;k<4;k++){
                    int dn = syn_info[sname - DD].data_qubits[k];
                    if (dn!=-1 && is_disabled[dn] != -1) c++;
                }
                return c;
            };

            // Handle a defect that is NOW on/near the *dynamic* boundary.
            // Returns true if it deletes something (i.e., made progress).
            auto handle_defect_on_new_boundary = [&](int qname)->bool {
                auto is_corner_dq = [&](int dn)->bool {
                    if (dynamic_boundaries[0].empty() || dynamic_boundaries[1].empty() ||
                        dynamic_boundaries[2].empty() || dynamic_boundaries[3].empty()) return false;
                    return dn==dynamic_boundaries[0].front() || dn==dynamic_boundaries[0].back()
                        || dn==dynamic_boundaries[1].front() || dn==dynamic_boundaries[1].back();
                };

                // data qubit case
                if (qname < DD) {
                    // Corner special case on dynamic boundary
                    if (is_corner_dq(qname)) {
                        std::vector<int> active_syn;
                        for (int s : data_info[qname]) if (is_disabled[s]!=-1) active_syn.push_back(s);
                        if ((int)active_syn.size()!=2) return false; // nothing to do
                        // If either is tentatively disabled, postpone (Python does this)
                        if (is_disabled[active_syn[0]]==1 || is_disabled[active_syn[1]]==1) return false;

                        if (verbose) std::cerr<<"Handle (new) corner defect at data "<<qname<<"\n";
                        // delete the lower-weight neighbor syndrome: choose the one leaving fewer remaining data
                        if (count_remaining_dataq(active_syn[0]) < count_remaining_dataq(active_syn[1])) {
                            is_disabled[active_syn[0]] = -1;
                            if (verbose) std::cerr<<"\tAlso delete syn "<<active_syn[0]<<"\n";
                        } else {
                            is_disabled[active_syn[1]] = -1;
                            if (verbose) std::cerr<<"\tAlso delete syn "<<active_syn[1]<<"\n";
                        }
                        is_disabled[qname] = -1;
                        return true;
                    }

                    auto in_vec = [&](const std::vector<int> &v, int x)->bool{
                        return std::find(v.begin(), v.end(), x)!=v.end();
                    };
                    bool on_red  = in_vec(dynamic_boundaries[0], qname) || in_vec(dynamic_boundaries[1], qname);
                    bool on_blue = in_vec(dynamic_boundaries[2], qname) || in_vec(dynamic_boundaries[3], qname);

                    if (on_red || on_blue){
                        // delete defective data and the adjacent opposite-color syndrome
                        is_disabled[qname] = -1;
                        if (verbose) std::cerr<<"Delete defective data "<<qname<<" on new "<<(on_red?"red":"blue")<<" boundary\n";
                        for (int s : data_info[qname]){
                            if (is_disabled[s]==-1) continue;
                            char b = syn_info[s - DD].basis;
                            if ((on_red  && b=='Z') || (on_blue && b=='X')) {
                                is_disabled[s] = -1;
                                if (verbose) std::cerr<<"\tAlso delete syn "<<s<<"\n";
                                return true;
                            }
                        }
                        return true;
                    }

                    // interior data near boundary: if any adjacent syn's superstabilizer would require a deleted qubit
                    char delete_basis = 'N';
                    for (int s : data_info[qname]){
                        if (count_remaining_dataq(s) < 4){
                            delete_basis = (syn_info[s - DD].basis=='X') ? 'Z' : 'X';
                            break;
                        }
                    }
                    if (delete_basis=='N') return false;

                    is_disabled[qname] = -1;
                    if (verbose) std::cerr<<"Delete defective data "<<qname<<" (superstabilizer requires deleted neighbors)\n";
                    for (int s : data_info[qname]){
                        if (syn_info[s - DD].basis == delete_basis){
                            is_disabled[s] = -1;
                            if (verbose) std::cerr<<"\tDelete syn "<<s<<"\n";
                        }
                    }
                    return true;
                }

                // syndrome case (on dynamic boundary iff any adjacent data is on a dynamic boundary)
                int sname = qname;
                char boundary_color = 'N';
                for (int k=0;k<4;k++){
                    int dn = syn_info[sname - DD].data_qubits[k];
                    if (dn==-1) continue;
                    if (dynamic_boundaries[0].empty()) continue;
                    if (dn==dynamic_boundaries[0].front() || dn==dynamic_boundaries[0].back()
                        || dn==dynamic_boundaries[1].front() || dn==dynamic_boundaries[1].back()){
                        // corner: delete both
                        is_disabled[sname] = -1;
                        is_disabled[dn] = -1;
                        if (verbose) std::cerr<<"Delete (new) corner syn "<<sname<<" and data "<<dn<<"\n";
                        return true;
                    }
                    auto in_vec = [&](const std::vector<int> &v, int x)->bool{
                        return std::find(v.begin(), v.end(), x)!=v.end();
                    };
                    if (in_vec(dynamic_boundaries[0], dn) || in_vec(dynamic_boundaries[1], dn)) { boundary_color = 'X'; break; }
                    if (in_vec(dynamic_boundaries[2], dn) || in_vec(dynamic_boundaries[3], dn)) { boundary_color = 'Z'; break; }
                }
                if (boundary_color=='N') return false;

                is_disabled[sname] = -1;
                if (syn_info[sname - DD].basis != boundary_color) {
                    if (verbose) std::cerr<<"Delete (new) boundary syn "<<sname<<"\n";
                } else {
                    // same color as boundary: delete neighboring opposite-color syns at (±2,0) and (0,±2)
                    auto [x,y] = syn_info[sname - DD].coords;
                    std::vector<std::pair<int,int>> neigh = {{x,y+2},{x,y-2},{x+2,y},{x-2,y}};
                    std::vector<int> to_delete;
                    for (auto [cx,cy] : neigh){
                        auto key = ((int64_t)cx<<32) ^ ( (int64_t)cy & 0xffffffffLL );
                        auto it = syndrome_matching.find(key);
                        if (it!=syndrome_matching.end()){
                            int sn = it->second;
                            if (is_disabled[sn] != -1){
                                is_disabled[sn] = -1;
                                to_delete.push_back(sn);
                            }
                        }
                    }
                    if (verbose){
                        std::cerr<<"Handle defective (new) boundary syn "<<sname<<" also delete [";
                        for (size_t i=0;i<to_delete.size();++i){ if (i) std::cerr<<","; std::cerr<<to_delete[i]; }
                        std::cerr<<"]\n";
                    }
                }
                return true;
            };

            // ========== New-boundary pass: keep handling until stable ========== PASS 2
            {
                bool change = true;
                while (change){
                    change = false;
                    for (int i=0;i<(int)is_disabled.size();++i){
                        if (is_disabled[i]==1){ // defective but not yet deleted
                            if (verbose) std::cerr<<"Inspect potential new-boundary defect "<<i<<"\n";
                            if (handle_defect_on_new_boundary(i)){
                                change = true;

                                // boundary-only cleanup, count, update; bail on failure/percolation
                                handle_disconnected_qubits_boundary();
                                if (check_remaining_qubits()) return;
                                if (update_dynamic_boundary()) { percolated = true; return; }
                                break; // restart scanning after boundary update
                            }
                        }
                    }
                }
            }

            // ========== Interior phase with tentative disables ==========
            std::vector<int> interior_defects;          // current interior disabled(==1) qubits
            std::vector<int> tentative_disable;         // items set to 1 tentatively (we may reset back to 0)

            // General-pass helpers (not boundary-only):
            auto handle_broken_syn = [&]()->bool{
                bool any=false;
                for (size_t si=0; si<syn_info.size(); ++si){
                    int sname = (int)si + DD;
                    if (is_disabled[sname]==1){ // syndrome disabled -> disable its data neighbors
                        for (int k=0;k<4;k++){
                            int dn = syn_info[si].data_qubits[k];
                            if (dn==-1) continue;
                            if (is_disabled[dn]==0){
                                is_disabled[dn]=1;
                                tentative_disable.push_back(dn);
                                any=true;
                                if (verbose) std::cerr<<"Interior: disable data "<<dn<<" due to broken syn "<<sname<<"\n";
                            }
                        }
                    }
                }
                return any;
            };

            auto handle_disconnected_syn_general = [&]()->bool{
                bool any=false;
                for (size_t si=0; si<syn_info.size(); ++si){
                    int sname = (int)si + DD;
                    if (is_disabled[sname]==0){
                        int active_data=0;
                        for (int k=0;k<4;k++){
                            int dn = syn_info[si].data_qubits[k];
                            if (dn!=-1 && is_disabled[dn]==0) active_data++;
                        }
                        if (active_data<=1){
                            is_disabled[sname] = 1;
                            tentative_disable.push_back(sname);
                            any=true;
                            if (verbose) std::cerr<<"Interior: disable syn "<<sname<<" (<=1 active data)\n";
                        }
                    }
                }
                return any;
            };

            auto handle_disconnected_data_general = [&]()->bool{
                bool any=false;
                for (int dn=0; dn<DD; ++dn){
                    if (is_disabled[dn]==0){
                        int active_x=0, active_z=0;
                        for (int s : data_info[dn]){
                            if (is_disabled[s]==0){
                                char b = syn_info[s - DD].basis;
                                if (b=='X') active_x++; else active_z++;
                            }
                        }
                        if (active_x==0 || active_z==0){
                            is_disabled[dn] = 1;
                            tentative_disable.push_back(dn);
                            any=true;
                            if (verbose) std::cerr<<"Interior: disable data "<<dn<<" (missing X/Z)\n";
                        }
                    }
                }
                return any;
            };

            auto remove_disconnected_component_general = [&]()->bool{
                // Components over *active* data (is_disabled==0). Keep largest, disable others.
                std::vector<int> active;
                for (int i=0;i<DD;i++) if (is_disabled[i]==0) active.push_back(i);
                int n=(int)active.size();
                if (n<=1) return false;
                DSU dsu(n);
                for (int i=0;i<n;i++){
                    int di=active[i];
                    std::unordered_set<int> syn_i;
                    for (int s : data_info[di]) if (is_disabled[s]==0) syn_i.insert(s);
                    for (int j=i+1;j<n;j++){
                        int dj=active[j];
                        bool share=false;
                        for (int s : data_info[dj]){
                            if (is_disabled[s]==0 && syn_i.count(s)){ share=true; break; }
                        }
                        if (share) dsu.unite(i,j);
                    }
                }
                auto comps = dsu.components();
                if (comps.size()<=1) return false;
                std::sort(comps.begin(), comps.end(),
                        [](const std::vector<int>&a,const std::vector<int>&b){return a.size()>b.size();});
                if (verbose) std::cerr<<"Interior: multiple components; disable all but largest\n";
                for (size_t k=1;k<comps.size();++k){
                    for (int idx : comps[k]){
                        int dn = active[idx];
                        if (is_disabled[dn]==0){
                            is_disabled[dn]=1;
                            tentative_disable.push_back(dn);
                            if (verbose) std::cerr<<"\tDisable data "<<dn<<"\n";
                        }
                    }
                }
                return true;
            };

            auto handle_disconnected_qubits_general = [&]()->bool{
                bool any=false, changed=true;
                while (changed){
                    changed=false;
                    if (handle_disconnected_syn_general()){ changed=true; any=true; }
                    if (handle_disconnected_data_general()){ changed=true; any=true; }
                    if (remove_disconnected_component_general()){ changed=true; any=true; }
                }
                return any;
            };

            // Prepare interior_defects and a helper to reset tentative disables (like Python)
            auto reset_tentative_disable = [&](){
                // keep disabled-but-not-deleted as interior list
                interior_defects.clear();
                for (int i=0;i<(int)is_disabled.size();++i){
                    if (is_disabled[i]==1) interior_defects.push_back(i);
                }
                // revert tentative disables that weren't deleted
                for (int q : tentative_disable){
                    if (is_disabled[q] != -1) is_disabled[q] = 0;
                }
                tentative_disable.clear();

                // refill tentatives by propagating from broken syns and disconnections
                bool changed=true;
                while (changed){
                    changed=false;
                    if (handle_broken_syn()) changed=true;
                    if (handle_disconnected_qubits_general()) changed=true;
                }
            };

            // Initialize interior state
            reset_tentative_disable();

            // Keep deforming boundary driven by interior defects
            // pass 3
            {
                bool change = true;
                while (change){
                    change=false;
                    // Scan interior defects + tentatives
                    std::vector<int> scan = interior_defects;
                    scan.insert(scan.end(), tentative_disable.begin(), tentative_disable.end());
                    for (int q : scan){
                        if (is_disabled[q]==1){
                            if (handle_defect_on_new_boundary(q)){
                                change = true;
                                handle_disconnected_qubits_boundary();
                                if (check_remaining_qubits()) return;
                                if (update_dynamic_boundary()) { percolated = true; return; }
                                break; // restart after boundary update
                            }
                        }
                    }
                    if (change) reset_tentative_disable();
                }
            }

            if (check_remaining_qubits()) return;

            // ========== Edge deformation metrics (left,right,top,bottom) ==========
            auto accumulate_deformation = [&](int edge_idx)->double {
                // edge_idx: 0=L,1=R use y regular; 2=Top,3=Bottom use x regular
                double acc=0.0;
                int last_regular = -2;
                bool across=false;
                for (int dn : dynamic_boundaries[edge_idx]){
                    auto [x,y]=data_list[dn].coords;
                    bool regular = (edge_idx==0 ? (y==0)
                                    : edge_idx==1 ? (y==2*(d-1))
                                    : edge_idx==2 ? (x==0)
                                                : (x==2*(d-1)));
                    int axis = (edge_idx<=1)? x : y;
                    if (regular){
                        if (across){
                            assert(axis - last_regular > 0);
                            acc += (axis - last_regular - 2)/2.0;
                            across=false;
                            if (verbose){
                                static const char* names[4]={"Left","Right","Top","Bottom"};
                                std::cerr<<names[edge_idx]<<" edge, deformed between "<<last_regular<<" and "<<axis<<"\n";
                            }
                        }
                        last_regular = axis;
                    } else {
                        if (last_regular!=-2) across=true;
                    }
                }
                return acc;
            };
            boundary_deformation[0]=accumulate_deformation(0);
            boundary_deformation[1]=accumulate_deformation(1);
            boundary_deformation[2]=accumulate_deformation(2);
            boundary_deformation[3]=accumulate_deformation(3);

            if (check_remaining_qubits()) return;

            // ========== Finalize sets: active data, ancilla vs gauges, clusters ==========

            // Active data list
            data.clear();
            for (int i=0;i<DD;i++){
                if (is_disabled[i]==0) data.push_back(data_list[i]);
            }

            // Disabled data (not deleted) for clustering
            std::vector<DataQubit> disabled_data_q;
            for (int i=0;i<DD;i++) if (is_disabled[i]==1) disabled_data_q.push_back(data_list[i]);

            std::vector<std::vector<int>> connected_defects;
            std::vector<std::vector<int>> cluster_x_idx; // indices into syn_info (we'll copy to MeasureQubit later)
            std::vector<std::vector<int>> cluster_z_idx;

            if (!disabled_data_q.empty()){
                DSU dsu((int)disabled_data_q.size());
                auto coord_of = [&](int idx)->std::pair<int,int>{ return disabled_data_q[idx].coords; };
                for (int i=0;i<(int)disabled_data_q.size();++i){
                    auto ci = coord_of(i);
                    for (int j=i+1;j<(int)disabled_data_q.size();++j){
                        auto cj = coord_of(j);
                        if (std::abs(ci.first - cj.first) <= 2 && std::abs(ci.second - cj.second) <= 2){
                            dsu.unite(i,j);
                        }
                    }
                }
                auto comps = dsu.components();
                connected_defects = comps;
                cluster_x_idx.assign(comps.size(), {});
                cluster_z_idx.assign(comps.size(), {});
            }

            // Split syns into stabilizers (ancilla) vs gauges (broken-neighbor) and bucket gauges by cluster
            x_ancilla.clear(); z_ancilla.clear();
            x_gauges.clear();  z_gauges.clear();

            auto name_in_disabled = [&](int dn)->int {
                // return cluster index containing dn (by name), or -1 if not found
                for (size_t ci=0; ci<connected_defects.size(); ++ci){
                    for (int local_idx : connected_defects[ci]){
                        if (disabled_data_q[local_idx].name == dn) return (int)ci;
                    }
                }
                return -1;
            };

            for (size_t si=0; si<syn_info.size(); ++si){
                int sname = (int)si + DD;
                if (is_disabled[sname]!=0) continue; // must be active to be measured now

                // Copy and mask data neighbors: -1 for deleted or disabled data
                MeasureQubit mq = syn_info[si];
                bool is_gauge=false;
                int any_disabled_dn = -1;

                for (int k=0;k<4;k++){
                    int dn = mq.data_qubits[k];
                    if (dn==-1) continue;
                    if (is_disabled[dn]==-1) { mq.data_qubits[k] = -1; continue; }
                    if (is_disabled[dn]==1)  { mq.data_qubits[k] = -1; is_gauge=true; any_disabled_dn=dn; }
                    // else keep as-is (active data)
                }

                if (is_gauge){
                    int ci = name_in_disabled(any_disabled_dn);
                    if (mq.basis=='X'){
                        x_gauges.push_back(mq);
                        if (ci!=-1) cluster_x_idx[ci].push_back((int)si);
                    } else {
                        z_gauges.push_back(mq);
                        if (ci!=-1) cluster_z_idx[ci].push_back((int)si);
                    }
                } else {
                    if (mq.basis=='X') x_ancilla.push_back(mq);
                    else               z_ancilla.push_back(mq);
                }
            }

            // Corrected defect construction block:
            // Build defect objects from clusters
            defect.clear();
            for (size_t ci=0; ci<connected_defects.size(); ++ci){
                std::vector<int> names; names.reserve(connected_defects[ci].size());
                std::vector<std::pair<int,int>> coords; coords.reserve(connected_defects[ci].size());
                for (int local_idx : connected_defects[ci]){
                    names.push_back(disabled_data_q[local_idx].name);
                    coords.push_back(disabled_data_q[local_idx].coords);
                }
                
                // Create maps to link the original syn_info index (si) to the actual filtered Gauge object.
                std::unordered_map<int, MeasureQubit> filtered_gauges;
                for (const auto& mq : x_gauges) filtered_gauges[mq.name - DD] = mq;
                for (const auto& mq : z_gauges) filtered_gauges[mq.name - DD] = mq;

                std::vector<MeasureQubit> xgs, zgs;
                
                // FIX: Retrieve the filtered MeasureQubit from the map.
                // The index si is relative to syn_info, so we need to adjust by DD for the name.
                for (int si : cluster_x_idx[ci]) xgs.push_back(filtered_gauges.at(si));
                for (int si : cluster_z_idx[ci]) zgs.push_back(filtered_gauges.at(si));
                
                defect.emplace_back(names, coords, xgs, zgs);
            }

            // All active qubits (data + syndrome)
            all_qubits.clear();
            for (int i=0;i<(int)is_disabled.size();++i) if (is_disabled[i]==0) all_qubits.push_back(i);

            if (verbose){
                std::vector<int> deleted;
                for (int i=0;i<(int)is_disabled.size();++i) if (is_disabled[i]==-1) deleted.push_back(i);
                std::cerr<<"Deleted due to boundary/defects: [";
                for (size_t i=0;i<deleted.size();++i){ if (i) std::cerr<<","; std::cerr<<deleted[i]; }
                std::cerr<<"]\n";
            }

            // ========== Observable & distances (Graph BFS) ==========
            auto add_complete_pairs = [&](Graph &G, const std::vector<int> &v){
                // Add edges for all pairs in v
                for (size_t i=0;i<v.size();++i){
                    for (size_t j=i+1;j<v.size();++j){
                        G.add_edge(v[i], v[j]);
                    }
                }
            };

            // One-shortest-path BFS (returns path of node ids)
            auto bfs_one_path = [&](const Graph &G, int64_t s, int64_t t)->std::vector<int64_t>{
                std::unordered_map<int64_t,int> dist;
                std::unordered_map<int64_t,int64_t> prev;
                std::deque<int64_t> q;
                dist[s]=0; prev[s]=-9; q.push_back(s);
                while(!q.empty()){
                    auto u=q.front(); q.pop_front();
                    if (u==t) break;
                    auto it=G.adj.find(u);
                    if (it==G.adj.end()) continue;
                    for (auto v: it->second){
                        if(!dist.count(v)){
                            dist[v]=dist[u]+1; prev[v]=u; q.push_back(v);
                        }
                    }
                }
                if (!dist.count(t)) throw std::runtime_error("Failed to find path.");
                std::vector<int64_t> path;
                for (int64_t cur=t; cur!=-9; cur=prev[cur]) path.push_back(cur);
                std::reverse(path.begin(), path.end());
                return path;
            };

            // bool get_metrics = true; // you can wire this to constructor arg if you want parity
            update_dynamic_boundary();

            if (get_metrics){
                // Vertical (top->bottom) graph using X-type connectivity
                Graph Gx;
                // add nodes (data names)
                for (auto &dq : data) Gx.add_node(dq.name);
                Gx.add_node(-1); Gx.add_node(-2); // top source, bottom target

                // edges from X stabilizers
                for (auto &mx : x_ancilla){
                    std::vector<int> active;
                    for (int k=0;k<4;k++){ int dn=mx.data_qubits[k]; if (dn!=-1) active.push_back(dn); }
                    add_complete_pairs(Gx, active);
                }
                // edges from X gauges inside defects
                for (auto &cl : defect){
                    std::unordered_set<int> acc;
                    for (auto &xg : cl.x_gauges){
                        for (int k=0;k<4;k++){ int dn=xg.data_qubits[k]; if (dn!=-1) acc.insert(dn); }
                    }
                    std::vector<int> v(acc.begin(), acc.end());
                    add_complete_pairs(Gx, v);
                }
                for (int q : dynamic_boundaries[2]) Gx.add_edge(-1, q); // top
                for (int q : dynamic_boundaries[3]) Gx.add_edge(-2, q); // bottom

                auto [nsp_v, len_v] = Gx.bfs_shortest_paths(-1, -2);
                vertical_distance = (len_v>0)? (len_v - 1) : 0;

                // Horizontal (left->right) graph using Z-type connectivity
                Graph Gz;
                for (auto &dq : data) Gz.add_node(dq.name);
                Gz.add_node(-1); Gz.add_node(-2); // left source, right target

                for (auto &mz : z_ancilla){
                    std::vector<int> active;
                    for (int k=0;k<4;k++){ int dn=mz.data_qubits[k]; if (dn!=-1) active.push_back(dn); }
                    add_complete_pairs(Gz, active);
                }
                for (auto &cl : defect){
                    std::unordered_set<int> acc;
                    for (auto &zg : cl.z_gauges){
                        for (int k=0;k<4;k++){ int dn=zg.data_qubits[k]; if (dn!=-1) acc.insert(dn); }
                    }
                    std::vector<int> v(acc.begin(), acc.end());
                    add_complete_pairs(Gz, v);
                }
                for (int q : dynamic_boundaries[0]) Gz.add_edge(-1, q); // left
                for (int q : dynamic_boundaries[1]) Gz.add_edge(-2, q); // right

                auto [nsp_h, len_h] = Gz.bfs_shortest_paths(-1, -2);
                horizontal_distance = (len_h>0)? (len_h - 1) : 0;

            }

            {
                // Observable path: top->bottom using X-type edges but NOT crossing red stabilizers
                Graph Gobs;
                for (auto &dq : data) Gobs.add_node(dq.name);
                Gobs.add_node(-1);
                Gobs.add_node(-2);
                for (auto &mx : x_ancilla) {
                    std::vector<int> active;
                    for (int k = 0; k < 4; k++) {
                        int dn = mx.data_qubits[k];
                        if (dn != -1) active.push_back(dn);
                    }
                    add_complete_pairs(Gobs, active);
                }
                for (auto &cl : defect) {
                    for (auto &xg : cl.x_gauges) {
                        std::vector<int> active;
                        for (int k = 0; k < 4; k++) {
                            int dn = xg.data_qubits[k];
                            if (dn != -1) active.push_back(dn);
                        }
                        add_complete_pairs(Gobs, active);
                    }
                }
                for (int q : dynamic_boundaries[2]) Gobs.add_edge(-1, q);
                for (int q : dynamic_boundaries[3]) Gobs.add_edge(-2, q);

                try {
                    auto path = bfs_one_path(Gobs, -1, -2);
                    observable.clear();
                    for (size_t i = 1; i + 1 < path.size(); ++i) {
                        if (path[i] >= 0) observable.push_back((int)path[i]);
                    }
                    if (verbose) {
                        std::cerr << "Observable path length=" << observable.size() << "\n";
                    }
                } catch (...) {
                    throw std::runtime_error("Failed to find observable path");
                }
            }



            (void)missing_coords; (void)data_loss_tolerance; (void)verbose; (void)get_metrics;
        }

        // --- simple counters (match Python API names) ---
        int num_disabled_qubits() const { return num_inactive_data + num_inactive_syn; }
        int num_disabled_data()   const { return num_inactive_data; }
        int num_data_in_superstabilizer() const { return num_data_superstabilizer; }
        int num_disabled_syndromes() const { return num_inactive_syn; }

        int actual_distance_vertical() const { return vertical_distance; }
        int actual_distance_horizontal() const { return horizontal_distance; }
        std::array<double,4> edge_deformation() const { return boundary_deformation; }

        bool terminated_due_to_qubit_loss() const { return too_many_qubits_lost; }
        bool is_percolated() const { return percolated; }
        int  num_clusters() const { return (int)defect.size(); }

        void change_shell_diameter(int new_size){
            for (auto &c: defect) c.change_diameter(new_size);
        }
        void reset_err_rates(double r, double g1, double g2){
            readout_err=r; gate1_err=g1; gate2_err=g2;
        }

        // --- Stim integration helpers ---
        void apply_1gate(stim::Circuit &circ, const char *gate, const std::vector<int> &qs){
            if (!qs.empty()) {
                circ.safe_append_u(gate, qubit_targets(qs), {});
            }
            // DEPOLARIZE1 on all active qubits
            if (!all_qubits.empty()) {
                circ.safe_append_u("DEPOLARIZE1", qubit_targets(all_qubits), {gate1_err});
            }
            circ.safe_append_u("TICK", {}, {});
        }

        void apply_2gate(stim::Circuit &circ, const char *gate, const std::vector<int> &pairwise_qubits){
            // pairwise_qubits is [c,t,c,t,...]
            if (!pairwise_qubits.empty()){
                circ.safe_append_u(gate, qubit_targets(pairwise_qubits), {});
                circ.safe_append_u("DEPOLARIZE2", qubit_targets(pairwise_qubits), {gate2_err});
            }
            // idle single-qubit noise on others
            if (pairwise_qubits.size() < all_qubits.size()){
                // build idle set = all_qubits \ used
                std::unordered_set<int> used(pairwise_qubits.begin(), pairwise_qubits.end());
                std::vector<int> idle;
                idle.reserve(all_qubits.size());
                for (int q: all_qubits) if (!used.count(q)) idle.push_back(q);
                if (!idle.empty()){
                    circ.safe_append_u("DEPOLARIZE1", qubit_targets(idle), {gate1_err});
                }
            }
            circ.safe_append_u("TICK", {}, {});
        }

        // Reset/Measure bookkeeping to emulate Python’s rec indices
        void reset_meas_qubits(stim::Circuit &circ, const char *op, const std::vector<int> &qs, bool last=false){
            if (std::string(op) == "R") {
                if (!qs.empty()) circ.safe_append_u("R", qubit_targets(qs), {});
            }
            // readout error before measurement
            if (!qs.empty()) circ.safe_append_u("X_ERROR", qubit_targets(qs), {readout_err});

            if (std::string(op) == "M" || std::string(op) == "MR") {
                if (!qs.empty()) {
                    circ.safe_append_u("M", qubit_targets(qs), {});
                    // Update measurement record offsets:
                    // newest round: map q -> -1, -2, ...
                    std::unordered_map<int,int64_t> round;
                    for (int i=(int)qs.size()-1, k=1; i>=0; --i, ++k) {
                        round[qs[i]] = -k;  // -1, -2, ...
                    }
                    // shift older rounds further back by |qs|
                    const int64_t shift = (int64_t)qs.size();
                    if (shift != 0) {
                        for (auto &r : meas_record) {
                            for (auto &kv : r) {
                                kv.second -= shift;
                            }
                        }
                    }
                    meas_record.push_back(std::move(round));
                    total_meas_count += shift;
                }
            }
            if (!last && qs.size() < all_qubits.size()){
                std::unordered_set<int> used(qs.begin(), qs.end());
                std::vector<int> idle;
                idle.reserve(all_qubits.size());
                for (int q: all_qubits) if (!used.count(q)) idle.push_back(q);
                if (!idle.empty()) circ.safe_append_u("DEPOLARIZE1", qubit_targets(idle), {gate1_err});
            }
        }

        stim::GateTarget get_meas_rec(int round_idx, int qubit_name) const {
            // round_idx is negative in Python; interpret relative to end
            int idx = (int)meas_record.size() + round_idx;

            if (idx < 0 || idx >= (int)meas_record.size()) {
                // This is where the error is thrown.
                throw std::runtime_error("get_meas_rec: round_idx out of range.");
            }
            auto it = meas_record[idx].find(qubit_name);
            if (it == meas_record[idx].end()){
                throw std::runtime_error("get_meas_rec: qubit not found in that round.");
            }
            return stim::GateTarget::rec((int32_t)it->second);
        }

        // This mirrors your Python syndrome_round. Most of the complexity
        // is in building the correct target lists; the flow is the same.
        stim::Circuit& syndrome_round(stim::Circuit &circ, bool first=false, bool double_half=false){
            // Build syn_except_* like Python (names only):
            std::vector<int> syn_except_xgauge; // all syn except X gauges
            {
                for (const auto &g: z_gauges) syn_except_xgauge.push_back(g.name);
                for (const auto &m: x_ancilla) syn_except_xgauge.push_back(m.name);
                for (const auto &m: z_ancilla) syn_except_xgauge.push_back(m.name);
            }
            std::vector<int> syn_except_zgauge; // all syn except Z gauges
            {
                for (const auto &g: x_gauges) syn_except_zgauge.push_back(g.name);
                for (const auto &m: x_ancilla) syn_except_zgauge.push_back(m.name);
                for (const auto &m: z_ancilla) syn_except_zgauge.push_back(m.name);
            }

            if (first) reset_meas_qubits(circ, "R", all_qubits);
            else       reset_meas_qubits(circ, "R", syn_except_xgauge);

            circ.safe_append_u("TICK", {}, {});

            // H on X-ancilla and X-gauges
            {
                std::vector<int> hs; hs.reserve(x_ancilla.size()+x_gauges.size());
                for (auto &m: x_ancilla) hs.push_back(m.name);
                for (auto &g: x_gauges)  hs.push_back(g.name);
                apply_1gate(circ, "H", hs);
            }

            // 4 CX layers as in Python
            for (int i=0;i<4;i++){
                std::vector<int> pairwise;
                // X ancilla/gauges: control is measure, target is data[i]
                auto emit_pair = [&](const MeasureQubit &mq, bool measure_is_ctrl, int di){
                    int dq = mq.data_qubits[di];
                    if (dq == -1) return;
                    if (measure_is_ctrl) {
                        pairwise.push_back(mq.name);
                        pairwise.push_back(dq);
                    } else {
                        pairwise.push_back(dq);
                        pairwise.push_back(mq.name);
                    }
                };
                for (auto &mx: x_ancilla) emit_pair(mx,true,i);
                for (auto &mx: x_gauges)  emit_pair(mx,true,i);
                for (auto &mz: z_ancilla) emit_pair(mz,false,i);
                apply_2gate(circ, "CX", pairwise);
            }

            // H again on X set
            {
                std::vector<int> hs; hs.reserve(x_ancilla.size()+x_gauges.size());
                for (auto &m: x_ancilla) hs.push_back(m.name);
                for (auto &g: x_gauges)  hs.push_back(g.name);
                apply_1gate(circ, "H", hs);
            }

            reset_meas_qubits(circ, "M", syn_except_zgauge); // R1: X-M, R2+: X-M

            if (!first){
                circ.safe_append_u("SHIFT_COORDS", {}, {0.0,0.0,1.0});
                
                // DETECTORs for stabilizers (X and Z)
                for (auto &a: x_ancilla){
                    auto tgts = std::vector<uint32_t>{ get_meas_rec(-1, a.name).data, get_meas_rec(-2, a.name).data };
                    circ.safe_append_u("DETECTOR", tgts, {(double)a.coords.first, (double)a.coords.second, 0.0});
                }
                for (auto &a: z_ancilla){
                    auto tgts = std::vector<uint32_t>{ get_meas_rec(-1, a.name).data, get_meas_rec(-2, a.name).data };
                    circ.safe_append_u("DETECTOR", tgts, {(double)a.coords.first, (double)a.coords.second, 0.0});
                }
                
                // cluster x-gauges superstabilizer detector
                for (auto &cl: defect){
                    if (!cl.x_gauges.empty()){
                        std::vector<uint32_t> tgts;
                        // X-Gauge: Correct Python indices for R2 X-M (size=3)
                        for (auto &xg: cl.x_gauges) tgts.push_back(get_meas_rec(-1, xg.name).data);
                        for (auto &xg: cl.x_gauges) tgts.push_back(get_meas_rec(-3, xg.name).data);
                        
                        auto c0 = cl.coords.empty()?std::pair<int,int>{0,0}:cl.coords[0];
                        circ.safe_append_u("DETECTOR", tgts, {(double)c0.first,(double)c0.second,0.0});
                    }
                }
            } else {
                // This part is the R1 X-M detector (correctly uses REC[-1] vs reset)
                for (auto &a: z_ancilla){
                    auto tgts = std::vector<uint32_t>{ get_meas_rec(-1, a.name).data };
                    circ.safe_append_u("DETECTOR", tgts, {(double)a.coords.first, (double)a.coords.second, 0.0});
                }
            }
            
            circ.safe_append_u("TICK", {}, {});

            if (!double_half) return circ;

            // second half (mirrors Python)
            reset_meas_qubits(circ, "R", syn_except_zgauge);
            circ.safe_append_u("TICK", {}, {});
            {
                std::vector<int> hs;
                for (auto &m: x_ancilla) hs.push_back(m.name);
                apply_1gate(circ, "H", hs);
            }
            for (int i=0;i<4;i++){
                std::vector<int> pairwise;
                auto emit_pair = [&](const MeasureQubit &mq, bool measure_is_ctrl, int di){
                    int dq = mq.data_qubits[di];
                    if (dq == -1) return;
                    if (measure_is_ctrl) {
                        pairwise.push_back(mq.name);
                        pairwise.push_back(dq);
                    } else {
                        pairwise.push_back(dq);
                        pairwise.push_back(mq.name);
                    }
                };
                for (auto &mx: x_ancilla) emit_pair(mx,true,i);
                for (auto &mz: z_ancilla) emit_pair(mz,false,i);
                for (auto &mz: z_gauges)  emit_pair(mz,false,i);
                apply_2gate(circ, "CX", pairwise);
            }
            {
                std::vector<int> hs;
                for (auto &m: x_ancilla) hs.push_back(m.name);
                apply_1gate(circ, "H", hs);
            }
            reset_meas_qubits(circ, "M", syn_except_xgauge); // R1: Z-M, R2+: Z-M
            circ.safe_append_u("SHIFT_COORDS", {}, {0.0,0.0,1.0});
            
            for (auto &a: x_ancilla){
                auto tgts = std::vector<uint32_t>{ get_meas_rec(-1, a.name).data, get_meas_rec(-2, a.name).data };
                circ.safe_append_u("DETECTOR", tgts, {(double)a.coords.first,(double)a.coords.second,0.0});
            }
            for (auto &a: z_ancilla){
                auto tgts = std::vector<uint32_t>{ get_meas_rec(-1, a.name).data, get_meas_rec(-2, a.name).data };
                circ.safe_append_u("DETECTOR", tgts, {(double)a.coords.first,(double)a.coords.second,0.0});
            }
            
            // cluster z-gauges superstabilizer detector (depends on first/not-first)
            for (auto &cl: defect){
                if (!cl.z_gauges.empty()){
                    std::vector<uint32_t> tgts;
                    // Z-Gauge: Correct history check for REC[-3]
                    if (!first) {
                        // Non-first: Compares R-Z[t] (at -1) with R-Z[t-1] (at -3)
                        for (auto &zg: cl.z_gauges) tgts.push_back(get_meas_rec(-1, zg.name).data);
                        for (auto &zg: cl.z_gauges) tgts.push_back(get_meas_rec(-3, zg.name).data);
                    } else {
                        // First round: only one Z-gauge measurement has occurred
                        for (auto &zg: cl.z_gauges) tgts.push_back(get_meas_rec(-1, zg.name).data);
                    }
                    auto c0 = cl.coords.empty()?std::pair<int,int>{0,0}:cl.coords[0];
                    circ.safe_append_u("DETECTOR", tgts, {(double)c0.first,(double)c0.second,0.0});
                }
            }
            circ.safe_append_u("TICK", {}, {});
            return circ;
        }
        stim::Circuit generate_stim(int rounds){
            stim::Circuit circ;

            // QUBIT_COORDS for data and ancilla/gauges
            for (auto &d: data) {
                circ.safe_append_u("QUBIT_COORDS", {static_cast<uint32_t>(d.name)}, {(double)d.coords.first, (double)d.coords.second});
            }
            for (auto &m: x_ancilla) circ.safe_append_u("QUBIT_COORDS", {static_cast<uint32_t>(m.name)}, {(double)m.coords.first, (double)m.coords.second});
            for (auto &m: x_gauges)  circ.safe_append_u("QUBIT_COORDS", {static_cast<uint32_t>(m.name)}, {(double)m.coords.first, (double)m.coords.second});
            for (auto &m: z_ancilla) circ.safe_append_u("QUBIT_COORDS", {static_cast<uint32_t>(m.name)}, {(double)m.coords.first, (double)m.coords.second});
            for (auto &m: z_gauges)  circ.safe_append_u("QUBIT_COORDS", {static_cast<uint32_t>(m.name)}, {(double)m.coords.first, (double)m.coords.second});

            if (!x_gauges.empty()){
                syndrome_round(circ, /*first=*/true, /*double_half=*/true);
                stim::Circuit body;
                syndrome_round(body, /*first=*/false, /*double_half=*/true);
                circ.append_repeat_block(rounds - 1, std::move(body), ""sv);
            } else {
                syndrome_round(circ, /*first=*/true, /*double_half=*/false);
                stim::Circuit body;
                syndrome_round(body, /*first=*/false, /*double_half=*/false);
                circ.append_repeat_block(rounds - 1, std::move(body), ""sv);
            }

            // Final data measurements
            std::vector<int> all_data_names; all_data_names.reserve(data.size());
            for (auto &d: data) all_data_names.push_back(d.name);
            reset_meas_qubits(circ, "M", all_data_names, /*last=*/true);
            circ.safe_append_u("SHIFT_COORDS", {}, {0.0,0.0,1.0});

            // DETECTORs that combine final data recs with ancilla/gauge recs (mirror Python finale)
            for (auto &a: z_ancilla){
                std::vector<uint32_t> tgts;
                for (int i=0;i<4;i++){
                    int dq = a.data_qubits[i];
                    if (dq!=-1) tgts.push_back(get_meas_rec(-1, dq).data);
                }
                tgts.push_back(get_meas_rec(-2, a.name).data);
                circ.safe_append_u("DETECTOR", tgts, {(double)a.coords.first,(double)a.coords.second,0.0});
            }
            
            for (auto &cl: defect){
                if (!cl.z_gauges.empty()){
                    std::vector<int> z_gauge_data_names;
                    for (auto &zg: cl.z_gauges){
                        for (int i=0;i<4;i++){
                            if (zg.data_qubits[i]!=-1) {
                                z_gauge_data_names.push_back(zg.data_qubits[i]);
                            }
                        }
                    }
                    std::vector<uint32_t> tgts;
                    for (int q: z_gauge_data_names) {
                        tgts.push_back(get_meas_rec(-1, q).data);
                    }
                    for (auto &zg: cl.z_gauges)   tgts.push_back(get_meas_rec(-2, zg.name).data);
                    auto c0 = cl.coords.empty()?std::pair<int,int>{0,0}:cl.coords[0];
                    circ.safe_append_u("DETECTOR", tgts, {(double)c0.first,(double)c0.second,0.0});
                }
            }

            // OBSERVABLE_INCLUDE for the logical (vertical) observable
            

            {
                std::vector<stim::GateTarget> recs;
                recs.reserve(this->observable.size());
                for (int q : this->observable) recs.push_back(get_meas_rec(-1, q));

                if (!recs.empty()) {
                    const auto &gate = stim::GATE_DATA.at("OBSERVABLE_INCLUDE");
                    std::vector<double> args_vec{0.0};
                    stim::SpanRef<const double> args(args_vec);
                    stim::SpanRef<const stim::GateTarget> tgts(recs);
                    circ.safe_append(stim::CircuitInstruction(gate.id, args, tgts, ""sv));
                } else {
                    std::cerr << "⚠️ No targets for OBSERVABLE_INCLUDE — observable list is empty.\n";
                }
            }

            return circ;
        }
    };
}

namespace dqec {

stim::Circuit build_superstabilized_circuit(const ExperimentConfig& cfg) {
    uint32_t effective_rounds = cfg.rounds == 0 ? cfg.code_dist : cfg.rounds;
    superstabilizer::LogicalQubit logi(
        static_cast<int>(cfg.code_dist),
        cfg.measurement_error,
        cfg.physical_error,
        cfg.physical_error,
        cfg.superstabilizer_qubits,
        /*data_loss_tolerance=*/1.0,
        /*verbose=*/false);

    return logi.generate_stim(static_cast<int>(effective_rounds));
}

std::string polygon_pragmas_str(const ExperimentConfig& cfg) {
    superstabilizer::LogicalQubit logi(
        static_cast<int>(cfg.code_dist),
        cfg.measurement_error,
        cfg.physical_error,
        cfg.physical_error,
        cfg.superstabilizer_qubits,
        /*data_loss_tolerance=*/1.0,
        /*verbose=*/false);

    // Build data qubit coord lookup: name -> (x, y)
    std::unordered_map<int, std::pair<int,int>> data_coords;
    for (auto& d : logi.data) {
        data_coords[d.name] = d.coords;
    }

    // Generate one "#!pragma POLYGON(r,g,b,a) q0 q1 ..." line per stabilizer.
    // Data qubits are sorted by ascending atan2(dy, dx) from the ancilla centre,
    // which gives a consistent counter-clockwise vertex order.
    auto make_polygon = [&](const superstabilizer::MeasureQubit& a, bool is_merge) -> std::string {
        bool is_x = (a.basis == 'X');
        const char* color;
        if (is_merge) {
            color = is_x ? "0,0.7,0,0.3" : "1,0.5,0,0.3";
        } else {
            color = is_x ? "0,0,1,0.25" : "1,0,0,0.25";
        }

        double ax = static_cast<double>(a.coords.first);
        double ay = static_cast<double>(a.coords.second);

        std::vector<std::pair<double, int>> angle_qubit;
        for (int i = 0; i < 4; i++) {
            int q = a.data_qubits[i];
            if (q == -1) continue;
            auto it = data_coords.find(q);
            if (it == data_coords.end()) continue;
            double dx = it->second.first  - ax;
            double dy = it->second.second - ay;
            angle_qubit.push_back({std::atan2(dy, dx), q});
        }

        if (angle_qubit.empty()) return "";

        std::sort(angle_qubit.begin(), angle_qubit.end());

        std::ostringstream oss;
        oss << "#!pragma POLYGON(" << color << ")";
        for (auto& [angle, q] : angle_qubit) {
            oss << " " << q;
        }
        return oss.str();
    };

    std::ostringstream out;
    for (auto& a : logi.x_ancilla) { auto l = make_polygon(a, false); if (!l.empty()) out << l << "\n"; }
    for (auto& a : logi.z_ancilla) { auto l = make_polygon(a, false); if (!l.empty()) out << l << "\n"; }
    for (auto& a : logi.x_gauges)  { auto l = make_polygon(a, true);  if (!l.empty()) out << l << "\n"; }
    for (auto& a : logi.z_gauges)  { auto l = make_polygon(a, true);  if (!l.empty()) out << l << "\n"; }
    return out.str();
}

}  // namespace dqec
