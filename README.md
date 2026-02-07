# Bucket Simulator - Surface Code Simulator

A high-performance quantum error correction simulator for surface codes using C++, Stim, PyMatching's Sparse Blossom decoder, and MPI for parallelization.

## Features

- Pure C++ implementation for maximum performance
- Stim library for fast quantum circuit simulation
- PyMatching Sparse Blossom for MWPM decoding
- MPI parallelization for scalable shot distribution
- Text-based configuration files
- Detailed output with logical error rates and runtime statistics

## Dependencies

- C++20 compatible compiler (GCC 10+, Clang 12+, or equivalent)
- CMake 3.20 or later
- MPI implementation (OpenMPI, MPICH, or Intel MPI)
- Git (for submodule management)

## Building

### 1. Clone the repository with submodules

```bash
git clone <repository-url> bucket_simulator
cd bucket_simulator
git submodule update --init --recursive
```

### 2. Build using CMake

```bash
mkdir build
cd build
cmake ..
make -j8
```

This will create the `bucket_simulator` executable in the `build` directory.

### Build Types

For optimized release build:
```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```

For debug build with symbols:
```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j8
```

## Configuration

Create a text configuration file (see [configs/example.txt](configs/example.txt)) with the following parameters:

```
code_distance 5          # Surface code distance
rounds 5                 # Number of QEC rounds
code_type rotated_memory_x  # Code type
physical_error 0.001     # Physical error rate
measurement_error 0.001  # Measurement error rate
reset_error 0.001        # Reset error rate
total_shots 10M          # Number of shots (supports K, M, B, G suffixes)
```

### Configuration Parameters

- `code_distance`: Distance of the surface code (odd integer, typically 3, 5, 7, ...)
- `rounds`: Number of syndrome extraction rounds per shot
- `code_type`: Type of surface code (e.g., `rotated_memory_x`, `rotated_memory_z`)
- `physical_error`: Depolarization error rate for data qubits (0.0 to 1.0)
- `measurement_error`: Bit flip probability for measurements (0.0 to 1.0)
- `reset_error`: Bit flip probability for qubit resets (0.0 to 1.0)
- `total_shots`: Total number of Monte Carlo shots to simulate

Magnitude suffixes:
- K = 1,000 (thousands)
- M = 1,000,000 (millions)
- B = 1,000,000,000 (billions)
- G = 1,000,000,000 (billions)

## Usage

### Single Process

Run the simulator with a single process (from build directory):

```bash
cd build
./bucket_simulator -config ../configs/example.txt
```

Or from project root:

```bash
./build/bucket_simulator -config configs/example.txt
```

### Multiple MPI Ranks

Run with MPI for parallel execution:

```bash
cd build
mpirun -np 4 ./bucket_simulator -config ../configs/example.txt
```

This will distribute the shots across 4 MPI processes.

### Command-Line Options

- `-config <file>`: Path to configuration file (required)
- `-output <dir>`: Output directory (default: `<project_root>/output`)
- `-h`, `--help`: Display usage information

### Custom Output Directory

Specify a custom location for results:

```bash
./bucket_simulator -config configs/example.txt -output /path/to/results
```

## Output

Results are written to the `output/` directory with a timestamp:

```
output/results_YYYYMMDD_HHMMSS.txt
```

### Output Format

```
Configuration:
  Code Distance: 5
  Rounds: 5
  Physical Error Rate: 0.001000
  Measurement Error Rate: 0.001000
  Reset Error Rate: 0.001000
  Total Shots: 10000000
  Code Type: rotated_memory_x

Results:
  Logical Error Rate: 0.001234 (12340/10000000)

Runtime Statistics:
  Total Runtime: 0 hours, 15 minutes, 32 seconds
  Shots per Second: 10752

Per-Rank Statistics:
  Rank 0: 2500000 shots, 18.5s, 135135 shots/s
  Rank 1: 2500000 shots, 18.2s, 137363 shots/s
  Rank 2: 2500000 shots, 18.7s, 133690 shots/s
  Rank 3: 2500000 shots, 18.4s, 135870 shots/s
```

## Performance

The simulator is designed for high performance:

- **Stim**: Extremely fast circuit sampling (>1 kHz for distance 100 codes)
- **Sparse Blossom**: 100-1000x faster than previous MWPM implementations
- **MPI Parallelization**: Linear scaling with number of processes (embarrassingly parallel)
- **Batch Processing**: Processes 1M shots per batch for memory efficiency

Typical performance (distance 5, 10M shots):
- Single core: ~30-60 seconds
- 4 cores: ~8-15 seconds (4x speedup)
- 16 cores: ~2-4 seconds (16x speedup)

## Project Structure

```
bucket_simulator/
├── CMakeLists.txt           # Build configuration
├── README.md                # This file
├── .gitignore              # Git ignore rules
├── Stim/                    # Stim library source (git submodule)
├── Pymatching/              # PyMatching library source (git submodule)
├── src/                     # Source files
│   ├── main.cpp            # MPI orchestration & output
│   ├── simulator.cpp       # Main simulation loop
│   ├── config.cpp          # Configuration parser
│   └── decoder.cpp         # MWPM decoder interface
├── include/                 # Header files
│   ├── simulator.hpp       # Simulator class
│   ├── config.hpp          # Configuration structures
│   └── decoder.hpp         # Decoder interface
├── configs/                 # Configuration files
│   ├── example.txt         # 10M shots example
│   └── test.txt            # Quick test config
├── output/                  # Simulation results (default location)
│   └── results_*.txt       # Timestamped result files
└── build/                   # Build artifacts (gitignored)
    ├── bucket_simulator    # Compiled executable
    ├── Stim/               # Stim build outputs (.o, .a files)
    ├── Pymatching/         # PyMatching build outputs
    └── ...                 # Other CMake artifacts
```

**Note**: The `build/Stim/` and `build/Pymatching/` directories contain only compiled build artifacts (object files, libraries), not source code. The actual source code is in the root-level `Stim/` and `Pymatching/` subdirectories.

## Troubleshooting

### Build Issues

**Problem**: Cannot find MPI
```
Solution: Ensure MPI is installed and in your PATH
  - Ubuntu/Debian: sudo apt-get install libopenmpi-dev
  - macOS: brew install open-mpi
  - Set MPI_HOME if needed: export MPI_HOME=/path/to/mpi
```

**Problem**: Submodules not initialized
```
Solution: Run git submodule update --init --recursive
```

**Problem**: C++20 not supported
```
Solution: Update your compiler to GCC 10+, Clang 12+, or equivalent
```

### Runtime Issues

**Problem**: Output file not created
```
Solution: Ensure you have write permissions in the current directory
The output/ directory is created automatically if it doesn't exist
```

**Problem**: MPI errors
```
Solution: Check your MPI installation and ensure mpirun is in PATH
Test with: mpirun -np 2 hostname
```

## Contributing

Contributions are welcome! Please ensure:
- Code follows existing style conventions
- All tests pass before submitting
- Documentation is updated for new features

## License

[Add your license here]

## References

- [Stim](https://github.com/quantumlib/Stim) - Quantum stabilizer circuit simulator
- [PyMatching](https://github.com/oscarhiggott/PyMatching) - MWPM decoder for QEC
- [Surface Codes](https://arxiv.org/abs/1208.0928) - Introduction to surface codes

## Acknowledgments

Based on the architecture of the DistributedQEC project.
