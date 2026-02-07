cmake_minimum_required(VERSION 3.20)
project(bucket_simulator VERSION 1.0 LANGUAGES CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add Stim as a subdirectory
set(STIM_BUILD_TESTS OFF CACHE BOOL "Disable Stim test build")
add_subdirectory(Stim)

# Add Pymatching as a subdirectory
add_subdirectory(Pymatching)

# Source files
set(SIMULATOR_FILES
    src/main.cpp
    src/simulator.cpp
    src/config.cpp
    src/decoder.cpp
    src/dem_analyzer.cpp
    src/bucket_sampler.cpp
    src/lattice_surgery.cpp
)

# Create executable
add_executable(bucket_simulator ${SIMULATOR_FILES})

# Optimization flags
if (CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(bucket_simulator PUBLIC -O3)
else()
    target_compile_options(bucket_simulator PUBLIC -ggdb3)
endif()

# Link libraries
target_link_libraries(bucket_simulator PRIVATE libstim libpymatching)

# Include directories
target_include_directories(bucket_simulator PRIVATE include)

# Define project source directory for default output path
target_compile_definitions(bucket_simulator PRIVATE PROJECT_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# MPI configuration
find_package(MPI REQUIRED)
target_include_directories(bucket_simulator PRIVATE
    ${MPI_C_INCLUDE_PATH}
    ${MPI_CXX_INCLUDE_PATH}
)
target_link_libraries(bucket_simulator PRIVATE
    ${MPI_C_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)
