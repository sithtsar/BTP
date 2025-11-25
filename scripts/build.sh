#!/bin/bash

# Create build directory
mkdir -p build
cd build

# Run CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build (use sysctl for macOS, nproc for Linux)
if command -v nproc &> /dev/null; then
    make -j$(nproc)
else
    make -j$(sysctl -n hw.ncpu)
fi

echo "Build complete. Executable: build/elbm_solver"
