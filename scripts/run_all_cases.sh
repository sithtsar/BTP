#!/bin/bash

# Script to run all rectangular pipe flow test cases
# Reproduces the Yu thesis comparison between BGK and ELBM

echo "==================================================="
echo "Running Rectangular Pipe Flow Test Cases"
echo "Comparing BGK vs ELBM at different Reynolds numbers"
echo "==================================================="
echo ""

# Create output directories
mkdir -p output/case1_re10_bgk
mkdir -p output/case1_re10_elbm
mkdir -p output/case2_re100_bgk
mkdir -p output/case2_re100_elbm
mkdir -p output/case3_re1000_bgk
mkdir -p output/case3_re1000_elbm

# Case 1: Re ~ 10 (ν = 0.1 m²/s, Δp = 5.0 Pa)
echo "=== Case 1: Re ~ 10 (ν = 0.1 m²/s) ==="
echo "Running BGK solver..."
./build/elbm_solver --solver BGK --viscosity 0.1 --nx 100 --ny 40 --steps 2500 --output-dir output/case1_re10_bgk

echo "Running ELBM solver..."
./build/elbm_solver --solver ELBM --viscosity 0.1 --nx 100 --ny 40 --steps 2500 --output-dir output/case1_re10_elbm

echo ""

# Case 2: Re ~ 100 (ν = 0.01 m²/s, Δp = 5.0 Pa)
echo "=== Case 2: Re ~ 100 (ν = 0.01 m²/s) ==="
echo "Running BGK solver..."
./build/elbm_solver --solver BGK --viscosity 0.01 --nx 100 --ny 40 --steps 2500 --output-dir output/case2_re100_bgk

echo "Running ELBM solver..."
./build/elbm_solver --solver ELBM --viscosity 0.01 --nx 100 --ny 40 --steps 2500 --output-dir output/case2_re100_elbm

echo ""

# Case 3: Re ~ 1000 (ν = 0.001 m²/s, Δp = 5.0 Pa)
echo "=== Case 3: Re ~ 1000 (ν = 0.001 m²/s) ==="
echo "Running BGK solver..."
./build/elbm_solver --solver BGK --viscosity 0.001 --nx 100 --ny 40 --steps 2500 --output-dir output/case3_re1000_bgk

echo "Running ELBM solver..."
./build/elbm_solver --solver ELBM --viscosity 0.001 --nx 100 --ny 40 --steps 2500 --output-dir output/case3_re1000_elbm

echo ""
echo "==================================================="
echo "All simulations complete!"
echo "Results saved in output/ directory"
echo "==================================================="
