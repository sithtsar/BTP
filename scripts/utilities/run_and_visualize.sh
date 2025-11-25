#!/bin/bash
# Complete workflow: Build → Run → Visualize

set -e  # Exit on error

echo "=========================================="
echo "ELBM Complete Workflow"
echo "=========================================="

# Step 1: Build
echo -e "\n[1/4] Building..."
./scripts/build.sh

# Step 2: Run analytical tests
echo -e "\n[2/4] Running analytical validation..."
./build/test_basic_flows all

# Step 3: Run cylinder benchmarks
echo -e "\n[3/4] Running cylinder flow benchmarks..."
./build/test_benchmark

# Step 4: Generate all visualizations
echo -e "\n[4/4] Generating visualizations..."
python plotting/generate_all_visualizations.py

echo -e "\n=========================================="
echo "Workflow Complete!"
echo "=========================================="
echo "View results:"
echo "  - Geometry: figures/geometry_diagrams/"
echo "  - Analytical: figures/analytical_validation/"
echo "  - Cylinder: figures/channel_cylinder/"
echo "=========================================="
