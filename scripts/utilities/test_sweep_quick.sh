#!/bin/bash

# Quick test of the sweep infrastructure (for validation purposes)
# Runs a single cylinder test at Re=10 for both solvers

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$REPO_ROOT/build"
OUTPUT_DIR="$REPO_ROOT/output"

echo "Quick Sweep Test - Single Re=10 case"
echo "===================================="

mkdir -p "$OUTPUT_DIR/sweep_cylinder_re10_bgk"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re10_elbm"

cd "$BUILD_DIR"

echo "Running BGK at Re=10..."
./test_benchmark cylinder 10 0 > "$OUTPUT_DIR/sweep_cylinder_re10_bgk/output.log" 2>&1

echo "Running ELBM at Re=10..."
./test_benchmark cylinder 10 1 > "$OUTPUT_DIR/sweep_cylinder_re10_elbm/output.log" 2>&1

echo "Aggregating results..."
cd "$REPO_ROOT"
python3 scripts/aggregate_results.py "$OUTPUT_DIR" "$OUTPUT_DIR/sweep_summary.csv"

echo "Generating visualizations..."
python3 plotting/plot_performance_analysis.py "$OUTPUT_DIR" "$REPO_ROOT/figures/performance"
python3 plotting/plot_cylinder_spatial_analysis.py "$OUTPUT_DIR" "$REPO_ROOT/figures/performance"

echo ""
echo "Quick test complete!"
echo "Results: $OUTPUT_DIR/sweep_summary.csv"
echo "Figures: $REPO_ROOT/figures/performance/"
