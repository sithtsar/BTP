#!/bin/bash

# Full Performance Sweep - All Cases
# Runs analytical tests + cylinder flows at Re=10, 100, 1000 for both BGK and ELBM

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$REPO_ROOT/build"
OUTPUT_DIR="$REPO_ROOT/output"

echo "=========================================="
echo "Full Performance Sweep - All Cases"
echo "=========================================="

mkdir -p "$OUTPUT_DIR/sweep_analytical"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re10_bgk"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re10_elbm"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re100_bgk"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re100_elbm"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re1000_bgk"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re1000_elbm"

cd "$BUILD_DIR"

echo ""
echo "Step 1/4: Analytical Tests (Couette, Poiseuille, Taylor-Green for both solvers)"
echo "====="
timeout 1200 ./test_basic_flows all > "$OUTPUT_DIR/sweep_analytical/output.log" 2>&1 &
ANALYTICAL_PID=$!

echo ""
echo "Step 2-4: Cylinder Flow Tests (Re=10, 100, 1000 for both solvers - parallel)"
echo "====="

# Re=10 - BGK
echo "  [Re=10] Running BGK..."
timeout 600 ./test_benchmark cylinder 10 0 > "$OUTPUT_DIR/sweep_cylinder_re10_bgk/output.log" 2>&1 &
PID_10_BGK=$!

# Re=10 - ELBM
echo "  [Re=10] Running ELBM..."
timeout 600 ./test_benchmark cylinder 10 1 > "$OUTPUT_DIR/sweep_cylinder_re10_elbm/output.log" 2>&1 &
PID_10_ELBM=$!

# Re=100 - BGK
echo "  [Re=100] Running BGK..."
timeout 600 ./test_benchmark cylinder 100 0 > "$OUTPUT_DIR/sweep_cylinder_re100_bgk/output.log" 2>&1 &
PID_100_BGK=$!

# Re=100 - ELBM
echo "  [Re=100] Running ELBM..."
timeout 600 ./test_benchmark cylinder 100 1 > "$OUTPUT_DIR/sweep_cylinder_re100_elbm/output.log" 2>&1 &
PID_100_ELBM=$!

# Re=1000 - BGK
echo "  [Re=1000] Running BGK..."
timeout 600 ./test_benchmark cylinder 1000 0 > "$OUTPUT_DIR/sweep_cylinder_re1000_bgk/output.log" 2>&1 &
PID_1000_BGK=$!

# Re=1000 - ELBM
echo "  [Re=1000] Running ELBM..."
timeout 600 ./test_benchmark cylinder 1000 1 > "$OUTPUT_DIR/sweep_cylinder_re1000_elbm/output.log" 2>&1 &
PID_1000_ELBM=$!

# Wait for all cylinder tests
echo ""
echo "Waiting for simulations..."
wait $PID_10_BGK && echo "  ✓ Re=10 BGK completed" || echo "  ✗ Re=10 BGK failed"
wait $PID_10_ELBM && echo "  ✓ Re=10 ELBM completed" || echo "  ✗ Re=10 ELBM failed"
wait $PID_100_BGK && echo "  ✓ Re=100 BGK completed" || echo "  ✗ Re=100 BGK failed"
wait $PID_100_ELBM && echo "  ✓ Re=100 ELBM completed" || echo "  ✗ Re=100 ELBM failed"
wait $PID_1000_BGK && echo "  ✓ Re=1000 BGK completed" || echo "  ✗ Re=1000 BGK failed"
wait $PID_1000_ELBM && echo "  ✓ Re=1000 ELBM completed" || echo "  ✗ Re=1000 ELBM failed"

# Wait for analytical tests
wait $ANALYTICAL_PID && echo "  ✓ Analytical tests completed" || echo "  ✗ Analytical tests failed"

echo ""
echo "Post-processing..."

# Move CSV outputs to proper locations
cd "$REPO_ROOT"
for Re in 10 100 1000; do
    for SOLVER in BGK ELBM; do
        CSV_FILE="$BUILD_DIR/LBM_output_Re${Re}_${SOLVER}.csv"
        if [ -f "$CSV_FILE" ]; then
            mkdir -p "$OUTPUT_DIR/sweep_cylinder_re${Re}_$(echo $SOLVER | tr '[:upper:]' '[:lower:]')"
            cp "$CSV_FILE" "$OUTPUT_DIR/sweep_cylinder_re${Re}_$(echo $SOLVER | tr '[:upper:]' '[:lower:]')/"
        fi
    done
done

echo "Aggregating results..."
python3 scripts/aggregate_results.py "$OUTPUT_DIR" "$OUTPUT_DIR/sweep_summary.csv"

echo ""
echo "Generating performance visualizations..."
source .venv/bin/activate 2>/dev/null || true
uv run plotting/plot_performance_analysis.py "$OUTPUT_DIR" "$REPO_ROOT/figures/performance"
uv run plotting/plot_cylinder_spatial_analysis.py "$OUTPUT_DIR" "$REPO_ROOT/figures/performance"

echo ""
echo "=========================================="
echo "Full Sweep Complete!"
echo "=========================================="
echo ""
echo "Results:"
echo "  Summary CSV: $OUTPUT_DIR/sweep_summary.csv"
echo "  Figures: $REPO_ROOT/figures/performance/"
echo ""
echo "Performance Profiles Generated:"
echo "  - Spatial cost distribution for each Re case (2_breakdown/)"
echo "  - Timing and throughput metrics (3_tables/)"
echo "  - Velocity field spatial analysis"
echo ""
