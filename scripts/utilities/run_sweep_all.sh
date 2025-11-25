#!/bin/bash

# Comprehensive Sweep Check & Performance Analysis
# Runs all validation cases (analytical + benchmark) for both BGK and ELBM solvers
# Collects timing metrics and prepares data for visualization

set -e  # Exit on error

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="$REPO_ROOT/build"
OUTPUT_DIR="$REPO_ROOT/output"
SCRIPTS_DIR="$REPO_ROOT/scripts"
PLOTTING_DIR="$REPO_ROOT/plotting"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Comprehensive Performance Sweep${NC}"
echo -e "${BLUE}========================================${NC}"

# Phase 0: Build
echo -e "\n${YELLOW}[Phase 0] Building executables...${NC}"
cd "$REPO_ROOT"
./scripts/build.sh > /dev/null 2>&1

if [ ! -f "$BUILD_DIR/test_analytical" ] || [ ! -f "$BUILD_DIR/test_benchmark" ]; then
    echo -e "${RED}Build failed! Missing executables.${NC}"
    exit 1
fi
echo -e "${GREEN}✓ Build successful${NC}"

# Phase 1: Setup directories
echo -e "\n${YELLOW}[Phase 1] Setting up output directories...${NC}"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/sweep_analytical_bgk"
mkdir -p "$OUTPUT_DIR/sweep_analytical_elbm"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re10_bgk"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re10_elbm"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re100_bgk"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re100_elbm"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re1000_bgk"
mkdir -p "$OUTPUT_DIR/sweep_cylinder_re1000_elbm"
mkdir -p "$REPO_ROOT/figures/performance/01_colorgradients"
mkdir -p "$REPO_ROOT/figures/performance/02_breakdown"
mkdir -p "$REPO_ROOT/figures/performance/03_tables"
echo -e "${GREEN}✓ Directories created${NC}"

# Phase 2: Run tests in parallel
echo -e "\n${YELLOW}[Phase 2] Running test suite (parallel execution)...${NC}"
echo "This may take several minutes depending on your system..."

# Store PIDs for all background jobs
declare -a PIDS
declare -a LABELS

# Analytical tests run all cases for both BGK and ELBM
echo "  Starting analytical tests..."
cd "$BUILD_DIR"

# Run all analytical tests (covers both BGK and ELBM internally)
# test_basic_flows runs Couette, Poiseuille, and Taylor-Green for both solvers
timeout 1200 ./test_basic_flows all > "$OUTPUT_DIR/sweep_analytical_bgk/output.log" 2>&1 &
PIDS+=($!)
LABELS+=("Analytical (all cases, both solvers)")

# Cylinder flow tests at Re=10, 100, 1000
for Re in 10 100 1000; do
    echo "  Starting cylinder flow tests (Re=$Re)..."

    # BGK
    timeout 300 ./test_benchmark "cylinder" "$Re" 0 > "$OUTPUT_DIR/sweep_cylinder_re${Re}_bgk/output.log" 2>&1 &
    PIDS+=($!)
    LABELS+=("Cylinder Re=$Re (BGK)")

    # ELBM
    timeout 300 ./test_benchmark "cylinder" "$Re" 1 > "$OUTPUT_DIR/sweep_cylinder_re${Re}_elbm/output.log" 2>&1 &
    PIDS+=($!)
    LABELS+=("Cylinder Re=$Re (ELBM)")
done

# Wait for all jobs and track results
echo -e "\n${YELLOW}Waiting for simulations to complete...${NC}"
FAILED_JOBS=()
for i in "${!PIDS[@]}"; do
    PID=${PIDS[$i]}
    LABEL=${LABELS[$i]}
    if wait $PID 2>/dev/null; then
        echo -e "  ${GREEN}✓${NC} ${LABEL}"
    else
        echo -e "  ${RED}✗${NC} ${LABEL}"
        FAILED_JOBS+=("$LABEL")
    fi
done

if [ ${#FAILED_JOBS[@]} -gt 0 ]; then
    echo -e "\n${YELLOW}Warning: Some tests failed:${NC}"
    for job in "${FAILED_JOBS[@]}"; do
        echo "  - $job"
    done
else
    echo -e "\n${GREEN}All tests completed successfully!${NC}"
fi

# Phase 3: Aggregate results
echo -e "\n${YELLOW}[Phase 3] Aggregating results...${NC}"
cd "$REPO_ROOT"
python3 "$SCRIPTS_DIR/aggregate_results.py" "$OUTPUT_DIR" 2>&1 | tail -20

if [ -f "$OUTPUT_DIR/sweep_summary.csv" ]; then
    echo -e "${GREEN}✓ Results aggregated${NC}"
else
    echo -e "${YELLOW}⚠ Results aggregation may have issues - check output${NC}"
fi

# Phase 4: Generate visualizations
echo -e "\n${YELLOW}[Phase 4] Generating visualizations...${NC}"
python3 "$PLOTTING_DIR/plot_performance_analysis.py" "$OUTPUT_DIR" "$REPO_ROOT/figures/performance" 2>&1 | tail -20

if [ -d "$REPO_ROOT/figures/performance/01_colorgradients" ] && [ $(ls -1 "$REPO_ROOT/figures/performance/01_colorgradients" | wc -l) -gt 0 ]; then
    echo -e "${GREEN}✓ Visualizations generated${NC}"
else
    echo -e "${YELLOW}⚠ Visualization generation may have issues${NC}"
fi

# Phase 5: Summary
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}Performance Sweep Complete!${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "\nResults locations:"
echo -e "  ${BLUE}Summary CSV:${NC} $OUTPUT_DIR/sweep_summary.csv"
echo -e "  ${BLUE}Figures:${NC} $REPO_ROOT/figures/performance/"
echo -e "  ${BLUE}Raw Output:${NC} $OUTPUT_DIR/sweep_*/"

# Print summary table if it exists
if [ -f "$OUTPUT_DIR/sweep_summary.csv" ]; then
    echo -e "\n${YELLOW}Summary Table:${NC}"
    column -t -s',' "$OUTPUT_DIR/sweep_summary.csv" | head -20
fi

echo -e "\n${GREEN}Next steps:${NC}"
echo "  1. Review figures in figures/performance/"
echo "  2. Check sweep_summary.csv for metrics"
echo "  3. Verify stability indicators (BGK at Re=1000 should diverge)"
echo "  4. Analyze ELBM performance scaling"
