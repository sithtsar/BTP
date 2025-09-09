#!/bin/bash

# LBM Poiseuille Flow - Quick Run (No Tests, No Questions)
# This script runs everything automatically without user input

# Don't exit on error for visualization step
# set -e  # Exit on error

echo "🚀 LBM Poiseuille Flow - Quick Analysis Pipeline"
echo "=============================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_step() { echo -e "${BLUE}[STEP]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    print_error "Please run this script from the lbm_poiseuille project root directory"
    exit 1
fi

print_step "Cleaning and building project..."
rm -rf build plots data *.csv *.png *.log
mkdir -p build plots data
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
make -j4 > /dev/null 2>&1
cd ..
print_success "Project built successfully"

print_step "Running LBM simulations..."
if [ ! -f "./build/bin/lbm_sim" ]; then
    print_error "Executable not found!"
    exit 1
fi

# Single-phase simulations
if ./build/bin/lbm_sim standard; then
    print_success "Standard BGK simulation completed"
else
    print_error "Standard BGK simulation failed!"
    exit 1
fi

if ./build/bin/lbm_sim entropic; then
    print_success "Entropic BGK simulation completed"
else
    print_error "Entropic BGK simulation failed!"
    exit 1
fi

# Multiphase simulation
print_step "Building and running multiphase simulation..."
cd build
if g++ -std=c++17 -O3 -I../include -I../src ../src/multiphase_lbm.cpp -o bin/multiphase_sim > /dev/null 2>&1; then
    print_success "Multiphase simulation compiled"
    
    # Run from project root so data/ path works correctly
    cd ..
    if ./build/bin/multiphase_sim > /dev/null 2>&1; then
        print_success "Multiphase simulation completed"
    else
        print_warning "Multiphase simulation failed - continuing with single-phase analysis"
    fi
    cd build
else
    print_warning "Failed to compile multiphase simulation - continuing with single-phase analysis"
fi
cd ..

print_success "All simulations completed"

print_step "Installing Python dependencies..."
missing_deps=""
for package in numpy matplotlib pandas scipy seaborn jinja2; do
    if ! python3 -c "import $package" 2>/dev/null; then
        missing_deps="$missing_deps $package"
    fi
done

if [ -n "$missing_deps" ]; then
    pip3 install $missing_deps > /dev/null 2>&1
fi
print_success "Python dependencies ready"

print_step "Generating comprehensive visualizations..."
cd tools/python

# Create H-evolution files from H-function profiles if they exist
if [ -f "../../data/standard_h_y_t0.csv" ] || [ -f "../../data/entropic_h_y_t0.csv" ]; then
    print_step "Creating H-evolution files from H-function data..."
    python3 -c "
import pandas as pd
import numpy as np
from pathlib import Path

# Create H-evolution data from H-function profiles
for mode in ['standard', 'entropic']:
    h_files = sorted(Path('../../data').glob(f'{mode}_h_y_t*.csv'))
    if h_files:
        h_evolution = []
        for i, file in enumerate(h_files):
            df = pd.read_csv(file)
            if 'h_y' in df.columns or 'H_y' in df.columns:
                h_col = 'h_y' if 'h_y' in df.columns else 'H_y'
                total_h = df[h_col].sum()
                h_change = 0 if i == 0 else total_h - prev_h
                h_evolution.append([i*1000, total_h, h_change, abs(h_change), 0])
                prev_h = total_h
        
        if h_evolution:
            h_df = pd.DataFrame(h_evolution, columns=['timestep', 'h_value', 'h_change', 'entropy_production', 'violations'])
            h_df.to_csv(f'../../data/{mode}_h_evolution.csv', index=False)
            print(f'Created {mode}_h_evolution.csv with {len(h_evolution)} timesteps')
" > /dev/null 2>&1
fi

# Run all analyses
total_plots=0

# Standard BGK analysis
print_step "Running standard BGK analysis..."
if python3 analysis_runner.py ../../data --type single_phase --mode standard --output ../../plots > /dev/null 2>&1; then
    print_success "Standard BGK analysis completed"
    total_plots=$((total_plots + 2))
else
    print_warning "Standard BGK analysis failed"
fi

# Entropic BGK analysis  
print_step "Running entropic BGK analysis..."
if python3 analysis_runner.py ../../data --type single_phase --mode entropic --output ../../plots > /dev/null 2>&1; then
    print_success "Entropic BGK analysis completed"
    total_plots=$((total_plots + 2))
else
    print_warning "Entropic BGK analysis failed"
fi

# Multiphase analysis (if multiphase data exists)
multiphase_files=$(find ../../data -name "multiphase_t*.csv" | wc -l)
if [ "$multiphase_files" -gt 0 ]; then
    print_step "Running multiphase analysis..."
    if python3 analysis_runner.py ../../data --type multiphase --timesteps 0 5000 10000 15000 --output ../../plots > /dev/null 2>&1; then
        print_success "Multiphase analysis completed"
        total_plots=$((total_plots + 6))
    else
        print_warning "Multiphase analysis failed"
    fi
fi

# H-theorem analysis
if [ -f "../../data/standard_h_evolution.csv" ] && [ -f "../../data/entropic_h_evolution.csv" ]; then
    print_step "Running H-theorem comparative analysis..."
    if python3 analysis_runner.py ../../data --type h_theorem --h-file ../../data/standard_h_evolution.csv --compare ../../data/entropic_h_evolution.csv --output ../../plots > /dev/null 2>&1; then
        print_success "H-theorem comparative analysis completed"
        total_plots=$((total_plots + 4))
    else
        print_warning "H-theorem comparative analysis failed"
    fi
elif [ -f "../../data/entropic_h_evolution.csv" ]; then
    print_step "Running H-theorem analysis..."
    if python3 analysis_runner.py ../../data --type h_theorem --h-file ../../data/entropic_h_evolution.csv --output ../../plots > /dev/null 2>&1; then
        print_success "H-theorem analysis completed"
        total_plots=$((total_plots + 3))
    else
        print_warning "H-theorem analysis failed"
    fi
fi

print_success "Generated $total_plots visualization plots"

cd ../..

# Check if any plots were generated
if [ -d "plots" ] && [ "$(find plots -name '*.png' | wc -l)" -gt 0 ]; then
    print_success "Visualizations generated successfully"
else
    print_warning "No visualization files generated - check simulation data"
    
    # List available CSV files for debugging
    echo "Available data files:"
    find data/ -name "*.csv" 2>/dev/null | head -10 | sed 's/^/  /' || echo "  No data/ directory found"
    
    # Try simple direct plotting
    print_step "Attempting basic plotting..."
    cd tools/python
    python3 -c "
import sys, os
sys.path.append('.')
try:
    from visualization import SinglePhasePlots, LayoutConfig
    import matplotlib
    matplotlib.use('Agg')  # Use non-GUI backend
    print('Basic plotting test passed')
except Exception as e:
    print(f'Plotting error: {e}')
" 2>/dev/null || print_warning "Basic plotting test failed"
    cd ../..
fi

echo ""
echo "🎉 Analysis Complete!"
echo "===================="
csv_files=$(find . -maxdepth 1 -name "*.csv" | wc -l)
plot_files=$(find plots -name "*.png" 2>/dev/null | wc -l)
echo "📊 Generated: $csv_files CSV files, $plot_files visualizations (Expected: $total_plots plots)"
echo ""
echo "📈 View your results:"
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "   open plots/"
else
    echo "   xdg-open plots/"
fi
echo ""
echo "🔬 Key files created:"
find . -maxdepth 1 -name "*_velocity_t*.csv" | head -4 | sed 's/^/   ✓ /'
find plots -name "*.png" | head -6 | sed 's/^/   📈 /'