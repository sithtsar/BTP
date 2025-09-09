#!/bin/bash

# LBM Poiseuille Flow - Complete Simulation and Analysis Pipeline
# This script runs the complete workflow from building to visualization

# Exit on error only for critical sections
# Visualization errors shouldn't stop the script

echo "🚀 LBM Poiseuille Flow - Complete Analysis Pipeline"
echo "=================================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if we're in the right directory
if [ ! -f "CMakeLists.txt" ]; then
    print_error "Please run this script from the lbm_poiseuille project root directory"
    exit 1
fi

# Step 1: Clean previous builds and outputs
print_step "Cleaning previous builds and outputs..."
rm -rf build plots data *.csv *.png *.log
mkdir -p plots data
print_success "Cleaned previous builds and outputs"

# Step 2: Build the project
print_step "Building the project..."
mkdir -p build
cd build

# Configure with CMake
print_step "Configuring with CMake..."
if ! cmake .. -DCMAKE_BUILD_TYPE=Release; then
    print_error "CMake configuration failed!"
    exit 1
fi

# Build the project
print_step "Compiling the project..."
if command -v nproc >/dev/null 2>&1; then
    BUILD_JOBS=$(nproc)
elif command -v sysctl >/dev/null 2>&1; then
    BUILD_JOBS=$(sysctl -n hw.ncpu)
else
    BUILD_JOBS=4
fi

if ! make -j${BUILD_JOBS}; then
    print_error "Build failed!"
    exit 1
fi

cd ..
print_success "Project built successfully"

# Step 3: Run tests (optional, but recommended)
read -p "🧪 Do you want to run the test suite? (y/N): " run_tests
if [[ $run_tests =~ ^[Yy]$ ]]; then
    print_step "Running test suite..."
    cd build
    
    # Try to fix GoogleTest library path issue on macOS
    if [[ "$OSTYPE" == "darwin"* ]]; then
        export DYLD_LIBRARY_PATH="/opt/homebrew/lib:/usr/local/lib:$DYLD_LIBRARY_PATH"
        print_step "Setting library paths for macOS..."
    fi
    
    # Run tests with timeout and error handling
    if command -v timeout >/dev/null 2>&1; then
        if timeout 300 ctest --output-on-failure 2>/dev/null; then
            print_success "All tests passed"
        else
            print_warning "Some tests failed or timed out - this won't affect the main simulation"
            print_step "Continuing with simulation..."
        fi
    else
        # No timeout command available
        if ctest --output-on-failure 2>/dev/null; then
            print_success "All tests passed"
        else
            print_warning "Some tests failed - this won't affect the main simulation"
            print_step "Continuing with simulation..."
        fi
    fi
    cd ..
else
    print_step "Skipping tests - proceeding to simulation"
fi

# Step 4: Run simulations
print_step "Running LBM simulations..."

# Check if executable exists
if [ ! -f "./build/bin/lbm_sim" ]; then
    print_error "Executable ./build/bin/lbm_sim not found!"
    exit 1
fi

# Run standard BGK simulation
print_step "Running standard BGK simulation..."
./build/bin/lbm_sim standard
if [ $? -eq 0 ]; then
    print_success "Standard BGK simulation completed"
else
    print_error "Standard BGK simulation failed"
    exit 1
fi

# Run entropic BGK simulation
print_step "Running entropic BGK simulation..."
./build/bin/lbm_sim entropic
if [ $? -eq 0 ]; then
    print_success "Entropic BGK simulation completed"
else
    print_error "Entropic BGK simulation failed"
    exit 1
fi

# Check if simulation outputs exist
print_step "Verifying simulation outputs..."
standard_files=$(find . -maxdepth 1 -name "standard_velocity_t*.csv" | wc -l)
entropic_files=$(find . -maxdepth 1 -name "entropic_velocity_t*.csv" | wc -l)

if [ "$standard_files" -eq 0 ]; then
    print_error "No standard BGK simulation output files found"
    echo "Expected files like: standard_velocity_t*.csv"
    exit 1
fi

if [ "$entropic_files" -eq 0 ]; then
    print_error "No entropic BGK simulation output files found"
    echo "Expected files like: entropic_velocity_t*.csv"
    exit 1
fi

print_success "Simulation output files verified ($standard_files standard, $entropic_files entropic)"

print_success "All simulations completed successfully"

# Step 5: Install Python dependencies if needed
print_step "Checking Python dependencies..."

# Check if python3 exists
if ! command -v python3 &> /dev/null; then
    print_error "Python 3 is not installed. Please install Python 3 and try again."
    exit 1
fi

# Check and install Python dependencies
missing_deps=""
for package in numpy matplotlib pandas scipy seaborn jinja2; do
    if ! python3 -c "import $package" 2>/dev/null; then
        missing_deps="$missing_deps $package"
    fi
done

if [ -n "$missing_deps" ]; then
    print_warning "Missing Python packages:$missing_deps"
    print_step "Installing missing packages..."
    
    # Try pip3 first, then pip
    if command -v pip3 &> /dev/null; then
        pip3 install $missing_deps
    elif command -v pip &> /dev/null; then
        pip install $missing_deps
    else
        print_error "pip is not installed. Please install pip and the following packages manually:"
        echo "$missing_deps"
        exit 1
    fi
    
    # Verify installation
    for package in $missing_deps; do
        if ! python3 -c "import $package" 2>/dev/null; then
            print_error "Failed to install $package. Please install manually."
            exit 1
        fi
    done
fi

print_success "Python dependencies verified"

# Step 6: Run comprehensive analysis
print_step "Running comprehensive analysis and visualization..."

# Create plots directory
mkdir -p plots

cd tools/python

# Auto-detect and analyze all data
print_step "Auto-detecting simulation data and generating visualizations..."
if python3 analysis_runner.py ../../ --output ../../plots --type auto; then
    print_success "Auto-analysis completed"
else
    print_warning "Auto-analysis failed, trying specific analyses..."
fi

# Generate specific analyses
print_step "Generating single-phase analysis..."
if python3 analysis_runner.py ../../ --type single_phase --mode standard --output ../../plots 2>/dev/null; then
    print_success "Standard BGK analysis completed"
else
    print_warning "Standard BGK analysis failed - checking data files..."
fi

if python3 analysis_runner.py ../../ --type single_phase --mode entropic --output ../../plots 2>/dev/null; then
    print_success "Entropic BGK analysis completed"
else
    print_warning "Entropic BGK analysis failed - checking data files..."
fi

# Generate H-theorem analysis if H-evolution files exist
if [ -f "../../standard_h_evolution.csv" ] && [ -f "../../entropic_h_evolution.csv" ]; then
    print_step "Generating H-theorem comparative analysis..."
    if python3 analysis_runner.py ../../ --type h_theorem --h-file ../../standard_h_evolution.csv --compare ../../entropic_h_evolution.csv --output ../../plots; then
        print_success "H-theorem comparative analysis completed"
    else
        print_warning "H-theorem comparative analysis failed"
    fi
elif [ -f "../../entropic_h_evolution.csv" ]; then
    print_step "Generating H-theorem analysis for entropic BGK..."
    if python3 analysis_runner.py ../../ --type h_theorem --h-file ../../entropic_h_evolution.csv --output ../../plots; then
        print_success "H-theorem analysis completed"
    else
        print_warning "H-theorem analysis failed"
    fi
elif [ -f "../../standard_h_evolution.csv" ]; then
    print_step "Generating H-theorem analysis for standard BGK..."
    if python3 analysis_runner.py ../../ --type h_theorem --h-file ../../standard_h_evolution.csv --output ../../plots; then
        print_success "H-theorem analysis completed"
    else
        print_warning "H-theorem analysis failed"
    fi
else
    print_warning "No H-evolution files found - skipping H-theorem analysis"
fi

cd ../..
print_success "Analysis and visualization completed"

# Step 7: Display results
print_step "Analysis Results Summary"
echo "========================="

# Count generated files
csv_files=$(find . -maxdepth 1 -name "*.csv" | wc -l)
plot_files=$(find plots -name "*.png" 2>/dev/null | wc -l)

echo "📊 Generated data files: $csv_files CSV files"
echo "📈 Generated visualizations: $plot_files plot files"

# List key output files
echo ""
echo "Key Output Files:"
echo "----------------"
echo "Simulation Data:"
find . -maxdepth 1 -name "*_velocity_t*.csv" | head -6 | sed 's/^/  ✓ /'
find . -maxdepth 1 -name "*_h_evolution.csv" | sed 's/^/  ✓ /'

echo ""
echo "Generated Visualizations:"
if [ -d "plots" ]; then
    find plots -name "*.png" | sed 's/^/  📈 /'
else
    echo "  No plots directory found"
fi

# Step 8: Performance summary
if [ -f "build/tests/test_performance_benchmarks" ]; then
    read -p "🏃 Do you want to run performance benchmarks? (y/N): " run_benchmarks
    if [[ $run_benchmarks =~ ^[Yy]$ ]]; then
        print_step "Running performance benchmarks (this may take a few minutes)..."
        cd build
        ./tests/test_performance_benchmarks
        cd ..
        print_success "Performance benchmarks completed"
    fi
fi

# Step 9: Final instructions
echo ""
print_success "🎉 Complete analysis pipeline finished!"
echo ""
echo "Next Steps:"
echo "----------"
echo "1. 📊 View visualizations in the 'plots/' directory"
echo "2. 📈 Open plot files with your preferred image viewer"
echo "3. 📋 Review simulation data in CSV files"
echo "4. 🔬 Examine H-theorem diagnostic reports (if generated)"
echo ""
echo "Visualization Commands:"
echo "----------------------"
# macOS
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "📱 Open plots folder:    open plots/"
    echo "🖼️  View all plots:      open plots/*.png"
# Linux
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "📁 Open plots folder:    xdg-open plots/"
    echo "🖼️  View all plots:      eog plots/*.png"
fi
echo ""
echo "Interactive Analysis:"
echo "--------------------"
echo "🎯 Interactive plots:   cd tools/python && python3 analysis_runner.py ../../ --show"
echo "🔧 Custom analysis:     cd tools/python && python3 -c \"from visualization import *; # Your code here\""
echo ""
echo "Testing and Validation:"
echo "----------------------"
echo "🧪 Run all tests:       cd build && ctest --verbose"
echo "⚡ Run benchmarks:      cd build && ./tests/test_performance_benchmarks"
echo ""
echo "📚 For more information, see README.md and CLAUDE.md"

# Check if plots were actually generated
if [ -d "plots" ] && [ "$(find plots -name '*.png' | wc -l)" -gt 0 ]; then
    print_success "Analysis pipeline completed successfully! 🎉"
    echo ""
    echo "📊 Your visualizations are ready in the 'plots/' directory"
    echo "🔬 Review the generated plots to analyze your LBM simulation results"
else
    print_warning "Analysis completed but no visualization files were generated"
    echo "📋 Check the console output above for any error messages"
    echo "🔧 You can run the visualization manually with:"
    echo "   cd tools/python && python3 analysis_runner.py ../../ --output ../../plots --type auto"
fi