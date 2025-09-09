# LBM Poiseuille Flow Simulation

A comprehensive Lattice Boltzmann Method (LBM) implementation for single-phase and multiphase fluid flow simulations with advanced H-theorem analysis and professional visualization capabilities.

## 🌟 Features

### Core Simulation Capabilities
- **Single-phase LBM**: Poiseuille flow with standard BGK and entropic BGK collision operators
- **Multiphase LBM**: Conservative phase-field model for immiscible fluids at high density ratios
- **H-theorem Analysis**: Automated thermodynamic consistency verification and entropy monitoring
- **Numerical Stability**: Real-time stability monitoring with comprehensive error detection

### Advanced Visualization System
- **Professional Plot Layouts**: Scientifically formatted multi-panel visualizations
- **Interactive Analysis**: Comprehensive error analysis and convergence studies
- **Multiphase Visualization**: 2D field plots, streamlines, and interface tracking
- **H-theorem Plots**: Evolution analysis, statistical validation, and comparative studies

### Comprehensive Testing Framework
- **Unit Tests**: Mathematical validation of collision operators and core components
- **Integration Tests**: End-to-end simulation workflow validation
- **Performance Benchmarks**: Scalability analysis and regression prevention
- **H-theorem Verification**: Automated thermodynamic consistency checking

## 🚀 Quick Start

### Prerequisites
```bash
# C++ dependencies
sudo apt-get install build-essential cmake libgtest-dev  # Ubuntu/Debian
# or
brew install cmake googletest  # macOS

# Python dependencies
pip install numpy matplotlib pandas scipy seaborn jinja2
```

### Building the Project
```bash
# Clone and build
git clone <repository-url>
cd lbm_poiseuille
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### Running Simulations

#### Single-phase Simulations
```bash
# From project root directory
cd build

# Run standard BGK simulation
echo "standard" | ./bin/lbm_sim

# Run entropic BGK simulation  
echo "entropic" | ./bin/lbm_sim
```

#### Multiphase Simulations
```bash
# Compile and run multiphase (if available)
# ./bin/multiphase_app
```

### Visualization and Analysis

#### Automatic Analysis with All Visualizations
```bash
# Run complete analysis pipeline (recommended)
cd tools/python
python analysis_runner.py ../../ --output ../../plots --type auto

# This will:
# 1. Auto-detect simulation type (single-phase, multiphase, or H-theorem)
# 2. Generate all appropriate visualizations
# 3. Save plots to the plots/ directory
# 4. Create diagnostic reports
```

#### Manual Visualization Commands
```bash
cd tools/python

# Single-phase analysis
python analysis_runner.py ../../ --type single_phase --mode standard --output ../../plots
python analysis_runner.py ../../ --type single_phase --mode entropic --output ../../plots

# Multiphase analysis (if data available)
python analysis_runner.py ../../ --type multiphase --output ../../plots --timesteps 0 5000 10000

# H-theorem analysis with BGK comparison
python analysis_runner.py ../../ --type h_theorem --h-file standard_h_evolution.csv --compare entropic_h_evolution.csv --output ../../plots

# Interactive plots (display instead of saving)
python analysis_runner.py ../../ --type auto --show
```

#### Individual Visualization Modules
```bash
# For advanced users - direct module usage
cd tools/python

# Single-phase plots
python -c "
from visualization import SinglePhasePlots, LayoutConfig
config = LayoutConfig(dpi=300, font_size=12)
plotter = SinglePhasePlots(config)
data = plotter.load_data_from_files('../../', 'standard')
figures = plotter.create_comprehensive_analysis(data, '../../plots')
print(f'Generated {len(figures)} plots')
"

# Multiphase plots (if data available)
python -c "
from visualization import MultiphasePlots
plotter = MultiphasePlots()
data = plotter.load_data_from_files('../../')
figures = plotter.create_comprehensive_analysis(data, timesteps_to_plot=[0, 1000, 5000], output_directory='../../plots')
print(f'Generated {len(figures)} plots')
"

# H-theorem analysis
python -c "
from visualization import HTheoremPlots
plotter = HTheoremPlots()
data = plotter.load_h_theorem_data('../../standard_h_evolution.csv')
figures = plotter.create_comprehensive_h_analysis(data, '../../plots')
print(f'Generated {len(figures)} plots')
"
```

## 📊 Generated Visualizations

After running the analysis, you'll find these visualizations in the `plots/` directory:

### Single-phase Analysis
- `single_phase_velocity_profiles.png` - Velocity evolution and analytical comparison
- `single_phase_error_analysis.png` - Convergence analysis and error metrics
- `single_phase_h_function_evolution.png` - H-theorem analysis (entropic BGK only)

### Multiphase Analysis
- `multiphase_fields_t*.png` - 2D field visualizations (phase, density, velocity)
- `multiphase_streamlines_t*.png` - Velocity streamlines with phase field overlay
- `multiphase_interface_evolution.png` - Interface tracking and stability analysis
- `multiphase_profiles_t*.png` - Y-averaged profiles with analytical comparison

### H-theorem Analysis
- `h_theorem_h_evolution.png` - H-function evolution with anomaly detection
- `h_theorem_statistical_analysis.png` - Spatial variance and energy evolution
- `h_theorem_entropy_production.png` - Entropy production analysis and distributions
- `h_theorem_bgk_comparison.png` - Standard vs entropic BGK comparison
- `h_theorem_diagnostic_report.txt` - Comprehensive analysis report

## 🧪 Testing

### Run All Tests
```bash
cd build
ctest --verbose
```

### Run Specific Test Categories
```bash
# Unit tests only
ctest -R "test_.*" --verbose

# Integration tests
ctest -R "test_.*_flow" --verbose

# Performance benchmarks (not run automatically)
./tests/test_performance_benchmarks
```

### Individual Test Execution
```bash
# Mathematical validation
./tests/test_collision_operators

# H-theorem verification
./tests/test_h_theorem_verification

# Complete workflow testing
./tests/test_poiseuille_flow
./tests/test_multiphase_flow
```

## 📁 Project Structure

```
lbm_poiseuille/
├── build/                          # Build directory (created by cmake)
├── include/                         # Header files
│   ├── lbm/                        # Core LBM solvers
│   ├── analysis/                   # H-theorem analysis
│   └── utils/                      # Utilities and validation
├── src/                            # Source implementations
├── tests/                          # Comprehensive test suite
│   ├── unit/                       # Unit tests
│   ├── integration/                # Integration tests
│   └── benchmarks/                 # Performance benchmarks
├── tools/python/                   # Python visualization toolkit
│   ├── visualization/              # Visualization modules
│   └── analysis_runner.py          # Main analysis script
├── plots/                          # Generated visualizations (created)
└── reports/                        # Report generation system
```

## 🔬 Simulation Parameters

### Single-phase Configuration
- **Grid**: 20×21 (Nx × Ny)
- **Time steps**: 10,000 with output every 1,000 steps
- **Relaxation time**: τ = 1.0
- **Body force**: 1×10⁻⁶
- **BGK methods**: Standard and entropic collision operators

### Multiphase Configuration
- **Grid**: 256×64 (Nx × Ny)
- **Time steps**: 20,000 with output every 1,000 steps
- **Density ratio**: 1000:1 (ρ_H/ρ_L)
- **Viscosity ratio**: 100:1 (μ_H/μ_L)
- **Interface thickness**: 4 lattice units
- **Surface tension**: Included via chemical potential

## 📈 Performance Metrics

The testing framework provides comprehensive performance analysis:

- **Single-phase baseline**: >50 timesteps/second on modern hardware
- **Entropic BGK overhead**: <50% performance impact
- **Multiphase performance**: >5 timesteps/second for typical problems 
- **Memory scaling**: ~0.2 MB per grid point (single-phase)
- **Grid scaling**: Near-linear performance scaling with problem size

## 🤝 Contributing

1. **Code Style**: Follow existing C++17 conventions
2. **Testing**: Add tests for new features
3. **Documentation**: Update README and inline documentation
4. **Validation**: Ensure H-theorem compliance for new collision operators

## 📝 Citation

If you use this code in your research, please cite:

```bibtex
@software{lbm_poiseuille_2025,
  title = {LBM Poiseuille Flow Simulation with H-theorem Analysis},
  author = {[Your Name]},
  year = {2025},
  url = {[Repository URL]}
}
```

## 📧 Support

- **Issues**: Report bugs and request features via GitHub issues
- **Documentation**: See `CLAUDE.md` for development details
- **Testing**: Run `ctest --verbose` for validation

---

**⚠️ Note**: This project implements advanced LBM methods with rigorous mathematical validation. The H-theorem analysis provides thermodynamic consistency verification essential for accurate fluid dynamics simulations.