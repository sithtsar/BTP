# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This project implements a comprehensive Lattice Boltzmann Method (LBM) simulation suite with advanced analysis capabilities:

1. **Single-phase LBM**: Poiseuille flow with standard BGK and entropic BGK collision operators
2. **Multiphase LBM**: Conservative phase-field model for immiscible fluids at high density ratios
3. **H-theorem Analysis**: Automated thermodynamic consistency verification and entropy monitoring
4. **Professional Visualization**: Comprehensive plotting system with scientific formatting
5. **Testing Framework**: Unit tests, integration tests, and performance benchmarks

## Architecture Overview

The codebase follows a modern three-layer architecture:

### **Core Layer (C++)**
- **LBM Solvers** (`include/lbm/`, `src/lbm/`): High-performance simulation engines
- **Analysis System** (`include/analysis/`, `src/analysis/`): H-theorem verification and entropy analysis
- **Utilities** (`include/utils/`, `src/utils/`): Mathematical functions, I/O, and validation
- **Build System**: CMake-based with comprehensive dependency management

### **Analysis Layer (Python)**
- **Visualization Toolkit** (`tools/python/visualization/`): Professional scientific plotting
- **Analysis Runner** (`tools/python/analysis_runner.py`): Automated workflow management
- **Layout System**: Consistent formatting and multi-panel arrangements

### **Testing Layer**
- **Unit Tests** (`tests/unit/`): Mathematical validation and component testing
- **Integration Tests** (`tests/integration/`): End-to-end workflow validation
- **Benchmarks** (`tests/benchmarks/`): Performance analysis and regression prevention

## Key Components

### Single-phase LBM (`include/lbm/single_phase.hpp`)
- **Implementation**: D2Q9 lattice with BGK collision operators
- **Features**: Standard and entropic BGK, H-theorem analysis, stability monitoring
- **Validation**: Poiseuille flow with analytical solution comparison
- **Output**: Velocity profiles, H-function evolution, diagnostic data

### Multiphase LBM (`include/lbm/multiphase.hpp`)
- **Implementation**: Phase-field method for two-phase flows
- **Features**: Conservative interface tracking, high density ratios (1000:1), surface tension
- **Validation**: Interface stability, mass conservation, analytical velocity profiles
- **Output**: Full field data, averaged profiles, interface metrics

### H-theorem Analysis (`include/analysis/h_theorem.hpp`)
- **Purpose**: Thermodynamic consistency verification for LBM simulations
- **Features**: Monotonicity checking, entropy production analysis, anomaly detection
- **Methods**: Both standard and entropic BGK H-function calculations
- **Output**: Evolution data, violation detection, statistical analysis

### Visualization System (`tools/python/visualization/`)
- **LayoutManager**: Professional plot arrangements and scientific formatting
- **SinglePhasePlots**: Velocity profiles, error analysis, H-function evolution
- **MultiphasePlots**: 2D fields, streamlines, interface tracking
- **HTheoremPlots**: Evolution analysis, comparative studies, anomaly detection

## Development Workflow

### Building and Testing
```bash
# Build project
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)

# Run tests
ctest --verbose

# Run benchmarks
./tests/test_performance_benchmarks
```

### Running Simulations
```bash
# Single-phase simulations
echo "standard" | ./bin/lbm_sim    # Standard BGK
echo "entropic" | ./bin/lbm_sim    # Entropic BGK

# Output files: *_velocity_t*.csv, *_h_y_t*.csv, *_h_evolution.csv
```

### Analysis and Visualization
```bash
# Automatic analysis (recommended)
cd tools/python
python analysis_runner.py ../../ --output ../../plots --type auto

# Manual analysis
python analysis_runner.py ../../ --type single_phase --mode standard --output ../../plots
python analysis_runner.py ../../ --type h_theorem --h-file standard_h_evolution.csv --output ../../plots

# Interactive plotting
python analysis_runner.py ../../ --show
```

## Data Structure and File Conventions

### Single-phase Output Files
- **Velocity data**: `{mode}_velocity_t{time}.csv` (y, avg_ux, analytic_ux)
- **H-function profiles**: `{mode}_h_y_t{time}.csv` (y, H_y) [entropic BGK only]
- **H-evolution**: `{mode}_h_evolution.csv` (timestep, h_value, h_change, entropy_production, violations)

### Multiphase Output Files
- **Full fields**: `multiphase_t{time}.csv` (x, y, phi, rho, ux, uy, p, analytical_ux)
- **Averaged profiles**: `multiphase_avg_t{time}.csv` (y, avg_phi, avg_rho, avg_ux, analytical_ux)

### Generated Visualizations
- **Single-phase**: `single_phase_*.png` (velocity profiles, error analysis, H-function evolution)
- **Multiphase**: `multiphase_*.png` (field plots, streamlines, interface evolution)
- **H-theorem**: `h_theorem_*.png` (evolution, statistics, BGK comparison)

## Simulation Parameters

### Single-phase Configuration
```cpp
config.Nx = 20;                     // Grid width
config.Ny = 21;                     // Grid height  
config.tau = 1.0;                   // Relaxation time
config.gravity = 1e-6;              // Body force
config.max_timesteps = 10000;       // Simulation length
config.use_entropic_bgk = false;    // BGK type
config.enable_h_theorem_analysis = true;  // H-theorem monitoring
```

### Multiphase Configuration  
```cpp
config.Nx = 256; config.Ny = 64;    // Grid size
config.rho_L = 1.0; config.rho_H = 1000.0;  // Density ratio
config.mu_L = 0.01; config.mu_H = 1.0;      // Viscosity ratio
config.sigma = 0.01;                 // Surface tension
config.xi = 4.0;                     // Interface thickness
```

## Key Features and Capabilities

### H-theorem Verification
- **Monotonicity**: Automated checking of H-function decrease
- **Entropy Production**: Positive entropy production verification
- **Anomaly Detection**: Statistical outliers and violation identification
- **Comparative Analysis**: Standard vs entropic BGK performance

### Numerical Stability
- **Real-time Monitoring**: Density, velocity, and distribution function bounds
- **Stability Indicators**: NaN/infinity detection, convergence analysis
- **Error Recovery**: Graceful handling of numerical instabilities

### Professional Visualization
- **Scientific Formatting**: Publication-ready plots with proper units and labels
- **Multi-panel Layouts**: Organized subplot arrangements for comprehensive analysis  
- **Interactive Features**: Zoom, pan, and data exploration capabilities
- **Export Formats**: High-resolution PNG, PDF, and vector graphics

### Performance Optimization
- **Memory Management**: Efficient allocation and deallocation strategies
- **Computational Kernels**: Optimized collision and streaming operations
- **Scalability**: Near-linear performance scaling with grid size
- **Benchmarking**: Automated performance regression detection

## Testing and Validation

### Mathematical Validation
- **Mass Conservation**: Machine precision conservation verification
- **Momentum Conservation**: Newton's laws compliance testing
- **Viscosity Scaling**: Kinematic viscosity accuracy validation
- **H-theorem Compliance**: Thermodynamic consistency verification

### Integration Testing
- **Workflow Validation**: End-to-end simulation pipeline testing
- **Parameter Sweeps**: Robustness across parameter ranges
- **Convergence Analysis**: Analytical solution comparison
- **Long-term Stability**: Extended simulation reliability

### Performance Benchmarking
- **Baseline Performance**: Standard and entropic BGK throughput
- **Scaling Analysis**: Grid size and parameter impact assessment
- **Memory Usage**: Allocation efficiency and memory scaling
- **Regression Prevention**: Automated performance monitoring

## Important Development Notes

### Code Style and Conventions
- **C++ Standard**: C++17 with modern features and best practices
- **Memory Safety**: RAII, smart pointers, and bounds checking
- **Error Handling**: Comprehensive exception handling and validation
- **Documentation**: Doxygen-compatible inline documentation

### File Organization
- **Headers**: Well-organized interface definitions with clear dependencies
- **Source**: Implementation files with optimized computational kernels
- **Tests**: Comprehensive coverage with mathematical validation
- **Tools**: Python analysis toolkit with modular design

### Performance Considerations
- **Hot Paths**: Collision and streaming operations are performance-critical
- **Memory Layout**: Optimize for cache efficiency and SIMD operations
- **Parallel Processing**: Design supports future OpenMP/CUDA acceleration
- **Profiling**: Built-in timing and memory usage measurement

## Troubleshooting

### Common Issues
- **Build Failures**: Check GoogleTest installation and CMake version (≥3.15)
- **Test Failures**: Verify numerical precision and parameter ranges
- **Visualization Errors**: Ensure Python dependencies are installed
- **Performance Issues**: Check compiler optimization flags and system resources

### Debug Information
- **Verbose Testing**: Use `ctest --verbose` for detailed test output  
- **Logging**: Enable debug builds for comprehensive diagnostic information
- **Memory Debugging**: Use valgrind or AddressSanitizer for memory issues
- **Performance Profiling**: Use perf or similar tools for hotspot analysis

## Future Extensions

### Planned Features
- **3D Simulations**: Extension to three-dimensional problems
- **GPU Acceleration**: CUDA implementation for high-performance computing
- **Advanced Boundary Conditions**: Complex geometries and moving boundaries
- **Parallel Processing**: MPI support for distributed computing

### Research Applications
- **Microfluidics**: Small-scale fluid flow analysis
- **Multiphase Systems**: Oil-water separation and emulsion dynamics  
- **Porous Media**: Flow through complex geometries
- **Heat Transfer**: Thermal LBM implementation

---

**Development Philosophy**: This codebase emphasizes mathematical rigor, computational efficiency, and scientific accuracy. All numerical methods are validated against analytical solutions and thermodynamic principles. The H-theorem analysis ensures physical consistency, while the comprehensive testing framework prevents regressions and validates correctness across parameter ranges.