#!/usr/bin/env python3
"""
Complete LBM Analysis Workflow
Runs simulations, generates plots, and creates presentations in one command.
"""

import subprocess
import sys
import os
from pathlib import Path
import time
import argparse

def run_command(command, description, cwd=None, timeout=None):
    """Run a command with error handling and progress reporting."""
    print(f"\n{'='*60}")
    print(f"STEP: {description}")
    print(f"Command: {' '.join(command) if isinstance(command, list) else command}")
    print(f"{'='*60}")
    
    try:
        start_time = time.time()
        result = subprocess.run(
            command, 
            cwd=cwd, 
            check=True, 
            capture_output=True, 
            text=True,
            timeout=timeout
        )
        
        elapsed = time.time() - start_time
        print(f"✓ SUCCESS ({elapsed:.1f}s)")
        
        if result.stdout:
            print("Output:", result.stdout[:500] + "..." if len(result.stdout) > 500 else result.stdout)
        
        return True, result.stdout
        
    except subprocess.CalledProcessError as e:
        print(f"✗ FAILED (exit code {e.returncode})")
        print(f"Error: {e.stderr}")
        return False, e.stderr
        
    except subprocess.TimeoutExpired:
        print(f"✗ TIMEOUT after {timeout}s")
        return False, "Command timed out"
        
    except Exception as e:
        print(f"✗ ERROR: {e}")
        return False, str(e)

def check_dependencies():
    """Check if all required dependencies are available."""
    print("Checking dependencies...")
    
    issues = []
    
    # Check C++ compiler
    try:
        subprocess.run(['g++', '--version'], capture_output=True, check=True)
        print("  ✓ g++ compiler available")
    except:
        issues.append("g++ compiler not found")
    
    # Check Python packages
    required_packages = ['numpy', 'matplotlib', 'pandas']
    for package in required_packages:
        try:
            __import__(package)
            print(f"  ✓ {package} available")
        except ImportError:
            issues.append(f"Python package '{package}' not installed")
    
    # Check Node.js for presentations
    try:
        subprocess.run(['node', '--version'], capture_output=True, check=True)
        print("  ✓ Node.js available")
    except:
        print("  ! Node.js not found (presentations will be skipped)")
    
    if issues:
        print(f"\n⚠️  Issues found:")
        for issue in issues:
            print(f"    - {issue}")
        return False
    
    print("  ✓ All dependencies available")
    return True

def compile_simulations():
    """Compile both single-phase and multiphase simulations."""
    project_root = Path.cwd()
    
    # Compile single-phase
    success, _ = run_command(
        ['g++', '-O3', '-o', 'lbm_sim', 'src/main.cpp'],
        "Compiling single-phase simulation",
        cwd=project_root,
        timeout=60
    )
    
    if not success:
        return False
    
    # Compile multiphase
    success, _ = run_command(
        ['g++', '-O3', '-o', 'lbm_multiphase', 'src/multiphase_lbm.cpp'],
        "Compiling multiphase simulation", 
        cwd=project_root,
        timeout=60
    )
    
    return success

def run_simulations(skip_existing=True):
    """Run both simulations if needed."""
    project_root = Path.cwd()
    data_dir = project_root / 'data'
    data_dir.mkdir(exist_ok=True)
    
    # Check if we should skip existing data
    if skip_existing:
        single_data = list(data_dir.glob('standard_velocity_t*.csv'))
        multi_data = list(data_dir.glob('multiphase_avg_t*.csv'))
        
        if single_data and multi_data:
            print("  ℹ️  Simulation data already exists, skipping simulations")
            print(f"     Found {len(single_data)} single-phase files")
            print(f"     Found {len(multi_data)} multiphase files")
            return True
    
    # Run single-phase simulations
    for mode in ['standard', 'entropic']:
        success, _ = run_command(
            ['./lbm_sim', mode],
            f"Running {mode} BGK simulation",
            cwd=project_root,
            timeout=300  # 5 minutes
        )
        if not success:
            print(f"Single-phase {mode} simulation failed")
            return False
    
    # Run multiphase simulation
    success, _ = run_command(
        ['./lbm_multiphase'],
        "Running multiphase simulation",
        cwd=project_root,
        timeout=600  # 10 minutes
    )
    
    if not success:
        print("Multiphase simulation failed")
        return False
    
    return True

def generate_plots():
    """Generate standardized plots using the new plotting system."""
    project_root = Path.cwd()
    
    # Use the new modular plotting system
    success, output = run_command(
        [sys.executable, 'src/plot_generators.py'],
        "Generating standardized plots",
        cwd=project_root,
        timeout=120
    )
    
    if success:
        print("Generated plots:")
        print(output)
    
    return success

def generate_presentation():
    """Generate presentation using the report system."""
    project_root = Path.cwd()
    
    # Check if Node.js is available
    try:
        subprocess.run(['node', '--version'], capture_output=True, check=True)
    except:
        print("  ℹ️  Node.js not available, skipping presentation generation")
        return True
    
    # Install Node dependencies if needed
    package_json = project_root / 'package.json'
    if package_json.exists():
        success, _ = run_command(
            ['npm', 'install'],
            "Installing Node.js dependencies",
            cwd=project_root,
            timeout=120
        )
        if not success:
            print("  ⚠️  Could not install Node.js dependencies")
    
    # Generate presentation
    success, _ = run_command(
        [sys.executable, 'generate_report.py', '--format', 'all'],
        "Generating presentation (PDF, HTML, PPTX)",
        cwd=project_root,
        timeout=180
    )
    
    return success

def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(description='Complete LBM Analysis Workflow')
    parser.add_argument('--skip-sim', action='store_true', 
                       help='Skip simulations if data already exists')
    parser.add_argument('--plots-only', action='store_true',
                       help='Only generate plots from existing data')
    parser.add_argument('--presentation-only', action='store_true',
                       help='Only generate presentation from existing plots')
    
    args = parser.parse_args()
    
    print("🚀 Starting Complete LBM Analysis Workflow")
    print(f"Working directory: {Path.cwd()}")
    
    # Check dependencies
    if not check_dependencies():
        print("\n❌ Dependency check failed. Please install missing components.")
        return 1
    
    # Step 1: Compile simulations (unless we're only doing plots/presentation)
    if not args.plots_only and not args.presentation_only:
        if not compile_simulations():
            print("\n❌ Compilation failed")
            return 1
    
    # Step 2: Run simulations (unless we're only doing plots/presentation)
    if not args.plots_only and not args.presentation_only:
        if not run_simulations(skip_existing=args.skip_sim):
            print("\n❌ Simulations failed")
            return 1
    
    # Step 3: Generate plots (unless we're only doing presentation)
    if not args.presentation_only:
        if not generate_plots():
            print("\n❌ Plot generation failed")
            return 1
    
    # Step 4: Generate presentation
    if not generate_presentation():
        print("\n⚠️  Presentation generation failed (but plots are available)")
    
    print("\n🎉 Workflow completed successfully!")
    print("\nGenerated files:")
    
    # List generated files
    plots_dir = Path('plots')
    if plots_dir.exists():
        for plot_file in plots_dir.rglob('*.png'):
            print(f"  📊 {plot_file}")
    
    reports_dir = Path('reports/generated')
    if reports_dir.exists():
        for report_file in reports_dir.glob('*'):
            if report_file.is_file():
                print(f"  📑 {report_file}")
    
    print("\nNext steps:")
    print("  • Review plots in the plots/ directory")
    print("  • Check presentations in reports/generated/")
    print("  • Use --plots-only for quick plot updates")
    print("  • Use --presentation-only to regenerate presentations")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())