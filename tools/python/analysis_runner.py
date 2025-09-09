#!/usr/bin/env python3
"""
LBM Analysis Runner

Main script for running comprehensive LBM simulation analysis including
single-phase, multiphase, and H-theorem visualizations.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List
import matplotlib.pyplot as plt

# Add visualization modules to path
sys.path.append(str(Path(__file__).parent))

from visualization import (
    LayoutManager, LayoutConfig,
    SinglePhasePlots, MultiphasePlots, HTheoremPlots
)


def run_single_phase_analysis(data_directory: Path, mode: str = "standard",
                             output_directory: Optional[Path] = None,
                             show_plots: bool = False) -> None:
    """Run comprehensive single-phase analysis."""
    print(f"Running single-phase analysis for {mode} BGK...")
    
    # Create visualizer
    layout_config = LayoutConfig(dpi=300, font_size=11)
    visualizer = SinglePhasePlots(layout_config)
    
    try:
        # Load data
        data = visualizer.load_data_from_files(data_directory, mode)
        print(f"Loaded data for {len(data.timesteps)} timesteps")
        
        # Generate comprehensive analysis
        figures = visualizer.create_comprehensive_analysis(data, output_directory)
        print(f"Generated {len(figures)} single-phase visualization figures")
        
        if show_plots:
            plt.show()
            
    except Exception as e:
        print(f"Error in single-phase analysis: {e}")


def run_multiphase_analysis(data_directory: Path,
                           output_directory: Optional[Path] = None,
                           timesteps: Optional[List[int]] = None,
                           show_plots: bool = False) -> None:
    """Run comprehensive multiphase analysis."""
    print("Running multiphase analysis...")
    
    # Create visualizer
    layout_config = LayoutConfig(dpi=300, font_size=11)
    visualizer = MultiphasePlots(layout_config)
    
    try:
        # Load data
        data = visualizer.load_data_from_files(data_directory)
        print(f"Loaded multiphase data for {len(data.timesteps)} timesteps")
        print(f"Grid size: {len(data.x_positions)}x{len(data.y_positions)}")
        
        # Generate comprehensive analysis
        figures = visualizer.create_comprehensive_analysis(
            data, timesteps, output_directory
        )
        print(f"Generated {len(figures)} multiphase visualization figures")
        
        if show_plots:
            plt.show()
            
    except Exception as e:
        print(f"Error in multiphase analysis: {e}")


def run_h_theorem_analysis(h_evolution_file: Path,
                          output_directory: Optional[Path] = None,
                          comparison_file: Optional[Path] = None,
                          show_plots: bool = False) -> None:
    """Run H-theorem analysis."""
    print("Running H-theorem analysis...")
    
    # Create visualizer
    layout_config = LayoutConfig(dpi=300, font_size=11)
    visualizer = HTheoremPlots(layout_config)
    
    try:
        # Load data
        data = visualizer.load_h_theorem_data(h_evolution_file)
        print(f"Loaded H-theorem data for {len(data.timesteps)} timesteps")
        print(f"Monotonicity violations: {data.violation_count}")
        print(f"Average entropy production: {data.average_entropy_production:.2e}")
        
        if comparison_file:
            # Comparison analysis
            data2 = visualizer.load_h_theorem_data(comparison_file)
            fig_comp, _ = visualizer.create_bgk_comparison_plot(data, data2)
            
            if output_directory:
                filename = output_directory / "h_theorem_bgk_comparison.png"
                visualizer.layout_manager.save_figure(fig_comp, str(filename), dpi=300)
                print("Generated BGK comparison plot")
        else:
            # Single analysis
            figures = visualizer.create_comprehensive_h_analysis(data, output_directory)
            print(f"Generated {len(figures)} H-theorem visualization figures")
        
        if show_plots:
            plt.show()
            
    except Exception as e:
        print(f"Error in H-theorem analysis: {e}")


def auto_detect_analysis_type(data_directory: Path) -> str:
    """Automatically detect the type of simulation data."""
    # Check for multiphase files
    if list(data_directory.glob("multiphase_t*.csv")):
        return "multiphase"
    
    # Check for single-phase files
    if (list(data_directory.glob("standard_velocity_t*.csv")) or 
        list(data_directory.glob("entropic_velocity_t*.csv"))):
        return "single_phase"
    
    # Check for H-theorem files
    if list(data_directory.glob("*_h_evolution.csv")):
        return "h_theorem"
    
    return "unknown"


def main():
    """Main analysis runner."""
    parser = argparse.ArgumentParser(
        description='Comprehensive LBM Simulation Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect and analyze all data in directory
  python analysis_runner.py /path/to/data --output /path/to/plots

  # Single-phase analysis for standard BGK
  python analysis_runner.py /path/to/data --type single_phase --mode standard

  # Multiphase analysis for specific timesteps
  python analysis_runner.py /path/to/data --type multiphase --timesteps 0 5000 10000

  # H-theorem analysis with BGK comparison
  python analysis_runner.py /path/to/data --type h_theorem --h-file standard_h_evolution.csv --compare entropic_h_evolution.csv

  # Show plots interactively instead of saving
  python analysis_runner.py /path/to/data --show
        """
    )
    
    parser.add_argument('data_directory', type=Path,
                        help='Directory containing simulation data')
    parser.add_argument('--type', choices=['single_phase', 'multiphase', 'h_theorem', 'auto'],
                        default='auto', help='Type of analysis to perform')
    parser.add_argument('--output', type=Path,
                        help='Output directory for plots and reports')
    parser.add_argument('--show', action='store_true',
                        help='Show plots interactively instead of saving')
    
    # Single-phase specific options
    parser.add_argument('--mode', choices=['standard', 'entropic'], default='standard',
                        help='BGK collision operator mode for single-phase')
    
    # Multiphase specific options
    parser.add_argument('--timesteps', nargs='+', type=int,
                        help='Specific timesteps to analyze for multiphase')
    
    # H-theorem specific options
    parser.add_argument('--h-file', type=Path,
                        help='H-theorem evolution CSV file')
    parser.add_argument('--compare', type=Path,
                        help='Second H-evolution file for BGK comparison')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.data_directory.exists():
        print(f"Error: Data directory {args.data_directory} does not exist")
        sys.exit(1)
    
    # Create output directory if specified
    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {args.output}")
    
    # Auto-detect analysis type if requested
    if args.type == 'auto':
        detected_type = auto_detect_analysis_type(args.data_directory)
        if detected_type == 'unknown':
            print("Could not auto-detect analysis type. Please specify --type")
            sys.exit(1)
        args.type = detected_type
        print(f"Auto-detected analysis type: {args.type}")
    
    # Run appropriate analysis
    try:
        if args.type == 'single_phase':
            run_single_phase_analysis(
                args.data_directory, args.mode, args.output, args.show
            )
        
        elif args.type == 'multiphase':
            run_multiphase_analysis(
                args.data_directory, args.output, args.timesteps, args.show
            )
        
        elif args.type == 'h_theorem':
            h_file = args.h_file
            if not h_file:
                # Try to find H-evolution file automatically
                h_files = list(args.data_directory.glob("*_h_evolution.csv"))
                if h_files:
                    h_file = h_files[0]
                    print(f"Auto-detected H-evolution file: {h_file}")
                else:
                    print("Error: No H-evolution file found. Use --h-file to specify.")
                    sys.exit(1)
            
            run_h_theorem_analysis(h_file, args.output, args.compare, args.show)
        
        print("\nAnalysis completed successfully!")
        
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nAnalysis failed with error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()