"""
Single-Phase LBM Visualization Module

This module provides comprehensive visualization capabilities for single-phase
LBM simulations including velocity profiles, H-function evolution, and error analysis.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from dataclasses import dataclass
from pathlib import Path
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

from .layout_manager import LayoutManager, PlotConfig, PlotType, LayoutConfig


@dataclass
class SinglePhaseData:
    """Data structure for single-phase simulation results."""
    timesteps: List[int]
    y_positions: np.ndarray
    velocities: Dict[int, np.ndarray]  # timestep -> velocity profile
    analytical_velocities: Optional[Dict[int, np.ndarray]] = None
    h_function_data: Optional[Dict[int, np.ndarray]] = None  # timestep -> H-function profile
    h_evolution_data: Optional[pd.DataFrame] = None  # Full H-theorem evolution
    simulation_params: Optional[Dict[str, Any]] = None
    mode: str = "standard"  # BGK mode (standard or entropic)


@dataclass
class ErrorMetrics:
    """Error analysis metrics."""
    l1_error: float
    l2_error: float
    max_error: float
    relative_error: float
    timestep: int


class SinglePhasePlots:
    """Comprehensive single-phase visualization with error analysis and H-function evolution."""
    
    def __init__(self, layout_config: Optional[LayoutConfig] = None):
        self.layout_manager = LayoutManager(layout_config)
        self.default_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
        
    def load_data_from_files(self, data_directory: Union[str, Path], 
                           mode: str = "standard") -> SinglePhaseData:
        """Load single-phase data from CSV files."""
        data_dir = Path(data_directory)
        
        # Find all velocity files
        velocity_files = list(data_dir.glob(f"{mode}_velocity_t*.csv"))
        velocity_files.sort(key=lambda x: int(x.stem.split('_t')[1]))
        
        timesteps = []
        velocities = {}
        analytical_velocities = {}
        
        for file_path in velocity_files:
            # Extract timestep from filename
            timestep = int(file_path.stem.split('_t')[1])
            timesteps.append(timestep)
            
            # Load velocity data
            df = pd.read_csv(file_path)
            y_positions = df['y'].values
            velocities[timestep] = df['avg_ux'].values
            
            if 'analytic_ux' in df.columns:
                analytical_velocities[timestep] = df['analytic_ux'].values
        
        # Load H-function data if available
        h_function_data = {}
        h_files = list(data_dir.glob(f"{mode}_h_y_t*.csv"))
        h_files.sort(key=lambda x: int(x.stem.split('_t')[1]))
        
        for file_path in h_files:
            timestep = int(file_path.stem.split('_t')[1])
            df = pd.read_csv(file_path)
            h_function_data[timestep] = df['H_y'].values
        
        # Load H-theorem evolution data if available
        h_evolution_data = None
        h_evolution_file = data_dir / f"{mode}_h_evolution.csv"
        if h_evolution_file.exists():
            h_evolution_data = pd.read_csv(h_evolution_file)
        
        return SinglePhaseData(
            timesteps=timesteps,
            y_positions=y_positions,
            velocities=velocities,
            analytical_velocities=analytical_velocities if analytical_velocities else None,
            h_function_data=h_function_data if h_function_data else None,
            h_evolution_data=h_evolution_data,
            mode=mode
        )
    
    def calculate_error_metrics(self, numerical: np.ndarray, 
                              analytical: np.ndarray, timestep: int) -> ErrorMetrics:
        """Calculate comprehensive error metrics."""
        error = numerical - analytical
        
        l1_error = np.mean(np.abs(error))
        l2_error = np.sqrt(np.mean(error**2))
        max_error = np.max(np.abs(error))
        
        # Relative error (avoid division by zero)
        analytical_max = np.max(np.abs(analytical))
        relative_error = l2_error / analytical_max if analytical_max > 1e-12 else 0.0
        
        return ErrorMetrics(
            l1_error=l1_error,
            l2_error=l2_error,
            max_error=max_error,
            relative_error=relative_error,
            timestep=timestep
        )
    
    def create_velocity_profile_plot(self, data: SinglePhaseData, 
                                   timesteps_to_plot: Optional[List[int]] = None,
                                   show_analytical: bool = True) -> Tuple[plt.Figure, plt.Axes]:
        """Create velocity profile visualization with analytical comparison."""
        if timesteps_to_plot is None:
            timesteps_to_plot = data.timesteps[-3:]  # Plot last 3 timesteps
        
        config = PlotConfig(
            plot_type=PlotType.VELOCITY_PROFILE,
            title="Velocity Profile Evolution",
            xlabel="Channel Height (y)",
            ylabel="Velocity (ux)",
            grid=True
        )
        
        fig, axes = self.layout_manager.create_multi_panel_layout([config])
        ax = axes[0]
        
        # Plot numerical results
        for i, timestep in enumerate(timesteps_to_plot):
            if timestep in data.velocities:
                color = self.default_colors[i % len(self.default_colors)]
                ax.plot(data.y_positions, data.velocities[timestep], 
                       'o-', color=color, markersize=4, linewidth=1.5,
                       label=f'Numerical t={timestep}')
                
                # Plot analytical solution if available
                if (show_analytical and data.analytical_velocities and 
                    timestep in data.analytical_velocities):
                    ax.plot(data.y_positions, data.analytical_velocities[timestep],
                           '--', color=color, linewidth=2, alpha=0.8,
                           label=f'Analytical t={timestep}')
        
        ax.legend(frameon=True, fancybox=True, shadow=True)
        ax.set_xlim(0, len(data.y_positions) - 1)
        
        return fig, ax
    
    def create_error_analysis_plot(self, data: SinglePhaseData) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create comprehensive error analysis visualization."""
        if not data.analytical_velocities:
            raise ValueError("Analytical velocities required for error analysis")
        
        # Calculate error metrics for all timesteps
        error_metrics = []
        for timestep in data.timesteps:
            if timestep in data.velocities and timestep in data.analytical_velocities:
                metrics = self.calculate_error_metrics(
                    data.velocities[timestep],
                    data.analytical_velocities[timestep],
                    timestep
                )
                error_metrics.append(metrics)
        
        # Create multi-panel layout
        plot_configs = [
            PlotConfig(PlotType.ERROR_ANALYSIS, "Error Evolution", 
                      "Timestep", "L2 Error"),
            PlotConfig(PlotType.ERROR_ANALYSIS, "Relative Error", 
                      "Timestep", "Relative Error (%)"),
            PlotConfig(PlotType.ERROR_ANALYSIS, "Spatial Error Distribution", 
                      "Channel Height (y)", "Error")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title="Error Analysis"
        )
        
        # Error evolution plot
        timesteps = [m.timestep for m in error_metrics]
        l1_errors = [m.l1_error for m in error_metrics]
        l2_errors = [m.l2_error for m in error_metrics]
        max_errors = [m.max_error for m in error_metrics]
        
        axes[0].semilogy(timesteps, l1_errors, 'o-', label='L1 Error', linewidth=2)
        axes[0].semilogy(timesteps, l2_errors, 's-', label='L2 Error', linewidth=2)
        axes[0].semilogy(timesteps, max_errors, '^-', label='Max Error', linewidth=2)
        axes[0].legend()
        axes[0].set_ylim(bottom=1e-10)
        
        # Relative error plot
        relative_errors = [m.relative_error * 100 for m in error_metrics]
        axes[1].plot(timesteps, relative_errors, 'o-', color='red', linewidth=2)
        axes[1].set_ylabel("Relative Error (%)")
        
        # Spatial error distribution for last timestep
        if error_metrics:
            last_timestep = error_metrics[-1].timestep
            numerical = data.velocities[last_timestep]
            analytical = data.analytical_velocities[last_timestep]
            error = numerical - analytical
            
            axes[2].plot(data.y_positions, error, 'o-', color='green', linewidth=2)
            axes[2].axhline(y=0, color='black', linestyle='--', alpha=0.5)
            axes[2].fill_between(data.y_positions, error, alpha=0.3, color='green')
        
        return fig, axes
    
    def create_h_function_evolution_plot(self, data: SinglePhaseData) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create H-function evolution visualization."""
        if not data.h_evolution_data is not None:
            raise ValueError("H-theorem evolution data required")
        
        plot_configs = [
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "H-Function Evolution", 
                      "Timestep", "H-Function Value"),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "H-Function Change Rate", 
                      "Timestep", "dH/dt"),
            PlotConfig(PlotType.STABILITY_INDICATORS, "Monotonicity Violations", 
                      "Timestep", "Violation Count")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title="H-Theorem Analysis"
        )
        
        df = data.h_evolution_data
        
        # H-function evolution
        axes[0].plot(df['timestep'], df['h_value'], 'b-', linewidth=2, label='H-Function')
        axes[0].set_yscale('log')
        axes[0].legend()
        
        # H-function change rate
        if 'h_change' in df.columns:
            axes[1].plot(df['timestep'], df['h_change'], 'r-', linewidth=2)
            axes[1].axhline(y=0, color='black', linestyle='--', alpha=0.5)
        
        # Monotonicity violations
        if 'monotonic_violation' in df.columns:
            violations = df['monotonic_violation'].astype(int)
            violation_positions = df.loc[violations == 1, 'timestep']
            if len(violation_positions) > 0:
                axes[2].scatter(violation_positions, 
                               np.ones(len(violation_positions)), 
                               color='red', s=50, marker='x', 
                               label='Violations')
                axes[2].set_ylim(0, 2)
                axes[2].legend()
            else:
                axes[2].text(0.5, 0.5, 'No Violations', 
                           transform=axes[2].transAxes, 
                           ha='center', va='center', fontsize=12)
        
        return fig, axes
    
    def create_comprehensive_analysis(self, data: SinglePhaseData, 
                                    output_directory: Optional[Union[str, Path]] = None) -> Dict[str, plt.Figure]:
        """Create comprehensive analysis with all visualizations."""
        figures = {}
        
        # Velocity profile plot
        fig_velocity, _ = self.create_velocity_profile_plot(data)
        figures['velocity_profiles'] = fig_velocity
        
        # Error analysis (if analytical data available)
        if data.analytical_velocities:
            fig_error, _ = self.create_error_analysis_plot(data)
            figures['error_analysis'] = fig_error
        
        # H-function evolution (if H-theorem data available)
        if data.h_evolution_data is not None:
            fig_h, _ = self.create_h_function_evolution_plot(data)
            figures['h_function_evolution'] = fig_h
        
        # Save figures if output directory specified
        if output_directory:
            output_dir = Path(output_directory)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Get mode from data for unique filenames
            mode_suffix = f"_{data.mode}" if hasattr(data, 'mode') and data.mode != 'standard' else ""
            
            for name, fig in figures.items():
                filename = output_dir / f"single_phase_{name}{mode_suffix}.png"
                self.layout_manager.save_figure(fig, str(filename), dpi=300)
        
        return figures
    
    def create_comparison_plot(self, standard_data: SinglePhaseData, 
                             entropic_data: SinglePhaseData,
                             timestep: int) -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
        """Create side-by-side comparison of standard vs entropic BGK."""
        left_config = PlotConfig(
            PlotType.VELOCITY_PROFILE,
            title="Standard BGK",
            xlabel="Channel Height (y)",
            ylabel="Velocity (ux)"
        )
        
        right_config = PlotConfig(
            PlotType.VELOCITY_PROFILE,
            title="Entropic BGK", 
            xlabel="Channel Height (y)",
            ylabel="Velocity (ux)"
        )
        
        fig, (ax_left, ax_right) = self.layout_manager.create_comparison_layout(
            left_config, right_config, 
            figure_title=f"BGK Comparison at t={timestep}"
        )
        
        # Plot standard BGK
        if timestep in standard_data.velocities:
            ax_left.plot(standard_data.y_positions, standard_data.velocities[timestep],
                        'o-', color='blue', label='Numerical')
            if (standard_data.analytical_velocities and 
                timestep in standard_data.analytical_velocities):
                ax_left.plot(standard_data.y_positions, 
                           standard_data.analytical_velocities[timestep],
                           '--', color='blue', alpha=0.7, label='Analytical')
            ax_left.legend()
        
        # Plot entropic BGK
        if timestep in entropic_data.velocities:
            ax_right.plot(entropic_data.y_positions, entropic_data.velocities[timestep],
                         'o-', color='red', label='Numerical')
            if (entropic_data.analytical_velocities and 
                timestep in entropic_data.analytical_velocities):
                ax_right.plot(entropic_data.y_positions,
                            entropic_data.analytical_velocities[timestep],
                            '--', color='red', alpha=0.7, label='Analytical')
            ax_right.legend()
        
        return fig, (ax_left, ax_right)


def main():
    """Example usage of single-phase visualization."""
    # This would typically be called from a separate script
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate single-phase LBM visualizations')
    parser.add_argument('data_directory', help='Directory containing simulation data')
    parser.add_argument('--mode', default='standard', choices=['standard', 'entropic'],
                       help='BGK collision operator mode')
    parser.add_argument('--output', help='Output directory for plots')
    
    args = parser.parse_args()
    
    # Create visualizer
    visualizer = SinglePhasePlots()
    
    # Load data
    data = visualizer.load_data_from_files(args.data_directory, args.mode)
    
    # Generate comprehensive analysis
    figures = visualizer.create_comprehensive_analysis(data, args.output)
    
    print(f"Generated {len(figures)} visualization figures")
    
    if not args.output:
        plt.show()


if __name__ == "__main__":
    main()