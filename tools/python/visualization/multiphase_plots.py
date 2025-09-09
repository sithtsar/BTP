"""
Multiphase LBM Visualization Module

This module provides advanced visualization capabilities for multiphase LBM simulations
including 2D field visualizations, interface tracking, and phase distribution analysis.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from dataclasses import dataclass
from pathlib import Path
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from scipy import ndimage
from scipy.ndimage import binary_erosion, binary_dilation

from .layout_manager import LayoutManager, PlotConfig, PlotType, LayoutConfig


@dataclass
class MultiphaseData:
    """Data structure for multiphase simulation results."""
    timesteps: List[int]
    x_positions: np.ndarray
    y_positions: np.ndarray
    # Full field data: timestep -> field_data[x, y]
    phi_fields: Dict[int, np.ndarray]  # Phase field
    rho_fields: Dict[int, np.ndarray]  # Density field
    ux_fields: Dict[int, np.ndarray]   # X-velocity field
    uy_fields: Dict[int, np.ndarray]   # Y-velocity field
    pressure_fields: Optional[Dict[int, np.ndarray]] = None
    # Averaged profile data
    averaged_profiles: Optional[Dict[int, pd.DataFrame]] = None
    analytical_velocities: Optional[Dict[int, np.ndarray]] = None
    simulation_params: Optional[Dict[str, Any]] = None


@dataclass
class InterfaceMetrics:
    """Interface analysis metrics."""
    interface_thickness: float
    interface_position: float
    surface_area: float
    mass_light_phase: float
    mass_heavy_phase: float
    mass_conservation_error: float
    timestep: int


class MultiphasePlots:
    """Advanced multiphase visualization with 2D fields and interface tracking."""
    
    def __init__(self, layout_config: Optional[LayoutConfig] = None):
        self.layout_manager = LayoutManager(layout_config)
        self.phase_colormap = 'RdYlBu_r'  # Good for phase field visualization
        self.velocity_colormap = 'viridis'
        self.density_colormap = 'plasma'
        
    def load_data_from_files(self, data_directory: Union[str, Path]) -> MultiphaseData:
        """Load multiphase data from CSV files."""
        data_dir = Path(data_directory)
        
        # Find all full field data files
        field_files = list(data_dir.glob("multiphase_t*.csv"))
        field_files.sort(key=lambda x: int(x.stem.split('_t')[1]))
        
        timesteps = []
        phi_fields = {}
        rho_fields = {}
        ux_fields = {}
        uy_fields = {}
        pressure_fields = {}
        
        x_positions = None
        y_positions = None
        
        for file_path in field_files:
            # Extract timestep from filename
            timestep = int(file_path.stem.split('_t')[1])
            timesteps.append(timestep)
            
            # Load field data
            df = pd.read_csv(file_path)
            
            # Get grid dimensions
            if x_positions is None:
                x_unique = df['x'].unique()
                y_unique = df['y'].unique()
                x_positions = np.sort(x_unique)
                y_positions = np.sort(y_unique)
                nx, ny = len(x_unique), len(y_unique)
            
            # Reshape data into 2D fields
            phi_fields[timestep] = df['phi'].values.reshape(nx, ny)
            rho_fields[timestep] = df['rho'].values.reshape(nx, ny)
            ux_fields[timestep] = df['ux'].values.reshape(nx, ny)
            uy_fields[timestep] = df['uy'].values.reshape(nx, ny)
            
            if 'p' in df.columns:
                pressure_fields[timestep] = df['p'].values.reshape(nx, ny)
        
        # Load averaged profile data if available
        averaged_profiles = {}
        avg_files = list(data_dir.glob("multiphase_avg_t*.csv"))
        avg_files.sort(key=lambda x: int(x.stem.split('_t')[1]))
        
        for file_path in avg_files:
            timestep = int(file_path.stem.split('_t')[1])
            averaged_profiles[timestep] = pd.read_csv(file_path)
        
        return MultiphaseData(
            timesteps=timesteps,
            x_positions=x_positions,
            y_positions=y_positions,
            phi_fields=phi_fields,
            rho_fields=rho_fields,
            ux_fields=ux_fields,
            uy_fields=uy_fields,
            pressure_fields=pressure_fields if pressure_fields else None,
            averaged_profiles=averaged_profiles if averaged_profiles else None
        )
    
    def detect_interface(self, phi_field: np.ndarray, 
                        threshold: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
        """Detect interface location from phase field."""
        # Find contour at phi = threshold
        contours = plt.contour(phi_field, levels=[threshold])
        plt.close()  # Close the temporary figure
        
        if len(contours.collections) > 0 and len(contours.collections[0].get_paths()) > 0:
            path = contours.collections[0].get_paths()[0]
            vertices = path.vertices
            return vertices[:, 0], vertices[:, 1]
        else:
            return np.array([]), np.array([])
    
    def calculate_interface_metrics(self, phi_field: np.ndarray, 
                                  rho_field: np.ndarray, timestep: int) -> InterfaceMetrics:
        """Calculate comprehensive interface analysis metrics."""
        # Interface thickness (gradient magnitude)
        grad_x, grad_y = np.gradient(phi_field)
        grad_magnitude = np.sqrt(grad_x**2 + grad_y**2)
        interface_thickness = np.mean(grad_magnitude[grad_magnitude > 0.01])
        
        # Interface position (center of mass of gradient)
        interface_mask = grad_magnitude > 0.1 * np.max(grad_magnitude)
        if np.any(interface_mask):
            y_coords, x_coords = np.where(interface_mask)
            interface_position = np.mean(y_coords)
        else:
            interface_position = phi_field.shape[1] / 2
        
        # Surface area (length of interface in 2D)
        interface_x, interface_y = self.detect_interface(phi_field)
        if len(interface_x) > 1:
            dx = np.diff(interface_x)
            dy = np.diff(interface_y) 
            surface_area = np.sum(np.sqrt(dx**2 + dy**2))
        else:
            surface_area = 0.0
        
        # Mass conservation
        light_mask = phi_field < 0.5
        heavy_mask = phi_field >= 0.5
        
        mass_light_phase = np.sum(rho_field[light_mask])
        mass_heavy_phase = np.sum(rho_field[heavy_mask])
        total_mass = mass_light_phase + mass_heavy_phase
        
        # Assuming initial total mass is known (would need to be passed in practice)
        expected_total_mass = total_mass  # Placeholder
        mass_conservation_error = abs(total_mass - expected_total_mass) / expected_total_mass
        
        return InterfaceMetrics(
            interface_thickness=interface_thickness,
            interface_position=interface_position,
            surface_area=surface_area,
            mass_light_phase=mass_light_phase,
            mass_heavy_phase=mass_heavy_phase,
            mass_conservation_error=mass_conservation_error,
            timestep=timestep
        )
    
    def create_field_visualization(self, data: MultiphaseData, timestep: int,
                                 show_interface: bool = True) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create 2D field visualization with proper colormaps and contours."""
        if timestep not in data.phi_fields:
            raise ValueError(f"No data available for timestep {timestep}")
        
        plot_configs = [
            PlotConfig(PlotType.PHASE_DISTRIBUTION, "Phase Field (φ)", 
                      "X Position", "Y Position", colorbar=True),
            PlotConfig(PlotType.FIELD_2D, "Density Field (ρ)", 
                      "X Position", "Y Position", colorbar=True),
            PlotConfig(PlotType.FIELD_2D, "X-Velocity Field (ux)", 
                      "X Position", "Y Position", colorbar=True),
            PlotConfig(PlotType.FIELD_2D, "Y-Velocity Field (uy)", 
                      "X Position", "Y Position", colorbar=True)
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, layout_shape=(2, 2),
            figure_title=f"Multiphase Fields at t={timestep}"
        )
        
        # Get field data
        phi = data.phi_fields[timestep]
        rho = data.rho_fields[timestep]
        ux = data.ux_fields[timestep]
        uy = data.uy_fields[timestep]
        
        # Create meshgrid for plotting
        X, Y = np.meshgrid(data.x_positions, data.y_positions, indexing='ij')
        
        # Phase field plot
        im1 = axes[0].contourf(X, Y, phi, levels=20, cmap=self.phase_colormap)
        self.layout_manager.style_manager.create_colorbar(axes[0], im1, "φ")
        if show_interface:
            axes[0].contour(X, Y, phi, levels=[0.5], colors='black', linewidths=2)
        
        # Density field plot
        im2 = axes[1].contourf(X, Y, rho, levels=20, cmap=self.density_colormap)
        self.layout_manager.style_manager.create_colorbar(axes[1], im2, "ρ")
        
        # X-velocity field plot
        im3 = axes[2].contourf(X, Y, ux, levels=20, cmap=self.velocity_colormap)
        self.layout_manager.style_manager.create_colorbar(axes[2], im3, "ux")
        
        # Y-velocity field plot
        im4 = axes[3].contourf(X, Y, uy, levels=20, cmap=self.velocity_colormap)
        self.layout_manager.style_manager.create_colorbar(axes[3], im4, "uy")
        
        return fig, axes
    
    def create_streamline_visualization(self, data: MultiphaseData, timestep: int,
                                      density: int = 2) -> Tuple[plt.Figure, plt.Axes]:
        """Create velocity field visualization with streamlines."""
        if timestep not in data.ux_fields:
            raise ValueError(f"No velocity data available for timestep {timestep}")
        
        config = PlotConfig(
            plot_type=PlotType.STREAMLINES,
            title=f"Velocity Streamlines at t={timestep}",
            xlabel="X Position",
            ylabel="Y Position",
            colorbar=True
        )
        
        fig, axes = self.layout_manager.create_multi_panel_layout([config])
        ax = axes[0]
        
        # Get velocity fields
        ux = data.ux_fields[timestep]
        uy = data.uy_fields[timestep]
        phi = data.phi_fields[timestep]
        
        # Create meshgrid
        X, Y = np.meshgrid(data.x_positions, data.y_positions, indexing='ij')
        
        # Plot phase field as background
        im = ax.contourf(X, Y, phi, levels=20, cmap=self.phase_colormap, alpha=0.6)
        
        # Calculate velocity magnitude for coloring streamlines
        velocity_magnitude = np.sqrt(ux**2 + uy**2)
        
        # Create streamlines
        strm = ax.streamplot(X, Y, ux, uy, 
                           color=velocity_magnitude, 
                           cmap='viridis',
                           density=density,
                           linewidth=1.5,
                           arrowsize=1.5)
        
        # Add interface contour
        ax.contour(X, Y, phi, levels=[0.5], colors='black', linewidths=2)
        
        # Create colorbar for velocity magnitude
        self.layout_manager.style_manager.create_colorbar(ax, strm.lines, "|V|")
        
        return fig, ax
    
    def create_interface_evolution_plot(self, data: MultiphaseData,
                                      timesteps_to_analyze: Optional[List[int]] = None) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create interface tracking and evolution analysis."""
        if timesteps_to_analyze is None:
            timesteps_to_analyze = data.timesteps[::max(1, len(data.timesteps)//10)]
        
        # Calculate interface metrics
        interface_metrics = []
        for timestep in timesteps_to_analyze:
            if timestep in data.phi_fields:
                metrics = self.calculate_interface_metrics(
                    data.phi_fields[timestep],
                    data.rho_fields[timestep],
                    timestep
                )
                interface_metrics.append(metrics)
        
        plot_configs = [
            PlotConfig(PlotType.PHASE_DISTRIBUTION, "Interface Position Evolution", 
                      "Timestep", "Interface Position"),
            PlotConfig(PlotType.STABILITY_INDICATORS, "Interface Thickness", 
                      "Timestep", "Thickness"),
            PlotConfig(PlotType.ERROR_ANALYSIS, "Mass Conservation", 
                      "Timestep", "Mass Conservation Error")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title="Interface Evolution Analysis"
        )
        
        # Extract data for plotting
        timesteps = [m.timestep for m in interface_metrics]
        positions = [m.interface_position for m in interface_metrics]
        thicknesses = [m.interface_thickness for m in interface_metrics]
        mass_errors = [m.mass_conservation_error for m in interface_metrics]
        
        # Interface position evolution
        axes[0].plot(timesteps, positions, 'o-', linewidth=2, markersize=5)
        axes[0].set_ylabel("Interface Position")
        
        # Interface thickness evolution
        axes[1].plot(timesteps, thicknesses, 's-', color='red', linewidth=2, markersize=5)
        axes[1].set_ylabel("Interface Thickness")
        
        # Mass conservation
        axes[2].semilogy(timesteps, mass_errors, '^-', color='green', linewidth=2, markersize=5)
        axes[2].set_ylabel("Mass Conservation Error")
        axes[2].set_ylim(bottom=1e-12)
        
        return fig, axes
    
    def create_profile_comparison_plot(self, data: MultiphaseData, timestep: int) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create y-averaged profile comparison with analytical solution."""
        if not data.averaged_profiles or timestep not in data.averaged_profiles:
            raise ValueError("Averaged profile data not available")
        
        df = data.averaged_profiles[timestep]
        
        plot_configs = [
            PlotConfig(PlotType.PHASE_DISTRIBUTION, "Phase Distribution Profile", 
                      "Channel Height (y)", "Average φ"),
            PlotConfig(PlotType.VELOCITY_PROFILE, "Velocity Profile", 
                      "Channel Height (y)", "Average ux"),
            PlotConfig(PlotType.FIELD_2D, "Density Profile", 
                      "Channel Height (y)", "Average ρ")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title=f"Averaged Profiles at t={timestep}"
        )
        
        # Phase distribution profile
        axes[0].plot(df['y'], df['avg_phi'], 'o-', linewidth=2, markersize=4, label='Numerical')
        axes[0].axhline(y=0.5, color='red', linestyle='--', alpha=0.7, label='Interface')
        axes[0].legend()
        axes[0].set_ylim(0, 1)
        
        # Velocity profile
        axes[1].plot(df['y'], df['avg_ux'], 'o-', linewidth=2, markersize=4, label='Numerical')
        if 'analytical_ux' in df.columns:
            axes[1].plot(df['y'], df['analytical_ux'], '--', linewidth=2, 
                        alpha=0.8, label='Analytical')
        axes[1].legend()
        
        # Density profile
        axes[2].plot(df['y'], df['avg_rho'], 'o-', linewidth=2, markersize=4, color='purple')
        
        return fig, axes
    
    def create_comprehensive_analysis(self, data: MultiphaseData,
                                    timesteps_to_plot: Optional[List[int]] = None,
                                    output_directory: Optional[Union[str, Path]] = None) -> Dict[str, plt.Figure]:
        """Create comprehensive multiphase analysis with all visualizations."""
        if timesteps_to_plot is None:
            # Select representative timesteps
            n_timesteps = len(data.timesteps)
            indices = [0, n_timesteps//4, n_timesteps//2, 3*n_timesteps//4, -1]
            timesteps_to_plot = [data.timesteps[i] for i in indices if i < n_timesteps]
        
        figures = {}
        
        # Field visualizations for each timestep
        for timestep in timesteps_to_plot:
            fig_fields, _ = self.create_field_visualization(data, timestep)
            figures[f'fields_t{timestep}'] = fig_fields
            
            # Streamline visualization
            fig_stream, _ = self.create_streamline_visualization(data, timestep)
            figures[f'streamlines_t{timestep}'] = fig_stream
            
            # Profile comparison (if averaged data available)
            if data.averaged_profiles and timestep in data.averaged_profiles:
                fig_profiles, _ = self.create_profile_comparison_plot(data, timestep)
                figures[f'profiles_t{timestep}'] = fig_profiles
        
        # Interface evolution analysis
        fig_interface, _ = self.create_interface_evolution_plot(data)
        figures['interface_evolution'] = fig_interface
        
        # Save figures if output directory specified
        if output_directory:
            output_dir = Path(output_directory)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            for name, fig in figures.items():
                filename = output_dir / f"multiphase_{name}.png"
                self.layout_manager.save_figure(fig, str(filename), dpi=300)
        
        return figures
    
    def create_phase_evolution_animation_data(self, data: MultiphaseData,
                                            timesteps: Optional[List[int]] = None) -> Dict[str, Any]:
        """Prepare data for phase field evolution animation."""
        if timesteps is None:
            timesteps = data.timesteps
        
        animation_data = {
            'timesteps': timesteps,
            'x_positions': data.x_positions,
            'y_positions': data.y_positions,
            'phi_fields': [data.phi_fields[t] for t in timesteps if t in data.phi_fields],
            'velocity_fields': [(data.ux_fields[t], data.uy_fields[t]) 
                              for t in timesteps if t in data.ux_fields]
        }
        
        return animation_data


def main():
    """Example usage of multiphase visualization."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate multiphase LBM visualizations')
    parser.add_argument('data_directory', help='Directory containing simulation data')
    parser.add_argument('--timesteps', nargs='+', type=int, 
                       help='Specific timesteps to visualize')
    parser.add_argument('--output', help='Output directory for plots')
    
    args = parser.parse_args()
    
    # Create visualizer
    visualizer = MultiphasePlots()
    
    # Load data
    data = visualizer.load_data_from_files(args.data_directory)
    
    # Generate comprehensive analysis
    figures = visualizer.create_comprehensive_analysis(
        data, args.timesteps, args.output
    )
    
    print(f"Generated {len(figures)} visualization figures")
    
    if not args.output:
        plt.show()


if __name__ == "__main__":
    main()