#!/usr/bin/env python3
"""
Modular LBM plotting system with standardized outputs for presentations.
Creates exactly 2 plots: velocity_analysis.png and multiphase_analysis.png
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configure matplotlib for high-quality output
plt.rcParams.update({
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'font.size': 12,
    'axes.linewidth': 1.2,
    'grid.alpha': 0.3,
    'lines.linewidth': 2,
    'font.family': 'sans-serif'
})

class DataLoader:
    """Handles safe loading and validation of simulation data."""
    
    def __init__(self, data_dir="data"):
        self.data_dir = Path(data_dir)
        
    def load_csv_safe(self, filepath):
        """Safely load CSV with error handling and NaN replacement."""
        try:
            data = pd.read_csv(filepath)
            # Replace infinite values with NaN
            data = data.replace([np.inf, -np.inf], np.nan)
            return data
        except Exception as e:
            print(f"Warning: Could not load {filepath}: {e}")
            return None
    
    def get_time_series(self, pattern):
        """Get time series data from files matching pattern."""
        files = sorted(glob.glob(str(self.data_dir / pattern)))
        times = []
        data_list = []
        
        for file in files:
            try:
                # Extract time from filename
                t = int(Path(file).stem.split('_t')[1])
                data = self.load_csv_safe(file)
                if data is not None and len(data) > 0:
                    times.append(t)
                    data_list.append(data)
            except Exception as e:
                print(f"Warning: Could not process {file}: {e}")
                continue
        
        return times, data_list
    
    def validate_data_quality(self, data, required_columns):
        """Check if data has required columns and valid values."""
        if data is None:
            return False, "Data is None"
        
        missing_cols = [col for col in required_columns if col not in data.columns]
        if missing_cols:
            return False, f"Missing columns: {missing_cols}"
        
        # Check for all-NaN columns
        nan_cols = []
        for col in required_columns:
            if data[col].isna().all():
                nan_cols.append(col)
        
        if nan_cols:
            return False, f"All-NaN columns: {nan_cols}"
        
        return True, "Data is valid"

class SinglePhaseAnalyzer:
    """Analyzes single-phase BGK simulation results."""
    
    def __init__(self, data_loader):
        self.loader = data_loader
        
    def analyze_convergence(self, mode):
        """Analyze convergence for a single BGK mode."""
        times, data_list = self.loader.get_time_series(f"{mode}_velocity_t*.csv")
        
        if not times:
            return None
        
        errors = []
        profiles = []
        
        for data in data_list:
            is_valid, msg = self.loader.validate_data_quality(
                data, ['y', 'avg_ux', 'analytic_ux']
            )
            
            if not is_valid:
                continue
            
            # Handle different column name variations
            if 'analytical_ux' in data.columns:
                analytical = data['analytical_ux'].values
            else:
                analytical = data['analytic_ux'].values
            
            simulated = data['avg_ux'].values
            y = data['y'].values
            
            # Calculate RMSE
            valid_idx = (~np.isnan(simulated)) & (~np.isnan(analytical))
            if np.any(valid_idx):
                error = np.sqrt(np.mean((simulated[valid_idx] - analytical[valid_idx])**2))
                errors.append(error)
                profiles.append((y, simulated, analytical))
        
        if not errors:
            return None
        
        return {
            'times': times[:len(errors)],
            'errors': errors,
            'profiles': profiles,
            'final_rmse': errors[-1] if errors else np.nan,
            'convergence_rate': (errors[0] - errors[-1]) / len(errors) if len(errors) > 1 else 0
        }

class MultiphaseAnalyzer:
    """Analyzes multiphase simulation results."""
    
    def __init__(self, data_loader):
        self.loader = data_loader
        
    def analyze_quality(self):
        """Analyze multiphase data quality and physics."""
        times, data_list = self.loader.get_time_series("multiphase_avg_t*.csv")
        
        if not times:
            return None
        
        valid_count = 0
        interface_positions = []
        mass_history = []
        
        for t, data in zip(times, data_list):
            is_valid, msg = self.loader.validate_data_quality(
                data, ['y', 'avg_phi', 'avg_rho', 'avg_ux']
            )
            
            if is_valid:
                valid_count += 1
                
                # Interface analysis
                phi = data['avg_phi'].values
                y = data['y'].values
                
                valid_phi = ~np.isnan(phi)
                if np.any(valid_phi):
                    # Find interface (phi ≈ 0.5)
                    phi_centered = np.abs(phi[valid_phi] - 0.5)
                    if len(phi_centered) > 0:
                        interface_idx = np.argmin(phi_centered)
                        interface_positions.append(y[valid_phi][interface_idx])
                
                # Mass conservation
                rho = data['avg_rho'].values
                valid_rho = ~np.isnan(rho)
                if np.any(valid_rho):
                    total_mass = np.mean(rho[valid_rho])  # Simplified mass
                    mass_history.append(total_mass)
        
        # Calculate metrics
        interface_stability = np.std(interface_positions) if len(interface_positions) > 1 else 0
        interface_drift = (interface_positions[-1] - interface_positions[0]) if len(interface_positions) > 1 else 0
        
        mass_conservation = 0
        if len(mass_history) > 1:
            mass_conservation = abs((mass_history[-1] - mass_history[0]) / mass_history[0] * 100)
        
        return {
            'total_files': len(data_list),
            'valid_files': valid_count,
            'interface_stability': interface_stability,
            'interface_drift': interface_drift,
            'mass_conservation': mass_conservation,
            'data_quality': valid_count / len(data_list) * 100 if data_list else 0
        }

class PlotGenerator:
    """Generates standardized plots for presentations."""
    
    def __init__(self, output_dir="plots"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
    def create_velocity_analysis(self, single_phase_results):
        """Create standardized velocity analysis plot."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Single-Phase LBM Velocity Analysis', fontsize=16, fontweight='bold')
        
        if not single_phase_results:
            # Create error plot if no data
            for ax in axes.flat:
                ax.text(0.5, 0.5, 'No Valid Data\nCheck Simulation', 
                       ha='center', va='center', transform=ax.transAxes,
                       bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
                ax.axis('off')
        else:
            # Plot 1: Final velocity profiles
            ax = axes[0, 0]
            colors = ['#1f77b4', '#ff7f0e']
            for i, (mode, result) in enumerate(single_phase_results.items()):
                if result and result['profiles']:
                    y, sim, ana = result['profiles'][-1]
                    ax.plot(sim, y, 'o-', color=colors[i], label=f'{mode.title()} BGK', markersize=4)
                    if i == 0:  # Show analytical only once
                        ax.plot(ana, y, 'k:', linewidth=3, label='Analytical', alpha=0.8)
            
            ax.set_xlabel('Velocity $u_x$')
            ax.set_ylabel('Channel height $y$')
            ax.set_title('Final Velocity Profiles')
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Plot 2: Error convergence
            ax = axes[0, 1]
            for i, (mode, result) in enumerate(single_phase_results.items()):
                if result and result['errors']:
                    ax.semilogy(result['times'], result['errors'], 'o-', 
                               color=colors[i], label=f'{mode.title()} BGK', markersize=4)
            
            ax.set_xlabel('Time Steps')
            ax.set_ylabel('RMSE (log scale)')
            ax.set_title('Error Convergence')
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Plot 3: Performance comparison
            ax = axes[1, 0]
            if len(single_phase_results) == 2:
                methods = list(single_phase_results.keys())
                final_errors = [single_phase_results[m]['final_rmse'] for m in methods]
                
                bars = ax.bar(methods, final_errors, color=colors[:len(methods)], alpha=0.7)
                ax.set_ylabel('Final RMSE')
                ax.set_title('Method Comparison')
                ax.grid(True, axis='y', alpha=0.3)
                
                # Add value labels
                for bar, error in zip(bars, final_errors):
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                           f'{error:.2e}', ha='center', va='bottom', fontsize=10)
            
            # Plot 4: Summary metrics
            ax = axes[1, 1]
            ax.axis('off')
            
            summary_text = "Analysis Summary\n\n"
            for mode, result in single_phase_results.items():
                if result:
                    summary_text += f"{mode.title()} BGK:\n"
                    summary_text += f"• Final RMSE: {result['final_rmse']:.2e}\n"
                    summary_text += f"• Convergence rate: {result['convergence_rate']:.2e}\n\n"
            
            ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
                   verticalalignment='top', fontsize=11,
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        output_file = self.output_dir / "velocity_analysis.png"
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        
        return str(output_file)
    
    def create_multiphase_analysis(self, multiphase_result):
        """Create standardized multiphase analysis plot."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Multiphase LBM Analysis', fontsize=16, fontweight='bold')
        
        if not multiphase_result:
            # Create error plot if no data
            for ax in axes.flat:
                ax.text(0.5, 0.5, 'No Valid Multiphase Data\nRun Simulation First', 
                       ha='center', va='center', transform=ax.transAxes,
                       bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
                ax.axis('off')
        else:
            # Plot 1: Data quality metrics
            ax = axes[0, 0]
            metrics = ['Data Quality', 'Interface Stability', 'Mass Conservation']
            values = [
                multiphase_result['data_quality'],
                100 - multiphase_result['interface_stability'] * 100,  # Convert to percentage
                100 - multiphase_result['mass_conservation']  # Good if low change
            ]
            
            colors = ['green' if v > 80 else 'orange' if v > 50 else 'red' for v in values]
            bars = ax.bar(metrics, values, color=colors, alpha=0.7)
            
            ax.set_ylabel('Quality Score (%)')
            ax.set_title('Simulation Quality Metrics')
            ax.set_ylim(0, 100)
            ax.grid(True, axis='y', alpha=0.3)
            
            for bar, value in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                       f'{value:.1f}%', ha='center', va='bottom', fontsize=10)
            
            # Plot 2: Interface analysis
            ax = axes[0, 1]
            ax.text(0.5, 0.5, f'Interface Analysis\n\n'
                              f'Stability: {multiphase_result["interface_stability"]:.3f}\n'
                              f'Drift: {multiphase_result["interface_drift"]:.3f}\n'
                              f'Mass Conservation: {multiphase_result["mass_conservation"]:.2f}%',
                   ha='center', va='center', transform=ax.transAxes,
                   bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
            ax.axis('off')
            
            # Plot 3: System parameters
            ax = axes[1, 0]
            ax.text(0.05, 0.95, 'Simulation Parameters\n\n'
                                '• Grid: 256×64\n'
                                '• Density ratio: 1000:1\n'
                                '• Viscosity ratio: 100:1\n'
                                '• Interface thickness: 4\n'
                                '• Time steps: 20,000\n'
                                '• Method: Conservative phase-field',
                   transform=ax.transAxes, verticalalignment='top', fontsize=11,
                   bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
            ax.axis('off')
            
            # Plot 4: Status summary
            ax = axes[1, 1]
            status_color = 'lightgreen' if multiphase_result['data_quality'] > 80 else 'lightcoral'
            status_text = 'GOOD' if multiphase_result['data_quality'] > 80 else 'NEEDS ATTENTION'
            
            ax.text(0.5, 0.7, f'Overall Status: {status_text}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=14, fontweight='bold')
            
            ax.text(0.5, 0.3, f'Files processed: {multiphase_result["total_files"]}\n'
                              f'Valid files: {multiphase_result["valid_files"]}\n'
                              f'Data quality: {multiphase_result["data_quality"]:.1f}%',
                   ha='center', va='center', transform=ax.transAxes, fontsize=12,
                   bbox=dict(boxstyle='round', facecolor=status_color, alpha=0.8))
            ax.axis('off')
        
        plt.tight_layout()
        output_file = self.output_dir / "multiphase_analysis.png"
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.close()
        
        return str(output_file)

def main():
    """Main function to generate standardized plots."""
    print("Generating standardized LBM analysis plots...")
    
    # Initialize components
    loader = DataLoader()
    single_analyzer = SinglePhaseAnalyzer(loader)
    multi_analyzer = MultiphaseAnalyzer(loader)
    plotter = PlotGenerator()
    
    # Analyze single-phase data
    print("Analyzing single-phase data...")
    single_results = {}
    for mode in ['standard', 'entropic']:
        result = single_analyzer.analyze_convergence(mode)
        if result:
            single_results[mode] = result
            print(f"  {mode}: Final RMSE = {result['final_rmse']:.2e}")
        else:
            print(f"  {mode}: No valid data found")
    
    # Analyze multiphase data
    print("Analyzing multiphase data...")
    multi_result = multi_analyzer.analyze_quality()
    if multi_result:
        print(f"  Data quality: {multi_result['data_quality']:.1f}%")
        print(f"  Interface stability: {multi_result['interface_stability']:.3f}")
    else:
        print("  No valid multiphase data found")
    
    # Generate plots
    print("Generating plots...")
    velocity_plot = plotter.create_velocity_analysis(single_results)
    multiphase_plot = plotter.create_multiphase_analysis(multi_result)
    
    print(f"\nGenerated plots:")
    print(f"  - {velocity_plot}")
    print(f"  - {multiphase_plot}")
    
    return [velocity_plot, multiphase_plot]

if __name__ == "__main__":
    main()