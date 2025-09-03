#!/usr/bin/env python3
"""
Comprehensive LBM plotting system with proper time-series analysis.
Addresses all identified issues:
- Proper interface stability calculation
- Time-dependent error evolution  
- BGK method differentiation
- Structured plot organization
- Verbose multi-panel plots
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set high-quality plotting defaults
plt.style.use('default')
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['lines.linewidth'] = 2

# Create structured directory layout
plots_dir = Path("plots")
plots_dir.mkdir(exist_ok=True)
(plots_dir / "single_phase").mkdir(exist_ok=True)
(plots_dir / "multiphase").mkdir(exist_ok=True)
(plots_dir / "analysis").mkdir(exist_ok=True)

def safe_load_csv(filename):
    """Safely load CSV with error handling."""
    try:
        data = pd.read_csv(filename)
        return data.replace([np.inf, -np.inf], np.nan)
    except Exception as e:
        print(f"Warning: Could not load {filename}: {e}")
        return None

def get_time_series_data(pattern, data_dir="data"):
    """Extract time series data from files matching pattern."""
    files = sorted(glob.glob(str(Path(data_dir) / pattern)))
    times = []
    data_list = []
    
    for file in files:
        try:
            t = int(Path(file).stem.split('_t')[1])
            data = safe_load_csv(file)
            if data is not None and len(data) > 0:
                times.append(t)
                data_list.append(data)
        except Exception as e:
            print(f"Warning: Could not process {file}: {e}")
            continue
    
    return times, data_list

def create_physical_system_diagram():
    """Create a diagram showing the physical system setup."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Single-phase Poiseuille flow diagram
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 5)
    
    # Channel walls
    ax1.fill_between([0, 10], [0, 0], [0.5, 0.5], color='gray', alpha=0.8, label='Wall')
    ax1.fill_between([0, 10], [4.5, 4.5], [5, 5], color='gray', alpha=0.8)
    
    # Flow visualization
    y_flow = np.linspace(0.5, 4.5, 20)
    x_start = 1
    for y in y_flow:
        # Parabolic velocity profile
        u = 4 * (y - 0.5) * (4.5 - y) / (4**2) * 2
        ax1.arrow(x_start, y, u, 0, head_width=0.1, head_length=0.2, fc='blue', ec='blue', alpha=0.7)
    
    ax1.text(5, 2.5, 'Poiseuille Flow\nu(y) = u_max[1-(2y/H-1)²]', ha='center', va='center', 
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    ax1.set_title('Single-Phase Channel Flow')
    ax1.set_xlabel('Flow Direction (x)')
    ax1.set_ylabel('Channel Height (y)')
    ax1.legend()
    
    # Multiphase system diagram
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 5)
    
    # Two fluid regions
    ax2.fill_between([0, 10], [0, 0], [2.5, 2.5], color='lightblue', alpha=0.7, label='Heavy Fluid (ρ_H)')
    ax2.fill_between([0, 10], [2.5, 2.5], [5, 5], color='lightcoral', alpha=0.7, label='Light Fluid (ρ_L)')
    
    # Interface
    ax2.plot([0, 10], [2.5, 2.5], 'k--', linewidth=3, label='Interface')
    
    # Surface tension visualization
    for x in np.linspace(1, 9, 5):
        ax2.annotate('', xy=(x, 2.7), xytext=(x, 2.3), 
                    arrowprops=dict(arrowstyle='<->', color='red', lw=2))
    
    ax2.text(5, 1.2, 'ρ_H = 1000ρ_L\nμ_H = 100μ_L', ha='center', va='center',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    ax2.text(5, 3.8, 'Surface Tension\nInterface Tracking', ha='center', va='center',
             bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
    
    ax2.set_title('Multiphase Interface Flow')
    ax2.set_xlabel('Domain Width (x)')
    ax2.set_ylabel('Domain Height (y)')
    ax2.legend()
    
    plt.tight_layout()
    output_file = plots_dir / "analysis" / "physical_system_diagram.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.show()
    
    print(f"✓ Physical system diagram saved: {output_file}")
    return str(output_file)

def plot_single_phase_evolution():
    """Create comprehensive single-phase evolution plots with proper data validation."""
    fig = plt.figure(figsize=(15, 12))
    fig.suptitle('Single-Phase LBM Flow Evolution: Standard vs Entropic BGK Analysis', fontsize=16, fontweight='bold')
    
    methods = ['standard', 'entropic']
    method_names = ['Standard BGK', 'Entropic BGK']
    colors = ['#1f77b4', '#ff7f0e']
    markers = ['o', 's']
    
    # Collect and validate time evolution data
    all_times = {}
    all_errors = {}
    all_profiles = {}
    
    print("Loading single-phase data...")
    for method in methods:
        times, data_list = get_time_series_data(f"{method}_velocity_t*.csv")
        all_times[method] = times
        all_errors[method] = []
        all_profiles[method] = []
        
        print(f"{method.title()} BGK: {len(times)} time points found")
        
        for data in data_list:
            if data is not None and len(data.columns) >= 3:
                # Handle different possible column names
                if 'y' in data.columns:
                    y = data['y'].values
                else:
                    y = data.iloc[:, 0].values
                    
                if 'avg_ux' in data.columns:
                    ux_sim = data['avg_ux'].values
                else:
                    ux_sim = data.iloc[:, 1].values
                    
                if 'analytic_ux' in data.columns:
                    ux_analytical = data['analytic_ux'].values
                elif 'analytical_ux' in data.columns:
                    ux_analytical = data['analytical_ux'].values
                else:
                    ux_analytical = data.iloc[:, 2].values
                
                # Validate data quality
                if not np.any(np.isnan(ux_sim)) and not np.any(np.isnan(ux_analytical)):
                    error = np.abs(ux_sim - ux_analytical)
                    rmse = np.sqrt(np.mean(error**2))
                    
                    all_errors[method].append(rmse)
                    all_profiles[method].append((y, ux_sim, ux_analytical))
                else:
                    print(f"Warning: NaN values found in {method} data")
    
    # Plot 1: Final velocity profiles comparison (most important)
    ax1 = plt.subplot(3, 3, 1)
    
    for method, name, color, marker in zip(methods, method_names, colors, markers):
        times = all_times[method]
        profiles = all_profiles[method]
        
        if len(profiles) > 0:
            # Use final time step
            y, ux_sim, ux_analytical = profiles[-1]
            final_time = times[-1]
            
            ax1.plot(ux_sim, y, color=color, marker=marker, markersize=3, 
                    label=f'{name} (t={final_time})', linewidth=2, alpha=0.8)
            
            # Show analytical solution only once
            if method == 'standard':
                ax1.plot(ux_analytical, y, 'k:', linewidth=3, label='Analytical', alpha=0.9)
    
    ax1.set_xlabel('Velocity $u_x$')
    ax1.set_ylabel('Channel height $y$')
    ax1.set_title('Final Velocity Profiles Comparison')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)
    
    # Plot 2: RMSE evolution over time with detailed analysis
    ax2 = plt.subplot(3, 3, 2)
    convergence_summary = {}
    
    for method, name, color, marker in zip(methods, method_names, colors, markers):
        times = all_times[method]
        errors = all_errors[method]
        
        if len(times) > 1 and len(errors) > 1:
            ax2.semilogy(times, errors, marker=marker, color=color, label=name, 
                        markersize=4, linewidth=2, alpha=0.8)
            
            # Calculate convergence metrics
            initial_rmse = errors[0]
            final_rmse = errors[-1]
            improvement_factor = initial_rmse / final_rmse
            rate = np.log(final_rmse/initial_rmse) / (times[-1] - times[0])
            
            convergence_summary[method] = {
                'initial': initial_rmse,
                'final': final_rmse,
                'improvement': improvement_factor,
                'rate': rate
            }
            
            print(f"{name}: RMSE {initial_rmse:.2e} → {final_rmse:.2e} (improvement: {improvement_factor:.1f}×)")
    
    ax2.set_xlabel('Time Steps')
    ax2.set_ylabel('RMSE (log scale)')
    ax2.set_title('Error Convergence')
    ax2.grid(True)
    ax2.legend()
    
    # Plot 3: Error profiles at final time
    ax3 = plt.subplot(3, 3, 3)
    for method, name, color in zip(methods, method_names, colors):
        times = all_times[method]
        profiles = all_profiles[method]
        
        if len(times) > 0:
            y, ux_sim, ux_analytical = profiles[-1]
            error = np.abs(ux_sim - ux_analytical)
            
            ax3.plot(error, y, color=color, label=name, marker='o', markersize=3)
            
            max_error = np.max(error)
            ax3.axvline(x=max_error, color=color, linestyle=':', alpha=0.7)
    
    ax3.set_xlabel('Absolute Error')
    ax3.set_ylabel('Channel height $y$')
    ax3.set_title('Final Error Profiles')
    ax3.grid(True)
    ax3.legend()
    
    # Plot 4: Error difference analysis between methods
    ax4 = plt.subplot(3, 3, 4)
    
    if len(methods) == 2 and all(method in convergence_summary for method in methods):
        # Compare convergence rates and final errors
        method_labels = [name.replace(' BGK', '') for name in method_names]
        final_errors = [convergence_summary[method]['final'] for method in methods]
        improvement_factors = [convergence_summary[method]['improvement'] for method in methods]
        
        x = np.arange(len(method_labels))
        width = 0.35
        
        bars1 = ax4.bar(x - width/2, final_errors, width, label='Final RMSE', 
                       alpha=0.8, color=colors)
        ax4_twin = ax4.twinx()
        bars2 = ax4_twin.bar(x + width/2, improvement_factors, width, label='Improvement Factor', 
                            alpha=0.8, color=['lightgreen', 'lightcoral'])
        
        # Add value labels
        for i, (err, imp) in enumerate(zip(final_errors, improvement_factors)):
            ax4.text(i - width/2, err, f'{err:.1e}', ha='center', va='bottom', fontsize=8)
            ax4_twin.text(i + width/2, imp, f'{imp:.1f}×', ha='center', va='bottom', fontsize=8)
        
        ax4.set_ylabel('Final RMSE', color='blue')
        ax4_twin.set_ylabel('Improvement Factor', color='green')
        ax4.set_title('Method Performance Comparison')
        ax4.set_xticks(x)
        ax4.set_xticklabels(method_labels)
        ax4.grid(True, axis='y', alpha=0.3)
        ax4.set_yscale('log')
        
        # Determine better method
        better_idx = np.argmin(final_errors)
        better_method = method_labels[better_idx]
        performance_diff = (max(final_errors) - min(final_errors)) / max(final_errors) * 100
        
        ax4.text(0.5, 0.95, f'Better: {better_method}\n{performance_diff:.1f}% improvement', 
                transform=ax4.transAxes, ha='center', va='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    else:
        ax4.text(0.5, 0.5, 'Insufficient data\nfor comparison', ha='center', va='center',
                transform=ax4.transAxes, fontsize=12)
    
    # Plot 5: Error statistics comparison
    ax5 = plt.subplot(3, 3, 5)
    methods_stats = []
    rmse_values = []
    max_error_values = []
    
    for method, name in zip(methods, method_names):
        if method in all_times and len(all_times[method]) > 0:
            y, ux_sim, ux_analytical = all_profiles[method][-1]
            error = np.abs(ux_sim - ux_analytical)
            rmse = np.sqrt(np.mean(error**2))
            max_error = np.max(error)
            
            methods_stats.append(name)
            rmse_values.append(rmse)
            max_error_values.append(max_error)
    
    if methods_stats:
        x = np.arange(len(methods_stats))
        width = 0.35
        
        bars1 = ax5.bar(x - width/2, rmse_values, width, label='RMSE', alpha=0.8, color=colors[:len(methods_stats)])
        bars2 = ax5.bar(x + width/2, max_error_values, width, label='Max Error', alpha=0.8, color=colors[:len(methods_stats)])
        
        for i, (rmse, max_err) in enumerate(zip(rmse_values, max_error_values)):
            ax5.text(i - width/2, rmse + rmse*0.1, f'{rmse:.1e}', ha='center', va='bottom', fontsize=8)
            ax5.text(i + width/2, max_err + max_err*0.1, f'{max_err:.1e}', ha='center', va='bottom', fontsize=8)
        
        ax5.set_ylabel('Error Magnitude')
        ax5.set_title('Final Error Metrics')
        ax5.set_xticks(x)
        ax5.set_xticklabels(methods_stats)
        ax5.legend()
        ax5.grid(True, axis='y')
    
    # Plot 6: Detailed error analysis summary
    ax6 = plt.subplot(3, 3, 6)
    ax6.axis('off')
    
    summary_text = "Analysis Summary\n\n"
    
    if convergence_summary:
        for method, name in zip(methods, method_names):
            if method in convergence_summary:
                data = convergence_summary[method]
                summary_text += f"{name}:\n"
                summary_text += f"• Initial RMSE: {data['initial']:.2e}\n"
                summary_text += f"• Final RMSE: {data['final']:.2e}\n"
                summary_text += f"• Improvement: {data['improvement']:.1f}×\n"
                summary_text += f"• Rate: {data['rate']:.2e}/step\n\n"
        
        # Overall assessment
        if len(convergence_summary) == 2:
            std_final = convergence_summary['standard']['final']
            ent_final = convergence_summary['entropic']['final']
            
            if std_final < ent_final:
                summary_text += "Result: Standard BGK performs better\n"
                summary_text += f"Advantage: {(ent_final/std_final - 1)*100:.1f}% lower RMSE"
            elif ent_final < std_final:
                summary_text += "Result: Entropic BGK performs better\n"
                summary_text += f"Advantage: {(std_final/ent_final - 1)*100:.1f}% lower RMSE"
            else:
                summary_text += "Result: Both methods equivalent\n"
                summary_text += "Performance difference < 1%"
    
    ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes, 
             verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    # Plot 7: Physical system representation
    ax7 = plt.subplot(3, 3, 7)
    ax7.set_xlim(0, 10)
    ax7.set_ylim(0, 5)
    
    # Channel walls
    ax7.fill_between([0, 10], [0, 0], [0.3, 0.3], color='gray', alpha=0.8)
    ax7.fill_between([0, 10], [4.7, 4.7], [5, 5], color='gray', alpha=0.8)
    
    # Velocity profile visualization
    if len(all_profiles['standard']) > 0:
        y, ux_sim, ux_analytical = all_profiles['standard'][-1]
        # Scale velocities for visualization
        ux_scaled = ux_sim * 20000  # Scale factor for arrows
        y_positions = np.linspace(0.5, 4.5, min(len(y), 15))
        
        for i, y_pos in enumerate(y_positions):
            if i < len(ux_scaled):
                arrow_length = max(0.1, ux_scaled[i])
                ax7.arrow(2, y_pos, arrow_length, 0, head_width=0.1, head_length=0.1, 
                         fc='blue', ec='blue', alpha=0.7)
    
    ax7.text(6, 2.5, 'Poiseuille Flow\nChannel\n\nGrid: 20×21\nτ = 1.0\nF = 1×10⁻⁶', 
             ha='center', va='center', 
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    ax7.set_title('Physical System')
    ax7.set_xlabel('Flow Direction')
    ax7.set_ylabel('Channel Height')
    ax7.grid(True, alpha=0.3)
    
    # Plot 8: Computational details
    ax8 = plt.subplot(3, 3, 8)
    ax8.axis('off')
    
    # Count total data points analyzed
    total_points = sum(len(all_times[method]) for method in methods)
    total_profiles = sum(len(all_profiles[method]) for method in methods)
    
    comp_text = "Computational Details\n\n"
    comp_text += f"Grid Resolution: 20×21\n"
    comp_text += f"Time Steps: 0 → 10,000\n"
    comp_text += f"Output Interval: 1,000 steps\n"
    comp_text += f"Data Points: {total_points}\n"
    comp_text += f"Profiles Analyzed: {total_profiles}\n\n"
    
    comp_text += "Parameters:\n"
    comp_text += "• Relaxation time: τ = 1.0\n"
    comp_text += "• Body force: 1×10⁻⁶\n"
    comp_text += "• Boundary: No-slip walls\n"
    comp_text += "• Collision: BGK operators\n\n"
    
    comp_text += "Validation: Poiseuille solution\n"
    if convergence_summary:
        avg_final_rmse = np.mean([convergence_summary[method]['final'] for method in convergence_summary])
        comp_text += f"Average RMSE: {avg_final_rmse:.2e}"
    
    ax8.text(0.05, 0.95, comp_text, transform=ax8.transAxes, 
             verticalalignment='top', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # Plot 9: Reserved for future analysis
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    ax9.text(0.5, 0.5, 'Additional Analysis\nSpace Reserved', 
             ha='center', va='center', transform=ax9.transAxes, fontsize=12,
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    
    plt.tight_layout()
    output_file = plots_dir / "single_phase" / "comprehensive_evolution.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.show()
    
    print(f"✓ Single-phase evolution plot saved: {output_file}")
    return str(output_file)

def create_2d_field_visualization():
    """Create 2D field visualization of multiphase flow."""
    # Look for full field data files
    field_files = sorted(glob.glob("data/multiphase_t*.csv"))
    if not field_files:
        print("Warning: No 2D field data found")
        return None
    
    # Select a representative time step
    target_time = 10000
    selected_file = None
    for file in field_files:
        try:
            t = int(Path(file).stem.split('_t')[1])
            if t == target_time:
                selected_file = file
                break
        except:
            continue
    
    if not selected_file:
        selected_file = field_files[-1]  # Use last available
        target_time = int(Path(selected_file).stem.split('_t')[1])
    
    data = safe_load_csv(selected_file)
    if data is None:
        print(f"Warning: Could not load {selected_file}")
        return None
    
    # Create 2D field visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'2D Multiphase Flow Fields at t={target_time}', fontsize=16, fontweight='bold')
    
    # Reshape data to 2D grid (assuming structured data)
    if 'x' in data.columns and 'y' in data.columns:
        x_unique = sorted(data['x'].unique())
        y_unique = sorted(data['y'].unique())
        
        if len(x_unique) * len(y_unique) == len(data):
            nx, ny = len(x_unique), len(y_unique)
            
            # Phase field visualization
            if 'phi' in data.columns:
                phi_2d = data['phi'].values.reshape(ny, nx)
                im1 = axes[0,0].imshow(phi_2d, origin='lower', cmap='RdBu_r', aspect='auto')
                axes[0,0].set_title('Phase Field φ')
                axes[0,0].set_xlabel('x')
                axes[0,0].set_ylabel('y')
                plt.colorbar(im1, ax=axes[0,0])
            
            # Density field
            if 'rho' in data.columns:
                rho_2d = data['rho'].values.reshape(ny, nx)
                im2 = axes[0,1].imshow(rho_2d, origin='lower', cmap='viridis', aspect='auto')
                axes[0,1].set_title('Density ρ')
                axes[0,1].set_xlabel('x')
                axes[0,1].set_ylabel('y')
                plt.colorbar(im2, ax=axes[0,1])
            
            # Velocity magnitude
            if 'ux' in data.columns and 'uy' in data.columns:
                u_mag = np.sqrt(data['ux']**2 + data['uy']**2)
                u_mag_2d = u_mag.values.reshape(ny, nx)
                im3 = axes[1,0].imshow(u_mag_2d, origin='lower', cmap='plasma', aspect='auto')
                axes[1,0].set_title('Velocity Magnitude |u|')
                axes[1,0].set_xlabel('x')
                axes[1,0].set_ylabel('y')
                plt.colorbar(im3, ax=axes[1,0])
                
                # Velocity vectors (subsampled)
                skip = max(1, nx//20)  # Show every 20th vector
                x_2d = data['x'].values.reshape(ny, nx)
                y_2d = data['y'].values.reshape(ny, nx)
                ux_2d = data['ux'].values.reshape(ny, nx)
                uy_2d = data['uy'].values.reshape(ny, nx)
                
                axes[1,1].quiver(x_2d[::skip, ::skip], y_2d[::skip, ::skip], 
                               ux_2d[::skip, ::skip], uy_2d[::skip, ::skip],
                               u_mag_2d[::skip, ::skip], cmap='plasma', scale=0.01)
                axes[1,1].set_title('Velocity Field')
                axes[1,1].set_xlabel('x')
                axes[1,1].set_ylabel('y')
                axes[1,1].set_aspect('equal')
    
    plt.tight_layout()
    output_file = plots_dir / "multiphase" / "2d_field_visualization.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.show()
    
    print(f"✓ 2D field visualization saved: {output_file}")
    return str(output_file)

def plot_multiphase_analysis():
    """Create multiphase data quality analysis and status report."""
    fig = plt.figure(figsize=(15, 10))
    fig.suptitle('Multiphase LBM Data Status Report', fontsize=16, fontweight='bold')
    
    # Get time series data
    times, data_list = get_time_series_data("multiphase_avg_t*.csv")
    
    if len(times) == 0:
        print("Warning: No multiphase data found")
        # Create empty report
        ax = plt.subplot(1, 1, 1)
        ax.text(0.5, 0.5, 'No Multiphase Data Found\n\nPlease run multiphase simulation first', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16,
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        ax.axis('off')
        plt.tight_layout()
        output_file = plots_dir / "multiphase" / "data_status_report.png"
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        plt.show()
        return str(output_file)
    
    print(f"Found {len(times)} multiphase time points: {times}")
    
    # Analyze data quality
    nan_analysis = []
    valid_data_count = 0
    
    for i, (t, data) in enumerate(zip(times, data_list)):
        if data is not None:
            total_cells = data.size
            nan_count = data.isna().sum().sum()
            nan_percentage = (nan_count / total_cells) * 100 if total_cells > 0 else 100
            
            # Check if any columns have valid data
            valid_columns = []
            for col in data.columns:
                if col != 'y':  # y-coordinate should always be valid
                    valid_count = data[col].notna().sum()
                    if valid_count > 0:
                        valid_columns.append(col)
            
            nan_analysis.append({
                'time': t,
                'nan_percentage': nan_percentage,
                'valid_columns': valid_columns,
                'total_cells': total_cells
            })
            
            if len(valid_columns) > 1:  # More than just coordinate columns
                valid_data_count += 1
        else:
            nan_analysis.append({
                'time': t,
                'nan_percentage': 100,
                'valid_columns': [],
                'total_cells': 0
            })
    
    # Create comprehensive data quality report
    print(f"Data quality analysis: {valid_data_count}/{len(times)} files have valid data")
    
    # Plot 1: Data quality timeline
    ax1 = plt.subplot(2, 3, 1)
    
    if nan_analysis:
        plot_times = [item['time'] for item in nan_analysis]
        nan_percentages = [item['nan_percentage'] for item in nan_analysis]
        
        ax1.plot(plot_times, nan_percentages, 'ro-', linewidth=2, markersize=6)
        ax1.axhline(y=100, color='red', linestyle='--', alpha=0.7, label='Complete NaN')
        ax1.axhline(y=50, color='orange', linestyle='--', alpha=0.7, label='50% NaN')
        ax1.axhline(y=0, color='green', linestyle='--', alpha=0.7, label='No NaN')
        
        ax1.set_xlabel('Time Steps')
        ax1.set_ylabel('NaN Percentage (%)')
        ax1.set_title('Data Quality Over Time')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        ax1.set_ylim(-5, 105)
        
        # Add status annotations
        if np.mean(nan_percentages) > 90:
            ax1.text(0.5, 0.5, 'CRITICAL\nData Corruption', ha='center', va='center',
                    transform=ax1.transAxes, fontsize=12, color='red',
                    bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
    
    # Plot 2: Column validity analysis
    ax2 = plt.subplot(2, 3, 2)
    
    if nan_analysis:
        # Count valid columns across all time steps
        all_columns = set()
        for item in nan_analysis:
            all_columns.update(item['valid_columns'])
        
        if all_columns:
            column_validity = {}
            for col in all_columns:
                valid_count = sum(1 for item in nan_analysis if col in item['valid_columns'])
                column_validity[col] = (valid_count / len(nan_analysis)) * 100
            
            cols = list(column_validity.keys())
            validity_pct = list(column_validity.values())
            
            bars = ax2.bar(cols, validity_pct, alpha=0.7, 
                          color=['green' if v > 50 else 'orange' if v > 0 else 'red' for v in validity_pct])
            
            for bar, pct in zip(bars, validity_pct):
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, 
                        f'{pct:.0f}%', ha='center', va='bottom', fontsize=8)
            
            ax2.set_ylabel('Validity Percentage (%)')
            ax2.set_title('Column Data Validity')
            ax2.set_ylim(0, 110)
            ax2.grid(True, axis='y', alpha=0.3)
            plt.setp(ax2.get_xticklabels(), rotation=45, ha='right')
        else:
            ax2.text(0.5, 0.5, 'No Valid Columns\nFound', ha='center', va='center',
                    transform=ax2.transAxes, fontsize=12,
                    bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
            ax2.axis('off')
    
    # Plot 3: Sample data structure
    ax3 = plt.subplot(2, 3, 3)
    ax3.axis('off')
    
    if data_list and data_list[0] is not None:
        sample_data = data_list[0]
        structure_text = "Sample Data Structure\n\n"
        structure_text += f"Columns: {list(sample_data.columns)}\n"
        structure_text += f"Shape: {sample_data.shape}\n"
        structure_text += f"Memory: {sample_data.memory_usage(deep=True).sum()/1024:.1f} KB\n\n"
        
        structure_text += "Data Types:\n"
        for col, dtype in sample_data.dtypes.items():
            structure_text += f"• {col}: {dtype}\n"
        
        ax3.text(0.05, 0.95, structure_text, transform=ax3.transAxes, 
                verticalalignment='top', fontsize=9, family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    # Phase field profiles at different times
    plt.subplot(2, 3, 2)
    selected_times = [0, 5000, 10000, 15000, 20000]
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_times)))
    
    for t, color in zip(selected_times, colors):
        filename = f"data/multiphase_avg_t{t}.csv"
        data = safe_load_csv(filename)
        if data is not None and 'avg_phi' in data.columns and len(data) > 0:
            y = data['y'].values
            phi = data['avg_phi'].values
            
            # Remove NaN values
            valid_idx = ~np.isnan(phi)
            if np.any(valid_idx):
                plt.plot(phi[valid_idx], y[valid_idx], color=color, linewidth=2, 
                        label=f't={t}', alpha=0.8)
    
    plt.xlabel('Phase field φ')
    plt.ylabel('Height y')
    plt.title('Phase Field Evolution')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xlim(0, 1)
    
    # Density profiles
    plt.subplot(2, 3, 3)
    for t, color in zip(selected_times, colors):
        filename = f"data/multiphase_avg_t{t}.csv"
        data = safe_load_csv(filename)
        if data is not None and 'avg_rho' in data.columns and len(data) > 0:
            y = data['y'].values
            rho = data['avg_rho'].values
            
            # Remove NaN values
            valid_idx = ~np.isnan(rho)
            if np.any(valid_idx):
                plt.plot(rho[valid_idx], y[valid_idx], color=color, linewidth=2, 
                        label=f't={t}', alpha=0.8)
    
    plt.xlabel('Density ρ')
    plt.ylabel('Height y')
    plt.title('Density Evolution')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Velocity profiles (if available and valid)
    plt.subplot(2, 3, 4)
    for t, color in zip(selected_times, colors):
        filename = f"data/multiphase_avg_t{t}.csv"
        data = safe_load_csv(filename)
        if data is not None and 'avg_ux' in data.columns and len(data) > 0:
            y = data['y'].values
            ux = data['avg_ux'].values
            
            # Remove NaN values
            valid_idx = ~np.isnan(ux)
            if np.any(valid_idx) and len(ux[valid_idx]) > 1:
                plt.plot(ux[valid_idx], y[valid_idx], color=color, linewidth=2, 
                        label=f't={t}', alpha=0.8)
    
    plt.xlabel('Velocity $u_x$')
    plt.ylabel('Height y')
    plt.title('Velocity Evolution')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Mass conservation check
    plt.subplot(2, 3, 5)
    
    # Get multiphase files for mass conservation analysis
    times, data_list = get_time_series_data("multiphase_avg_t*.csv")
    mass_history = []
    valid_times = []
    
    for t, data in zip(times, data_list):
        if data is not None and 'avg_rho' in data.columns and len(data) > 0:
            rho = data['avg_rho'].values
            y = data['y'].values
            
            # Calculate total mass (simple integration)
            valid_idx = ~np.isnan(rho)
            if np.any(valid_idx) and len(rho[valid_idx]) > 1:
                dy = y[1] - y[0] if len(y) > 1 else 1.0
                total_mass = np.sum(rho[valid_idx]) * dy
                mass_history.append(total_mass)
                valid_times.append(t)
    
    if len(mass_history) > 1:
        mass_history = np.array(mass_history)
        mass_change = (mass_history - mass_history[0]) / mass_history[0] * 100
        plt.plot(valid_times, mass_change, 'ro-', linewidth=2, markersize=5)
        plt.xlabel('Time Steps')
        plt.ylabel('Mass Change (%)')
        plt.title('Mass Conservation')
        plt.grid(True, alpha=0.3)
        plt.axhline(y=0, color='black', linestyle='--', alpha=0.7)
    
    # Summary statistics
    plt.subplot(2, 3, 6)
    plt.axis('off')
    
    # Calculate interface positions for stability analysis
    interface_positions = []
    times, data_list = get_time_series_data("multiphase_avg_t*.csv")
    
    for t, data in zip(times, data_list):
        if data is not None and 'avg_phi' in data.columns and len(data) > 0:
            phi = data['avg_phi'].values
            y = data['y'].values
            
            # Find interface position (where phi = 0.5)
            valid_idx = ~np.isnan(phi)
            if np.any(valid_idx) and len(phi[valid_idx]) > 1:
                # Find closest position to phi = 0.5
                phi_centered = phi[valid_idx] - 0.5
                interface_idx = np.argmin(np.abs(phi_centered))
                if interface_idx < len(y[valid_idx]):
                    interface_positions.append(y[valid_idx][interface_idx])
    
    # Calculate interface drift
    drift = 0.0
    if len(interface_positions) > 1:
        drift = interface_positions[-1] - interface_positions[0]
    
    # Create summary text
    summary_text = "Multiphase Simulation Summary\n\n"
    summary_text += f"Total files processed: {len(data_list)}\n"
    summary_text += f"Time steps analyzed: {len(valid_times)}\n"
    
    if len(interface_positions) > 1:
        summary_text += f"Interface drift: {drift:.3f}\n"
        summary_text += f"Interface stability: {np.std(interface_positions):.3f}\n"
    
    if len(mass_history) > 1:
        max_mass_change = np.max(np.abs(mass_change))
        summary_text += f"Max mass change: {max_mass_change:.2f}%\n"
    
    summary_text += "\nDensity ratio: 1000:1\n"
    summary_text += "Viscosity ratio: 100:1\n"
    summary_text += "Grid: 256×64"
    
    plt.text(0.1, 0.9, summary_text, transform=plt.gca().transAxes, 
            verticalalignment='top', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    output_file = plots_dir / "multiphase" / "comprehensive_analysis.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.show()
    
    print(f"✓ Multiphase analysis plot saved: {output_file}")
    return str(output_file)

def main():
    """Generate all comprehensive plots with physical system visualizations."""
    print("Starting comprehensive LBM analysis plotting...")
    
    generated_plots = []
    
    # Generate physical system diagrams
    try:
        system_diagram = create_physical_system_diagram()
        if system_diagram:
            generated_plots.append(system_diagram)
    except Exception as e:
        print(f"Error generating system diagram: {e}")
    
    # Generate 2D field visualization
    try:
        field_viz = create_2d_field_visualization()
        if field_viz:
            generated_plots.append(field_viz)
    except Exception as e:
        print(f"Error generating 2D field visualization: {e}")
    
    # Generate single-phase analysis
    try:
        single_phase_plot = plot_single_phase_evolution()
        if single_phase_plot:
            generated_plots.append(single_phase_plot)
    except Exception as e:
        print(f"Error generating single-phase plots: {e}")
    
    # Generate multiphase analysis
    try:
        multiphase_plot = plot_multiphase_analysis()
        if multiphase_plot:
            generated_plots.append(multiphase_plot)
    except Exception as e:
        print(f"Error generating multiphase plots: {e}")
    
    print(f"\nPlot generation completed!")
    print(f"Generated plots: {len(generated_plots)}")
    for plot in generated_plots:
        print(f"  - {plot}")
    
    return generated_plots

if __name__ == "__main__":
    main()