"""
H-Theorem Specific Visualization Module

This module provides specialized visualization for H-theorem analysis including
evolution plots, statistical analysis, and comparative visualizations between
standard and entropic BGK methods.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from dataclasses import dataclass
from pathlib import Path
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import seaborn as sns
from scipy import stats
from scipy.signal import savgol_filter

from .layout_manager import LayoutManager, PlotConfig, PlotType, LayoutConfig


@dataclass
class HTheoremData:
    """Data structure for H-theorem analysis results."""
    evolution_data: pd.DataFrame  # Full evolution data
    timesteps: np.ndarray
    h_values: np.ndarray
    h_changes: np.ndarray
    entropy_production_rates: np.ndarray
    monotonicity_violations: np.ndarray
    # Optional statistical data
    h_variances: Optional[np.ndarray] = None
    kinetic_energies: Optional[np.ndarray] = None
    total_masses: Optional[np.ndarray] = None
    # Derived metrics
    violation_count: int = 0
    overall_monotonic: bool = True
    average_entropy_production: float = 0.0
    decay_rate: float = 0.0


@dataclass
class AnomalyDetection:
    """Anomaly detection results for H-theorem analysis."""
    anomaly_timesteps: List[int]
    anomaly_types: List[str]  # 'violation', 'spike', 'plateau'
    anomaly_magnitudes: List[float]
    statistical_outliers: List[int]


class HTheoremPlots:
    """Specialized H-theorem visualization with statistical analysis and anomaly detection."""
    
    def __init__(self, layout_config: Optional[LayoutConfig] = None):
        self.layout_manager = LayoutManager(layout_config)
        self.violation_color = '#ff4444'
        self.normal_color = '#4444ff'
        self.entropy_color = '#44ff44'
        
    def load_h_theorem_data(self, filepath: Union[str, Path]) -> HTheoremData:
        """Load H-theorem evolution data from CSV file."""
        df = pd.read_csv(filepath)
        
        # Extract basic arrays
        timesteps = df['timestep'].values
        h_values = df['h_value'].values
        h_changes = df['h_change'].values if 'h_change' in df.columns else np.gradient(h_values)
        entropy_production = df['entropy_production'].values if 'entropy_production' in df.columns else -h_changes
        violations = df['monotonic_violation'].values.astype(bool) if 'monotonic_violation' in df.columns else np.zeros_like(timesteps, dtype=bool)
        
        # Optional statistical data
        h_variances = df['h_variance'].values if 'h_variance' in df.columns else None
        kinetic_energies = df['kinetic_energy'].values if 'kinetic_energy' in df.columns else None
        total_masses = df['total_mass'].values if 'total_mass' in df.columns else None
        
        # Calculate derived metrics
        violation_count = int(np.sum(violations))
        overall_monotonic = violation_count == 0
        average_entropy_production = np.mean(entropy_production[~np.isnan(entropy_production)])
        
        # Calculate decay rate (exponential fit)
        decay_rate = 0.0
        if len(h_values) > 10:
            try:
                # Fit exponential decay: H(t) = H0 * exp(-r*t)
                valid_mask = (h_values > 0) & ~np.isnan(h_values)
                if np.sum(valid_mask) > 5:
                    log_h = np.log(h_values[valid_mask])
                    t_valid = timesteps[valid_mask]
                    slope, _, _, _, _ = stats.linregress(t_valid, log_h)
                    decay_rate = -slope
            except:
                decay_rate = 0.0
        
        return HTheoremData(
            evolution_data=df,
            timesteps=timesteps,
            h_values=h_values,
            h_changes=h_changes,
            entropy_production_rates=entropy_production,
            monotonicity_violations=violations,
            h_variances=h_variances,
            kinetic_energies=kinetic_energies,
            total_masses=total_masses,
            violation_count=violation_count,
            overall_monotonic=overall_monotonic,
            average_entropy_production=average_entropy_production,
            decay_rate=decay_rate
        )
    
    def detect_anomalies(self, data: HTheoremData, 
                        violation_threshold: float = 1e-12,
                        spike_threshold: float = 3.0) -> AnomalyDetection:
        """Detect various types of anomalies in H-theorem evolution."""
        anomaly_timesteps = []
        anomaly_types = []
        anomaly_magnitudes = []
        
        # Monotonicity violations
        violation_indices = np.where(data.monotonicity_violations)[0]
        for idx in violation_indices:
            if idx < len(data.timesteps):
                anomaly_timesteps.append(data.timesteps[idx])
                anomaly_types.append('violation')
                anomaly_magnitudes.append(data.h_changes[idx])
        
        # Statistical outliers in H-function changes
        if len(data.h_changes) > 10:
            # Use z-score for outlier detection
            z_scores = np.abs(stats.zscore(data.h_changes, nan_policy='omit'))
            outlier_indices = np.where(z_scores > spike_threshold)[0]
            
            for idx in outlier_indices:
                if idx < len(data.timesteps):
                    anomaly_timesteps.append(data.timesteps[idx])
                    anomaly_types.append('spike')
                    anomaly_magnitudes.append(z_scores[idx])
        
        # Detect plateaus (regions with very small changes)
        if len(data.h_changes) > 20:
            # Smooth the data and look for flat regions
            window_size = min(11, len(data.h_changes) // 5)
            if window_size >= 3:
                smoothed_changes = savgol_filter(np.abs(data.h_changes), window_size, 3)
                plateau_threshold = np.percentile(smoothed_changes, 5)  # Bottom 5%
                plateau_indices = np.where(smoothed_changes < plateau_threshold)[0]
                
                # Group consecutive indices
                if len(plateau_indices) > 0:
                    groups = np.split(plateau_indices, np.where(np.diff(plateau_indices) != 1)[0] + 1)
                    for group in groups:
                        if len(group) > 5:  # At least 5 consecutive points
                            mid_idx = group[len(group)//2]
                            if mid_idx < len(data.timesteps):
                                anomaly_timesteps.append(data.timesteps[mid_idx])
                                anomaly_types.append('plateau')
                                anomaly_magnitudes.append(len(group))
        
        # Statistical outliers using IQR method
        statistical_outliers = []
        if len(data.h_values) > 10:
            q1, q3 = np.percentile(data.h_values, [25, 75])
            iqr = q3 - q1
            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr
            outlier_mask = (data.h_values < lower_bound) | (data.h_values > upper_bound)
            statistical_outliers = data.timesteps[outlier_mask].tolist()
        
        return AnomalyDetection(
            anomaly_timesteps=anomaly_timesteps,
            anomaly_types=anomaly_types,
            anomaly_magnitudes=anomaly_magnitudes,
            statistical_outliers=statistical_outliers
        )
    
    def create_h_evolution_plot(self, data: HTheoremData, 
                               show_anomalies: bool = True) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create comprehensive H-function evolution plot with anomaly detection."""
        plot_configs = [
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "H-Function Evolution", 
                      "Timestep", "H-Function Value"),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "H-Function Change Rate", 
                      "Timestep", "dH/dt"),
            PlotConfig(PlotType.STABILITY_INDICATORS, "Cumulative Violations", 
                      "Timestep", "Violation Count")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title="H-Theorem Evolution Analysis"
        )
        
        # H-function evolution with log scale
        axes[0].semilogy(data.timesteps, np.abs(data.h_values), 
                        linewidth=2, color=self.normal_color, label='|H(t)|')
        axes[0].set_ylabel("log |H-Function|")
        
        # Add exponential fit if decay rate is available
        if data.decay_rate > 0:
            t_fit = np.linspace(data.timesteps[0], data.timesteps[-1], 100)
            h0 = np.abs(data.h_values[0])
            h_fit = h0 * np.exp(-data.decay_rate * (t_fit - data.timesteps[0]))
            axes[0].plot(t_fit, h_fit, '--', color='red', alpha=0.7, 
                        label=f'Exp fit (r={data.decay_rate:.2e})')
        
        axes[0].legend()
        
        # H-function change rate
        axes[1].plot(data.timesteps, data.h_changes, linewidth=1.5, color=self.entropy_color)
        axes[1].axhline(y=0, color='black', linestyle='--', alpha=0.5)
        axes[1].set_ylabel("dH/dt")
        
        # Highlight negative regions (entropy increase - violations)
        violation_mask = data.h_changes > 0
        if np.any(violation_mask):
            axes[1].fill_between(data.timesteps, data.h_changes, 0, 
                                where=violation_mask, alpha=0.3, color=self.violation_color,
                                label='Violations')
            axes[1].legend()
        
        # Cumulative violations
        cumulative_violations = np.cumsum(data.monotonicity_violations)
        axes[2].step(data.timesteps, cumulative_violations, where='post', 
                    linewidth=2, color=self.violation_color)
        axes[2].set_ylabel("Cumulative Violations")
        
        # Add anomaly markers if requested
        if show_anomalies:
            anomalies = self.detect_anomalies(data)
            for timestep, anomaly_type in zip(anomalies.anomaly_timesteps, anomalies.anomaly_types):
                if anomaly_type == 'violation':
                    for ax in axes:
                        ax.axvline(x=timestep, color=self.violation_color, 
                                 linestyle=':', alpha=0.7, linewidth=1)
                elif anomaly_type == 'spike':
                    for ax in axes:
                        ax.axvline(x=timestep, color='orange', 
                                 linestyle=':', alpha=0.7, linewidth=1)
        
        return fig, axes
    
    def create_statistical_analysis_plot(self, data: HTheoremData) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create statistical analysis of H-theorem behavior."""
        if data.h_variances is None:
            raise ValueError("Statistical data (h_variances) required for this plot")
        
        plot_configs = [
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "H-Function Variance", 
                      "Timestep", "Spatial Variance"),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "Energy Evolution", 
                      "Timestep", "Kinetic Energy"),
            PlotConfig(PlotType.STABILITY_INDICATORS, "Mass Conservation", 
                      "Timestep", "Total Mass")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title="Statistical Analysis"
        )
        
        # H-function spatial variance
        axes[0].semilogy(data.timesteps, data.h_variances, 
                        'o-', markersize=3, linewidth=1.5, color='purple')
        axes[0].set_ylabel("log(Variance)")
        
        # Kinetic energy evolution
        if data.kinetic_energies is not None:
            axes[1].plot(data.timesteps, data.kinetic_energies, 
                        'o-', markersize=3, linewidth=1.5, color='green')
            axes[1].set_ylabel("Kinetic Energy")
        
        # Mass conservation
        if data.total_masses is not None:
            # Calculate relative mass change
            initial_mass = data.total_masses[0]
            relative_mass = (data.total_masses - initial_mass) / initial_mass
            axes[2].plot(data.timesteps, relative_mass, 
                        'o-', markersize=3, linewidth=1.5, color='red')
            axes[2].set_ylabel("Relative Mass Change")
            axes[2].ticklabel_format(style='scientific', axis='y', scilimits=(-3, 3))
        
        return fig, axes
    
    def create_entropy_production_analysis(self, data: HTheoremData) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create detailed entropy production analysis."""
        plot_configs = [
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "Entropy Production Rate", 
                      "Timestep", "Entropy Production"),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "Cumulative Entropy Production", 
                      "Timestep", "Cumulative Entropy"),
            PlotConfig(PlotType.STABILITY_INDICATORS, "Production Rate Distribution", 
                      "Entropy Production Rate", "Frequency")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title="Entropy Production Analysis"
        )
        
        # Entropy production rate
        axes[0].plot(data.timesteps, data.entropy_production_rates, 
                    linewidth=1.5, color=self.entropy_color)
        axes[0].axhline(y=0, color='black', linestyle='--', alpha=0.5)
        axes[0].set_ylabel("dS/dt")
        
        # Smooth trend line
        if len(data.entropy_production_rates) > 20:
            window_size = min(21, len(data.entropy_production_rates) // 10)
            if window_size >= 3:
                smoothed = savgol_filter(data.entropy_production_rates, window_size, 3)
                axes[0].plot(data.timesteps, smoothed, '--', linewidth=2, 
                           alpha=0.8, color='red', label='Trend')
                axes[0].legend()
        
        # Cumulative entropy production
        cumulative_entropy = np.cumsum(data.entropy_production_rates)
        axes[1].plot(data.timesteps, cumulative_entropy, 
                    linewidth=2, color='darkgreen')
        axes[1].set_ylabel("Cumulative Entropy")
        
        # Distribution of production rates
        valid_rates = data.entropy_production_rates[~np.isnan(data.entropy_production_rates)]
        if len(valid_rates) > 10:
            axes[2].hist(valid_rates, bins=30, alpha=0.7, color=self.entropy_color, 
                        edgecolor='black', linewidth=0.5)
            
            # Add statistical markers
            mean_rate = np.mean(valid_rates)
            median_rate = np.median(valid_rates)
            axes[2].axvline(x=mean_rate, color='red', linestyle='-', 
                          linewidth=2, label=f'Mean: {mean_rate:.2e}')
            axes[2].axvline(x=median_rate, color='blue', linestyle='--', 
                          linewidth=2, label=f'Median: {median_rate:.2e}')
            axes[2].legend()
            axes[2].set_xlabel("Entropy Production Rate")
            axes[2].set_ylabel("Frequency")
        
        return fig, axes
    
    def create_bgk_comparison_plot(self, standard_data: HTheoremData, 
                                 entropic_data: HTheoremData) -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create comparative analysis between standard and entropic BGK methods."""
        plot_configs = [
            PlotConfig(PlotType.COMPARATIVE_ANALYSIS, "H-Function Comparison", 
                      "Timestep", "H-Function Value"),
            PlotConfig(PlotType.COMPARATIVE_ANALYSIS, "Entropy Production Comparison", 
                      "Timestep", "Entropy Production Rate"),
            PlotConfig(PlotType.COMPARATIVE_ANALYSIS, "Violation Analysis", 
                      "BGK Method", "Total Violations")
        ]
        
        fig, axes = self.layout_manager.create_multi_panel_layout(
            plot_configs, figure_title="Standard vs Entropic BGK Comparison"
        )
        
        # H-function comparison
        axes[0].semilogy(standard_data.timesteps, np.abs(standard_data.h_values), 
                        linewidth=2, label='Standard BGK', color='blue')
        axes[0].semilogy(entropic_data.timesteps, np.abs(entropic_data.h_values), 
                        linewidth=2, label='Entropic BGK', color='red')
        axes[0].legend()
        axes[0].set_ylabel("log |H-Function|")
        
        # Entropy production comparison
        axes[1].plot(standard_data.timesteps, standard_data.entropy_production_rates, 
                    linewidth=2, label='Standard BGK', color='blue', alpha=0.7)
        axes[1].plot(entropic_data.timesteps, entropic_data.entropy_production_rates, 
                    linewidth=2, label='Entropic BGK', color='red', alpha=0.7)
        axes[1].axhline(y=0, color='black', linestyle='--', alpha=0.5)
        axes[1].legend()
        axes[1].set_ylabel("Entropy Production Rate")
        
        # Violation comparison (bar plot)
        methods = ['Standard BGK', 'Entropic BGK']
        violations = [standard_data.violation_count, entropic_data.violation_count]
        bars = axes[2].bar(methods, violations, color=['blue', 'red'], alpha=0.7)
        axes[2].set_ylabel("Total Violations")
        
        # Add value labels on bars
        for bar, value in zip(bars, violations):
            height = bar.get_height()
            axes[2].text(bar.get_x() + bar.get_width()/2., height,
                        f'{value}', ha='center', va='bottom')
        
        return fig, axes
    
    def create_diagnostic_report(self, data: HTheoremData, 
                               output_path: Optional[Union[str, Path]] = None) -> str:
        """Generate comprehensive diagnostic report."""
        anomalies = self.detect_anomalies(data)
        
        report = f"""
=== H-Theorem Analysis Diagnostic Report ===

Basic Statistics:
- Total timesteps analyzed: {len(data.timesteps)}
- Timestep range: {data.timesteps[0]} to {data.timesteps[-1]}
- Initial H-value: {data.h_values[0]:.6e}
- Final H-value: {data.h_values[-1]:.6e}
- Total H-function change: {data.h_values[-1] - data.h_values[0]:.6e}

Monotonicity Analysis:
- Total violations: {data.violation_count}
- Violation rate: {data.violation_count/len(data.timesteps)*100:.2f}%
- Overall monotonic: {'YES' if data.overall_monotonic else 'NO'}

Entropy Production:
- Average entropy production rate: {data.average_entropy_production:.6e}
- Decay rate (exponential fit): {data.decay_rate:.6e}

Anomaly Detection:
- Total anomalies detected: {len(anomalies.anomaly_timesteps)}
- Monotonicity violations: {sum(1 for t in anomalies.anomaly_types if t == 'violation')}
- Statistical spikes: {sum(1 for t in anomalies.anomaly_types if t == 'spike')}
- Plateau regions: {sum(1 for t in anomalies.anomaly_types if t == 'plateau')}
- Statistical outliers: {len(anomalies.statistical_outliers)}

Recommendations:
"""
        
        if data.violation_count > 0:
            report += "- CAUTION: Monotonicity violations detected. Check numerical stability.\n"
        else:
            report += "- GOOD: No monotonicity violations detected.\n"
        
        if data.decay_rate > 0:
            report += f"- H-function decays exponentially with rate {data.decay_rate:.2e}\n"
        else:
            report += "- WARNING: No clear exponential decay pattern detected.\n"
        
        if len(anomalies.statistical_outliers) > len(data.timesteps) * 0.05:
            report += "- WARNING: High number of statistical outliers detected.\n"
        
        if output_path:
            with open(output_path, 'w') as f:
                f.write(report)
        
        return report
    
    def create_comprehensive_h_analysis(self, data: HTheoremData,
                                      output_directory: Optional[Union[str, Path]] = None) -> Dict[str, plt.Figure]:
        """Create comprehensive H-theorem analysis with all visualizations."""
        figures = {}
        
        # H-function evolution
        fig_evolution, _ = self.create_h_evolution_plot(data)
        figures['h_evolution'] = fig_evolution
        
        # Statistical analysis (if data available)
        if data.h_variances is not None:
            fig_stats, _ = self.create_statistical_analysis_plot(data)
            figures['statistical_analysis'] = fig_stats
        
        # Entropy production analysis
        fig_entropy, _ = self.create_entropy_production_analysis(data)
        figures['entropy_production'] = fig_entropy
        
        # Generate diagnostic report
        if output_directory:
            output_dir = Path(output_directory)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            report_path = output_dir / "h_theorem_diagnostic_report.txt"
            self.create_diagnostic_report(data, report_path)
            
            # Save figures
            for name, fig in figures.items():
                filename = output_dir / f"h_theorem_{name}.png"
                self.layout_manager.save_figure(fig, str(filename), dpi=300)
        
        return figures


def main():
    """Example usage of H-theorem visualization."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate H-theorem visualizations')
    parser.add_argument('h_evolution_file', help='H-theorem evolution CSV file')
    parser.add_argument('--output', help='Output directory for plots and report')
    parser.add_argument('--compare', help='Second H-evolution file for comparison')
    
    args = parser.parse_args()
    
    # Create visualizer
    visualizer = HTheoremPlots()
    
    # Load data
    data = visualizer.load_h_theorem_data(args.h_evolution_file)
    
    if args.compare:
        # Comparison mode
        data2 = visualizer.load_h_theorem_data(args.compare)
        fig_comp, _ = visualizer.create_bgk_comparison_plot(data, data2)
        
        if args.output:
            output_dir = Path(args.output)
            output_dir.mkdir(parents=True, exist_ok=True)
            visualizer.layout_manager.save_figure(
                fig_comp, str(output_dir / "bgk_comparison.png"), dpi=300
            )
        else:
            plt.show()
    else:
        # Single analysis mode
        figures = visualizer.create_comprehensive_h_analysis(data, args.output)
        print(f"Generated {len(figures)} H-theorem visualization figures")
        
        if not args.output:
            plt.show()


if __name__ == "__main__":
    main()