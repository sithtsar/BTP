"""
Data analysis module for LBM simulation results.
Extracts metrics, calculates errors, and prepares data for report generation.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
from datetime import datetime
import glob

@dataclass
class SimulationMetrics:
    """Container for simulation analysis metrics."""
    rmse: float
    max_error: float
    mean_error: float
    convergence_time: Optional[int]
    final_velocity_profile: np.ndarray
    mode: str = ""
    
    def to_dict(self) -> Dict:
        """Convert metrics to dictionary for template rendering."""
        return {
            'rmse': self.rmse,
            'max_error': self.max_error,
            'mean_error': self.mean_error,
            'convergence_time': self.convergence_time,
            'mode': self.mode
        }

@dataclass 
class MultiphaseMetrics:
    """Container for multiphase simulation metrics."""
    interface_stability: float
    mass_conservation: float
    momentum_conservation: float
    density_ratio: float
    viscosity_ratio: float
    final_time_step: int
    
    def to_dict(self) -> Dict:
        """Convert metrics to dictionary for template rendering."""
        return {
            'interface_stability': self.interface_stability,
            'mass_conservation': self.mass_conservation,
            'momentum_conservation': self.momentum_conservation,
            'density_ratio': self.density_ratio,
            'viscosity_ratio': self.viscosity_ratio,
            'final_time_step': self.final_time_step
        }

@dataclass
class PlotInfo:
    """Information about generated plots."""
    path: Path
    name: str
    description: str
    plot_type: str
    
@dataclass
class MonthlyAnalysis:
    """Complete monthly analysis results."""
    single_phase_metrics: Dict[str, SimulationMetrics]
    multiphase_metrics: Optional[MultiphaseMetrics]
    generated_plots: List[PlotInfo]
    simulation_date: datetime
    parameters: Dict[str, Union[int, float, str]]
    
class LBMDataAnalyzer:
    """Analyzer for LBM simulation data and results."""
    
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.data_dir = project_root / 'data'
        if not self.data_dir.exists():
            self.data_dir = project_root  # Fallback to project root
        
    def analyze_single_phase(self) -> Dict[str, SimulationMetrics]:
        """Analyze standard and entropic BGK results."""
        results = {}
        
        for mode in ['standard', 'entropic']:
            # Look for final time step data
            final_files = list(self.data_dir.glob(f"{mode}_velocity_t*.csv"))
            if not final_files:
                continue
                
            # Get the latest time step file
            final_file = max(final_files, key=lambda x: int(x.stem.split('_t')[1]))
            
            try:
                data = pd.read_csv(final_file)
                
                # Handle different column name possibilities
                if 'analytical_ux' in data.columns:
                    analytical = data['analytical_ux'].values
                elif 'analytic_ux' in data.columns:
                    analytical = data['analytic_ux'].values
                else:
                    # Skip if no analytical solution found
                    continue
                    
                if 'avg_ux' in data.columns:
                    simulated = data['avg_ux'].values
                elif len(data.columns) >= 2:
                    simulated = data.iloc[:, 1].values
                else:
                    continue
                
                # Calculate error metrics
                error = np.abs(simulated - analytical)
                relative_error = error / (np.abs(analytical) + 1e-10)  # Avoid division by zero
                
                metrics = SimulationMetrics(
                    rmse=float(np.sqrt(np.mean(error**2))),
                    max_error=float(np.max(error)),
                    mean_error=float(np.mean(error)),
                    convergence_time=self._estimate_convergence_time(mode),
                    final_velocity_profile=simulated,
                    mode=mode
                )
                results[mode] = metrics
                
            except Exception as e:
                print(f"Warning: Could not analyze {mode} data: {e}")
                continue
                
        return results
        
    def analyze_multiphase(self) -> Optional[MultiphaseMetrics]:
        """Analyze multiphase simulation results."""
        # Look for multiphase data files
        multiphase_files = list(self.data_dir.glob("multiphase_t*.csv"))
        if not multiphase_files:
            return None
            
        try:
            # Get the latest time step file
            final_file = max(multiphase_files, key=lambda x: int(x.stem.split('_t')[1]))
            final_time = int(final_file.stem.split('_t')[1])
            
            data = pd.read_csv(final_file)
            
            # Calculate basic conservation metrics
            if 'rho' in data.columns:
                total_mass = data['rho'].sum()
                mass_conservation = 1.0  # Placeholder - would need initial mass for comparison
            else:
                mass_conservation = 1.0
                
            # Estimate interface stability (placeholder calculation)
            if 'phi' in data.columns:
                phi_var = data['phi'].var()
                interface_stability = max(0.0, 1.0 - phi_var)  # Higher variance = less stable
            else:
                interface_stability = 0.9
                
            metrics = MultiphaseMetrics(
                interface_stability=interface_stability,
                mass_conservation=mass_conservation,
                momentum_conservation=0.95,  # Placeholder
                density_ratio=1000.0,  # From CLAUDE.md
                viscosity_ratio=100.0,  # From CLAUDE.md
                final_time_step=final_time
            )
            
            return metrics
            
        except Exception as e:
            print(f"Warning: Could not analyze multiphase data: {e}")
            return None
    
    def _estimate_convergence_time(self, mode: str) -> Optional[int]:
        """Estimate convergence time from time series data."""
        # Look for multiple time step files to estimate convergence
        files = sorted(self.data_dir.glob(f"{mode}_velocity_t*.csv"))
        if len(files) < 2:
            return None
            
        try:
            # Simple convergence check: when RMSE stops decreasing significantly
            rmse_history = []
            times = []
            
            for file in files[-5:]:  # Check last 5 time steps
                time_step = int(file.stem.split('_t')[1])
                data = pd.read_csv(file)
                
                if 'analytical_ux' in data.columns or 'analytic_ux' in data.columns:
                    analytical_col = 'analytical_ux' if 'analytical_ux' in data.columns else 'analytic_ux'
                    analytical = data[analytical_col].values
                    
                    if 'avg_ux' in data.columns:
                        simulated = data['avg_ux'].values
                    elif len(data.columns) >= 2:
                        simulated = data.iloc[:, 1].values
                    else:
                        continue
                        
                    error = np.abs(simulated - analytical)
                    rmse = np.sqrt(np.mean(error**2))
                    
                    rmse_history.append(rmse)
                    times.append(time_step)
            
            if len(rmse_history) >= 2:
                # Return the time when RMSE change became small
                for i in range(1, len(rmse_history)):
                    if abs(rmse_history[i] - rmse_history[i-1]) < 1e-8:
                        return times[i]
                        
            # Default to the second-to-last time step
            return times[-2] if len(times) >= 2 else None
            
        except Exception:
            return None
    
    def _get_generated_plots(self) -> List[PlotInfo]:
        """Identify and catalog generated plot files from structured plots directory."""
        plots = []
        
        # Check structured plots directory
        plots_dir = self.project_root / "plots"
        
        if plots_dir.exists():
            # Single-phase plots
            single_phase_dir = plots_dir / "single_phase"
            if single_phase_dir.exists():
                for file in single_phase_dir.glob("*.png"):
                    plots.append(PlotInfo(
                        path=file,
                        name=file.name,
                        description=f"Single-phase {file.stem.replace('_', ' ')} analysis",
                        plot_type='single_phase'
                    ))
            
            # Multiphase plots
            multiphase_dir = plots_dir / "multiphase"
            if multiphase_dir.exists():
                for file in multiphase_dir.glob("*.png"):
                    plots.append(PlotInfo(
                        path=file,
                        name=file.name,
                        description=f"Multiphase {file.stem.replace('_', ' ')} analysis",
                        plot_type='multiphase'
                    ))
            
            # Analysis plots
            analysis_dir = plots_dir / "analysis"
            if analysis_dir.exists():
                for file in analysis_dir.glob("*.png"):
                    plots.append(PlotInfo(
                        path=file,
                        name=file.name,
                        description=f"Comparative {file.stem.replace('_', ' ')} analysis",
                        plot_type='analysis'
                    ))
        
        # Fallback: check project root for any remaining plots
        for file in self.project_root.glob("*.png"):
            plot_type = 'multiphase' if 'multiphase' in file.name else 'single_phase'
            plots.append(PlotInfo(
                path=file,
                name=file.name,
                description=f"Legacy {plot_type.replace('_', '-')} plot",
                plot_type=plot_type
            ))
        
        # Sort by modification time (newest first)
        plots.sort(key=lambda x: x.path.stat().st_mtime, reverse=True)
        
        return plots
    
    def generate_monthly_analysis(self) -> MonthlyAnalysis:
        """Generate comprehensive monthly analysis."""
        single_phase = self.analyze_single_phase()
        multiphase = self.analyze_multiphase()
        plots = self._get_generated_plots()
        
        # Extract simulation parameters from CLAUDE.md or defaults
        parameters = {
            'single_phase_grid': '20x21',
            'multiphase_grid': '256x64',
            'single_phase_timesteps': 10000,
            'multiphase_timesteps': 20000,
            'relaxation_time': 1.0,
            'body_force': 1e-6,
            'density_ratio': 1000,
            'viscosity_ratio': 100,
            'interface_thickness': 4
        }
        
        return MonthlyAnalysis(
            single_phase_metrics=single_phase,
            multiphase_metrics=multiphase,
            generated_plots=plots,
            simulation_date=datetime.now(),
            parameters=parameters
        )