"""
Visualization package for LBM simulations.

This package provides comprehensive visualization capabilities for both
single-phase and multiphase LBM simulations, including H-theorem analysis.
"""

from .layout_manager import LayoutManager, PlotConfig, PlotType, LayoutConfig
from .single_phase_plots import SinglePhasePlots, SinglePhaseData
from .multiphase_plots import MultiphasePlots, MultiphaseData
from .h_theorem_plots import HTheoremPlots, HTheoremData

__all__ = [
    'LayoutManager', 'PlotConfig', 'PlotType', 'LayoutConfig',
    'SinglePhasePlots', 'SinglePhaseData',
    'MultiphasePlots', 'MultiphaseData', 
    'HTheoremPlots', 'HTheoremData'
]