"""
Layout Manager for LBM Visualization System

This module provides organized plot arrangements with consistent styling
and adaptive layout algorithms for different plot combinations.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from enum import Enum


class PlotType(Enum):
    """Enumeration of different plot types for layout organization."""
    VELOCITY_PROFILE = "velocity_profile"
    H_FUNCTION_EVOLUTION = "h_function_evolution"
    ERROR_ANALYSIS = "error_analysis"
    FIELD_2D = "field_2d"
    STREAMLINES = "streamlines"
    PHASE_DISTRIBUTION = "phase_distribution"
    STABILITY_INDICATORS = "stability_indicators"
    COMPARATIVE_ANALYSIS = "comparative_analysis"


@dataclass
class PlotConfig:
    """Configuration for individual plots."""
    plot_type: PlotType
    title: str
    xlabel: str = ""
    ylabel: str = ""
    colorbar: bool = False
    aspect_ratio: Optional[float] = None
    grid: bool = True


@dataclass
class LayoutConfig:
    """Configuration for overall layout."""
    figure_width: float = 12.0
    figure_height: float = 8.0
    dpi: int = 300
    font_size: int = 10
    title_font_size: int = 12
    label_font_size: int = 10
    legend_font_size: int = 9
    tight_layout: bool = True
    subplot_spacing: float = 0.3


class FigureSizeCalculator:
    """Calculates optimal figure sizes based on plot content."""
    
    def __init__(self):
        self.base_width = 4.0
        self.base_height = 3.0
        self.aspect_ratios = {
            PlotType.VELOCITY_PROFILE: 1.2,
            PlotType.H_FUNCTION_EVOLUTION: 1.5,
            PlotType.ERROR_ANALYSIS: 1.3,
            PlotType.FIELD_2D: 1.0,
            PlotType.STREAMLINES: 1.0,
            PlotType.PHASE_DISTRIBUTION: 1.0,
            PlotType.STABILITY_INDICATORS: 1.4,
            PlotType.COMPARATIVE_ANALYSIS: 1.6
        }
    
    def calculate_subplot_size(self, plot_type: PlotType) -> Tuple[float, float]:
        """Calculate optimal size for a single subplot."""
        aspect = self.aspect_ratios.get(plot_type, 1.2)
        width = self.base_width
        height = self.base_height / aspect
        return width, height
    
    def calculate_figure_size(self, plot_configs: List[PlotConfig], 
                            layout_shape: Tuple[int, int]) -> Tuple[float, float]:
        """Calculate optimal figure size for multiple subplots."""
        rows, cols = layout_shape
        
        # Calculate average subplot dimensions
        total_width = 0
        total_height = 0
        
        for config in plot_configs:
            w, h = self.calculate_subplot_size(config.plot_type)
            total_width += w
            total_height += h
        
        avg_width = total_width / len(plot_configs)
        avg_height = total_height / len(plot_configs)
        
        # Account for spacing and margins
        figure_width = cols * avg_width + (cols - 1) * 0.5 + 2.0
        figure_height = rows * avg_height + (rows - 1) * 0.5 + 1.5
        
        return figure_width, figure_height


class SubplotOrganizer:
    """Organizes subplots in optimal arrangements."""
    
    def __init__(self):
        self.preferred_arrangements = {
            1: [(1, 1)],
            2: [(1, 2), (2, 1)],
            3: [(1, 3), (3, 1)],
            4: [(2, 2), (1, 4), (4, 1)],
            5: [(2, 3), (3, 2)],
            6: [(2, 3), (3, 2)],
            7: [(3, 3)],
            8: [(3, 3)],
            9: [(3, 3)]
        }
    
    def determine_layout_shape(self, num_plots: int, 
                             preferred_aspect: Optional[float] = None) -> Tuple[int, int]:
        """Determine optimal grid shape for given number of plots."""
        if num_plots <= 0:
            return (1, 1)
        
        if num_plots in self.preferred_arrangements:
            arrangements = self.preferred_arrangements[num_plots]
            
            if preferred_aspect is not None:
                # Choose arrangement closest to preferred aspect ratio
                best_arrangement = arrangements[0]
                best_diff = float('inf')
                
                for rows, cols in arrangements:
                    aspect = cols / rows
                    diff = abs(aspect - preferred_aspect)
                    if diff < best_diff:
                        best_diff = diff
                        best_arrangement = (rows, cols)
                
                return best_arrangement
            else:
                return arrangements[0]
        
        # For larger numbers, calculate square-ish arrangement
        rows = int(np.ceil(np.sqrt(num_plots)))
        cols = int(np.ceil(num_plots / rows))
        return (rows, cols)
    
    def create_subplot_grid(self, fig: plt.Figure, plot_configs: List[PlotConfig],
                          layout_shape: Optional[Tuple[int, int]] = None) -> List[plt.Axes]:
        """Create organized subplot grid."""
        num_plots = len(plot_configs)
        
        if layout_shape is None:
            layout_shape = self.determine_layout_shape(num_plots)
        
        rows, cols = layout_shape
        axes = []
        
        for i, config in enumerate(plot_configs):
            if i >= rows * cols:
                break
                
            ax = fig.add_subplot(rows, cols, i + 1)
            
            # Apply plot-specific configurations
            if config.aspect_ratio is not None:
                ax.set_aspect(config.aspect_ratio)
            
            if config.grid:
                ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            
            axes.append(ax)
        
        return axes


class StyleManager:
    """Manages consistent styling across all plots."""
    
    def __init__(self, config: LayoutConfig):
        self.config = config
        self._setup_matplotlib_style()
    
    def _setup_matplotlib_style(self):
        """Configure matplotlib with scientific formatting standards."""
        # Set publication-ready defaults
        rcParams.update({
            'font.size': self.config.font_size,
            'axes.titlesize': self.config.title_font_size,
            'axes.labelsize': self.config.label_font_size,
            'xtick.labelsize': self.config.font_size - 1,
            'ytick.labelsize': self.config.font_size - 1,
            'legend.fontsize': self.config.legend_font_size,
            'figure.titlesize': self.config.title_font_size + 2,
            'font.family': 'serif',
            'font.serif': ['Times New Roman', 'DejaVu Serif'],
            'mathtext.fontset': 'stix',
            'axes.linewidth': 1.0,
            'axes.spines.top': False,
            'axes.spines.right': False,
            'axes.grid': True,
            'grid.alpha': 0.3,
            'grid.linewidth': 0.5,
            'lines.linewidth': 1.5,
            'lines.markersize': 4,
            'figure.dpi': self.config.dpi,
            'savefig.dpi': self.config.dpi,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1
        })
    
    def apply_scientific_formatting(self, ax: plt.Axes, config: PlotConfig):
        """Apply scientific formatting to a single axis."""
        # Set title and labels
        if config.title:
            ax.set_title(config.title, fontweight='bold', pad=10)
        
        if config.xlabel:
            ax.set_xlabel(config.xlabel, fontweight='normal')
        
        if config.ylabel:
            ax.set_ylabel(config.ylabel, fontweight='normal')
        
        # Format tick labels for scientific notation when appropriate
        ax.ticklabel_format(style='scientific', axis='both', scilimits=(-3, 3))
        
        # Ensure proper spacing
        ax.margins(0.02)
    
    def create_colorbar(self, ax: plt.Axes, mappable, label: str = "") -> plt.Axes:
        """Create properly formatted colorbar."""
        cbar = plt.colorbar(mappable, ax=ax, shrink=0.8, aspect=20)
        if label:
            cbar.set_label(label, rotation=270, labelpad=15)
        cbar.ax.tick_params(labelsize=self.config.font_size - 1)
        return cbar


class LayoutManager:
    """Main layout manager for organized plot arrangements."""
    
    def __init__(self, config: Optional[LayoutConfig] = None):
        self.config = config or LayoutConfig()
        self.figure_calculator = FigureSizeCalculator()
        self.subplot_organizer = SubplotOrganizer()
        self.style_manager = StyleManager(self.config)
    
    def create_multi_panel_layout(self, plot_configs: List[PlotConfig],
                                layout_shape: Optional[Tuple[int, int]] = None,
                                figure_title: str = "") -> Tuple[plt.Figure, List[plt.Axes]]:
        """Create well-organized multi-panel layout."""
        if not plot_configs:
            raise ValueError("At least one plot configuration is required")
        
        # Determine layout shape
        if layout_shape is None:
            layout_shape = self.subplot_organizer.determine_layout_shape(len(plot_configs))
        
        # Calculate figure size
        fig_width, fig_height = self.figure_calculator.calculate_figure_size(
            plot_configs, layout_shape)
        
        # Create figure
        fig = plt.figure(figsize=(fig_width, fig_height), dpi=self.config.dpi)
        
        if figure_title:
            fig.suptitle(figure_title, fontsize=self.config.title_font_size + 2,
                        fontweight='bold', y=0.95)
        
        # Create subplot grid
        axes = self.subplot_organizer.create_subplot_grid(fig, plot_configs, layout_shape)
        
        # Apply styling to each subplot
        for ax, config in zip(axes, plot_configs):
            self.style_manager.apply_scientific_formatting(ax, config)
        
        # Apply tight layout
        if self.config.tight_layout:
            plt.tight_layout(rect=[0, 0, 1, 0.93] if figure_title else [0, 0, 1, 1])
        
        return fig, axes
    
    def apply_consistent_styling(self, fig: plt.Figure) -> None:
        """Apply consistent fonts, colors, and formatting to existing figure."""
        for ax in fig.get_axes():
            # Apply grid styling
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            
            # Format spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(1.0)
            ax.spines['bottom'].set_linewidth(1.0)
            
            # Format ticks
            ax.tick_params(direction='out', length=4, width=1)
    
    def save_figure(self, fig: plt.Figure, filename: str, 
                   format: str = 'png', **kwargs) -> None:
        """Save figure with consistent formatting."""
        default_kwargs = {
            'dpi': self.config.dpi,
            'bbox_inches': 'tight',
            'pad_inches': 0.1,
            'facecolor': 'white',
            'edgecolor': 'none'
        }
        default_kwargs.update(kwargs)
        
        fig.savefig(filename, format=format, **default_kwargs)
    
    def create_comparison_layout(self, left_config: PlotConfig, right_config: PlotConfig,
                               figure_title: str = "") -> Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]:
        """Create side-by-side comparison layout."""
        fig_width = 10.0
        fig_height = 4.5
        
        fig = plt.figure(figsize=(fig_width, fig_height), dpi=self.config.dpi)
        
        if figure_title:
            fig.suptitle(figure_title, fontsize=self.config.title_font_size + 2,
                        fontweight='bold', y=0.95)
        
        # Create two subplots side by side
        ax_left = fig.add_subplot(1, 2, 1)
        ax_right = fig.add_subplot(1, 2, 2)
        
        # Apply styling
        self.style_manager.apply_scientific_formatting(ax_left, left_config)
        self.style_manager.apply_scientific_formatting(ax_right, right_config)
        
        if self.config.tight_layout:
            plt.tight_layout(rect=[0, 0, 1, 0.93] if figure_title else [0, 0, 1, 1])
        
        return fig, (ax_left, ax_right)