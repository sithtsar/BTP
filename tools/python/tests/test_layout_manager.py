"""
Unit tests for LayoutManager and related classes.

Tests layout calculation, styling functions, and adaptive layout algorithms.
"""

import pytest
import matplotlib.pyplot as plt
import numpy as np
from unittest.mock import patch, MagicMock
import tempfile
import os

# Add the parent directory to the path to import our modules
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from visualization.layout_manager import (
    LayoutManager, PlotConfig, LayoutConfig, PlotType,
    FigureSizeCalculator, SubplotOrganizer, StyleManager
)


class TestPlotConfig:
    """Test PlotConfig dataclass."""
    
    def test_plot_config_creation(self):
        """Test basic PlotConfig creation."""
        config = PlotConfig(
            plot_type=PlotType.VELOCITY_PROFILE,
            title="Test Plot",
            xlabel="X axis",
            ylabel="Y axis"
        )
        
        assert config.plot_type == PlotType.VELOCITY_PROFILE
        assert config.title == "Test Plot"
        assert config.xlabel == "X axis"
        assert config.ylabel == "Y axis"
        assert config.grid is True  # Default value
        assert config.colorbar is False  # Default value
    
    def test_plot_config_defaults(self):
        """Test PlotConfig with minimal parameters."""
        config = PlotConfig(
            plot_type=PlotType.H_FUNCTION_EVOLUTION,
            title="Minimal Config"
        )
        
        assert config.xlabel == ""
        assert config.ylabel == ""
        assert config.aspect_ratio is None
        assert config.grid is True


class TestLayoutConfig:
    """Test LayoutConfig dataclass."""
    
    def test_layout_config_defaults(self):
        """Test default LayoutConfig values."""
        config = LayoutConfig()
        
        assert config.figure_width == 12.0
        assert config.figure_height == 8.0
        assert config.dpi == 300
        assert config.font_size == 10
        assert config.tight_layout is True
    
    def test_layout_config_custom(self):
        """Test custom LayoutConfig values."""
        config = LayoutConfig(
            figure_width=15.0,
            dpi=150,
            font_size=12
        )
        
        assert config.figure_width == 15.0
        assert config.dpi == 150
        assert config.font_size == 12


class TestFigureSizeCalculator:
    """Test FigureSizeCalculator class."""
    
    def setUp(self):
        self.calculator = FigureSizeCalculator()
    
    def test_calculate_subplot_size(self):
        """Test subplot size calculation for different plot types."""
        calculator = FigureSizeCalculator()
        
        # Test velocity profile
        width, height = calculator.calculate_subplot_size(PlotType.VELOCITY_PROFILE)
        assert width == 4.0
        assert height == pytest.approx(2.5, rel=1e-2)  # 3.0 / 1.2
        
        # Test 2D field plot
        width, height = calculator.calculate_subplot_size(PlotType.FIELD_2D)
        assert width == 4.0
        assert height == pytest.approx(3.0, rel=1e-2)  # 3.0 / 1.0
    
    def test_calculate_figure_size(self):
        """Test figure size calculation for multiple subplots."""
        calculator = FigureSizeCalculator()
        
        configs = [
            PlotConfig(PlotType.VELOCITY_PROFILE, "Plot 1"),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "Plot 2")
        ]
        
        width, height = calculator.calculate_figure_size(configs, (1, 2))
        
        # Should be reasonable figure size
        assert width > 8.0  # At least 2 subplots wide plus margins
        assert height > 3.0  # At least 1 subplot high plus margins
        assert width < 20.0  # Not unreasonably large
        assert height < 15.0


class TestSubplotOrganizer:
    """Test SubplotOrganizer class."""
    
    def test_determine_layout_shape(self):
        """Test layout shape determination for different numbers of plots."""
        organizer = SubplotOrganizer()
        
        # Test specific cases
        assert organizer.determine_layout_shape(1) == (1, 1)
        assert organizer.determine_layout_shape(2) in [(1, 2), (2, 1)]
        assert organizer.determine_layout_shape(4) in [(2, 2), (1, 4), (4, 1)]
        
        # Test larger numbers
        rows, cols = organizer.determine_layout_shape(10)
        assert rows * cols >= 10
        assert abs(rows - cols) <= 1  # Should be roughly square
    
    def test_determine_layout_shape_with_aspect(self):
        """Test layout shape with preferred aspect ratio."""
        organizer = SubplotOrganizer()
        
        # Prefer wide layout (aspect > 1)
        rows, cols = organizer.determine_layout_shape(4, preferred_aspect=2.0)
        assert cols >= rows  # Should prefer wider layout
        
        # Prefer tall layout (aspect < 1)
        rows, cols = organizer.determine_layout_shape(4, preferred_aspect=0.5)
        assert rows >= cols  # Should prefer taller layout
    
    @patch('matplotlib.pyplot.Figure')
    def test_create_subplot_grid(self, mock_figure):
        """Test subplot grid creation."""
        organizer = SubplotOrganizer()
        
        # Create mock figure and axes
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_fig.add_subplot.side_effect = [mock_ax1, mock_ax2]
        
        configs = [
            PlotConfig(PlotType.VELOCITY_PROFILE, "Plot 1", grid=True),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "Plot 2", grid=False, aspect_ratio=1.5)
        ]
        
        axes = organizer.create_subplot_grid(mock_fig, configs, (1, 2))
        
        # Check that subplots were created
        assert len(axes) == 2
        assert mock_fig.add_subplot.call_count == 2
        
        # Check that aspect ratio was set for second plot
        mock_ax2.set_aspect.assert_called_once_with(1.5)


class TestStyleManager:
    """Test StyleManager class."""
    
    def test_style_manager_initialization(self):
        """Test StyleManager initialization with configuration."""
        config = LayoutConfig(font_size=12, dpi=150)
        style_manager = StyleManager(config)
        
        assert style_manager.config.font_size == 12
        assert style_manager.config.dpi == 150
    
    @patch('matplotlib.pyplot.Axes')
    def test_apply_scientific_formatting(self, mock_axes):
        """Test scientific formatting application."""
        config = LayoutConfig()
        style_manager = StyleManager(config)
        
        mock_ax = MagicMock()
        plot_config = PlotConfig(
            PlotType.VELOCITY_PROFILE,
            title="Test Title",
            xlabel="X Label",
            ylabel="Y Label"
        )
        
        style_manager.apply_scientific_formatting(mock_ax, plot_config)
        
        # Check that title and labels were set
        mock_ax.set_title.assert_called_once_with("Test Title", fontweight='bold', pad=10)
        mock_ax.set_xlabel.assert_called_once_with("X Label", fontweight='normal')
        mock_ax.set_ylabel.assert_called_once_with("Y Label", fontweight='normal')
        
        # Check that scientific notation was applied
        mock_ax.ticklabel_format.assert_called_once()
        mock_ax.margins.assert_called_once_with(0.02)
    
    @patch('matplotlib.pyplot.colorbar')
    def test_create_colorbar(self, mock_colorbar):
        """Test colorbar creation."""
        config = LayoutConfig()
        style_manager = StyleManager(config)
        
        mock_ax = MagicMock()
        mock_mappable = MagicMock()
        mock_cbar = MagicMock()
        mock_colorbar.return_value = mock_cbar
        
        result = style_manager.create_colorbar(mock_ax, mock_mappable, "Test Label")
        
        # Check colorbar was created with correct parameters
        mock_colorbar.assert_called_once_with(mock_mappable, ax=mock_ax, shrink=0.8, aspect=20)
        mock_cbar.set_label.assert_called_once_with("Test Label", rotation=270, labelpad=15)


class TestLayoutManager:
    """Test main LayoutManager class."""
    
    def test_layout_manager_initialization(self):
        """Test LayoutManager initialization."""
        # Test with default config
        manager = LayoutManager()
        assert manager.config is not None
        assert isinstance(manager.figure_calculator, FigureSizeCalculator)
        assert isinstance(manager.subplot_organizer, SubplotOrganizer)
        assert isinstance(manager.style_manager, StyleManager)
        
        # Test with custom config
        custom_config = LayoutConfig(font_size=14)
        manager = LayoutManager(custom_config)
        assert manager.config.font_size == 14
    
    def test_create_multi_panel_layout_empty_configs(self):
        """Test error handling for empty plot configurations."""
        manager = LayoutManager()
        
        with pytest.raises(ValueError, match="At least one plot configuration is required"):
            manager.create_multi_panel_layout([])
    
    @patch('matplotlib.pyplot.figure')
    @patch('matplotlib.pyplot.tight_layout')
    def test_create_multi_panel_layout(self, mock_tight_layout, mock_figure):
        """Test multi-panel layout creation."""
        manager = LayoutManager()
        
        # Create mock figure
        mock_fig = MagicMock()
        mock_figure.return_value = mock_fig
        mock_ax = MagicMock()
        mock_fig.add_subplot.return_value = mock_ax
        
        configs = [
            PlotConfig(PlotType.VELOCITY_PROFILE, "Plot 1"),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "Plot 2")
        ]
        
        fig, axes = manager.create_multi_panel_layout(configs, figure_title="Test Figure")
        
        # Check that figure was created
        mock_figure.assert_called_once()
        mock_fig.suptitle.assert_called_once()
        
        # Check that tight layout was applied
        mock_tight_layout.assert_called_once()
        
        assert fig == mock_fig
        assert len(axes) == 2
    
    @patch('matplotlib.pyplot.figure')
    def test_create_comparison_layout(self, mock_figure):
        """Test comparison layout creation."""
        manager = LayoutManager()
        
        # Create mock figure
        mock_fig = MagicMock()
        mock_figure.return_value = mock_fig
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_fig.add_subplot.side_effect = [mock_ax1, mock_ax2]
        
        left_config = PlotConfig(PlotType.VELOCITY_PROFILE, "Left Plot")
        right_config = PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "Right Plot")
        
        fig, (ax_left, ax_right) = manager.create_comparison_layout(
            left_config, right_config, "Comparison"
        )
        
        # Check that figure was created with correct size
        mock_figure.assert_called_once_with(figsize=(10.0, 4.5), dpi=300)
        mock_fig.suptitle.assert_called_once_with(
            "Comparison", fontsize=14, fontweight='bold', y=0.95
        )
        
        # Check that two subplots were created
        assert mock_fig.add_subplot.call_count == 2
        assert ax_left == mock_ax1
        assert ax_right == mock_ax2
    
    @patch('matplotlib.pyplot.Figure.savefig')
    def test_save_figure(self, mock_savefig):
        """Test figure saving with consistent formatting."""
        manager = LayoutManager()
        mock_fig = MagicMock()
        
        manager.save_figure(mock_fig, "test.png", format='png', custom_param='value')
        
        # Check that savefig was called with correct parameters
        mock_savefig.assert_called_once()
        args, kwargs = mock_savefig.call_args
        
        assert args[0] == "test.png"
        assert kwargs['format'] == 'png'
        assert kwargs['dpi'] == 300
        assert kwargs['bbox_inches'] == 'tight'
        assert kwargs['custom_param'] == 'value'
    
    def test_apply_consistent_styling(self):
        """Test consistent styling application to existing figure."""
        manager = LayoutManager()
        
        # Create mock figure with axes
        mock_fig = MagicMock()
        mock_ax1 = MagicMock()
        mock_ax2 = MagicMock()
        mock_fig.get_axes.return_value = [mock_ax1, mock_ax2]
        
        # Create mock spines
        for ax in [mock_ax1, mock_ax2]:
            ax.spines = {
                'top': MagicMock(),
                'right': MagicMock(),
                'left': MagicMock(),
                'bottom': MagicMock()
            }
        
        manager.apply_consistent_styling(mock_fig)
        
        # Check that styling was applied to all axes
        for ax in [mock_ax1, mock_ax2]:
            ax.grid.assert_called_once_with(True, alpha=0.3, linestyle='-', linewidth=0.5)
            ax.spines['top'].set_visible.assert_called_once_with(False)
            ax.spines['right'].set_visible.assert_called_once_with(False)
            ax.tick_params.assert_called_once_with(direction='out', length=4, width=1)


# Integration tests
class TestLayoutManagerIntegration:
    """Integration tests for LayoutManager with actual matplotlib."""
    
    def test_real_figure_creation(self):
        """Test creating actual matplotlib figures."""
        # Use Agg backend to avoid display issues in testing
        plt.switch_backend('Agg')
        
        manager = LayoutManager()
        configs = [
            PlotConfig(PlotType.VELOCITY_PROFILE, "Velocity Profile", "x", "u(x)"),
            PlotConfig(PlotType.H_FUNCTION_EVOLUTION, "H-function", "Time", "H(t)")
        ]
        
        fig, axes = manager.create_multi_panel_layout(configs)
        
        # Check that we got a real figure and axes
        assert isinstance(fig, plt.Figure)
        assert len(axes) == 2
        assert all(isinstance(ax, plt.Axes) for ax in axes)
        
        # Check that titles were set
        assert axes[0].get_title() == "Velocity Profile"
        assert axes[1].get_title() == "H-function"
        
        plt.close(fig)
    
    def test_save_real_figure(self):
        """Test saving actual figure to file."""
        plt.switch_backend('Agg')
        
        manager = LayoutManager()
        config = PlotConfig(PlotType.VELOCITY_PROFILE, "Test Plot")
        
        fig, axes = manager.create_multi_panel_layout([config])
        
        # Add some dummy data
        x = np.linspace(0, 1, 100)
        y = x**2
        axes[0].plot(x, y)
        
        # Save to temporary file
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            manager.save_figure(fig, tmp.name)
            
            # Check that file was created and has reasonable size
            assert os.path.exists(tmp.name)
            assert os.path.getsize(tmp.name) > 1000  # Should be at least 1KB
            
            # Clean up
            os.unlink(tmp.name)
        
        plt.close(fig)


if __name__ == '__main__':
    pytest.main([__file__])