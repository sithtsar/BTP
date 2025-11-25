
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys

# Set publication-quality plotting defaults
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

def smooth_2d(array, kernel_size=3):
    """Simple 2D mean smoothing filter"""
    kernel = np.ones((kernel_size, kernel_size)) / (kernel_size * kernel_size)
    return np.convolve(array.flatten(), kernel.flatten(), mode='same').reshape(array.shape)

def gaussian_kernel(size, sigma=1.0):
    """Create a 2D Gaussian kernel"""
    kernel = np.zeros((size, size))
    center = size // 2

    for i in range(size):
        for j in range(size):
            x = i - center
            y = j - center
            kernel[i, j] = np.exp(-(x**2 + y**2) / (2 * sigma**2))

    return kernel / kernel.sum()

def gaussian_filter_masked(array, sigma=1.0, solid_mask=None, kernel_size=5):
    """Apply Gaussian filter only to fluid regions, preserving solid values"""
    if solid_mask is None or np.all(solid_mask):
        return array.copy()

    # Create a copy to work with
    filtered = array.copy()

    # Get kernel
    kernel = gaussian_kernel(kernel_size, sigma)
    k_center = kernel_size // 2

    # Get array dimensions
    ny, nx = array.shape

    # Only process fluid regions (where solid_mask is False)
    fluid_mask = ~solid_mask

    # Create temporary array for convolution result
    temp_filtered = np.zeros_like(array)

    # Apply 2D convolution only to fluid cells
    for i in range(ny):
        for j in range(nx):
            if not fluid_mask[i, j]:
                continue  # Skip solid cells

            # Accumulate weighted sum from neighboring cells
            weighted_sum = 0.0
            weight_sum = 0.0

            for ki in range(kernel_size):
                for kj in range(kernel_size):
                    ni = i + ki - k_center
                    nj = j + kj - k_center

                    # Check bounds
                    if 0 <= ni < ny and 0 <= nj < nx:
                        # Only include fluid cells in the convolution
                        if fluid_mask[ni, nj]:
                            weight = kernel[ki, kj]
                            weighted_sum += array[ni, nj] * weight
                            weight_sum += weight

            # Normalize by actual weight sum (accounts for boundary effects)
            if weight_sum > 0:
                temp_filtered[i, j] = weighted_sum / weight_sum
            else:
                temp_filtered[i, j] = array[i, j]  # Fallback if no valid neighbors

    # Update only fluid regions
    filtered[fluid_mask] = temp_filtered[fluid_mask]

    return filtered

def load_data_safe(filepath):
    """Load CSV data file, return data even if NaN (for showing instability)"""
    try:
        data = np.loadtxt(filepath, delimiter=',', skiprows=1)  # Skip header
        if np.any(np.isnan(data)):
            print(f"  ⚠️  {filepath.name} contains NaN (expected for unstable cases)")
        return data
    except:
        print(f"  ⚠️  Skipping {filepath.name} (cannot load)")
        return None

def reshape_csv_data(data, nx, ny):
    """Helper function to reshape CSV data to 2D arrays"""
    x = data[:, 0].astype(int)
    y = data[:, 1].astype(int)
    ux = data[:, 2]
    uy = data[:, 3]
    speed = data[:, 4]
    vorticity = data[:, 5]

    # Create 2D arrays, initialize to NaN for missing data
    ux_2d = np.full((ny, nx), np.nan)
    uy_2d = np.full((ny, nx), np.nan)
    speed_2d = np.full((ny, nx), np.nan)
    vort_2d = np.full((ny, nx), np.nan)

    # Fill fluid data
    for i in range(len(data)):
        xi, yi = x[i], y[i]
        ux_2d[yi, xi] = ux[i]
        uy_2d[yi, xi] = uy[i]
        speed_2d[yi, xi] = speed[i]
        vort_2d[yi, xi] = vorticity[i]

    # Create solid mask based on geometry
    solid_2d = np.zeros((ny, nx), dtype=bool)
    cx, cy, radius = 50, 25, 5  # Fixed for this simulation
    for i in range(nx):
        for j in range(ny):
            dx = i - cx
            dy = j - cy
            if dx*dx + dy*dy <= radius*radius or j == 0 or j == ny-1:
                solid_2d[j, i] = True

    return ux_2d, uy_2d, speed_2d, vort_2d, solid_2d

def plot_channel_cylinder(Re):
    """Plot Channel Flow with Cylinder - BGK vs ELBM with vorticity"""
    print(f"\n📊 Plotting Channel Flow with Cylinder for Re = {Re}...")

    # Load CSV data (format: x,y,u,v,speed,vorticity)
    bgk_data = load_data_safe(Path(f'LBM_output_Re{Re}_BGK.csv'))
    elbm_data = load_data_safe(Path(f'LBM_output_Re{Re}_ELBM.csv'))

    if bgk_data is None and elbm_data is None:
        print(f"  ❌ No valid data for Re = {Re}")
        return

    # Fixed grid dimensions
    nx, ny = 200, 50

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True, sharey=True)
    fig.suptitle(f'Channel Flow with Cylinder - Re = {Re}', fontsize=16)

    # BGK Plots
    if bgk_data is not None:
        ux, uy, speed, vort, solid = reshape_csv_data(bgk_data, nx, ny)

        # Velocity Magnitude
        ax = axes[0, 0]
        speed_masked = np.ma.masked_where(solid, speed)
        im = ax.imshow(speed_masked, extent=[0, nx, 0, ny], origin='lower', cmap='viridis', aspect='equal')
        ax.set_title('BGK - Velocity Magnitude')
        ax.set_ylabel('y')
        fig.colorbar(im, ax=ax)

        # Vorticity
        ax = axes[0, 1]
        vort_masked = np.ma.masked_where(solid, vort)
        vmax = max(abs(np.min(vort_masked)), abs(np.max(vort_masked)))
        im = ax.imshow(vort_masked, extent=[0, nx, 0, ny], origin='lower',
                      cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
        ax.set_title('BGK - Vorticity')
        fig.colorbar(im, ax=ax)



    # ELBM Plots
    if elbm_data is not None:
        ux, uy, speed, vort, solid = reshape_csv_data(elbm_data, nx, ny)

        # Velocity Magnitude
        ax = axes[1, 0]
        speed_masked = np.ma.masked_where(solid, speed)
        im = ax.imshow(speed_masked, extent=[0, nx, 0, ny], origin='lower', cmap='viridis', aspect='equal')
        ax.set_title('ELBM - Velocity Magnitude')
        ax.set_ylabel('y')
        ax.set_xlabel('x')
        fig.colorbar(im, ax=ax)

        # Vorticity
        ax = axes[1, 1]
        vort_masked = np.ma.masked_where(solid, vort)
        vmax = max(abs(np.min(vort_masked)), abs(np.max(vort_masked)))
        im = ax.imshow(vort_masked, extent=[0, nx, 0, ny], origin='lower',
                      cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
        ax.set_title('ELBM - Vorticity')
        ax.set_xlabel('x')
        fig.colorbar(im, ax=ax)



    plt.tight_layout()

    outdir = Path('figures/channel_cylinder')
    outdir.mkdir(parents=True, exist_ok=True)
    filename = outdir / f're_{Re}_comparison.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {filename}")
    plt.close()

def plot_streamlines_separate(Re):
    """Plot separate streamlines for BGK and ELBM with proper smoothing"""
    print(f"\n🌊 Plotting Separate Streamlines for Re = {Re}...")

    # Load CSV data
    bgk_data = load_data_safe(Path(f'LBM_output_Re{Re}_BGK.csv'))
    elbm_data = load_data_safe(Path(f'LBM_output_Re{Re}_ELBM.csv'))

    if bgk_data is None and elbm_data is None:
        print(f"  ❌ No valid data for Re = {Re}")
        return

    # Fixed grid dimensions
    nx, ny = 200, 50

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle(f'Streamlines - Channel Flow with Cylinder (Re = {Re})', fontsize=16)

    # BGK Streamlines
    if bgk_data is not None:
        ux, uy, speed, vort, solid = reshape_csv_data(bgk_data, nx, ny)

        # No smoothing applied to velocity fields
        ux_smooth = ux
        uy_smooth = uy
        speed_smooth = speed

        ax = axes[0]
        y_coords, x_coords = np.mgrid[0:ny, 0:nx]
        ux_masked = np.ma.masked_where(solid, ux_smooth)
        uy_masked = np.ma.masked_where(solid, uy_smooth)
        speed_masked = np.ma.masked_where(solid, speed_smooth)

        strm = ax.streamplot(x_coords, y_coords, ux_masked, uy_masked,
                           density=2.0, linewidth=1.2, color=speed_masked,
                           cmap='plasma', arrowsize=1.5, minlength=0.1)

        # Add cylinder
        cx, cy, radius = 50, 25, 5
        from matplotlib.patches import Circle
        circle = Circle((cx, cy), radius, color='black', fill=True, alpha=0.9, zorder=5)
        ax.add_patch(circle)

        ax.set_title('BGK - Streamlines', fontsize=14)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        fig.colorbar(strm.lines, ax=ax, label='Velocity Magnitude')

    # ELBM Streamlines
    if elbm_data is not None:
        ux, uy, speed, vort, solid = reshape_csv_data(elbm_data, nx, ny)

        # No smoothing applied to velocity fields
        ux_smooth = ux
        uy_smooth = uy
        speed_smooth = speed

        ax = axes[1]
        y_coords, x_coords = np.mgrid[0:ny, 0:nx]
        ux_masked = np.ma.masked_where(solid, ux_smooth)
        uy_masked = np.ma.masked_where(solid, uy_smooth)
        speed_masked = np.ma.masked_where(solid, speed_smooth)

        strm = ax.streamplot(x_coords, y_coords, ux_masked, uy_masked,
                           density=2.0, linewidth=1.2, color=speed_masked,
                           cmap='plasma', arrowsize=1.5, minlength=0.1)

        # Add cylinder
        cx, cy, radius = 50, 25, 5
        from matplotlib.patches import Circle
        circle = Circle((cx, cy), radius, color='black', fill=True, alpha=0.9, zorder=5)
        ax.add_patch(circle)

        ax.set_title('ELBM - Streamlines', fontsize=14)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        fig.colorbar(strm.lines, ax=ax, label='Velocity Magnitude')

    plt.tight_layout()

    outdir = Path('figures/channel_cylinder')
    outdir.mkdir(parents=True, exist_ok=True)
    filename = outdir / f're_{Re}_streamlines_separate.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"  ✅ Saved: {filename}")
    plt.close()




def test_gaussian_masking():
    """Test that Gaussian filtering preserves solid regions"""
    print("\n🧪 Testing Gaussian masking...")

    # Create test array with known pattern
    test_array = np.zeros((12, 12))
    # Create a step function - high values on left, low on right
    test_array[:, :6] = 10.0
    test_array[:, 6:] = 2.0
    # Solid region in center
    test_array[5:7, 5:7] = -999

    # Create solid mask
    solid_mask = np.zeros((12, 12), dtype=bool)
    solid_mask[5:7, 5:7] = True

    print("  Original array (showing smoothing effect):")
    print(test_array.astype(int))

    # Apply filtering
    filtered = gaussian_filter_masked(test_array, sigma=1.0, solid_mask=solid_mask)

    print("\n  Filtered array:")
    print(filtered.astype(int))

    # Check that solid regions are unchanged
    solid_preserved = np.allclose(filtered[solid_mask], test_array[solid_mask])

    # Check that fluid regions have been smoothed (step function should be blurred)
    fluid_original = test_array[~solid_mask]
    fluid_filtered = filtered[~solid_mask]
    fluid_changed = not np.allclose(fluid_original, fluid_filtered, rtol=1e-10)

    # Check smoothing effect - values near the step should be intermediate
    step_region = filtered[4:8, 4:8]  # Around the solid region
    has_intermediate_values = np.any((step_region > 3) & (step_region < 9))

    print(f"\n  Solid regions preserved: {solid_preserved}")
    print(f"  Fluid regions smoothed: {fluid_changed}")
    print(f"  Smoothing effect visible: {has_intermediate_values}")

    # Show values around the solid region
    print(f"  Values around solid region: {step_region.flatten()[:10]}")

    success = solid_preserved and fluid_changed and has_intermediate_values

    if success:
        print("  ✅ Gaussian masking working correctly")
    else:
        print("  ❌ Gaussian masking has issues")

    return success

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        test_gaussian_masking()
    else:
        if len(sys.argv) > 1:
            re_numbers = [int(re) for re in sys.argv[1:] if re != "test"]
        else:
            re_numbers = [10, 100, 1000]  # Match the test_benchmark.cpp Reynolds numbers

        for re in re_numbers:
            plot_channel_cylinder(re)
            plot_streamlines_separate(re)
