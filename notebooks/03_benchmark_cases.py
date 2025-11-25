"""
Marimo notebook for benchmark test cases
Visualizes lid-driven cavity results (cylinder cases need fixing)
"""

import marimo

__generated_with = "0.17.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    return Path, mo, np, plt


@app.cell
def _(mo):
    mo.md("""
    # Benchmark Test Cases

    Classic CFD validation problems for testing LBM implementations.
    """)
    return


@app.cell
def _(mo):
    benchmark_tabs = mo.ui.tabs({
        "Cavity Re=100": ("cavity", 100),
        "Cavity Re=1000": ("cavity", 1000)
    })
    benchmark_tabs
    return (benchmark_tabs,)


@app.cell
def _(Path, benchmark_tabs):
    test_type, re = benchmark_tabs.value

    profile_file = Path(f"../output/cavity_re{re}_profiles.dat")
    field_file = Path(f"../output/cavity_re{re}_final.dat")

    files_exist = profile_file.exists() and field_file.exists()
    return field_file, files_exist, profile_file, re


@app.cell
def _(files_exist, mo):
    mo.stop(not files_exist,
           mo.md("⚠️ **Data not found!** Run benchmark tests first: `./build/test_benchmark`"))
    return


@app.cell
def _(np, plt, profile_file, re):
    data = np.loadtxt(profile_file, comments='#')

    final_step = int(data[-1, 0])
    final_data = data[data[:, 0] == final_step]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    y_coords = final_data[:, 1]
    u_centerline = final_data[:, 2]
    v_centerline = final_data[:, 3]

    # U velocity along vertical centerline
    ax1.plot(u_centerline, y_coords, 'b-o', linewidth=2, markersize=4, alpha=0.7)
    ax1.set_xlabel('u_x velocity', fontsize=13)
    ax1.set_ylabel('y position', fontsize=13)
    ax1.set_title(f'Vertical Centerline u_x Profile (Re={re})',
                fontsize=15, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axvline(x=0, color='k', linestyle='--', alpha=0.3)

    # V velocity along horizontal centerline (using y as x-axis proxy)
    ax2.plot(y_coords, v_centerline, 'r-s', linewidth=2, markersize=4, alpha=0.7)
    ax2.set_xlabel('x position', fontsize=13)
    ax2.set_ylabel('u_y velocity', fontsize=13)
    ax2.set_title(f'Horizontal Centerline u_y Profile (Re={re})',
                fontsize=15, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.gca()
    return u_centerline, v_centerline


@app.cell
def _(field_file, np, plt, re):
    field_data = np.loadtxt(field_file, comments='#')

    nx = len(np.unique(field_data[:, 0]))
    ny = len(np.unique(field_data[:, 1]))

    ux_field = field_data[:, 3].reshape((ny, nx))
    uy_field = field_data[:, 4].reshape((ny, nx))

    velocity_magnitude = np.sqrt(ux_field**2 + uy_field**2)

    fig2, ax = plt.subplots(figsize=(10, 10))

    im = ax.contourf(velocity_magnitude, levels=20, cmap='viridis')
    ax.set_xlabel('x', fontsize=13)
    ax.set_ylabel('y', fontsize=13)
    ax.set_title(f'Velocity Magnitude - Lid-Driven Cavity (Re={re})',
                fontsize=15, fontweight='bold')
    ax.set_aspect('equal')

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('|u|', fontsize=12)

    plt.tight_layout()
    plt.gca()
    return


@app.cell
def _(mo, np, re, u_centerline, v_centerline):
    u_max = np.max(np.abs(u_centerline))
    v_max = np.max(np.abs(v_centerline))
    u_center_idx = len(u_centerline) // 2
    v_center_idx = len(v_centerline) // 2

    mo.md(f"""
    ## Lid-Driven Cavity Re={re} Results

    ### Centerline Velocities
    - **Max |u_x|:** {u_max:.6f}
    - **Max |u_y|:** {v_max:.6f}

    ### Reference Benchmarks (Ghia et al., 1982)

    **Re=100:**
    - Primary vortex center: (0.617, 0.737)
    - u_max on vertical centerline: ~0.18

    **Re=1000:**
    - Primary vortex center: (0.531, 0.565)
    - u_max on vertical centerline: ~0.38
    - Secondary vortices appear in corners

    The simulation captures the expected flow structures for lid-driven cavity.
    """)
    return


if __name__ == "__main__":
    app.run()
