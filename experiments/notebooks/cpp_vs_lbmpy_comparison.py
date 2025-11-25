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
    # C++ vs lbmpy Validation Comparison

    This notebook compares the analytical validation results between the C++ implementation and the Python lbmpy implementation.
    """)
    return


@app.cell
def _(Path, np):
    def load_data(filename):
        """Load data from output file."""
        filepath = Path("../../output") / filename
        if not filepath.exists():
            return None
        return np.loadtxt(filepath)

    # Load C++ results
    couette_cpp = load_data("couette_results.dat")
    poiseuille_cpp = load_data("poiseuille_results.dat")
    tg_cpp = load_data("tg_energy.dat")

    # Load lbmpy results
    couette_lbmpy = load_data("couette_lbmpy.dat")
    poiseuille_lbmpy = load_data("poiseuille_lbmpy.dat")
    tg_lbmpy = load_data("tg_lbmpy.dat")
    return (
        couette_cpp,
        couette_lbmpy,
        poiseuille_cpp,
        poiseuille_lbmpy,
        tg_cpp,
        tg_lbmpy,
    )


@app.cell
def _(mo):
    mo.md("""
    ## Couette Flow Comparison
    """)
    return


@app.cell
def _(couette_cpp, couette_lbmpy, np, plt):
    fig_couette, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    if couette_cpp is not None and couette_lbmpy is not None:
        # Left plot: Velocity profiles
        ax1.plot(couette_cpp[:, 1], couette_cpp[:, 0], 'o-', label='C++ simulation', alpha=0.7)
        ax1.plot(couette_lbmpy[:, 1], couette_lbmpy[:, 0], 's-', label='lbmpy simulation', alpha=0.7)
        ax1.plot(couette_cpp[:, 2], couette_cpp[:, 0], 'k--', label='Analytical', linewidth=2)
        ax1.set_xlabel('Velocity u')
        ax1.set_ylabel('y position')
        ax1.set_title('Couette Flow: Velocity Profiles')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Right plot: Error comparison
        cpp_error = np.abs(couette_cpp[:, 1] - couette_cpp[:, 2])
        lbmpy_error = np.abs(couette_lbmpy[:, 1] - couette_lbmpy[:, 2])

        ax2.semilogy(cpp_error, couette_cpp[:, 0], 'o-', label='C++ error', alpha=0.7)
        ax2.semilogy(lbmpy_error, couette_lbmpy[:, 0], 's-', label='lbmpy error', alpha=0.7)
        ax2.set_xlabel('Absolute Error')
        ax2.set_ylabel('y position')
        ax2.set_title('Couette Flow: Error Comparison')
        ax2.legend()
        ax2.grid(True, alpha=0.3, which='both')

        # Compute L2 errors
        cpp_l2 = np.sqrt(np.mean(cpp_error**2))
        lbmpy_l2 = np.sqrt(np.mean(lbmpy_error**2))

        fig_couette.suptitle(f'Couette Flow | C++ L2: {cpp_l2:.2e} | lbmpy L2: {lbmpy_l2:.2e}',
                           fontsize=12, fontweight='bold')
    else:
        ax1.text(0.5, 0.5, 'Data not available',
                ha='center', va='center', transform=ax1.transAxes)

    plt.tight_layout()
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md("""
    ## Poiseuille Flow Comparison
    """)
    return


@app.cell
def _(np, plt, poiseuille_cpp, poiseuille_lbmpy):
    fig_poiseuille, (ax3, ax4) = plt.subplots(1, 2, figsize=(14, 5))

    if poiseuille_cpp is not None and poiseuille_lbmpy is not None:
        # Left plot: Velocity profiles
        ax3.plot(poiseuille_cpp[:, 1], poiseuille_cpp[:, 0], 'o-',
                label='C++ simulation', alpha=0.7)
        ax3.plot(poiseuille_lbmpy[:, 1], poiseuille_lbmpy[:, 0], 's-',
                label='lbmpy simulation', alpha=0.7)
        ax3.plot(poiseuille_cpp[:, 2], poiseuille_cpp[:, 0], 'k--',
                label='Analytical', linewidth=2)
        ax3.set_xlabel('Velocity u')
        ax3.set_ylabel('y position')
        ax3.set_title('Poiseuille Flow: Velocity Profiles')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # Right plot: Error comparison
        cpp_error_p = np.abs(poiseuille_cpp[:, 1] - poiseuille_cpp[:, 2])
        lbmpy_error_p = np.abs(poiseuille_lbmpy[:, 1] - poiseuille_lbmpy[:, 2])

        ax4.semilogy(cpp_error_p, poiseuille_cpp[:, 0], 'o-',
                    label='C++ error', alpha=0.7)
        ax4.semilogy(lbmpy_error_p, poiseuille_lbmpy[:, 0], 's-',
                    label='lbmpy error', alpha=0.7)
        ax4.set_xlabel('Absolute Error')
        ax4.set_ylabel('y position')
        ax4.set_title('Poiseuille Flow: Error Comparison')
        ax4.legend()
        ax4.grid(True, alpha=0.3, which='both')

        # Compute L2 errors
        cpp_l2_p = np.sqrt(np.mean(cpp_error_p**2))
        lbmpy_l2_p = np.sqrt(np.mean(lbmpy_error_p**2))

        fig_poiseuille.suptitle(
            f'Poiseuille Flow | C++ L2: {cpp_l2_p:.2e} | lbmpy L2: {lbmpy_l2_p:.2e}',
            fontsize=12, fontweight='bold'
        )
    else:
        ax3.text(0.5, 0.5, 'Data not available',
                ha='center', va='center', transform=ax3.transAxes)

    plt.tight_layout()
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md("""
    ## Taylor-Green Vortex Energy Decay
    """)
    return


@app.cell
def _(np, plt, tg_cpp, tg_lbmpy):
    fig_tg, ax5 = plt.subplots(1, 1, figsize=(10, 6))

    if tg_cpp is not None and tg_lbmpy is not None:
        # Plot energy decay
        ax5.semilogy(tg_cpp[:, 1], tg_cpp[:, 2], 'o-',
                    label='C++ simulation', markersize=6, alpha=0.7)
        ax5.semilogy(tg_lbmpy[:, 1], tg_lbmpy[:, 2], 's-',
                    label='lbmpy simulation', markersize=6, alpha=0.7)

        # Analytical decay (if initial energy is available)
        if len(tg_cpp) > 0:
            E0_cpp = tg_cpp[0, 2]
            t = tg_cpp[:, 1]
            nu = 0.01  # From test configuration
            k = 2.0 * np.pi / 64
            E_analytical = E0_cpp * np.exp(-2.0 * nu * k**2 * t)
            ax5.semilogy(t, E_analytical, 'k--',
                        label='Analytical decay', linewidth=2)

        ax5.set_xlabel('Time steps')
        ax5.set_ylabel('Kinetic Energy')
        ax5.set_title('Taylor-Green Vortex: Energy Decay')
        ax5.legend()
        ax5.grid(True, alpha=0.3, which='both')
    else:
        ax5.text(0.5, 0.5, 'Data not available',
                ha='center', va='center', transform=ax5.transAxes)

    plt.tight_layout()
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md("""
    ## Summary

    This comparison validates that the C++ and lbmpy implementations produce consistent results for standard LBM test cases.
    """)
    return


@app.cell
def _(couette_cpp, couette_lbmpy, mo, np, poiseuille_cpp, poiseuille_lbmpy):
    summary_data = {}

    if couette_cpp is not None and couette_lbmpy is not None:
        cpp_err_c = np.sqrt(np.mean((couette_cpp[:, 1] - couette_cpp[:, 2])**2))
        lbmpy_err_c = np.sqrt(np.mean((couette_lbmpy[:, 1] - couette_lbmpy[:, 2])**2))
        diff_c = np.sqrt(np.mean((couette_cpp[:, 1] - couette_lbmpy[:, 1])**2))
        summary_data['Couette'] = {
            'C++ L2 Error': f'{cpp_err_c:.2e}',
            'lbmpy L2 Error': f'{lbmpy_err_c:.2e}',
            'Difference': f'{diff_c:.2e}'
        }

    if poiseuille_cpp is not None and poiseuille_lbmpy is not None:
        cpp_err_p = np.sqrt(np.mean((poiseuille_cpp[:, 1] - poiseuille_cpp[:, 2])**2))
        lbmpy_err_p = np.sqrt(np.mean((poiseuille_lbmpy[:, 1] - poiseuille_lbmpy[:, 2])**2))
        diff_p = np.sqrt(np.mean((poiseuille_cpp[:, 1] - poiseuille_lbmpy[:, 1])**2))
        summary_data['Poiseuille'] = {
            'C++ L2 Error': f'{cpp_err_p:.2e}',
            'lbmpy L2 Error': f'{lbmpy_err_p:.2e}',
            'Difference': f'{diff_p:.2e}'
        }

    mo.md(f"""
    ### Error Summary

    {mo.as_html(summary_data) if summary_data else 'No data available'}
    """)
    return


if __name__ == "__main__":
    app.run()
