"""
Marimo notebook for rectangular pipe flow pressure profile analysis
Reproduces Figures 3.2-3.19 from Yu thesis comparing BGK vs ELBM stability
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
    # ELBM Rectangular Pipe Flow Analysis

    Replication of Keming Yu's honor thesis results comparing BGK and ELBM solvers at different Reynolds numbers.

    ## Test Cases
    - **Case 1 (Re~10)**: Both solvers stable, similar behavior
    - **Case 2 (Re~100)**: BGK shows numerical diffusion, ELBM remains stable
    - **Case 3 (Re~1000)**: BGK complete breakdown, ELBM demonstrates unconditional stability
    """)
    return


@app.cell
def _(mo):
    case_selector = mo.ui.dropdown(
        options={
            "Case 1: Re~10 (ν=0.1 m²/s)": "case1_re10",
            "Case 2: Re~100 (ν=0.01 m²/s)": "case2_re100",
            "Case 3: Re~1000 (ν=0.001 m²/s)": "case3_re1000"
        },
        value="case1_re10",
        label="Test Case"
    )
    return (case_selector,)


@app.cell
def _(mo):
    step_selector = mo.ui.slider(
        start=0,
        stop=2500,
        step=250,
        value=2500,
        label="Timestep",
        show_value=True
    )
    return (step_selector,)


@app.cell
def _(case_selector, mo, step_selector):
    mo.hstack([case_selector, step_selector], justify="start")
    return


@app.cell
def _(Path, case_selector, np, step_selector):
    case = case_selector.value
    step = step_selector.value

    output_dir = Path("../output")
    bgk_file = output_dir / f"{case}_bgk" / f"bgk_centerline_{step}.dat"
    elbm_file = output_dir / f"{case}_elbm" / f"elbm_centerline_{step}.dat"

    data_exists = bgk_file.exists() and elbm_file.exists()

    if data_exists:
        bgk_data = np.loadtxt(bgk_file, comments='#')
        elbm_data = np.loadtxt(elbm_file, comments='#')

        bgk_x = bgk_data[:, 0]
        bgk_p = bgk_data[:, 1]
        elbm_x = elbm_data[:, 0]
        elbm_p = elbm_data[:, 1]
    else:
        bgk_x = bgk_p = elbm_x = elbm_p = None
    return bgk_data, bgk_p, bgk_x, data_exists, elbm_data, elbm_p, elbm_x


@app.cell
def _(data_exists, mo):
    mo.stop(not data_exists, mo.md("⚠️ **Data not found!** Run simulations first: `./run_all_cases.sh`"))
    return


@app.cell
def _(bgk_p, bgk_x, case_selector, elbm_p, elbm_x, plt, step_selector):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # BGK
    ax1.plot(bgk_x, bgk_p, 'b-', linewidth=2.5, label='BGK')
    ax1.set_xlabel('x (lattice units)', fontsize=13)
    ax1.set_ylabel('Pressure (lattice units)', fontsize=13)
    ax1.set_title('LBGK Solver', fontsize=15, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.legend(fontsize=11)

    # ELBM
    ax2.plot(elbm_x, elbm_p, 'r-', linewidth=2.5, label='ELBM')
    ax2.set_xlabel('x (lattice units)', fontsize=13)
    ax2.set_ylabel('Pressure (lattice units)', fontsize=13)
    ax2.set_title('ELBM Solver', fontsize=15, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.legend(fontsize=11)

    case_name = case_selector.options_labels[case_selector.value]
    fig.suptitle(f'{case_name} - Centerline Pressure at t={step_selector.value}',
                fontsize=16, fontweight='bold', y=1.02)

    plt.tight_layout()
    plt.gca()
    return


@app.cell
def _(bgk_data, elbm_data, mo, np):
    bgk_mean = np.mean(bgk_data[:, 1])
    bgk_std = np.std(bgk_data[:, 1])
    elbm_mean = np.mean(elbm_data[:, 1])
    elbm_std = np.std(elbm_data[:, 1])

    rel_diff = abs(bgk_mean - elbm_mean) / elbm_mean * 100

    mo.md(f"""
    ## Statistics

    | Solver | Mean Pressure | Std Deviation |
    |--------|---------------|---------------|
    | BGK    | {bgk_mean:.6f} | {bgk_std:.6f} |
    | ELBM   | {elbm_mean:.6f} | {elbm_std:.6f} |

    **Relative Difference:** {rel_diff:.2f}%
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    ## Expected Behavior

    ### Case 1 (Re~10)
    Both BGK and ELBM should show stable, similar behavior with pressure profile flattening due to high viscosity.

    ### Case 2 (Re~100)
    BGK exhibits excessive numerical diffusion while ELBM maintains stability with minimal diffusion.

    ### Case 3 (Re~1000)
    BGK completely breaks down showing instability, while ELBM maintains stable physical behavior demonstrating unconditional stability.
    """)
    return


if __name__ == "__main__":
    app.run()
