"""
Run All lbmpy Validations
==========================
Master script to run Couette, Poiseuille, and Taylor-Green validations
comparing BGK, ELBM, and analytical solutions.
"""

import sys
import subprocess
from pathlib import Path


def run_validation(script_name):
    """Run a validation script and report results."""
    script_path = Path(__file__).parent / script_name

    print(f"\n{'='*80}")
    print(f"Running {script_name}...")
    print(f"{'='*80}\n")

    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            capture_output=False,
            text=True,
            check=True
        )
        print(f"\n✓ {script_name} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n✗ {script_name} failed with exit code {e.returncode}")
        return False
    except Exception as e:
        print(f"\n✗ {script_name} failed with exception: {e}")
        return False


def main():
    """Run all validation scripts."""
    print("\n" + "="*80)
    print("LBMPY VALIDATION SUITE: BGK vs ELBM vs ANALYTICAL")
    print("="*80)
    print("\nThis will run three validation cases:")
    print("  1. Couette Flow (shear-driven)")
    print("  2. Poiseuille Flow (pressure-driven)")
    print("  3. Taylor-Green Vortex (decaying vortex)")
    print("\nEach case will compare BGK, ELBM, and analytical solutions.")
    print("="*80 + "\n")

    scripts = [
        'couette_lbmpy.py',
        'poiseuille_lbmpy.py',
        'taylor_green_lbmpy.py'
    ]

    results = {}
    for script in scripts:
        results[script] = run_validation(script)

    # Summary
    print("\n" + "="*80)
    print("VALIDATION SUITE SUMMARY")
    print("="*80)

    all_passed = True
    for script, passed in results.items():
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{script:<30} {status}")
        if not passed:
            all_passed = False

    print("="*80)

    if all_passed:
        print("\n✓ All validations PASSED!\n")
        print("Generated figures are available in:")
        print("  figures/lbmpy_validation/\n")
        sys.exit(0)
    else:
        print("\n✗ Some validations FAILED\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
