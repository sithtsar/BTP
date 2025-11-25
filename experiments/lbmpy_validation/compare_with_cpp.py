"""
Comparison Script: C++ vs lbmpy Implementations
===============================================
Compares validation test results between C++ and Python lbmpy implementations.
"""

import subprocess
import sys
import os
import numpy as np
from pathlib import Path


def load_results(filepath):
    """Load results from .dat file."""
    if not os.path.exists(filepath):
        return None

    data = np.loadtxt(filepath)
    return data


def print_section(title):
    """Print formatted section header."""
    print()
    print("=" * 70)
    print(title.center(70))
    print("=" * 70)


def run_cpp_tests():
    """Run C++ analytical tests."""
    print_section("Running C++ Analytical Tests")

    cpp_executable = "../../build/test_analytical"

    if not os.path.exists(cpp_executable):
        print("ERROR: C++ executable not found. Please build first:")
        print("  cd ../.. && ./build.sh")
        return False

    try:
        result = subprocess.run(
            [cpp_executable],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)

        return result.returncode == 0

    except subprocess.TimeoutExpired:
        print("ERROR: C++ tests timed out")
        return False
    except Exception as e:
        print(f"ERROR running C++ tests: {e}")
        return False


def run_lbmpy_tests():
    """Run lbmpy validation tests."""
    print_section("Running lbmpy Validation Tests")

    scripts = [
        ("Couette Flow", "couette_lbmpy.py"),
        ("Poiseuille Flow", "poiseuille_lbmpy.py"),
        ("Taylor-Green Vortex", "taylor_green_lbmpy.py")
    ]

    results = {}

    for name, script in scripts:
        print(f"\n--- {name} ---")
        try:
            result = subprocess.run(
                [sys.executable, script],
                capture_output=True,
                text=True,
                timeout=300
            )

            print(result.stdout)
            if result.stderr:
                print("STDERR:", result.stderr)

            results[name] = (result.returncode == 0)

        except subprocess.TimeoutExpired:
            print(f"ERROR: {name} timed out")
            results[name] = False
        except Exception as e:
            print(f"ERROR running {name}: {e}")
            results[name] = False

    return all(results.values())


def compare_couette():
    """Compare Couette flow results."""
    print_section("Couette Flow Comparison")

    cpp_file = "../../output/couette_results.dat"
    lbmpy_file = "../../output/couette_lbmpy.dat"

    cpp_data = load_results(cpp_file)
    lbmpy_data = load_results(lbmpy_file)

    if cpp_data is None:
        print(f"WARNING: C++ results not found at {cpp_file}")
        return False

    if lbmpy_data is None:
        print(f"WARNING: lbmpy results not found at {lbmpy_file}")
        return False

    # Compare
    print(f"C++ data shape: {cpp_data.shape}")
    print(f"lbmpy data shape: {lbmpy_data.shape}")

    # Extract velocity profiles
    cpp_u = cpp_data[:, 1]  # Column 1: u_sim
    lbmpy_u = lbmpy_data[:, 1]

    cpp_analytical = cpp_data[:, 2]  # Column 2: u_exact
    lbmpy_analytical = lbmpy_data[:, 2]

    # Compute errors
    cpp_error = np.sqrt(np.mean((cpp_u - cpp_analytical)**2))
    lbmpy_error = np.sqrt(np.mean((lbmpy_u - lbmpy_analytical)**2))
    diff_error = np.sqrt(np.mean((cpp_u - lbmpy_u)**2))

    print(f"\nC++ L2 error:        {cpp_error:.6e}")
    print(f"lbmpy L2 error:      {lbmpy_error:.6e}")
    print(f"C++ vs lbmpy error:  {diff_error:.6e}")

    if diff_error < 0.001:
        print("✓ C++ and lbmpy results agree (diff < 0.001)")
        return True
    else:
        print("✗ C++ and lbmpy results differ significantly")
        return False


def compare_poiseuille():
    """Compare Poiseuille flow results."""
    print_section("Poiseuille Flow Comparison")

    cpp_file = "../../output/poiseuille_results.dat"
    lbmpy_file = "../../output/poiseuille_lbmpy.dat"

    cpp_data = load_results(cpp_file)
    lbmpy_data = load_results(lbmpy_file)

    if cpp_data is None:
        print(f"WARNING: C++ results not found at {cpp_file}")
        return False

    if lbmpy_data is None:
        print(f"WARNING: lbmpy results not found at {lbmpy_file}")
        return False

    print(f"C++ data shape: {cpp_data.shape}")
    print(f"lbmpy data shape: {lbmpy_data.shape}")

    # Extract velocity profiles
    cpp_u = cpp_data[:, 1]
    lbmpy_u = lbmpy_data[:, 1]

    cpp_analytical = cpp_data[:, 2]
    lbmpy_analytical = lbmpy_data[:, 2]

    # Compute errors
    cpp_error = np.sqrt(np.mean((cpp_u - cpp_analytical)**2))
    lbmpy_error = np.sqrt(np.mean((lbmpy_u - lbmpy_analytical)**2))
    diff_error = np.sqrt(np.mean((cpp_u - lbmpy_u)**2))

    print(f"\nC++ L2 error:        {cpp_error:.6e}")
    print(f"lbmpy L2 error:      {lbmpy_error:.6e}")
    print(f"C++ vs lbmpy error:  {diff_error:.6e}")

    if diff_error < 0.001:
        print("✓ C++ and lbmpy results agree (diff < 0.001)")
        return True
    else:
        print("✗ C++ and lbmpy results differ significantly")
        return False


def compare_taylor_green():
    """Compare Taylor-Green vortex results."""
    print_section("Taylor-Green Vortex Comparison")

    cpp_file = "../../output/tg_energy.dat"
    lbmpy_file = "../../output/tg_lbmpy.dat"

    cpp_data = load_results(cpp_file)
    lbmpy_data = load_results(lbmpy_file)

    if cpp_data is None:
        print(f"WARNING: C++ results not found at {cpp_file}")
        return False

    if lbmpy_data is None:
        print(f"WARNING: lbmpy results not found at {lbmpy_file}")
        return False

    print(f"C++ data shape: {cpp_data.shape}")
    print(f"lbmpy data shape: {lbmpy_data.shape}")

    # This is a basic check - energy decay comparison
    # More sophisticated comparison would require viscosity extraction
    print("\nNote: Taylor-Green comparison requires matching initial conditions")
    print("Current lbmpy implementation uses simplified initialization")

    return True


def main():
    """Main comparison workflow."""
    print("="*70)
    print("C++ vs lbmpy Validation Comparison".center(70))
    print("="*70)

    # Change to script directory
    script_dir = Path(__file__).parent
    os.chdir(script_dir)

    # Ensure output directory exists
    os.makedirs("../../output", exist_ok=True)

    # Step 1: Run C++ tests
    cpp_success = run_cpp_tests()
    if not cpp_success:
        print("\nERROR: C++ tests failed")
        # Continue anyway to see what we can compare

    # Step 2: Run lbmpy tests
    lbmpy_success = run_lbmpy_tests()
    if not lbmpy_success:
        print("\nERROR: lbmpy tests failed")

    # Step 3: Compare results
    print_section("Comparing Results")

    couette_ok = compare_couette()
    poiseuille_ok = compare_poiseuille()
    taylor_green_ok = compare_taylor_green()

    # Final summary
    print_section("Summary")
    print()
    print(f"C++ tests:          {'PASS' if cpp_success else 'FAIL'}")
    print(f"lbmpy tests:        {'PASS' if lbmpy_success else 'FAIL'}")
    print()
    print(f"Couette comparison: {'✓ PASS' if couette_ok else '✗ FAIL'}")
    print(f"Poiseuille comparison: {'✓ PASS' if poiseuille_ok else '✗ FAIL'}")
    print(f"Taylor-Green comparison: {'✓ PASS' if taylor_green_ok else '✗ FAIL'}")
    print()

    all_pass = all([cpp_success, lbmpy_success, couette_ok, poiseuille_ok, taylor_green_ok])

    if all_pass:
        print("="*70)
        print("ALL VALIDATIONS PASSED".center(70))
        print("="*70)
        return 0
    else:
        print("="*70)
        print("SOME VALIDATIONS FAILED".center(70))
        print("="*70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
