#!/usr/bin/env python3
"""
Master script to generate ALL enhanced visualizations
Runs all plotting modules in correct order
"""

import subprocess
import sys
from pathlib import Path

def run_script(script_path, description):
    """Run a plotting script and report status"""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")

    try:
        result = subprocess.run([sys.executable, script_path],
                              check=True, capture_output=True, text=True)
        print(result.stdout)
        print(f"✓ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ {description} failed:")
        print(e.stderr)
        return False

def main():
    plotting_dir = Path(__file__).parent

    scripts = [
        (plotting_dir / 'geometry_diagrams.py',
         'Geometry Diagrams (all test cases)'),

        (plotting_dir / 'plot_analytical_entropy.py',
         'Analytical Validation - H-Theorem Plots'),

        (plotting_dir / 'plot_cylinder_entropy.py',
         'Cylinder Flow - Entropy Validation'),
    ]

    print("="*60)
    print("ELBM Visualization Suite")
    print("Generating all enhanced plots...")
    print("="*60)

    success_count = 0
    total_count = len(scripts)

    for script, description in scripts:
        if script.exists():
            if run_script(script, description):
                success_count += 1
        else:
            print(f"⚠ Skipping {description}: script not found at {script}")

    print(f"\n{'='*60}")
    print(f"Summary: {success_count}/{total_count} plotting modules succeeded")
    print(f"{'='*60}")

    if success_count == total_count:
        print("\n✓ All visualizations generated successfully!")
        print("\nOutput locations:")
        print("  - Geometry diagrams: figures/geometry_diagrams/")
        print("  - Analytical validation: figures/analytical_validation/")
        print("  - Cylinder flow: figures/channel_cylinder/")
    else:
        print(f"\n⚠ {total_count - success_count} plotting modules failed")
        print("Check output above for error details")
        sys.exit(1)

if __name__ == '__main__':
    main()
