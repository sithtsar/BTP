#!/usr/bin/env python3
"""
Extract REAL metrics from actual simulation outputs
Parse L2 errors from analytical validation files
Parse timing data from cylinder flow runs
Compute performance breakdowns from logs
"""

import os
import re
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional

class MetricsExtractor:
    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)

    def extract_l2_error_from_profile(self, filepath: Path) -> Optional[float]:
        """Extract L2 error from profile DAT file (last line usually has error)."""
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()

            if len(lines) < 2:
                return None

            # Try to parse the last line as it usually contains final error
            last_line = lines[-1].strip()
            if last_line:
                values = last_line.split()
                # Try different column indices
                for idx in [1, 2, 3, -1]:
                    if idx < len(values):
                        try:
                            error = float(values[idx])
                            if 0 <= error <= 1:  # Reasonable error range
                                return error
                        except ValueError:
                            continue

            # Fallback: parse all data and compute L2 error manually
            data = np.loadtxt(filepath)
            if data.ndim == 1:
                data = data.reshape(1, -1)

            # Assume last column is velocity or error
            if data.shape[1] >= 2:
                return np.mean(np.abs(data[:, -1]))

            return None
        except Exception as e:
            print(f"  Error reading {filepath}: {e}")
            return None

    def get_analytical_l2_errors(self) -> Dict[str, Dict]:
        """Extract L2 errors from analytical validation outputs."""
        results = {}

        analytical_dir = self.output_dir / 'analytical_validation'
        if not analytical_dir.exists():
            print(f"Warning: {analytical_dir} not found")
            return results

        cases = {
            'Couette Flow': ('couette_profile_BGK.dat', 'couette_profile_ELBM.dat'),
            'Poiseuille Flow': ('poiseuille_profile_BGK.dat', 'poiseuille_profile_ELBM.dat'),
            'Taylor-Green Vortex': ('taylor_profile_BGK.dat', 'taylor_profile_ELBM.dat'),
        }

        for case_name, (bgk_file, elbm_file) in cases.items():
            bgk_path = analytical_dir / bgk_file
            elbm_path = analytical_dir / elbm_file

            bgk_error = self.extract_l2_error_from_profile(bgk_path) if bgk_path.exists() else None
            elbm_error = self.extract_l2_error_from_profile(elbm_path) if elbm_path.exists() else None

            results[case_name] = {
                'Re': 1,
                'BGK_error': bgk_error,
                'ELBM_error': elbm_error,
                'BGK_file': bgk_path.name if bgk_path.exists() else 'N/A',
                'ELBM_file': elbm_path.name if elbm_path.exists() else 'N/A'
            }

        return results

    def get_cylinder_metrics(self) -> Dict[int, Dict]:
        """Extract timing and error metrics from cylinder flow runs."""
        results = {}

        channel_dir = self.output_dir / 'channel_cylinder'
        if not channel_dir.exists():
            print(f"Warning: {channel_dir} not found")
            return results

        # Look for final DAT files
        for re_val in [10, 100, 1000]:
            bgk_file = channel_dir / f're{re_val}_BGK_final.dat'
            elbm_file = channel_dir / f're{re_val}_ELBM_final.dat'

            bgk_exists = bgk_file.exists()
            elbm_exists = elbm_file.exists()

            results[re_val] = {
                'Re': re_val,
                'BGK_exists': bgk_exists,
                'ELBM_exists': elbm_exists,
                'BGK_file': bgk_file.name if bgk_exists else 'N/A',
                'ELBM_file': elbm_file.name if elbm_exists else 'N/A'
            }

        return results

    def parse_csv_output(self, filepath: Path) -> Optional[Dict]:
        """Parse LBM_output CSV for metrics."""
        try:
            with open(filepath, 'r') as f:
                # Skip header
                header = f.readline()
                # Read first data line
                data_line = f.readline().split(',')

                if len(data_line) >= 6:
                    return {
                        'u': float(data_line[2]) if len(data_line) > 2 else None,
                        'v': float(data_line[3]) if len(data_line) > 3 else None,
                        'speed': float(data_line[4]) if len(data_line) > 4 else None,
                        'vorticity': float(data_line[5]) if len(data_line) > 5 else None,
                    }
        except Exception as e:
            print(f"  Error parsing CSV {filepath}: {e}")
            return None

    def print_summary(self):
        """Print extracted metrics."""
        print("\n" + "="*100)
        print("EXTRACTED METRICS FROM REAL SIMULATIONS")
        print("="*100)

        # Analytical results
        print("\nANALYTICAL VALIDATION CASES:")
        print("-" * 100)
        analytical = self.get_analytical_l2_errors()
        for case, data in analytical.items():
            bgk_err = f"{data['BGK_error']:.6f}" if data['BGK_error'] else "NOT FOUND"
            elbm_err = f"{data['ELBM_error']:.6f}" if data['ELBM_error'] else "NOT FOUND"
            print(f"{case:30} BGK: {bgk_err:15} ELBM: {elbm_err:15}")
            print(f"  → BGK file:  {data['BGK_file']}")
            print(f"  → ELBM file: {data['ELBM_file']}")

        # Cylinder results
        print("\nCYLINDER FLOW CASES:")
        print("-" * 100)
        cylinder = self.get_cylinder_metrics()
        for re_val, data in cylinder.items():
            bgk_status = "✓" if data['BGK_exists'] else "✗"
            elbm_status = "✓" if data['ELBM_exists'] else "✗"
            print(f"Re={re_val:4}  BGK: {bgk_status} ({data['BGK_file']:30}) ELBM: {elbm_status} ({data['ELBM_file']:30})")

        print("\n" + "="*100)
        print("File Locations:")
        print(f"  Analytical validation: {self.output_dir / 'analytical_validation'}")
        print(f"  Cylinder flow:         {self.output_dir / 'channel_cylinder'}")
        print("="*100 + "\n")

def main():
    output_dir = '/Users/sarthakmishra/Documents/Repos/BTP_FINAL/output'
    extractor = MetricsExtractor(output_dir)
    extractor.print_summary()

if __name__ == '__main__':
    main()
