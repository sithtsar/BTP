#!/usr/bin/env python3
"""
Results Aggregation Script for Performance Sweep

Parses simulation output files and creates comprehensive results summary table.
Extracts metrics like L2 errors, drag/lift coefficients, wall-clock times.
"""

import os
import sys
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

class ResultsAggregator:
    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.results = []

    def parse_analytical_output(self, solver_type: str, case_name: str) -> Optional[Dict]:
        """Parse analytical test output for L2 errors and timing."""
        test_cases = {
            'couette': 'Couette Flow',
            'poiseuille': 'Poiseuille Flow',
            'taylor_green': 'Taylor-Green Vortex'
        }

        solver_name = 'BGK' if solver_type == '0' else 'ELBM'
        search_dir = self.output_dir / f'sweep_analytical_{solver_name.lower()}'

        if not search_dir.exists():
            return None

        # Look for output files
        output_log = search_dir / 'output.log'
        if not output_log.exists():
            return None

        try:
            with open(output_log, 'r') as f:
                content = f.read()

            # Extract L2 error from output
            l2_patterns = [
                rf'{case_name}.*?L2 error:\s*([\d.e+-]+)',
                rf'L2 error.*?([\d.e+-]+)',
                rf'Error:\s*([\d.e+-]+)'
            ]

            l2_error = None
            for pattern in l2_patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    try:
                        l2_error = float(match.group(1))
                        break
                    except (ValueError, IndexError):
                        continue

            # Extract timing info if available
            timing_match = re.search(r'Time.*?(\d+\.?\d*)\s*ms', content)
            wall_time = float(timing_match.group(1)) if timing_match else None

            # Check for stability (NaN or divergence)
            is_stable = 'nan' not in content.lower() and 'diverged' not in content.lower()
            is_stable = is_stable and 'inf' not in content.lower()

            return {
                'case': case_name,
                'solver': solver_name,
                'type': 'analytical',
                'l2_error': l2_error,
                'wall_time': wall_time,
                'stable': is_stable,
                'nodes': 128 * 128  # Analytical tests use 128x128
            }

        except Exception as e:
            print(f"Warning: Could not parse {output_log}: {e}", file=sys.stderr)
            return None

    def parse_cylinder_output(self, Re: int, solver_type: str) -> Optional[Dict]:
        """Parse cylinder flow output for drag/lift and timing."""
        solver_name = 'BGK' if solver_type == '0' else 'ELBM'
        search_dir = self.output_dir / f'sweep_cylinder_re{Re}_{solver_name.lower()}'

        if not search_dir.exists():
            return None

        try:
            output_log = search_dir / 'output.log'
            if not output_log.exists():
                return None

            with open(output_log, 'r') as f:
                content = f.read()

            # Extract final drag/lift coefficients
            cd_patterns = [
                r'Cd = ([-\d.e+-]+)',
                r'Cd.*?([+-]?[\d.e+-]+)',
            ]
            cl_patterns = [
                r'Cl = ([-\d.e+-]+)',
                r'Cl.*?([+-]?[\d.e+-]+)',
            ]

            cd = None
            for pattern in cd_patterns:
                matches = re.findall(pattern, content)
                if matches:
                    try:
                        # Get last match (final value)
                        cd = float(matches[-1])
                        break
                    except (ValueError, IndexError):
                        continue

            cl = None
            for pattern in cl_patterns:
                matches = re.findall(pattern, content)
                if matches:
                    try:
                        cl = float(matches[-1])
                        break
                    except (ValueError, IndexError):
                        continue

            # Extract timing
            timing_match = re.search(r'Time.*?(\d+\.?\d*)\s*s', content)
            wall_time = float(timing_match.group(1)) * 1000 if timing_match else None

            # Check for stability
            is_stable = 'nan' not in content.lower() and 'diverged' not in content.lower()

            return {
                'case': f'Cylinder Re={Re}',
                'solver': solver_name,
                'type': 'cylinder',
                'Re': Re,
                'Cd': cd,
                'Cl': cl,
                'wall_time': wall_time,
                'stable': is_stable,
                'nodes': 400 * 120
            }

        except Exception as e:
            print(f"Warning: Could not parse cylinder output for Re={Re}, {solver_name}: {e}", file=sys.stderr)
            return None

    def aggregate(self) -> List[Dict]:
        """Aggregate all results."""
        results = []

        # Parse analytical tests
        for case in ['couette', 'poiseuille', 'taylor_green']:
            for solver_type in ['0', '1']:
                result = self.parse_analytical_output(solver_type, case)
                if result:
                    results.append(result)

        # Parse cylinder tests
        for Re in [10, 100, 1000]:
            for solver_type in ['0', '1']:
                result = self.parse_cylinder_output(Re, solver_type)
                if result:
                    results.append(result)

        self.results = results
        return results

    def save_csv(self, output_file: str):
        """Save aggregated results to CSV."""
        if not self.results:
            print("No results to save!", file=sys.stderr)
            return

        # Determine which columns we have
        has_l2 = any('l2_error' in r for r in self.results)
        has_forces = any('Cd' in r for r in self.results)

        # Build header
        header = ['Case', 'Solver', 'Type', 'Stable']
        if has_l2:
            header.append('L2 Error')
        if has_forces:
            header.extend(['Re', 'Cd', 'Cl'])
        header.extend(['Wall Time (ms)', 'Nodes', 'Throughput (M nodes/s)'])

        # Write CSV
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)

            for result in self.results:
                row = [
                    result.get('case', ''),
                    result.get('solver', ''),
                    result.get('type', ''),
                    '✓' if result.get('stable', False) else '✗'
                ]

                if has_l2:
                    l2 = result.get('l2_error', '')
                    if l2 is not None:
                        row.append(f'{l2:.2e}')
                    else:
                        row.append('N/A')

                if has_forces:
                    row.append(result.get('Re', ''))
                    cd = result.get('Cd')
                    cl = result.get('Cl')
                    row.append(f'{cd:.4f}' if cd is not None else 'N/A')
                    row.append(f'{cl:.6f}' if cl is not None else 'N/A')

                # Timing and throughput
                wall_time = result.get('wall_time')
                nodes = result.get('nodes', 1)

                if wall_time:
                    row.append(f'{wall_time:.2f}')
                    throughput = (nodes / 1e6) / (wall_time / 1000) if wall_time > 0 else 0
                    row.append(f'{throughput:.2f}')
                else:
                    row.extend(['N/A', 'N/A'])

                row.append(nodes)

                writer.writerow(row)

        print(f"Results saved to: {output_file}")

    def print_summary(self):
        """Print summary statistics."""
        if not self.results:
            print("No results to summarize")
            return

        print("\n" + "=" * 80)
        print("PERFORMANCE SWEEP SUMMARY")
        print("=" * 80)

        # Group by solver
        bgk_results = [r for r in self.results if r.get('solver') == 'BGK']
        elbm_results = [r for r in self.results if r.get('solver') == 'ELBM']

        print(f"\nBGK Results: {len(bgk_results)} tests")
        bgk_stable = sum(1 for r in bgk_results if r.get('stable', False))
        print(f"  Stable: {bgk_stable}/{len(bgk_results)}")

        print(f"\nELBM Results: {len(elbm_results)} tests")
        elbm_stable = sum(1 for r in elbm_results if r.get('stable', False))
        print(f"  Stable: {elbm_stable}/{len(elbm_results)}")

        # Timing comparison
        bgk_times = [r.get('wall_time') for r in bgk_results if r.get('wall_time')]
        elbm_times = [r.get('wall_time') for r in elbm_results if r.get('wall_time')]

        if bgk_times:
            avg_bgk = sum(bgk_times) / len(bgk_times)
            print(f"\nAverage BGK Time: {avg_bgk:.2f} ms")
        if elbm_times:
            avg_elbm = sum(elbm_times) / len(elbm_times)
            print(f"Average ELBM Time: {avg_elbm:.2f} ms")

        if bgk_times and elbm_times:
            ratio = (sum(elbm_times) / len(elbm_times)) / (sum(bgk_times) / len(bgk_times))
            print(f"ELBM/BGK Time Ratio: {ratio:.2f}x")

        print("\n" + "=" * 80)


def main():
    if len(sys.argv) < 2:
        print("Usage: aggregate_results.py <output_dir> [output_csv]")
        sys.exit(1)

    output_dir = sys.argv[1]
    output_csv = sys.argv[2] if len(sys.argv) > 2 else os.path.join(output_dir, 'sweep_summary.csv')

    aggregator = ResultsAggregator(output_dir)
    results = aggregator.aggregate()

    print(f"Aggregated {len(results)} test results")
    aggregator.save_csv(output_csv)
    aggregator.print_summary()


if __name__ == '__main__':
    main()
