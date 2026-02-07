#!/usr/bin/env python3
"""
Analyze experiment results and generate summary tables/plots.

Usage:
    python analyze_results.py <run_dir>
    python analyze_results.py <run_dir> --csv
    python analyze_results.py <run_dir> --plot

Example:
    python analyze_results.py experiments/radar_sweep/runs/20260201_103000
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict
import argparse


def parse_result_file(filepath: Path) -> dict:
    """Parse a single result file and extract key metrics."""
    result = {
        'file': filepath.name,
        'config_name': filepath.stem.replace('_results', ''),
    }

    with open(filepath, 'r') as f:
        content = f.read()

    # Extract configuration values
    patterns = {
        'code_distance': r'Code Distance:\s*(\d+)',
        'rounds': r'Rounds:\s*(\d+)',
        'physical_error': r'Physical Error Rate:\s*([\d.e-]+)',
        'total_shots': r'Total Shots:\s*(\d+)',
        'distillation_protocol': r'Distillation Protocol:\s*(.+)',
        'distillation_rounds': r'Distillation Rounds:\s*(\d+)',
        'raw_epr_fidelity': r'Raw EPR Fidelity:\s*([\d.]+)',
        'entanglement_rate': r'Entanglement Rate:\s*([\d.]+)\s*MHz',
        'logical_error_rate': r'Logical Error Rate:\s*([\d.e-]+)',
        'runtime': r'Total Runtime:\s*(\d+)\s*hours,\s*(\d+)\s*minutes,\s*(\d+)\s*seconds',
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            if key == 'runtime':
                hours, mins, secs = match.groups()
                result[key] = int(hours) * 3600 + int(mins) * 60 + int(secs)
            elif key in ['code_distance', 'rounds', 'total_shots', 'distillation_rounds']:
                result[key] = int(match.group(1))
            elif key in ['physical_error', 'raw_epr_fidelity', 'logical_error_rate', 'entanglement_rate']:
                result[key] = float(match.group(1))
            else:
                result[key] = match.group(1).strip()

    return result


def load_all_results(run_dir: Path) -> list:
    """Load all result files from a run directory (including nested sub-experiments)."""
    results = []

    # Try direct results directory first
    results_dir = run_dir / 'results'

    if results_dir.exists():
        # Check if results_dir contains sub-experiment directories or direct results
        subdirs = [d for d in results_dir.iterdir() if d.is_dir()]

        if subdirs:
            # Nested structure: results/sub_experiment/*.txt
            for subdir in sorted(subdirs):
                sub_name = subdir.name
                for result_file in sorted(subdir.glob('*_results.txt')):
                    try:
                        result = parse_result_file(result_file)
                        result['sub_experiment'] = sub_name
                        results.append(result)
                    except Exception as e:
                        print(f"Warning: Failed to parse {result_file}: {e}")
        else:
            # Flat structure: results/*.txt
            for result_file in sorted(results_dir.glob('*_results.txt')):
                try:
                    result = parse_result_file(result_file)
                    results.append(result)
                except Exception as e:
                    print(f"Warning: Failed to parse {result_file}: {e}")
    else:
        print(f"Error: Results directory not found: {results_dir}")

    return results


def print_summary_table(results: list):
    """Print a formatted summary table."""
    if not results:
        print("No results to display.")
        return

    # Group by code distance
    by_distance = defaultdict(list)
    for r in results:
        d = r.get('code_distance', 'unknown')
        by_distance[d].append(r)

    print("\n" + "=" * 100)
    print("EXPERIMENT RESULTS SUMMARY")
    print("=" * 100)

    for distance in sorted(by_distance.keys()):
        print(f"\n--- Code Distance d={distance} ---")
        print(f"{'Protocol':<20} {'Rate (MHz)':<12} {'Raw F':<8} {'LER':<14} {'Shots':<12}")
        print("-" * 70)

        dist_results = sorted(by_distance[distance],
                              key=lambda x: (x.get('entanglement_rate', 0),
                                           x.get('distillation_protocol', '')))

        for r in dist_results:
            protocol = r.get('distillation_protocol', 'None')
            if protocol == 'None' or not protocol:
                protocol = 'none'
            rate = r.get('entanglement_rate', 0)
            fidelity = r.get('raw_epr_fidelity', 0)
            ler = r.get('logical_error_rate', 0)
            shots = r.get('total_shots', 0)

            print(f"{protocol:<20} {rate:<12.0f} {fidelity:<8.4f} {ler:<14.2e} {shots:<12,}")

    print("\n" + "=" * 100)


def export_csv(results: list, output_path: Path):
    """Export results to CSV format."""
    if not results:
        print("No results to export.")
        return

    # Determine all columns
    columns = ['config_name', 'code_distance', 'entanglement_rate', 'distillation_protocol',
               'distillation_rounds', 'raw_epr_fidelity', 'logical_error_rate', 'total_shots', 'runtime']

    with open(output_path, 'w') as f:
        # Header
        f.write(','.join(columns) + '\n')

        # Data rows
        for r in results:
            row = [str(r.get(col, '')) for col in columns]
            f.write(','.join(row) + '\n')

    print(f"CSV exported to: {output_path}")


def generate_plot(results: list, output_path: Path):
    """Generate a plot of LER vs entanglement rate by protocol."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not installed. Skipping plot generation.")
        print("Install with: pip install matplotlib")
        return

    if not results:
        print("No results to plot.")
        return

    # Group by (distance, protocol)
    grouped = defaultdict(lambda: {'rates': [], 'lers': []})
    for r in results:
        key = (r.get('code_distance', 0), r.get('distillation_protocol', 'none'))
        rate = r.get('entanglement_rate', 0)
        ler = r.get('logical_error_rate', 0)
        if rate > 0 and ler > 0:
            grouped[key]['rates'].append(rate)
            grouped[key]['lers'].append(ler)

    # Create plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    distances = sorted(set(r.get('code_distance', 0) for r in results))

    protocol_colors = {
        'none': 'black',
        'None': 'black',
        '2→1 Pumping': 'blue',
        '3→1 Pumping': 'green',
        '2→1 Recurrence': 'orange',
        '3→1 Recurrence': 'red',
    }

    for idx, d in enumerate(distances[:3]):
        ax = axes[idx] if len(distances) > 1 else axes
        ax.set_title(f'Code Distance d={d}')
        ax.set_xlabel('Entanglement Rate (MHz)')
        ax.set_ylabel('Logical Error Rate')
        ax.set_yscale('log')

        for (dist, protocol), data in grouped.items():
            if dist != d:
                continue
            if data['rates']:
                color = protocol_colors.get(protocol, 'gray')
                ax.scatter(data['rates'], data['lers'], label=protocol, color=color, s=50)
                # Sort for line plot
                sorted_pairs = sorted(zip(data['rates'], data['lers']))
                rates, lers = zip(*sorted_pairs)
                ax.plot(rates, lers, color=color, alpha=0.5)

        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Plot saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description='Analyze experiment results')
    parser.add_argument('run_dir', type=Path, help='Path to run directory')
    parser.add_argument('--csv', action='store_true', help='Export to CSV')
    parser.add_argument('--plot', action='store_true', help='Generate plot')
    args = parser.parse_args()

    run_dir = args.run_dir
    if not run_dir.exists():
        print(f"Error: Run directory not found: {run_dir}")
        sys.exit(1)

    print(f"Analyzing results in: {run_dir}")
    results = load_all_results(run_dir)
    print(f"Found {len(results)} result files")

    if not results:
        sys.exit(1)

    # Always print summary
    print_summary_table(results)

    # Export CSV if requested
    if args.csv:
        csv_path = run_dir / 'results_summary.csv'
        export_csv(results, csv_path)

    # Generate plot if requested
    if args.plot:
        plot_path = run_dir / 'results_plot.png'
        generate_plot(results, plot_path)


if __name__ == "__main__":
    main()
