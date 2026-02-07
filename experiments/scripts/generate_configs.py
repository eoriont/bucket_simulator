#!/usr/bin/env python3
"""
Config generator for bucket_simulator experiments.

Usage:
    python generate_configs.py <experiment_spec.json> <output_dir>

Experiment spec format (JSON):
{
    "template": "distributed_base.txt",
    "name": "radar_sweep",
    "defaults": {
        "physical_error": 0.001,
        "total_shots": "100K"
    },
    "sweep": [
        {"var": "code_distance", "values": [13, 15, 17]},
        {"var": "entanglement_rate", "values": ["100e6", "200e6"]}
    ]
}
"""

import sys
import os
import json
import itertools
from pathlib import Path


def parse_value(value):
    """Convert values to appropriate format for config files."""
    if isinstance(value, str):
        return value
    elif isinstance(value, bool):
        return "true" if value else "false"
    elif isinstance(value, (int, float)):
        # Format large numbers in scientific notation if needed
        if isinstance(value, float) and abs(value) >= 1e6:
            return f"{value:.0e}"
        return str(value)
    return str(value)


def generate_config(template_path: Path, output_path: Path, variables: dict):
    """Generate a config file from template with variable substitution."""
    with open(template_path, 'r') as f:
        template = f.read()

    # Substitute all variables
    config = template
    for var, value in variables.items():
        placeholder = "{{" + var.upper() + "}}"
        config = config.replace(placeholder, parse_value(value))

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(config)

    return output_path


def load_experiment_spec(spec_path: Path) -> dict:
    """Load experiment specification from JSON file."""
    with open(spec_path, 'r') as f:
        return json.load(f)


def generate_sweep_configs(spec: dict, template_dir: Path, output_dir: Path) -> list:
    """Generate all config files for a parameter sweep."""
    template_path = template_dir / spec['template']
    defaults = spec.get('defaults', {})
    sweep_vars = spec.get('sweep', [])

    # Build list of (var_name, values) for each sweep variable
    sweep_items = [(s['var'], s['values']) for s in sweep_vars]

    # Generate all combinations
    if sweep_items:
        var_names = [item[0] for item in sweep_items]
        var_values = [item[1] for item in sweep_items]
        combinations = list(itertools.product(*var_values))
    else:
        var_names = []
        combinations = [()]

    configs = []
    for combo in combinations:
        # Build variables dict: defaults + this combination
        variables = dict(defaults)
        for var_name, value in zip(var_names, combo):
            variables[var_name] = value

        # Auto-set rounds = code_distance if not explicitly set and distance is being swept
        if 'code_distance' in variables and 'rounds' not in variables:
            variables['rounds'] = variables['code_distance']

        # Generate config filename from varied parameters
        name_parts = []
        for var_name, value in zip(var_names, combo):
            # Shorten variable names for filename
            short_names = {
                'code_distance': 'd',
                'entanglement_rate': 'rate',
                'distillation_protocol': 'dist',
                'distillation_rounds': 'rnd',
                'raw_epr_fidelity': 'fid',
            }
            short_name = short_names.get(var_name, var_name[:4])
            short_val = str(value).replace('e6', 'M').replace('e9', 'G').replace('e-', 'em').replace('_', '').replace('to', '')
            name_parts.append(f"{short_name}{short_val}")

        config_name = "_".join(name_parts) if name_parts else "default"
        config_path = output_dir / f"{config_name}.txt"

        generate_config(template_path, config_path, variables)
        configs.append({
            'path': config_path,
            'variables': variables,
            'name': config_name
        })

    return configs


def main():
    if len(sys.argv) < 3:
        print("Usage: python generate_configs.py <experiment_spec.json> <output_dir>")
        sys.exit(1)

    spec_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    # Determine template directory (relative to script location)
    script_dir = Path(__file__).parent.parent
    template_dir = script_dir / "templates"

    spec = load_experiment_spec(spec_path)
    configs = generate_sweep_configs(spec, template_dir, output_dir)

    print(f"Generated {len(configs)} config files in {output_dir}")
    for cfg in configs:
        print(f"  - {cfg['path'].name}")

    # Write manifest
    manifest_path = output_dir / "manifest.json"
    manifest = {
        'experiment': spec.get('name', 'unnamed'),
        'template': spec['template'],
        'num_configs': len(configs),
        'configs': [{'name': c['name'], 'file': c['path'].name} for c in configs]
    }
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"Manifest written to {manifest_path}")


if __name__ == "__main__":
    main()
