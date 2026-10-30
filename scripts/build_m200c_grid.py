"""Rebuild the July 2024 M200m -> M200c table with the pinned modern CCL.

Run from the repository root with the optional `conversion` dependencies.
Released-catalog validation measures the difference from the historical CCL.
"""
import argparse
import hashlib
import json
from pathlib import Path

import camb
import numpy as np
import pyccl
import scipy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Output NPZ file')
    args = parser.parse_args()
    parameters = dict(Omega_c=0.2603, Omega_b=0.0486, h=0.6774,
                      sigma8=0.8159, n_s=0.9667, transfer_function='boltzmann_camb')
    cosmology = pyccl.Cosmology(**parameters)
    translate = pyccl.halos.mass_translator(
        mass_in='200m', mass_out='200c', concentration='Bhattacharya13')
    log_mass = np.linspace(12, 16.5, 200)
    log1pz = np.linspace(np.log(1), np.log(5.0), 30)
    # Preserve the original numeric convention: no factors of h at the CCL
    # boundary. Applying them would change the values in the released catalog.
    grid = np.column_stack([
        np.log10(translate(cosmology, 10**log_mass, np.exp(-u))) for u in log1pz])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, log_mass=log_mass, log1pz=log1pz, log_m200c=grid)
    root = Path(__file__).resolve().parents[1]
    provenance = {
        'method': 'historical_CCL_recipe_rebuilt_with_pinned_modern_CCL',
        'source_history': 'provenance/m200c-history.json',
        'history_sha256': hashlib.sha256((root/'provenance/m200c-history.json').read_bytes()).hexdigest(),
        'source_session': 749, 'source_line': 96, 'source_date': '2024-07-24',
        'cosmology': parameters, 'profile': 'NFW', 'concentration': 'Bhattacharya13',
        'mass_in': '200m', 'mass_out': '200c',
        'mass_convention': 'Released numeric Msun/h masses passed directly to CCL Msun interface; no h rescaling',
        'interpolation': 'linear in log10(M200m) and ln(1+z), returning log10(M200c)',
        'grid_shape': list(grid.shape), 'redshift_range': [0, 4], 'log10_mass_range': [12, 16.5],
        'versions': {name: module.__version__ for name, module in
                     [('pyccl', pyccl), ('camb', camb), ('numpy', np), ('scipy', scipy)]},
        'table_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
    }
    args.output.with_suffix('.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(args.output, flush=True)


if __name__ == '__main__':
    main()
