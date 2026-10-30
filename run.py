"""Command-line entry point; all output files go to explicit destinations."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from remapping import MassRemapping, apply_hdf5, fit_released_pairs, validate_pairs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest='command', required=True)
    fit = commands.add_parser('recover-grid', help='Fit a calibration spline to released mass pairs')
    fit.add_argument('--samples', type=Path, required=True)
    fit.add_argument('--output', type=Path, required=True)
    fit.add_argument('--regularization', type=float, default=1e-6)
    fit.add_argument('--tail-samples', type=Path)
    plot = commands.add_parser('plot', help='Reproduce the paper figure using released calibrated masses')
    plot.add_argument('--samples', type=Path, required=True)
    plot.add_argument('--output', type=Path, required=True)
    apply = commands.add_parser('apply', help='Write a new mass-only HDF5 sidecar')
    apply.add_argument('--grid', type=Path, required=True)
    apply.add_argument('--input', type=Path, required=True)
    apply.add_argument('--output', type=Path, required=True)
    apply.add_argument('--allow-extrapolation', action='store_true')
    apply.add_argument('--mass-definition', choices=['200m', '200c', 'both'], default='200m',
                       help='Output mass definition (default: 200m)')
    args = p.parse_args()
    if args.command == 'recover-grid':
        data = np.load(args.samples)
        training = [data['sample_100']]
        training_seeds = [100]
        if args.tail_samples:
            tail = np.load(args.tail_samples)
            for key in sorted(tail.files):
                if not key.startswith('sample_'):
                    continue
                seed = int(key.removeprefix('sample_'))
                if seed in (100, 102, 104) or not len(tail[key]):
                    continue
                training.append(tail[key])
                training_seeds.append(seed)
        training = np.concatenate(training)
        mapping = fit_released_pairs(training, regularization=args.regularization)
        mapping.provenance.update({
            'sample_sha256': hashlib.sha256(args.samples.read_bytes()).hexdigest(),
            'training_seeds': sorted(training_seeds),
            'development_validation_seed': 102,
            'training_validation': validate_pairs(mapping, training),
            'independent_validation': validate_pairs(mapping, data['sample_102']),
        })
        if args.tail_samples:
            mapping.provenance['tail_samples_sha256'] = hashlib.sha256(args.tail_samples.read_bytes()).hexdigest()
        mapping.save(args.output)
        print(json.dumps(mapping.provenance['independent_validation'], indent=2))
    elif args.command == 'plot':
        from paper_plot import plot_release
        print(plot_release(args.samples, args.output))
    elif args.command == 'apply':
        apply_hdf5(MassRemapping.load(args.grid), args.input, args.output,
                   allow_extrapolation=args.allow_extrapolation, mass_definition=args.mass_definition)
        print(args.output)


if __name__ == '__main__':
    main()
