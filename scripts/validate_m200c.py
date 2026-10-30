"""Check the conversion and complete RFoF -> M200c chain against release pairs."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mass_definition import m200c_provenance, released_log10_m200c
from remapping import MassRemapping


def metrics(predicted, actual, supported):
    error = np.abs(np.expm1(np.log(10) * (predicted - actual)))
    return {
        'rows': len(error), 'supported_rows': int(supported.sum()),
        'supported_fraction': float(supported.mean()),
        'median_absolute_relative_mass_error': float(np.median(error)),
        'p99_absolute_relative_mass_error': float(np.quantile(error, .99)),
        'p999_absolute_relative_mass_error': float(np.quantile(error, .999)),
        'max_absolute_relative_mass_error': float(error.max()),
        'supported_max_absolute_relative_mass_error': float(error[supported].max()),
    }


def validate(mapping, samples):
    # Samples store all columns as float64; recover the catalog's actual
    # float32 redshift arithmetic for the historical ln(1+z) step.
    z, raw, m200m, actual = samples.T
    z = z.astype(np.float32)
    supported = mapping.supported(z, raw)
    converted = released_log10_m200c(z, m200m)
    predicted = mapping.log10_m200c(z, raw, allow_extrapolation=True)
    return {
        'conversion_only': metrics(converted, actual, np.ones(len(z), dtype=bool)),
        'rfof_to_m200c': metrics(predicted, actual, supported),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--grid', type=Path, default=Path('data/recovered-grid'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mapping = MassRemapping.load(args.grid)
    with np.load(args.samples) as data:
        reports = {key: validate(mapping, data[key]) for key in data.files if key.startswith('sample_')}
    result = {
        'sample_sha256': hashlib.sha256(args.samples.read_bytes()).hexdigest(),
        'm200m_grid_sha256': hashlib.sha256((args.grid/'conversion_grid.npy').read_bytes()).hexdigest(),
        'conversion_provenance': m200c_provenance(),
        'checks': reports,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(reports, indent=2))


if __name__ == '__main__':
    main()
