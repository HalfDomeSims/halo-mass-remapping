from pathlib import Path

import numpy as np

from remapping import MassRemapping, validate_pairs


def test_final_independent_released_catalog():
    root = Path(__file__).resolve().parents[1]
    mapping = MassRemapping.load(root/'data/recovered-grid')
    samples = np.load(root/'data/release_final_validation.npz')['sample_104']
    assert 104 not in mapping.provenance['training_seeds']
    metrics = validate_pairs(mapping, samples)
    assert metrics['p99_absolute_relative_mass_error'] < 1e-8
    assert metrics['supported_fraction'] > 0.995
    assert metrics['supported_max_absolute_relative_mass_error'] < 1e-5
