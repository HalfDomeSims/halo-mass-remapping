import h5py
import numpy as np
import pytest

from scripts.sample_release import EDGES, Z_LEFT, descending_bound, histograms, stratified_sample


def test_small_catalog_slices_match_direct_selection(tmp_path):
    z = np.array([3.1, 3.05, 3, 2.03, 1.05, 1.05, 1, .05, 0], dtype='f8')
    lm = np.full(len(z), 13.0)
    with h5py.File(tmp_path/'input.h5', 'w') as f:
        f['redshift'], f['log10_m_rfof'], f['halo_mass_m200m'] = z, lm, 10**lm
        raw, calibrated, slices = histograms(f)
        for i, left in enumerate(Z_LEFT):
            indices = np.flatnonzero((z >= left) & (z <= left+.05))
            expected = np.histogram(lm[indices], bins=EDGES)[0]
            np.testing.assert_array_equal(raw[i], expected)
            np.testing.assert_array_equal(calibrated[i], expected)
        for threshold in [0, .05, 1.05, 4]:
            for strict in [False, True]:
                matching = np.flatnonzero(z < threshold if strict else z <= threshold)
                expected = matching[0] if len(matching) else len(z)
                assert descending_bound(f['redshift'], threshold, strict) == expected


@pytest.mark.parametrize('include_m200c', [False, True])
def test_stratified_reservoir_is_independent_of_chunk_boundaries(tmp_path, include_m200c):
    rng = np.random.default_rng(2)
    z = np.sort(rng.uniform(0, 4, 2000))[::-1]
    lm = rng.uniform(13, 13.1, len(z))
    with h5py.File(tmp_path/'input.h5', 'w') as f:
        f['redshift'], f['log10_m_rfof'], f['halo_mass_m200m'] = z, lm, 10**lm
        f['halo_mass_m200c'] = 10**(lm - .1 - .01*z)
        a = stratified_sample(f, per_bin=3, chunk_size=113, seed=71, include_m200c=include_m200c)
        b = stratified_sample(f, per_bin=3, chunk_size=2000, seed=71, include_m200c=include_m200c)
        np.testing.assert_array_equal(a, b)
        if include_m200c:
            np.testing.assert_allclose(a[:, 3], a[:, 2] - .1 - .01*a[:, 0], rtol=1e-14)
