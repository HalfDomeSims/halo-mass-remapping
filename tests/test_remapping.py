import h5py
import numpy as np
import pytest

from remapping import MassRemapping, apply_hdf5, historical_axes, fit_released_pairs, validate_pairs


def test_historical_formula_and_paired_not_cartesian_evaluation(tmp_path):
    z, m = historical_axes()
    corr = 0.001*z[:, None] + 0.002*(m[None, :]-13)
    mapping = MassRemapping(z, m, corr)
    zi, mi = np.array([0.4, 1.7]), np.array([13.2, 14.1])
    expected = mi*(1+0.001*zi+0.002*(mi-13))
    np.testing.assert_allclose(mapping.log10_m200m(zi, mi), expected, rtol=1e-14)
    mapping.save(tmp_path)
    np.testing.assert_allclose(MassRemapping.load(tmp_path).log10_m200m(zi, mi), expected, rtol=1e-14)


def test_historical_boundary_clamping():
    z, m = historical_axes()
    mapping = MassRemapping(z, m, np.broadcast_to(z[:, None]*0.001, (len(z), len(m))))
    np.testing.assert_allclose(mapping.log10_m200m([0, 4], 13), 13*(1+z[[0, -1]]*.001))
    with pytest.raises(ValueError):
        mapping.log10_m200m(1, 12)
    with pytest.raises(ValueError):
        mapping.log10_m200m(np.nan, 13)


def test_recovery_against_an_independently_known_mapping():
    rng = np.random.default_rng(903)
    z, m = historical_axes()
    zz, mm = np.meshgrid(z, m, indexing='ij')
    samples = np.column_stack((zz.ravel(), mm.ravel()))
    truth = lambda zi, mi: mi*(1-0.0004*zi+0.0002*(mi-13))
    training = np.column_stack((samples, truth(samples[:, 0], samples[:, 1])))
    recovered = fit_released_pairs(training)
    zi = rng.uniform(z[0], z[-1], 1000)
    mi = rng.uniform(m[0], m[-1], 1000)
    heldout = np.column_stack((zi, mi, truth(zi, mi)))
    assert validate_pairs(recovered, heldout)['max_absolute_relative_mass_error'] < 1e-8


def test_sidecar_preserves_source_and_row_order(tmp_path):
    z, m = historical_axes()
    mapping = MassRemapping(z, m, np.full((len(z), len(m)), -0.001))
    src, dst = tmp_path/'input.h5', tmp_path/'output.h5'
    redshift, mass = np.array([2, 0.5, 3.8]), np.array([14, 13, 15.2])
    with h5py.File(src, 'w') as f:
        f['redshift'], f['log10_m_rfof'] = redshift, mass
    before = src.read_bytes()
    apply_hdf5(mapping, src, dst, chunk_size=2)
    assert src.read_bytes() == before
    with h5py.File(dst) as f:
        np.testing.assert_allclose(f['halo_mass_m200m'][:], 10**(mass*0.999), rtol=1e-13)
    with pytest.raises(FileExistsError):
        apply_hdf5(mapping, src, dst)


def test_unsupported_mass_does_not_leave_partial_catalog(tmp_path):
    z, m = historical_axes()
    mapping = MassRemapping(z, m, np.zeros((len(z), len(m))))
    src, dst = tmp_path/'input.h5', tmp_path/'output.h5'
    with h5py.File(src, 'w') as f:
        f['redshift'], f['log10_m_rfof'] = [1, 1], [13, 17]
    with pytest.raises(ValueError):
        apply_hdf5(mapping, src, dst, chunk_size=1)
    assert not dst.exists()
    assert not list(tmp_path.glob('*.partial-*'))
