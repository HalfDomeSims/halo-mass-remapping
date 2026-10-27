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


@pytest.mark.parametrize('chunk_size', [-1, 0, 1.5, True])
def test_invalid_chunk_size_never_writes_a_catalog(tmp_path, chunk_size):
    z, m = historical_axes()
    mapping = MassRemapping(z, m, np.zeros((len(z), len(m))))
    src, dst = tmp_path/'input.h5', tmp_path/'output.h5'
    with h5py.File(src, 'w') as f:
        f['redshift'], f['log10_m_rfof'] = [1, 2], [13, 14]
    before = src.read_bytes()
    with pytest.raises(ValueError, match='positive integer'):
        apply_hdf5(mapping, src, dst, chunk_size=chunk_size)
    assert src.read_bytes() == before
    assert not dst.exists()
    assert not list(tmp_path.glob('*.partial-*'))


def test_report_support_agrees_with_conversion_domain():
    z, m = historical_axes()
    mapping = MassRemapping(z, m, np.zeros((len(z), len(m))), {
        'support': {'redshift': z.tolist(), 'upper_log_mass': [17.0]*len(z)}})
    samples = np.array([[1, 13, 13], [5, 13, 13], [1, 12, 12], [1, 17, 17]])
    report = validate_pairs(mapping, samples)
    assert report['supported_rows'] == 1
    np.testing.assert_array_equal(mapping.supported(samples[:, 0], samples[:, 1]), [True, False, False, False])
    for row in samples[1:]:
        with pytest.raises(ValueError):
            mapping.log10_m200m(row[0], row[1])


@pytest.mark.parametrize('bad_mass', [1000.0, -1000.0])
def test_unrepresentable_linear_mass_cannot_create_sidecar(tmp_path, bad_mass):
    z, m = historical_axes()
    mapping = MassRemapping(z, m, np.zeros((len(z), len(m))))
    src, dst = tmp_path/'input.h5', tmp_path/'output.h5'
    with h5py.File(src, 'w') as f:
        f['redshift'], f['log10_m_rfof'] = [1], [bad_mass]
    with pytest.raises(ValueError, match='finite and positive'):
        apply_hdf5(mapping, src, dst, allow_extrapolation=True)
    assert not dst.exists()
    assert not list(tmp_path.glob('*.partial-*'))


@pytest.mark.parametrize('samples,kwargs', [
    (np.empty((0, 3)), {}),
    ([[1, 13, 13]], {'regularization': 0}),
    ([[1, 13, 13]], {'regularization': -1}),
    ([[1, 13, 13]], {'regularization': np.nan}),
    ([[1, 13, 13]], {'support_margin': -0.1}),
    ([[5, 13, 13]], {}),
    ([[1, 0, 13]], {}),
])
def test_invalid_training_requests_fail_before_solving(samples, kwargs):
    with pytest.raises(ValueError):
        fit_released_pairs(samples, **kwargs)


def test_training_provenance_records_actual_seeds(tmp_path, monkeypatch):
    import json
    import run as cli
    import sys

    base, tail, output = tmp_path/'base.npz', tmp_path/'tail.npz', tmp_path/'grid'
    np.savez(base, sample_100=[[1, 13, 13]], sample_102=[[2, 13, 13]])
    np.savez(tail, sample_100=[[1, 14, 14]], sample_104=[[1, 15, 15]], sample_106=[[3, 13, 13]])
    z, m = historical_axes()
    def fit(training, **kwargs):
        np.testing.assert_array_equal(training, [[1, 13, 13], [3, 13, 13]])
        return MassRemapping(z, m, np.zeros((len(z), len(m))))
    monkeypatch.setattr(cli, 'fit_released_pairs', fit)
    monkeypatch.setattr(sys, 'argv', ['run.py', 'recover-grid', '--samples', str(base), '--tail-samples', str(tail), '--output', str(output)])
    cli.main()
    assert json.loads((output/'provenance.json').read_text())['training_seeds'] == [100, 106]
