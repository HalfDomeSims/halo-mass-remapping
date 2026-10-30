import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from mass_definition import released_log10_m200c
from remapping import MassRemapping, apply_hdf5
from scripts.validate_m200c import validate


ROOT = Path(__file__).resolve().parents[1]


def test_conversion_broadcasting_and_redshift_evolution():
    z = np.array([0., 1., 2., 4.])[:, None]
    mass = np.array([12., 14., 16.5])[None, :]
    result = released_log10_m200c(z, mass)
    assert result.shape == (4, 3)
    assert np.all(result < mass)
    assert np.all(np.diff(result - mass, axis=0) > 0)
    for i in range(4):
        for j in range(3):
            assert result[i, j] == released_log10_m200c(z[i, 0], mass[0, j])
    assert np.isfinite(released_log10_m200c(np.float32(4), 14))


@pytest.mark.parametrize('z,lm', [(np.nan, 13), (1, np.inf), (-.01, 13),
                                (4.01, 13), (1, 11.9), (1, 16.6)])
def test_conversion_rejects_invalid_or_unsupported_inputs(z, lm):
    with pytest.raises(ValueError):
        released_log10_m200c(z, lm)


def test_final_independent_m200c_release():
    mapping = MassRemapping.load(ROOT/'data/recovered-grid')
    assert 104 not in mapping.provenance['training_seeds']
    with np.load(ROOT/'data/m200c-samples/release_samples.npz') as data:
        report = validate(mapping, data['sample_104'])
    assert report['conversion_only']['max_absolute_relative_mass_error'] < 1e-5
    assert report['rfof_to_m200c']['supported_fraction'] > .995
    assert report['rfof_to_m200c']['supported_max_absolute_relative_mass_error'] < 2e-5


@pytest.mark.parametrize('definition', ['200c', 'both'])
def test_sidecar_columns_match_released_masses(tmp_path, definition):
    mapping = MassRemapping.load(ROOT/'data/recovered-grid')
    with np.load(ROOT/'data/m200c-samples/release_samples.npz') as data:
        sample = data['sample_104']
    sample = sample[mapping.supported(sample[:, 0], sample[:, 1])][::1000]
    src, dst = tmp_path/'input.h5', tmp_path/'output.h5'
    with h5py.File(src, 'w') as f:
        f['redshift'] = sample[:, 0].astype('f4')
        f['log10_m_rfof'] = sample[:, 1].astype('f4')
    before = src.read_bytes()
    apply_hdf5(mapping, src, dst, chunk_size=7, mass_definition=definition)
    assert src.read_bytes() == before
    with h5py.File(dst) as f:
        expected = {'halo_mass_m200c'}
        if definition == 'both':
            expected.add('halo_mass_m200m')
            np.testing.assert_allclose(f['halo_mass_m200m'][:], 10**sample[:, 2], rtol=1e-5)
        assert set(f) == expected
        np.testing.assert_allclose(f['halo_mass_m200c'][:], 10**sample[:, 3], rtol=2e-5)
        assert f.attrs['mass_definition'] == definition
        assert json.loads(f.attrs['m200c_conversion_provenance'])['concentration'] == 'Bhattacharya13'


def test_failed_m200c_conversion_leaves_no_partial_sidecar(tmp_path):
    mapping = MassRemapping.load(ROOT/'data/recovered-grid')
    src, dst = tmp_path/'input.h5', tmp_path/'output.h5'
    with h5py.File(src, 'w') as f:
        f['redshift'], f['log10_m_rfof'] = [1, 5], [13, 13]
    with pytest.raises(ValueError, match='M200c conversion table'):
        apply_hdf5(mapping, src, dst, chunk_size=1,
                   mass_definition='both', allow_extrapolation=True)
    assert not dst.exists()
    assert not list(tmp_path.glob('*.partial-*'))


def test_cli_writes_selected_definition(tmp_path, monkeypatch):
    import run
    import sys
    src, dst = tmp_path/'input.h5', tmp_path/'output.h5'
    with h5py.File(src, 'w') as f:
        f['redshift'], f['log10_m_rfof'] = [1], [13]
    monkeypatch.setattr(sys, 'argv', ['run.py', 'apply', '--grid', str(ROOT/'data/recovered-grid'),
                                     '--input', str(src), '--output', str(dst), '--mass-definition', '200c'])
    run.main()
    with h5py.File(dst) as f:
        assert set(f) == {'halo_mass_m200c'}
