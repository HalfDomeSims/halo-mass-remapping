"""HalfDome's historical log-mass correction, with explicit provenance.

The table stores delta(log10 M)/log10 M, not delta M/M.
Original grid files can be loaded directly. Recovered grids fitted to released
catalog pairs are a distinct product and carry their validation and support.
"""
from dataclasses import dataclass, field
import json
import os
from pathlib import Path
from uuid import uuid4

import h5py
import numpy as np
from scipy.interpolate import NdBSpline, RectBivariateSpline
from scipy.sparse import diags, eye, kron, vstack
from scipy.sparse.linalg import spsolve

PARTICLE_MASS = 1.9487254007799495e10
LOG_MASS_MIN = np.log10(320 * PARTICLE_MASS)


def historical_axes():
    return np.arange(0, 3.5, 0.05) + 0.025, np.arange(LOG_MASS_MIN, 16, 0.05)


@dataclass
class MassRemapping:
    redshift: np.ndarray
    log_mass: np.ndarray
    correction: np.ndarray
    provenance: dict = field(default_factory=dict)

    def __post_init__(self):
        self.redshift = np.asarray(self.redshift, dtype=float)
        self.log_mass = np.asarray(self.log_mass, dtype=float)
        self.correction = np.asarray(self.correction, dtype=float)
        if self.correction.shape != (len(self.redshift), len(self.log_mass)):
            raise ValueError('Correction shape must be (redshift, log_mass)')
        if min(len(self.redshift), len(self.log_mass)) < 4:
            raise ValueError('Cubic interpolation needs at least four points on each axis')
        if not (np.all(np.diff(self.redshift) > 0) and np.all(np.diff(self.log_mass) > 0)):
            raise ValueError('Grid axes must be strictly increasing')
        if not np.all(np.isfinite(self.correction)):
            raise ValueError('Correction grid must be finite')
        self._spline = RectBivariateSpline(self.redshift, self.log_mass, self.correction, s=0)

    def log10_m200m(self, redshift, log10_m_rfof, *, allow_extrapolation=False):
        z, lm = np.broadcast_arrays(np.asarray(redshift, dtype=float), np.asarray(log10_m_rfof, dtype=float))
        if not (np.all(np.isfinite(z)) and np.all(np.isfinite(lm))):
            raise ValueError('Masses and redshifts must be finite')
        if not allow_extrapolation:
            if np.any((z < 0) | (z > 4) | (lm < self.log_mass[0]) | (lm > self.log_mass[-1])):
                raise ValueError('Input outside the historical calibrated domain; pass allow_extrapolation=True to use historical boundary clamping')
            support = self.provenance.get('support')
            if support:
                upper = np.interp(z, support['redshift'], support['upper_log_mass'])
                if np.any(lm > upper):
                    raise ValueError('Mass above the released-catalog recovery support; original grid required for unvalidated extrapolation')
        # FITPACK clamps outside the tabulated axes, as in the original code.
        corr = self._spline.ev(z.ravel(), lm.ravel()).reshape(z.shape)
        return lm * (1.0 + corr)

    def save(self, directory):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        np.save(directory / 'conversion_grid.npy', self.correction)
        np.save(directory / 'conversion_grid_z.npy', self.redshift)
        np.save(directory / 'conversion_grid_m.npy', self.log_mass)
        (directory / 'provenance.json').write_text(json.dumps(self.provenance, indent=2) + '\n')

    @classmethod
    def load(cls, directory):
        directory = Path(directory)
        meta = directory / 'provenance.json'
        provenance = json.loads(meta.read_text()) if meta.exists() else {'method': 'provided_historical_grid'}
        return cls(np.load(directory / 'conversion_grid_z.npy'), np.load(directory / 'conversion_grid_m.npy'),
                   np.load(directory / 'conversion_grid.npy'), provenance)


def fit_released_pairs(samples, *, regularization=1e-6, support_margin=0.25):
    """Recover a spline on the historical axes from (z, logM_RF, logM_200m).

    Tiny curvature regularization completes regions with no released halos.
    Those regions are excluded by the default support check and must not be
    mistaken for recovery of the original high-mass extrapolation.
    """
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 2 or samples.shape[1] != 3 or not np.isfinite(samples).all():
        raise ValueError('Expected finite N x 3 paired samples')
    z_grid, m_grid = historical_axes()
    zero = RectBivariateSpline(z_grid, m_grid, np.zeros((len(z_grid), len(m_grid))))
    knots = zero.get_knots()
    coords = np.column_stack((np.clip(samples[:, 0], z_grid[0], z_grid[-1]),
                              np.clip(samples[:, 1], m_grid[0], m_grid[-1])))
    design = NdBSpline.design_matrix(coords, knots, (3, 3))
    target = samples[:, 2] / samples[:, 1] - 1
    nz, nm = len(z_grid), len(m_grid)
    design.resize((len(samples), nz*nm))
    dm = diags([np.ones(nm-2), -2*np.ones(nm-2), np.ones(nm-2)], [0, 1, 2], shape=(nm-2, nm))
    dz = diags([np.ones(nz-2), -2*np.ones(nz-2), np.ones(nz-2)], [0, 1, 2], shape=(nz-2, nz))
    penalty = vstack((kron(eye(nz), dm), kron(dz, eye(nm))), format='csr')
    normal = design.T @ design + regularization**2 * (penalty.T @ penalty + eye(nz*nm)*1e-4)
    rhs = design.T @ target
    coefficients = spsolve(normal.tocsc(), rhs)
    # Residual correction reduces normal-equation roundoff in sparse tail bins.
    for _ in range(2):
        gradient = design.T @ (target - design @ coefficients) - regularization**2 * (
            penalty.T @ (penalty @ coefficients) + 1e-4*coefficients)
        coefficients += spsolve(normal.tocsc(), gradient)
    spline = NdBSpline(knots, coefficients.reshape(nz, nm), (3, 3))
    zz, mm = np.meshgrid(z_grid, m_grid, indexing='ij')
    grid = spline(np.column_stack((zz.ravel(), mm.ravel()))).reshape(nz, nm)
    # A local ceiling identifies supported masses, rather than legitimizing the
    # arbitrary regularized continuation into unobserved high-mass regions.
    upper = []
    for z in z_grid:
        local = samples[np.abs(samples[:, 0]-z) < 0.075, 1]
        upper.append(max(float(m_grid[0]), float(local.max())-support_margin) if len(local) else float(m_grid[0]))
    return MassRemapping(z_grid, m_grid, grid, {
        'method': 'spline_recovered_from_released_mass_pairs',
        'original_grid_recovered': False,
        'training_rows': len(samples), 'curvature_regularization': regularization,
        'support_margin_dex': support_margin,
        'support': {'redshift': z_grid.tolist(), 'upper_log_mass': upper},
        'historical_formula': 'log10_M200m = log10_Mrfof * (1 + spline(z, log10_Mrfof))',
    })


def validate_pairs(mapping, samples):
    samples = np.asarray(samples)
    predicted = mapping.log10_m200m(samples[:, 0], samples[:, 1], allow_extrapolation=True)
    relative = np.expm1(np.log(10) * (predicted - samples[:, 2]))
    absolute = np.abs(relative)
    result = {'rows': len(samples), 'median_absolute_relative_mass_error': float(np.median(absolute)),
              'p99_absolute_relative_mass_error': float(np.quantile(absolute, 0.99)),
              'p999_absolute_relative_mass_error': float(np.quantile(absolute, 0.999)),
              'max_absolute_relative_mass_error': float(absolute.max())}
    worst = int(np.argmax(absolute))
    result['worst_input'] = samples[worst, :2].tolist()
    support = mapping.provenance.get('support')
    if support:
        upper = np.interp(samples[:, 0], support['redshift'], support['upper_log_mass'])
        inside = samples[:, 1] <= upper
        result['supported_rows'] = int(inside.sum())
        result['supported_fraction'] = float(inside.mean())
        result['supported_max_absolute_relative_mass_error'] = float(absolute[inside].max()) if inside.any() else None
    return result


def apply_hdf5(mapping, source, destination, *, chunk_size=100_000, allow_extrapolation=False):
    """Create a new mass-only sidecar with row correspondence to the input."""
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Destination must be a separate sidecar')
    if destination.exists():
        raise FileExistsError(destination)
    temporary = destination.with_name(destination.name + '.partial-' + uuid4().hex)
    try:
        _write_sidecar(mapping, source, temporary, chunk_size, allow_extrapolation)
        # Link atomically without ever replacing an existing destination.
        os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _write_sidecar(mapping, source, destination, chunk_size, allow_extrapolation):
    with h5py.File(source, 'r') as src:
        raw_name = 'log10_m_rfof' if 'log10_m_rfof' in src else 'log10_m'
        z, lm = src['redshift'], src[raw_name]
        if z.ndim != 1 or lm.ndim != 1 or len(z) != len(lm):
            raise ValueError('Expected aligned one-dimensional redshift and log mass columns')
        with h5py.File(destination, 'x') as dst:
            dst.attrs['source_catalog'] = str(source.resolve())
            dst.attrs['units'] = 'Msun/h'
            dst.attrs['mapping_provenance'] = json.dumps(mapping.provenance)
            dst.attrs['row_alignment'] = 'Same row order and row count as source_catalog'
            dst.attrs['allow_extrapolation'] = allow_extrapolation
            masses = dst.create_dataset('halo_mass_m200m', shape=lm.shape, dtype='f8')
            for start in range(0, len(z), chunk_size):
                end = min(start + chunk_size, len(z))
                values = mapping.log10_m200m(z[start:end], lm[start:end], allow_extrapolation=allow_extrapolation)
                masses[start:end] = 10**values
