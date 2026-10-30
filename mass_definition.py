"""The released HalfDome M200m-to-M200c interpolation convention.

This reproduces the historical catalog conversion, including its numerical mass
convention. It is not a general-purpose, unit-corrected NFW mass translator.
"""
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator


DATA = Path(__file__).resolve().parent / 'data'


@lru_cache(maxsize=1)
def _conversion():
    with np.load(DATA / 'm200c-conversion.npz') as data:
        return RegularGridInterpolator(
            (data['log_mass'], data['log1pz']), data['log_m200c'],
            method='linear', bounds_error=True)


def m200c_provenance():
    return json.loads((DATA / 'm200c-conversion.json').read_text())


def released_log10_m200c(redshift, log10_m200m):
    """Convert released M200m log masses to M200c log masses, both in Msun/h.

The July 2024 code interpolated log10(M200c) on a 200 x 30 grid in
log10(M200m) and ln(1+z). It passed the numeric Msun/h masses directly to CCL's
Msun interface. We preserve that release convention, not silently correct it.
The conversion table covers 0 <= z <= 4 and 12 <= log10(M200m) <= 16.5.
There is no extrapolation outside this table.
"""
    # Keep float32 arithmetic in ln(1+z) when supplied the stored redshift
    # column, matching the original notebook. Other inputs use float64.
    zin = np.asarray(redshift)
    if zin.dtype != np.dtype('float32'):
        zin = zin.astype(float)
    z, lm = np.broadcast_arrays(zin, np.asarray(log10_m200m, dtype=float))
    if not (np.isfinite(z).all() and np.isfinite(lm).all()):
        raise ValueError('Masses and redshifts must be finite')
    if np.any((z < 0) | (z > 4) | (lm < 12) | (lm > 16.5)):
        raise ValueError('Outside M200c conversion table: require 0 <= z <= 4 and 12 <= log10(M200m) <= 16.5')
    u = np.log(1 + z).astype(float)
    # Float32 ln(5) is slightly above the double-precision grid endpoint.
    u = np.clip(u, 0, np.log(5.0))
    coordinates = np.column_stack((lm.ravel(), u.ravel()))
    return _conversion()(coordinates).reshape(z.shape)
