# HalfDome halo mass remapping

Reproduce the purple RFoF / green calibrated halo mass-function figure from the
HalfDome N-body paper, and apply the July 2024 RFoF-to-Tinker08 M200m calibration
and the released M200m-to-M200c conversion.

[HalfDome paper](https://arxiv.org/abs/2407.17462) ·
[Public data](https://halfdomesims.github.io/data/)

The runnable notebook is **[HalfDome_mass_remapping.ipynb](HalfDome_mass_remapping.ipynb)**.
It runs entirely from the small files included here. No NERSC connection,
`nbodykit`, Aemulus emulator, or historical container is needed to run it.

**M200c support added October 29, 2026:** use `mapping.log10_m200c(...)` or
`run.py apply --mass-definition 200c` (choose `both` for both mass columns).
The restored historical conversion passes all 36 tests and the executed
notebook. On 76,980 independent seed-104 halos, the full RFoF→M200c chain agrees
with released masses to within 4.27 ppm inside the calibrated range. See
[M200c conversion](#m200c-conversion) for the method and historical unit
convention, and [validation results](VALIDATION.md#m200c-extension--october-29-2026)
for the measured errors and support limits.

## Method

This repository provides the calibration and figure used for the HalfDome halo
mass-function analysis:

* **Paper plot:** histogram counts from all 11 released HDF5 catalogs, including
  their existing calibrated masses. The grey theory curves use a pinned `hmf`
  installation.
* **Mass calibration:** a cubic spline fitted to released pairs of RFoF and
  calibrated masses, on a 70 × 65 grid in redshift and log mass. Independent
  validation and a conservative upper-mass support limit accompany the table.
* **M200c:** the archived NFW/Bhattacharya13 conversion of calibrated M200m,
  interpolated on the historical 200 × 30 grid. The small bundled table runs
  with NumPy/SciPy; CCL is needed only to rebuild it.

![Reproduced paper HMF figure](figures/reproduced-HMF.png)

## Install and run

Python 3.12 and [uv](https://docs.astral.sh/uv/) are recommended:

```bash
uv sync --frozen --extra notebook --extra test
uv run python run.py plot --samples data/release_samples.npz --output /tmp/halfdome-HMF
uv run python -m pytest -q
```

Open `HalfDome_mass_remapping.ipynb` with the project's Python kernel, starting in
the repository root. The notebook writes PNG/PDF figures under
`validation/` in the repository; change `OUTPUT` to
choose another destination. An executed copy is included.

To fit the calibration table from the included released-catalog pairs:

```bash
uv run python run.py recover-grid \
  --samples data/release_samples.npz \
  --tail-samples data/release_tail_samples.npz \
  --output /tmp/halfdome-recovered-grid
```

## Apply to a catalog

```python
import numpy as np
from remapping import MassRemapping

mapping = MassRemapping.load("data/recovered-grid")
# redshift and log10_m_rfof are aligned arrays from a halo catalog.
log10_m200m = mapping.log10_m200m(redshift, log10_m_rfof)
m200m = 10**log10_m200m  # Msun/h
log10_m200c = mapping.log10_m200c(redshift, log10_m_rfof)
m200c = 10**log10_m200c  # Msun/h, matching the released catalog convention
```

The historical correction is applied to **log mass**:

```text
log10(M200m) = log10(Mrfof) × [1 + spline(z, log10(Mrfof))]
```

It is not a multiplicative correction to linear mass. The released catalogs store
the RFoF log mass and redshift as float32; use those stored values when comparing
to the existing M200m fields. When preparing input from particle counts, the old
code used `log10(Length * M0 * 1e10)` and `1/Aemit - 1`, saved to float32. The
particle mass is approximately `1.9487254007799495e10 Msun/h`; the calibration
applies above 320 particles, about `6.236e12 Msun/h`.

For an HDF5 file containing `redshift` and `log10_m_rfof` (or `log10_m`), write a
new mass-only file with the same row order:

```bash
uv run python run.py apply --grid data/recovered-grid \
  --input halos.hdf5 --output calibrated-masses.hdf5

# Choose M200c, or write both aligned mass columns.
uv run python run.py apply --grid data/recovered-grid \
  --input halos.hdf5 --output calibrated-both.hdf5 --mass-definition both
```

The input is read-only. Existing output files are not replaced. The sidecar
contains `halo_mass_m200m` by default. `--mass-definition 200c` writes
`halo_mass_m200c`, and `both` writes both columns. Attributes identify the source,
units, mapping and conversion provenance, and row alignment. Failed conversions
leave no partial output.

## M200c conversion

The recovered July 24, 2024 inputs are archived in
[`provenance/m200c-history.json`](provenance/m200c-history.json). They first apply
the M200m calibration, then translate M200m to M200c with CCL's NFW mass
translator and the **Bhattacharya13** concentration relation for M200m.
This is a profile conversion, not a second abundance match to a M200c HMF.

The historical interpolator is linear in `log10(M200m)` and `ln(1+z)` and
returns `log10(M200c)`. Its grid covers log mass 12–16.5 and redshift 0–4.
Those limits are enforced even when `allow_extrapolation=True` permits the
RFoF spline's unvalidated continuation. For existing M200m masses alone:

```python
from mass_definition import released_log10_m200c

log10_m200c = released_log10_m200c(redshift, np.log10(m200m))
```

The release code passed numeric Msun/h masses directly to CCL's Msun interface,
without factors of h. This implementation preserves that historical convention
to reproduce released masses. It is not a general physical mass converter.
It also preserves the original float32 `ln(1+z)` arithmetic when given the
stored float32 redshift column.

The bundled table was rebuilt with a pinned modern CCL/CAMB environment; small
differences from the original environment are measured against released pairs.
The recipe, versions, and table hash are in
[`data/m200c-conversion.json`](data/m200c-conversion.json). To rebuild and validate:

```bash
uv sync --frozen --extra notebook --extra test --extra conversion
uv run python scripts/build_m200c_grid.py --output /tmp/m200c-conversion.npz
uv run python scripts/validate_m200c.py \
  --samples data/m200c-samples/release_samples.npz \
  --output /tmp/m200c-validation.json
```

The checked-in M200c samples contain `(redshift, log10_m_rfof, log10_m200m,
log10_m200c)`. They are validation inputs, not three-column `recover-grid`
training inputs. Re-export them on NERSC with
`scripts/sample_release.py --include-m200c --sample-seeds 102 104 --output NEW_DIRECTORY`.

## Calibration range and validation

See [VALIDATION.md](VALIDATION.md) and
[`data/recovered-grid/provenance.json`](data/recovered-grid/provenance.json) for
the measured errors. Seed 100 and high-mass pairs from seeds 106–120 train the
calibration. Seed 102 is development validation; seed 104 is reserved for a
separate final check. Samples are stratified in mass and redshift and therefore
overrepresent rare massive halos.

The default calibrated range ends 0.25 dex (five mass bins) below the local
maximum training mass. This reserves
room for cubic interpolation near the sparse tail. Calls above that limit fail
explicitly. Passing `allow_extrapolation=True`, or the equivalent CLI flag,
permits the regularized continuation and boundary clamping. That continuation
is outside the validated range.

The original interpolation grid covers redshift 0.025–3.475, although catalogs
extend to z=4. The historical code clamps to the grid boundaries; this is
preserved. New input below the original mass range is rejected by default.

`MassRemapping.load(directory)` accepts the three `conversion_grid*.npy` files
with the same grid convention. M200c adds the separately documented conversion
table. The fnl=20 production catalogs are outside the current validation.

## Provenance

* [`scripts/sample_release.py`](scripts/sample_release.py): reads the public
  catalogs under `/global/cfs/cdirs/cmb/halfdome/full_res/halos/`, exports the
  small stratified samples, and computes the four redshift-slice histograms.
* [`data/release_samples.json`](data/release_samples.json) and
  [`data/release_tail_samples.json`](data/release_tail_samples.json): source
  filenames, sizes, modification times, selected row intervals, and seeds.
* [`uv.lock`](uv.lock): package resolution; `requirements-lock.txt` additionally
  records the exact environment used for validation.
