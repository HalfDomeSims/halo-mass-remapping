"""Build the small runnable notebook; the separate history archive is unchanged."""
from pathlib import Path

import nbformat as nbf


root = Path(__file__).resolve().parents[1]
nb = nbf.v4.new_notebook()
nb.cells = [
    nbf.v4.new_markdown_cell('''# HalfDome halo mass calibration and paper figure

Runnable recovery of the July 2024 purple RFoF / green Tinker08-calibrated HMF figure.
The figure uses histogram counts from all 11 released HDF5 catalogs. It does not
depend on recovery of the missing original mapping tables.

The mapping below was **recovered from released pairs**, not copied from the
original `conversion_grid*.npy` files. Its provenance and independent validation
are shown explicitly. The archived notebook inputs are preserved separately.'''),
    nbf.v4.new_code_cell('''from pathlib import Path
import json
import tempfile
import h5py
import numpy as np
from IPython.display import Image, display
from paper_plot import plot_release
from remapping import MassRemapping, apply_hdf5, validate_pairs

ROOT = Path.cwd()
DATA = ROOT / "data"
OUTPUT = Path.home() / "projects/halfdome/recovered-remapping-20261027/validation"
OUTPUT.mkdir(parents=True, exist_ok=True)
samples = np.load(DATA / "release_samples.npz")
print("Released seeds:", samples["seeds"].tolist())'''),
    nbf.v4.new_markdown_cell('''## Reproduce the paper HMF figure

Redshift slices are 0–0.05, 1–1.05, 2–2.05, and 3–3.05. Error bars are the
standard deviation over the 11 realizations, as in the recovered inputs.
Grey curves use the installed HMF package with the HalfDome cosmology.'''),
    nbf.v4.new_code_cell('''figure = plot_release(DATA / "release_samples.npz", OUTPUT / "reproduced-HMF")
display(Image(filename=str(figure)))'''),
    nbf.v4.new_markdown_cell('''## Apply the recovered mapping

The historical formula is `log10_M200m = log10_Mrfof * (1 + correction(z, log10_Mrfof))`.
It is not a multiplicative correction to linear mass. Training excludes seed 102;
the following report uses its independent released mass pairs. Both samples are
stratified in mass and redshift and therefore emphasize rare massive halos.

The default call rejects masses above the recovery support. Setting
`allow_extrapolation=True` reproduces the historical boundary behavior, but the
unobserved high-mass continuation is not established by recovery.'''),
    nbf.v4.new_code_cell('''mapping = MassRemapping.load(DATA / "recovered-grid")
final_validation = np.load(DATA / "release_final_validation.npz")
report = validate_pairs(mapping, final_validation["sample_104"])
print("Final independent check: seed 104")
print(json.dumps(report, indent=2))
print("Original table recovered:", mapping.provenance["original_grid_recovered"])'''),
    nbf.v4.new_code_cell('''# A small independent example within the supported mass range.
example = samples["sample_102"]
example = example[example[:, 1] < 13.3][:20]
predicted = mapping.log10_m200m(example[:, 0], example[:, 1])
relative_error = np.expm1(np.log(10) * (predicted - example[:, 2]))
print("Maximum relative mass error in example:", np.max(np.abs(relative_error)))'''),
    nbf.v4.new_markdown_cell('''## Write an aligned HDF5 mass sidecar

This demonstration uses a temporary 20-row catalog. Real inputs need `redshift`
and `log10_m_rfof` (or the historical `log10_m`) in the same row order.
The output is a new HDF5 file containing `halo_mass_m200m` in solar masses/h,
plus source, units, mapping provenance, and row-alignment attributes. Existing
catalogs are never opened for writing.'''),
    nbf.v4.new_code_cell('''with tempfile.TemporaryDirectory(prefix="halfdome-example-") as temporary:
    source = Path(temporary) / "input.h5"
    destination = Path(temporary) / "m200m.h5"
    with h5py.File(source, "w") as f:
        f["redshift"] = example[:, 0]
        f["log10_m_rfof"] = example[:, 1]
    apply_hdf5(mapping, source, destination, chunk_size=7)
    with h5py.File(destination, "r") as f:
        np.testing.assert_allclose(f["halo_mass_m200m"][:], 10**predicted)
        print("Sidecar rows:", len(f["halo_mass_m200m"]))
        print("Units:", f.attrs["units"])'''),
    nbf.v4.new_markdown_cell('''The source README records commands, provenance, validation limits, and the
relationship to the earlier June Watson/Aemulus scripts. This recovery provides
M200m; it does not silently substitute a different concentration prescription
for the historical M200m→M200c conversion.'''),
]
nb.metadata.kernelspec = {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'}
nbf.write(nb, root / 'HalfDome_mass_remapping.ipynb')
