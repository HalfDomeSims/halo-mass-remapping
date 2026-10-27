# Validation of the recovery

Recovery date: October 27, 2026. This validates a reconstruction against the
released catalogs. It does not establish that the original missing table has
been recovered bit for bit.

## Figure

The plotted histogram counts were computed directly from all 11 released HDF5
catalogs, for seeds 100, 102, …, 120, with read-only access. The four slices are
z=0–0.05, 1–1.05, 2–2.05, and 3–3.05. The 0.1-dex mass bins, shell volumes,
standard deviation over realizations, colors, and two-row layout follow the
July 2024 inputs. Calibrated histogram counts use the existing released M200m
column, independently of the fitted reconstruction. Theory curves use the
current pinned `hmf`/CAMB environment; historical theory-library bitwise identity
is not asserted.

## Mapping

The spline uses the original 70 × 65 interpolation geometry and the historical
correction to log mass. Training uses a stratified sample from seed 100, plus
high-mass pairs (log10 RFoF mass ≥14) from seeds 106–120. Seed 102 was used during
development. Seed 104 was excluded from training and used for a separate final
check after selecting the recovery method and support margin.

The support limit is 0.25 dex below the maximum training mass within 0.075 in
redshift of each original grid point. It is deliberately conservative near
sparsely populated spline cells. This is an empirical domain, not a mathematical
guarantee for every new halo. Samples are stratified and do not reflect the
population frequency of rare halos.

| Final independent check: seed 104 | Result |
| --- | ---: |
| Paired halos checked | 76,980 |
| Pairs inside default recovery support | 76,810 (99.779%) |
| Median absolute relative mass error, all pairs | 4.09 × 10⁻¹⁵ |
| 99th percentile, all pairs | 1.79 × 10⁻¹¹ |
| 99.9th percentile, all pairs | 7.25 × 10⁻¹⁰ |
| Maximum inside default support | 4.74 × 10⁻⁶ (4.74 ppm) |
| Maximum including unsupported tail | 2.11 × 10⁻³ (0.211%) |

The unsupported tail is rejected by default. The regularized continuation into
that tail is not a substitute for the original high-mass extrapolation. The
original raw NPY files would still be preferable for strict historical identity.

Machine-readable metrics, source-sample SHA256 hashes, the support envelope,
and training identities are in
[`data/recovered-grid/provenance.json`](data/recovered-grid/provenance.json).

## Executed checks

* The complete runnable notebook executes successfully from the pinned local
  environment and exports the paper-style PNG and PDF.
* Numerical tests check the log-mass convention against a known analytic result,
  boundary behavior, recovery of a known independent mapping, HDF5 row alignment
  and source preservation, partial-output cleanup, and the final released seed.
* The CLI creates a new aligned M200m sidecar. Production files were not opened
  for writing. The fnl=20 catalog was not processed.

See the GitHub Actions test run for ongoing numerical checks. Historical
M200m→M200c code is retained in provenance but has not been ported or validated as
part of the runnable mapper.
