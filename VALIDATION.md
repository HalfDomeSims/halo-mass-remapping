# Validation

Validation date: October 27, 2026. The calibration is checked against independent
released HalfDome catalogs.

## Figure

The plotted histogram counts were computed directly from all 11 released HDF5
catalogs, for seeds 100, 102, …, 120, with read-only access. The four slices are
z=0–0.05, 1–1.05, 2–2.05, and 3–3.05. The 0.1-dex mass bins, shell volumes,
standard deviation over realizations, colors, and two-row layout follow the
July 2024 inputs. Calibrated histogram counts use the existing released M200m
column, independently of the fitted calibration. Theory curves use the pinned
`hmf`/CAMB environment.

## Mapping

The spline uses the original 70 × 65 interpolation geometry and the historical
correction to log mass. Training uses a stratified sample from seed 100, plus
high-mass pairs (log10 RFoF mass ≥14) from seeds 106–120. Seed 102 was used during
development. Seed 104 was excluded from training and used for a separate final
check after selecting the fitting method and support margin.

The support limit is 0.25 dex below the maximum training mass within 0.075 in
redshift of each original grid point. It is deliberately conservative near
sparsely populated spline cells. This is an empirical domain, not a mathematical
guarantee for every new halo. Samples are stratified and do not reflect the
population frequency of rare halos.

| Final independent check: seed 104 | Result |
| --- | ---: |
| Paired halos checked | 76,980 |
| Pairs inside default calibrated range | 76,810 (99.779%) |
| Median absolute relative mass error, all pairs | 4.09 × 10⁻¹⁵ |
| 99th percentile, all pairs | 1.79 × 10⁻¹¹ |
| 99.9th percentile, all pairs | 7.25 × 10⁻¹⁰ |
| Maximum inside default support | 4.74 × 10⁻⁶ (4.74 ppm) |
| Maximum including unsupported tail | 2.11 × 10⁻³ (0.211%) |

The unsupported tail is rejected by default. The optional regularized
continuation into that tail is outside the validated range.

Machine-readable metrics, source-sample SHA256 hashes, the support envelope,
and training identities are in
[`data/recovered-grid/provenance.json`](data/recovered-grid/provenance.json).

## Executed checks

* The complete runnable notebook executes successfully from the pinned local
  environment and exports the paper-style PNG and PDF.
* Numerical tests check the log-mass convention against a known analytic result,
  boundary behavior, recovery of a known independent mapping, HDF5 row alignment
  and source preservation, partial-output cleanup, and the final released seed.
* Input and reporting checks reject invalid chunk sizes and unrepresentable
  masses, use the same domain mask for conversion and validation, and record
  the training seeds actually loaded. Sampling tests check slice boundaries
  and independence from chunk boundaries.
* The CLI creates a new aligned M200m sidecar. Production files were not opened
  for writing. The fnl=20 catalog was not processed.

See the GitHub Actions test run for ongoing numerical checks.

## M200c extension — October 29, 2026

The July 24, 2024 IPython session 749, input 96, supplies the M200m→M200c
recipe: an NFW conversion using CCL's Bhattacharya13 concentration relation,
evaluated on 200 log-mass points from 12 to 16.5 and 30 ln(1+z) points from
0 to ln(5). The final catalog write is input 98. These inputs are preserved in
[`provenance/m200c-history.json`](provenance/m200c-history.json).

The original grid was rebuilt using CCL 3.3.6 and CAMB 2.0.5, retaining the
release's numerical convention of passing Msun/h values directly to CCL's
Msun interface without h rescaling. This is explicitly a reproduction of the
released masses. The rebuild is close to, but not bitwise identical to, the
historical conversion. No fit to M200c catalog values was performed.

Fresh four-column samples retain the released RFoF, M200m, M200c, and redshift
values from seeds 102 and 104. Each uses the same reservoir selection (32
halos per 0.05-redshift × 0.05-dex bin) as the existing sampler. The conversion
uses the stored float32 redshifts in the historical `ln(1+z)` calculation.
Only read access was used on the public catalogs.

| Independent check: seed 104 | Result |
| --- | ---: |
| Paired halos checked | 76,980 |
| Pairs inside default RFoF calibration range | 76,810 (99.779%) |
| M200m→M200c maximum relative mass error, all pairs | 2.47 × 10⁻⁶ (2.47 ppm) |
| RFoF→M200c median relative mass error, all pairs | 4.79 × 10⁻⁷ |
| RFoF→M200c 99th percentile, all pairs | 1.91 × 10⁻⁶ |
| RFoF→M200c maximum inside default support | 4.27 × 10⁻⁶ (4.27 ppm) |
| RFoF→M200c maximum including unsupported tail | 2.10 × 10⁻³ (0.210%) |

Seed 102 provides an additional 77,238 pairs: the conversion-only maximum is
2.60 ppm, and the full-chain supported maximum is 2.45 ppm. The unsupported
RFoF tail remains excluded by default. The M200c interpolation table never
extrapolates, including when RFoF extrapolation is explicitly requested.

[`data/m200c-validation.json`](data/m200c-validation.json) records both checks,
sample and grid hashes, and the complete conversion provenance. The M200m grid
and its original validation data have not changed. The original paper plot
continues to use M200m.

New regression checks cover the independent released M200c sample, scalar and
broadcast inputs, table limits, redshift evolution, M200c/both sidecar selection,
source preservation, row alignment, provenance, CLI selection, and cleanup when
a later chunk cannot be converted. Runtime conversion uses the bundled table;
the optional CCL dependency is required only for regeneration.

The extended suite passes all 36 tests, and the updated notebook executes
successfully, including its M200c validation and two-column HDF5 sidecar example.
