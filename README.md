# modelOCP

MATLAB processing of Arbin battery measurements, diagonal-average OCP/OCV
estimation, and multi-site multi-reaction (MSMR) electrode fitting. The methods
and their relation to the article are documented in the companion modelOCV
repository's `doc/REPRODUCIBILITY.md`.

## Requirements

MATLAB R2023b or newer; Optimization Toolbox (`fmincon`), Global Optimization
Toolbox (`ga`), and Signal Processing Toolbox (`findpeaks`, `xcorr`). The numerical
unit tests do not need the experimental workbooks; the small synthetic fitting
test requires the optimization toolboxes. Preparation is verified on R2024b.

## Run the workflows

```matlab
setupModelOCP()
run('example_modelOCPelectrode_graphite.m')
% or: run('example_modelOCPelectrode_NMC.m')
```

The half-cell examples require the original experimental XLSX workbooks at the
configured `input/neg/` and `input/pos/` paths. Workbooks are not part of the public
source distribution. The small capacity references are in `input/reference/`.
Paths and figure formatting are local to this repository; no ECS 2027 directory
or externally installed `thesisFormat` is required. Do not recursively add this
repository and modelOCV together: call the appropriate setup when switching.

Entry points:

- `modelElectrodeCapacities.m`: calculate capacity from the stated host-mass
  convention and coating masses. Values are in mAh before the active-material
  fraction is applied; results remain in the workspace.
- `example_modelOCPelectrode_graphite.m` and `example_modelOCPelectrode_NMC.m`:
  electrode-specific settings, processing and MSMR fitting.
- `fun/modelOCPelectrode.m`: shared half-cell workflow.
- `modelOCV_commercial.m`, `modelOCV_NMC_graphite.m`, and
  `modelOCV_NMC_graphite_aged.m`: full-cell processing, requiring their configured
  private workbooks.

```matlab
[model, ocvrel, supplementaryResults, ocpFigures] = modelOCPelectrode(settings);
```

The examples display figures and return models in memory. Model exports are
opt-in: `settings.outputFile = ""` disables half-cell saving; the full-cell scripts
use `saveModel = false` by default. `getDataArbin_XLSX` separately caches imports
under `.cache/` beside the workbooks, or in an explicit third-argument directory.
It checks source sizes/timestamps, schema version and the counter setting. Raw workbooks are never
rewritten. To force reimport, pass `true` as the fourth argument.

## Numerical conventions and diagnostics

Potentials are in V versus Li; temperature inputs are in degrees Celsius and are
converted to kelvin internally. Capacities in the processing/model schema are Ah;
capacity reference MAT files store `Q_r` in mAh. MSMR `w` is dimensionless in the
`F*(U-U0)/(R*T*w)` convention. `dZ` is the positive magnitude `-dtheta/dU` [1/V].
The legacy gallery output `dZj` excludes `F/(R*T)`; multiply by that factor for
physical gallery derivative magnitudes.

Current integration retains the original right-endpoint rectangular rule in Ah.
Decreasing test times are rejected rather than integrated across a reset. The
configured temperature is propagated to processed branches and the MSMR fit.

`fitMSMR` performs GA followed by `fmincon`, restoring the caller's random stream.
The parameter bounds, minimum gallery-potential separation, sum of site fractions,
and positive fitted interval are constrained. A locally converged feasible fit
can be accepted even when GA used its generation budget. Both stages' exit flags,
objectives and final constraint violation remain in `fitInfo`.

```matlab
settings.optimizer.seed = 42;
settings.optimizer.gaOptions = {'Display', 'off', 'MaxGenerations', 200};
settings.optimizer.fminconOptions = {'Display', 'off', 'MaxIterations', 1500};
```

Check `ocvrel.MSMR.(branch).fitInfo` and voltage residuals before interpreting a
calibration. Voltage agreement does not identify unobserved transitions or prove
absolute stoichiometry. Hidden-gallery and active-material corrections remain
explicit downstream assumptions; their physical interpretation has not been
changed by this preparation. Initial gallery fractions follow the same
ascending-potential order as their peak locations.

## Figures and supplementary snapshot

`ocpFigures` is 7-by-1, indexed in window order:

1. Measured half-cell voltage.
2. Relative processed branches and diagonal-average OCP.
3. Initial cutoff-based absolute coordinates.
4. Measured branches and separate MSMR fits on relative coordinates.
5. Branch fits and mean-parameter curve on absolute coordinates.
6. Mean-parameter composition, occupations and positive gallery derivatives.
7. Relative and initial cutoff-based measured OCP, plus the MSMR mean curve
   after the configured hidden-gallery transformation.

`supplementaryResults` retains values and figure handles 1--6 before the
hidden-gallery/active-material corrections. The returned measured LUT is remapped
later to the final interval, so figure 7 is not a plot of that later remapping.
The fit to the diagonal-average curve and the arithmetic mean of branch
parameters are distinct objects. Residual fields use model minus measurement.

## Public numerical references

`input/reference/ocpReferences.mat` is a numerical snapshot of the existing local
model outputs, without machine-specific function handles. It allows inspection
without raw workbooks; it does not replace the private inputs required to rerun
experimental fitting. Its positive absolute calibration differs from modelOCV's
frozen article inputs. See `input/reference/README.md`.

## Tests and publication

```matlab
runModelOCPTests()
```

Tests cover interpolation, capacity extraction, MSMR conservation/derivatives and
saturation, asymmetric gallery ordering, import/cache behavior, phase processing,
and a small synthetic optimizer case. Tiny synthetic workbook fixtures are kept
in the ignored `.build/` directory for inspection. See `VALIDATION.md` for actual
results and limitations. Large experimental refits are separate scientific runs.

Code contributions are CC BY-SA 4.0; numerical references and capacity values are
CC BY 4.0. Preserve inherited BMS3 attribution; see `LICENSE.md`,
`DATA_LICENSE.md`, and `THIRD_PARTY_NOTICES.md`. `PUBLISHING.md` explains the
existing repository history and private-data exclusions. No new release identifier
or archived DOI is asserted.
