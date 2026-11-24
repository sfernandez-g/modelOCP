# modelOCP

MATLAB tools for processing Arbin half-cell measurements, estimating electrode
open-circuit potential (OCP) by diagonal averaging, and fitting a multi-site
multi-reaction (MSMR) model. The examples configure graphite and NMC811 electrodes.

## Requirements

- MATLAB R2023b or later.
- Optimization Toolbox (`fmincon`).
- Global Optimization Toolbox (`ga`).
- Signal Processing Toolbox (`findpeaks`, `xcorr`).

## Usage

From the repository folder, run either example:

```matlab
run('example_modelOCPelectrode_graphite.m')
% Alternatively:
run('example_modelOCPelectrode_NMC.m')
```

Each script adds `fun/` to the MATLAB path and defines the experimental inputs,
temperature, fitting parameters, and figure settings. Update `settings.inputFile`
to point to your Arbin XLSX workbook. The original experimental workbooks are
not distributed; the examples require them or compatible replacement data.
Review the step/cycle configuration in `fun/modelOCPelectrode.m` when using a
different measurement protocol. Use this repository's helper functions when
switching between this project and the companion modelOCV repository.

| Entry point | Purpose |
| --- | --- |
| `example_modelOCPelectrode_graphite.m` | Process and fit the graphite half-cell. |
| `example_modelOCPelectrode_NMC.m` | Process and fit the NMC811 half-cell. |
| `modelElectrodeCapacities.m` | Calculate reference capacities from coating masses and the stated host-mass convention. |
| `fun/modelOCPelectrode.m` | Run the shared processing and fitting workflow with a settings structure. |

For direct function use, first add `fun/` to the path. Each function documents
its syntax, inputs, and outputs in its header.

```matlab
addpath(fullfile(pwd, 'fun'))
[model, ocvrel, supplementaryResults, ocpFigures] = modelOCPelectrode(settings);
```

The examples return results in memory and display seven figures. Set
`settings.outputFile` to a MAT file path to save the model; the default empty
string disables saving. The destination folder must already exist.
`getDataArbin_XLSX` caches imported data under `.cache/` beside the workbooks,
or in its optional third-argument directory. It checks source sizes and
timestamps, the cache schema, and the counter setting before reusing a cache.
Pass `true` as its fourth argument to force reimport. Raw workbooks are not
rewritten.

## Units and modeling conventions

- Electrode potentials are in V versus Li. Temperature inputs are in degrees
  Celsius and are converted to kelvin internally.
- Processing/model capacities are in Ah. The capacity reference files store
  `Q_r` in mAh, before multiplication by the active-material fraction.
- MSMR `w` is dimensionless in the `F*(U-U0)/(R*T*w)` convention.
  `compositionMSMR` returns `dZ = -dtheta/dU`, a positive magnitude in 1/V.
  Its legacy gallery output `dZj` excludes `F/(R*T)`; multiply by this factor
  to obtain physical gallery derivative magnitudes.
- Current integration uses a right-endpoint rectangular rule in Ah. Decreasing
  test times are rejected. The configured temperature is propagated to the
  processed branches and the MSMR fit.

The capacity script uses the documented host-mass convention, including the
omission of lithium from the NMC host molar mass. Active-material fractions,
endpoint anchors, and the NMC hidden-gallery correction are explicit assumptions
in the examples. They must be reviewed when interpreting absolute lithiation.

`fitMSMR` uses a genetic algorithm followed by `fmincon`. Both examples set the
random seed to 0, and the fitter restores the caller's random stream. To adjust
the optimizer, set:

```matlab
settings.optimizer.seed = 42;
settings.optimizer.gaOptions = {'Display', 'off', 'MaxGenerations', 200};
settings.optimizer.fminconOptions = {'Display', 'off', 'MaxIterations', 1500};
```

The fit constrains parameter bounds, minimum gallery-potential separation,
the sum of site fractions, and a positive fitted interval. Check
`ocvrel.MSMR.(branch).fitInfo` and voltage residuals before interpreting a fit.
The diagnostics include both optimizers' exit flags and objectives, the final
constraint violation, and the seed. A feasible locally converged fit can be
accepted even when the genetic algorithm exhausts its generation budget.
Residual fields use model minus measurement. Voltage agreement alone does not
establish parameter uniqueness, hidden transitions, or absolute stoichiometry.

## Returned figures and results

`ocpFigures` contains seven figure handles in window order:

1. Measured half-cell voltage.
2. Relative processed branches and diagonal-average OCP.
3. Initial cutoff-based absolute coordinates.
4. Measured branches and separate MSMR fits on relative coordinates.
5. Branch fits and the mean-parameter curve on absolute coordinates.
6. Mean-parameter composition, occupations, and positive gallery derivatives.
7. Relative and initial cutoff-based measured OCP, with the MSMR mean curve
   after the configured hidden-gallery transformation.

`supplementaryResults` retains values and figure handles 1-6 before the
hidden-gallery and active-material corrections. The measured lookup table
returned in `model` is remapped later to the final interval; figure 7 does not
show that final remapping. The fit to the diagonal-average curve and the
arithmetic mean of branch parameters are distinct results.

## Public reference data and reproducibility

The files in `input/reference/` support inspection of the numerical references
without the experimental workbooks:

| File | Contents |
| --- | --- |
| `theoPracCapNeg.mat` | Scalar graphite reference capacity `Q_r` [mAh]. |
| `theoPracCapPos.mat` | Scalar NMC811 reference capacity `Q_r` [mAh]. |
| `ocpReferences.mat` | Structures `negative`, `positive`, and `provenance`. Each electrode has matching column vectors `thetaAbsolute` [-], `thetaRelative` [-], and `voltage_V` [V]. |

The capacity files are unchanged copies of the existing references. Their
calculation and assumed coating masses are documented in
`modelElectrodeCapacities.m`. The OCP arrays were copied without resampling from
existing average lookup tables, without machine-specific function handles.
The snapshot provenance identifies its source model files by basename and SHA-256:

| Source model | SHA-256 |
| --- | --- |
| `Uocp_n.mat` | `B0761AE0D302650F217BD03196B457563F3DD11AA763473B85515F5A735EB7E3` |
| `Uocp_p.mat` | `356B4F2F5894C05E9E8E0683BCF292E190B64BD10EB054FB6378FA2D8D70601A` |

These are existing model outputs, not newly fitted results. The companion
modelOCV article uses an older frozen snapshot with matching average potential
arrays but a positive absolute coordinate differing by up to approximately
0.0041244. The snapshots do not have identical absolute calibration.
The public references do not replace the private experimental inputs needed to
repeat the fits. Record the source revision, MATLAB version, fitting settings,
and data provenance when reporting results.

The development validation record dated November 20, 2026 reports 18 tests passed
under MATLAB R2024b Update 6 on Windows. Coverage included interpolation,
capacity extraction, MSMR derivatives and saturation, gallery ordering,
synthetic Arbin import/cache behavior, phase processing, and a small synthetic
fit. The synthetic fit checks the optimizer interface with fixed reference
bounds; it does not demonstrate general experimental parameter recovery. These
checks do not establish the physical validity of endpoint or hidden-gallery
assumptions, optimizer uniqueness, or reproducibility of a full experimental
fit. Development tests and their generated fixtures are kept locally and are
excluded from the public distribution. The publication cleanup did not rerun
these tests or refit the experimental models.

## Attribution and citation

Use `CITATION.cff` to cite this implementation and report the source revision.
The inherited BMS3 functions accompany Gregory L. Plett and M. Scott Trimboli,
*Battery Management Systems, Volume III: Physics-Based Methods*, Artech House,
2024. OCP processing follows the methodology of Lu et al. (2021),
DOI: 10.1149/1945-7111/ac11a4. Cite the relevant methodology sources along with
the software. The companion article's bibliographic metadata can be added when
available; no release DOI is assigned here.

Code contributions are licensed under CC BY-SA 4.0. The specified numerical
reference data are licensed under CC BY 4.0. See `LICENSE.md` for scope and
third-party attribution.
