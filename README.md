# modelOCP

MATLAB scripts and functions for processing Arbin battery test data, estimating OCP/OCV curves, and fitting multi-site multi-reaction (MSMR) electrode models.

## Requirements

- MATLAB R2023b or newer (the scripts use `orderedcolors` and recent name/value syntax).
- Optimization Toolbox for `fmincon`; Global Optimization Toolbox for `ga`; Signal Processing Toolbox for `xcorr`.
- The Arbin XLSX files under `input/` and capacity MAT files under `output/` for the workflows that use them.

## Run a workflow

Run a top-level script from the MATLAB editor or with `run`, for example:

```matlab
run('modelOCP_n_graphite.m')
```

Scripts resolve input, output, and `fun/` paths relative to their own location, so MATLAB's current folder does not need to be the repository root. The workflows are interactive: they create figures and save model MAT files in `output/`.

The main entry scripts are:

- `modelElectrodeCapacities.m`: calculate practical electrode capacities.
- `modelOCP_n_graphite.m` and `modelOCP_p_NMC.m`: process half-cell data and fit MSMR models.
- `modelOCV_commercial.m` and `modelOCV_NMC_graphite.m`: process full-cell OCV datasets.

Reusable processing is implemented in `fun/`, including `getDataArbin_XLSX`, `processCU`, `estimateOCV`, and `fitMSMR`.

`getDataArbin_XLSX` caches processed input in a MAT file beside the source data (or in the optional save directory). It checks source file size/timestamp and the counter-reconstruction setting before reusing a cache. Force a fresh read with:

```matlab
test_data = getDataArbin_XLSX(file, 7, [], true);
```

## MSMR fitting settings

`fitMSMR` returns a ninth output, `fitInfo`, with the random seed, objective values, exit flags, and solver output structures. Existing calls requesting eight outputs remain valid. Set `config.optimizer.seed` to control GA initialization and pass solver name/value pairs through `config.optimizer.gaOptions` or `config.optimizer.fminconOptions`, for example:

```matlab
config.optimizer.seed = 42;
config.optimizer.gaOptions = {'Display', 'off', 'MaxGenerations', 200};
config.optimizer.fminconOptions = {'Display', 'off', 'MaxIterations', 1500};
```

The random number generator state is restored when fitting returns.

## Tests

Run the lightweight function tests from the repository root:

```matlab
results = runtests('tests');
assertSuccess(results);
```

These tests cover interpolation input handling and extrapolation, MSMR composition behavior, and capacity segment extraction. They do not require the experimental XLSX data or the optimization toolboxes.
