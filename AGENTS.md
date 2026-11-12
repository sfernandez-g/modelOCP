# Repository Guidelines

## Project Structure & Module Organization

This repository contains MATLAB workflows for battery OCP/OCV analysis. Top-level scripts (`modelOCP_*.m`, `modelOCV_*.m`, and `modelElectrodeCapacities.m`) orchestrate processing, plots, fitting, and saved models. Reusable functions live in `fun/`; MATLAB unit tests live in `tests/`. Experimental Arbin workbooks are under `input/`, and generated MAT files are stored in `output/`. Resolve repository files relative to the script using `fileparts(mfilename('fullpath'))` and `fullfile`.

## Build, Test, and Development Commands

There is no separate build step. From the MATLAB Command Window, run a workflow with `run('modelOCP_n_graphite.m')` (or another top-level script). Run the test suite from the repository root with:

```matlab
results = runtests('tests');
assertSuccess(results);
```

The workflows require MATLAB R2023b or newer and the toolboxes documented in `README.md`. They create figures and write model MAT files under `output/`.

## Coding Style & Naming Conventions

Use MATLAB `.m` files, four-space indentation in functions, descriptive camelCase function and variable names, and uppercase names for established scientific quantities only where already used. Keep reusable calculations in `fun/` as functions with explicit inputs and outputs; leave workflow-specific plotting and configuration in top-level scripts. Use `fullfile` for paths and include units in comments for physical quantities. Avoid changing experimental source files as part of routine code edits.

## Testing Guidelines

Tests use MATLAB's function-based `matlab.unittest` style. Name test files `test*.m` and test functions with a `test` prefix. Add focused tests for numerical edge cases and helper behavior; tests should not require the large experimental workbooks or Optimization Toolbox unless essential to the behavior under test. Run `runtests('tests')` before submitting changes.

## Commit & Pull Request Guidelines

No established commit convention is visible in the repository. Use short imperative commit subjects, for example `Validate Arbin capacity segments`. Pull requests should explain the scientific or workflow impact, list relevant MATLAB checks, identify any changed input/output assumptions, and include plot images only when visual output changes are central to the review.

## Data and Cache Notes

`getDataArbin_XLSX` caches processed data as MAT files and checks source workbook metadata and counter settings before reuse. Use its `forceReprocess` argument when intentionally rebuilding a cache after changing processing behavior.
