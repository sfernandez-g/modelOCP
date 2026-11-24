# Repository preparation checks

Checked on 2026-11-20 using MATLAB R2024b Update 6 on Windows.

**18 tests passed, 0 failures, 0 incomplete**, using `runModelOCPTests`:

| Area | Tests | Scope |
| --- | ---: | --- |
| Synthetic Arbin import | 3 | Ah integration, separate cache/reuse, time-reset rejection, missing-source error |
| MSMR composition | 4 | Bounds/monotonicity, finite-difference derivative, saturated galleries, sum of gallery derivatives |
| Synthetic fitting | 1 | Known one-gallery case, inverse potential, feasibility diagnostics, random-stream restoration |
| Capacity extraction | 2 | Segment reference subtraction and empty selection |
| MSMR initialization | 1 | Asymmetric two-gallery peak/fraction ordering |
| Interpolation | 5 | Unsorted/duplicate samples, interpolation/extrapolation and invalid inputs |
| Phase processing | 2 | Mixed row/column inputs, 1 Ah reference, temperature propagation, branch/estimator pairing |

The importer tests create only tiny synthetic XLSX files under the ignored
`.build/` directory; they never use or overwrite experimental workbooks. The
synthetic fitting case checks the optimizer interface and conservation constraints,
not recovery of a general experimental parameter set. Its intentionally fixed
bounds make the reference answer known.

Saturation tests exercise potentials far beyond the physical operating window
only as numerical limiting cases. `compositionMSMR` returns finite occupations,
derivatives and supplied gallery potentials there. The positive derivative
convention and the legacy thermal-factor convention of `dZj` are retained.

The existing experimental model MAT outputs were not re-fitted or rewritten.
Their hashes are recorded in `input/reference/README.md`. The portable numerical
snapshot is copied without resampling; the positive calibration difference from
the modelOCV article snapshot is documented rather than silently removed.

Hidden-gallery and active-material corrections were retained as downstream
assumptions. These tests do not establish their physical validity, absolute
stoichiometry accuracy, optimizer uniqueness, or reproducibility of a full private
experimental fitting run. No runtime speedup factor is claimed.
