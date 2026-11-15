# Validation

Local test suite: 10 tests passed across the two repositories (3 original-study,
7 ML). The source audit and descriptive reanalysis completed. All three ML feature
sets completed with eight stability repeats and optional t-SNE. Independent
nearest-center, inertia, and PCA projection checks passed on the real data.

The local Accelerate/NumPy runtime emits matrix-operation warnings despite finite,
consistent results. These are recorded as a reproducibility concern, not suppressed
or dismissed. Linux CI is provided to test a clean, separate environment and save
fresh outputs as an Actions artifact. Inspect its actual status before claiming a
cross-platform reproduction. The R/GSVA rerun has not been executed.

## Independent Linux reproduction

GitHub Actions run https://github.com/sahilnpatel/nsclc-circrna-analysis/actions/runs/43683057794 succeeded on November 15, 2026. Fresh dependency installation, all three tests, the source audit, and the descriptive reanalysis completed successfully. The run saved a reproduced-results artifact. This does not include the unexecuted R/GSVA or historical web-based analyses.
