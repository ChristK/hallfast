# Validation

These scripts reproduce the accuracy and speed results quoted in the README
and the vignette. They are kept in the repository but are not part of the
package build (see `.Rbuildignore`).

| Script | What it does | Run time |
|---|---|---|
| `validate_vs_bw.R` | Compares `hallfast` with the independent reference solver `hall_bw_reference()` (full four-compartment system, RK4, dt = 0.05 day), which is the gold standard, and with the `bw` package (dt = 1 and 0.1 day). Covers the realistic and stress populations, all initialisation paths, appetite feedback, the size of `bw`'s departures from the appendix, tolerance sensitivity and benchmarks. | ≈15 min |
| `check_bw_convergence.R` | Shows that `bw` converges to `hallfast(model = "bw")` as its time step shrinks, on the cases where the two differ most. | ≈2 min |

Requirements: the installed `hallfast` package, plus `bw` and `data.table`
from CRAN. Run from the repository root:

```sh
Rscript validation/validate_vs_bw.R
Rscript validation/check_bw_convergence.R
```

Results are written to `validation/results/`. The CSV files are also copied to
`inst/extdata/validation/`, where the vignette reads them.
`results/session_info.txt` records the environment used.
