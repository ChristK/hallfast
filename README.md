# hallfast

<!-- badges: start -->
[![R-CMD-check](https://github.com/ChristK/hallfast/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChristK/hallfast/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Fast solver for the Hall adult body weight model.** `hallfast` converts
changes in energy intake (kcal/day) into body weight and BMI at the end of each
year, or of any other interval. It uses the dynamic adult model of
[Hall et al. (Lancet 2011)](https://doi.org/10.1016/S0140-6736(11)60812-X),
with the web appendix as the reference specification. The package is built for
large microsimulation models such as IMPACTncd.

| | Daily RK4 (`bw` package) | `hallfast` |
|---|---|---|
| Time per person-year (1 thread) | 230–340 µs | **≈1.1–1.6 µs** (≈0.4 µs for sustained changes); 140–320× faster |
| 1 million people × 10 years | not feasible (> 250 GB of daily matrices) | **≈11 s** on 1 thread; scales with threads |
| Maximum error against a high-precision solution (3,000 people, 10 years) | 0.025 kg | **< 0.0001 kg** |

## Installation

```r
# install.packages("remotes")
remotes::install_github("ChristK/hallfast", build_vignettes = TRUE)
```

If the repository is private, set a GitHub personal access token first, for
example with `Sys.setenv(GITHUB_PAT = "...")` or `gitcreds::gitcreds_set()`.

## Usage

```r
library(hallfast)

bw  <- c(70, 85, 100, 120); ht <- c(1.65, 1.75, 1.70, 1.80)
age <- c(30, 45, 60, 50);   sex <- c("female", "male", "female", "male")
dEI <- matrix(c(-50, -100, -100, -100, -100), 4, 5, byrow = TRUE)  # kcal/day vs baseline

BW   <- hall_bw_annual(bw, ht, age, sex, dEI)   # n x years: weight at end of each year
dBMI <- (BW - bw) / ht^2

# Long format (one row per person-year), e.g. in data.table:
# dt[, bw_end := hall_bw_long(pid, bw, ht, age, sex, dEI)]

# Appetite feedback (Hall 2024): EIchange is then the intervention "push"
hall_bw_annual(bw, ht, age, sex, dEI, appetite_k = "baseline")   # k = 95 kcal/day per kg

# Reproduce the 'bw' package (v1.0.0)
hall_bw_annual(bw, ht, age, sex, dEI, PAL = 1.5, model = "bw")
```

See `vignette("hallfast")` for:

- the model and numerical method;
- why the default PAL is 1.6;
- appetite feedback: the published k values, when to use them, and their
  limitations;
- the differences from `bw`;
- validation results.

## How it works

Within each interval the change in intake is constant. This allows four
simplifications, none of which approximate the model:

1. **Exact closed forms.** Adaptive thermogenesis, extracellular fluid and
   glycogen each have one (exponential relaxation, and a Riccati equation for
   glycogen).
2. **One equation for fat and lean mass.** Forbes partitioning reduces them to
   a single scalar ODE.
3. **Exact integrals for the fast compartments.** The 1–2-day transients enter
   as exact time-integrals, and an exact change of variables removes most of
   the 14-day adaptive-thermogenesis transient.
4. **Adaptive integration.** The remaining equation is solved with an adaptive
   Dormand–Prince 5(4) method, at about 10–30 evaluations per person-year.

`hall_bw_reference()` provides an independent fixed-step solution of the full
four-compartment system, so you can check the accuracy yourself.

## Features

- **Two interfaces:** a matrix interface (`hall_bw_annual()`) and a
  long-format interface (`hall_bw_long()`) suited to data.table pipelines.
- **The full set of model inputs:** sodium and carbohydrate changes, and known
  baseline energy intake or fat mass.
- **Optional proportional appetite feedback** (Hall 2024). `appetite_k` takes
  a number, a preset from `hall_appetite_presets`, a value per person, or a
  person × year matrix.
- **`model = "bw"`** reproduces `bw::adult_weight()`.
- **OpenMP multi-threading** through `threads`, or
  `options(hallfast.threads = )`.
- **Standalone C++ header** (`src/hall_model.h`) with no R dependency.

## Validation

The scripts in [`validation/`](validation) reproduce all accuracy and speed
results; they need the `bw` and `data.table` packages. The unit tests use
[tinytest](https://cran.r-project.org/package=tinytest):

```r
tinytest::test_package("hallfast")
```

## Citation

```r
citation("hallfast")
```

Please also cite Hall et al. (2011), and Hall (2024) if you use appetite
feedback.

## Licence

GPL-3. See [LICENSE.md](LICENSE.md).

## Acknowledgements

The R package [`bw`](https://CRAN.R-project.org/package=bw) by Dalia
Camacho-García-Formentí and Rodrigo Zepeda-Tello (MIT licence) implements the
same model with daily time-stepping. It was used as a comparator during
development and validation.
