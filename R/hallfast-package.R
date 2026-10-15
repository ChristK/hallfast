# Package-level documentation and native routine registration.
# Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos. GPL-3.

#' hallfast: Fast Solver for the Hall Adult Body Weight Model
#'
#' Converts changes in energy intake into body weight and BMI at the end of
#' each year (or any other interval) using the dynamic adult body-weight model
#' of Hall et al. (2011). The model's web appendix is the reference
#' specification. The solver is several hundred times faster than daily
#' time-stepping, with negligible numerical error, and is designed for large
#' microsimulation models such as IMPACTncd.
#'
#' Main functions:
#' \itemize{
#'   \item [hall_bw_annual()]: individuals by intervals (matrix interface).
#'   \item [hall_bw_long()]: long-format interface (one row per person-year).
#'   \item [hall_bw_reference()]: slow high-precision reference solver for
#'     validation.
#'   \item [hall_appetite_presets]: published appetite feedback gains
#'     (Hall 2024).
#' }
#'
#' Options:
#' \itemize{
#'   \item `hallfast.threads`: default number of OpenMP threads (1).
#' }
#'
#' See `vignette("hallfast")` for the model, the numerical method, guidance
#' on physical activity level and appetite feedback, the differences from the
#' 'bw' package, and validation results.
#'
#' @references
#' Hall KD, Sacks G, Chandramohan D, Chow CC, Wang YC, Gortmaker SL, Swinburn
#' BA (2011). Quantification of the effect of energy imbalance on bodyweight.
#' *Lancet* 378(9793):826-837. \doi{10.1016/S0140-6736(11)60812-X}
#'
#' Hall KD (2024). Physiology of the weight-loss plateau in response to diet
#' restriction, GLP-1 receptor agonism, and bariatric surgery. *Obesity*
#' 32(6):1163-1168. \doi{10.1002/oby.24027}
#' @keywords internal
"_PACKAGE"

#' @useDynLib hallfast, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Was the package compiled with OpenMP?
#'
#' @return `TRUE` if the compiled code supports multi-threading through the
#'   `threads` argument of [hall_bw_annual()] and [hall_bw_long()], otherwise
#'   `FALSE`. On macOS, R's default toolchain usually lacks OpenMP.
#' @examples
#' hall_has_openmp()
#' @export
hall_has_openmp <- function() .hall_has_openmp()
