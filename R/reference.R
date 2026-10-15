# Reference (slow, high-precision) solver for validation.
# Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos. GPL-3.

#' Reference solver: full model integrated with a small fixed step
#'
#' Integrates the full four-compartment system of the Hall et al. (2011)
#' appendix (adaptive thermogenesis, extracellular fluid, glycogen and lean
#' mass, with fat mass from the Forbes relation) jointly, using classical
#' fourth-order Runge-Kutta with `steps_per_day` steps per day. It uses no
#' closed-form solutions, no quasi-steady-state reductions and no change of
#' variables. It is therefore an independent check of [hall_bw_annual()],
#' which should agree with it to within about 1e-4 kg. With
#' `appetite_k > 0`, the realised intake is computed from the full state at
#' every stage.
#'
#' This function is slow (about 30,000 derivative evaluations per
#' person-year at the default step) and is intended for validation, not for
#' production use.
#'
#' @inheritParams hall_bw_annual
#' @param appetite_k Appetite feedback gain(s), as in [hall_bw_annual()]: a
#'   scalar, a length-\eqn{n} vector, an \eqn{n \times Y}{n x Y} matrix or a
#'   preset name.
#' @param steps_per_day Number of RK4 steps per day. Default 20
#'   (\eqn{dt = 0.05} day).
#' @param output `"bw"` (default) returns body weight; `"all"` returns a list
#'   with `BW` and `tissue` (fat plus lean mass).
#' @return An \eqn{n \times Y}{n x Y} matrix of body weight (kg) at the end of
#'   each interval, or a list (see `output`).
#' @seealso [hall_bw_annual()]
#' @examples
#' dEI <- matrix(c(-100, -200, 50, 0), 2, 2)
#' fast <- hall_bw_annual(c(80, 95), c(1.75, 1.65), c(40, 60), c("male", "female"), dEI)
#' ref  <- hall_bw_reference(c(80, 95), c(1.75, 1.65), c(40, 60), c("male", "female"), dEI,
#'                           steps_per_day = 4)
#' max(abs(fast - ref))   # well below 0.001 kg
#' @export
hall_bw_reference <- function(bw, ht, age, sex, EIchange, NAchange = 0,
                              EI = NA_real_, fat = NA_real_, PAL = 1.6,
                              pcarb_base = 0.5, pcarb = pcarb_base, days = 365,
                              appetite_k = 0, model = c("hall2011", "bw"),
                              steps_per_day = 20L, output = c("bw", "all")) {
  model <- match.arg(model); output <- match.arg(output)
  if (!is.numeric(steps_per_day) || length(steps_per_day) != 1L || steps_per_day < 1)
    stop("'steps_per_day' must be a single integer >= 1.", call. = FALSE)
  .hall_check_controls(days, 1, 1L)
  n <- length(bw)
  if (is.null(dim(EIchange))) EIchange <- matrix(EIchange, nrow = n)
  Y <- ncol(EIchange)
  if (nrow(EIchange) != n) stop("'EIchange' must have one row per individual.", call. = FALSE)
  if (length(ht) != n || length(age) != n || length(sex) != n)
    stop("'bw', 'ht', 'age' and 'sex' must have the same length.", call. = FALSE)
  female     <- .hall_code_sex(sex)
  PAL        <- .hall_rep_n(PAL, n, "PAL")
  EI         <- .hall_rep_n(EI, n, "EI")
  fat        <- .hall_rep_n(fat, n, "fat")
  pcarb_base <- .hall_rep_n(pcarb_base, n, "pcarb_base")
  dEI <- .hall_as_matrix(EIchange, n, Y, "EIchange")
  dNa <- .hall_as_matrix(NAchange, n, Y, "NAchange")
  pc  <- .hall_as_matrix(pcarb, n, Y, "pcarb")
  km  <- .hall_resolve_k(appetite_k)
  km  <- if (length(km) == n * Y) matrix(km, n, Y) else .hall_as_matrix(km, n, Y, "appetite_k")
  EI0 <- ifelse(is.na(EI), .hall_ei0(bw, ht, age, female, PAL), EI)
  .hall_check_inputs(bw, ht, age, PAL, pcarb_base, pc, EI0, dEI, km)
  res <- .hall_reference_core(as.numeric(bw), as.numeric(ht), as.numeric(age), female, PAL,
                              EI, fat, pcarb_base, dEI, dNa, pc, km, days,
                              as.integer(steps_per_day), model == "bw")
  dn <- list(names(bw), paste0("year", seq_len(Y)))
  dimnames(res$BW) <- dn; dimnames(res$tissue) <- dn
  if (output == "bw") res$BW else res
}
