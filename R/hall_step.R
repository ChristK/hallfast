# One interval at a time: a carried state for microsimulations that advance
# their population year by year (hallfast 0.3.0).
# Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos. GPL-3.

# Columns of a hall_state, in order. Baseline inputs (fixed), the derived
# baseline intake and BMI, then the dynamic state carried between intervals.
.HALL_STATE_COLS <- c("bw0", "ht", "age0", "female", "PAL", "EI", "fat", "pcarb_base",
                      "EI0", "bmi0", "height_imputed", "bw_compat",
                      "L", "AT", "ECF", "G", "dEI", "h", "BW", "BMI", "intervals")

#' A carried model state for stepping one interval at a time
#'
#' Creates the baseline state of the Hall et al. (2011) adult model for
#' \eqn{n} individuals, to be advanced one interval at a time with
#' [hall_step()]. This is the interface for microsimulations that simulate
#' their population year by year and cannot hand the whole history of intake
#' changes to [hall_bw_annual()] at once: keep the state between years (one
#' row per individual), and step it with each year's change in intake.
#'
#' Give either `bw` (body weight, kg) with `ht`, or `bmi`. With `bmi`, height
#' is optional: missing heights are replaced by the sex-specific `ref_height`,
#' as in [hall_bmi_annual()], and such individuals can then only be stepped
#' with changes in intake given as a fraction of baseline intake
#' (`EIchange_unit = "fraction"` in [hall_step()]).
#'
#' @param age Baseline age (years); length \eqn{n}. As in the model, age
#'   enters only the baseline equations.
#' @param sex Sex; length \eqn{n}, coded as in [hall_bw_annual()].
#' @param bw Baseline body weight (kg), length \eqn{n}; or `NULL` when `bmi`
#'   is given.
#' @param bmi Baseline BMI (kg/m\eqn{^2}), length \eqn{n}; or `NULL` when `bw`
#'   is given.
#' @param ht Height (m): required with `bw`; with `bmi`, `NULL` (default), a
#'   scalar or a length-\eqn{n} vector that may contain `NA`.
#' @param PAL,EI,fat,pcarb_base Baseline physical activity level, energy
#'   intake (kcal/day; `NA` = Mifflin-St Jeor RMR \eqn{\times}{x} `PAL`), fat
#'   mass (kg; `NA` = Jackson et al.) and carbohydrate share, as in
#'   [hall_bw_annual()]. Scalars or length-\eqn{n} vectors. `EI` and `fat`
#'   cannot be given for individuals whose height is imputed.
#' @param ref_height Reference heights used for missing heights; see
#'   [hall_ref_height].
#' @param model `"hall2011"` (default) or `"bw"`, as in [hall_bw_annual()].
#'
#' @return A `hall_state`: a data frame with one row per individual and the
#'   columns
#'   * `bw0`, `ht`, `age0`, `female`, `PAL`, `EI`, `fat`, `pcarb_base`: the
#'     baseline inputs (`EI` and `fat` as given, possibly `NA`);
#'   * `EI0`: baseline energy intake (kcal/day), the reference for changes
#'     given as fractions;
#'   * `bmi0`: baseline BMI; `height_imputed`: whether `ht` is a reference
#'     height; `bw_compat`: `model == "bw"`;
#'   * `L`, `AT`, `ECF`, `G`, `h`: the dynamic state carried between
#'     intervals (lean mass, kg; adaptive thermogenesis, kcal/day;
#'     extracellular fluid, kg; glycogen, kg; the integrator's step-size
#'     hint, days, `NA` before the first interval);
#'   * `dEI`: the realised change in intake at the end of the last interval
#'     (kcal/day; with appetite feedback it differs from the push), an
#'     output;
#'   * `BW`, `BMI`: body weight (kg) and BMI at the end of the last interval
#'     (at baseline: `bw0` and `bmi0`); `intervals`: the number of intervals
#'     advanced.
#'
#'   Rows can be subset, reordered or bound together (for example to add
#'   individuals who enter a simulation later): each row is independent.
#'   Individuals with a missing baseline input have an all-`NA` state.
#'
#' @seealso [hall_step()] to advance the state; [hall_bw_annual()] for the
#'   model and its arguments.
#' @examples
#' s <- hall_state(age = c(30, 55), sex = c("female", "male"),
#'                 bw = c(70, 95), ht = c(1.65, 1.80))
#' s <- hall_step(s, c(-100, -200))   # year 1, kcal/day below baseline
#' s <- hall_step(s, c(-100, -200))   # year 2
#' s[, c("BW", "BMI", "intervals")]
#'
#' # The same as solving both years at once
#' hall_bw_annual(c(70, 95), c(1.65, 1.80), c(30, 55), c("female", "male"),
#'                matrix(c(-100, -200), 2, 2))[, 2]
#' @export
hall_state <- function(age, sex, bw = NULL, bmi = NULL, ht = NULL, PAL = 1.6,
                       EI = NA_real_, fat = NA_real_, pcarb_base = 0.5,
                       ref_height = hall_ref_height, model = c("hall2011", "bw")) {
  model <- match.arg(model)
  if (is.null(bw) == is.null(bmi))
    stop("Give exactly one of 'bw' (with 'ht') and 'bmi'.", call. = FALSE)
  n <- length(if (is.null(bw)) bmi else bw)
  if (n == 0L) stop("'bw' or 'bmi' has length 0.", call. = FALSE)
  if (length(age) != n || length(sex) != n)
    stop("'age', 'sex' and 'bw' (or 'bmi') must have the same length.", call. = FALSE)
  female     <- .hall_code_sex(sex)
  PAL        <- .hall_rep_n(PAL, n, "PAL")
  EI         <- .hall_rep_n(EI, n, "EI")
  fat        <- .hall_rep_n(fat, n, "fat")
  pcarb_base <- .hall_rep_n(pcarb_base, n, "pcarb_base")
  age        <- as.numeric(age)

  if (!is.null(bw)) {
    if (is.null(ht)) stop("'ht' is required with 'bw'.", call. = FALSE)
    ht  <- .hall_rep_n(ht, n, "ht")
    bw0 <- as.numeric(bw)
    if (any(bw0 <= 0, na.rm = TRUE)) stop("'bw' must be positive (kg).", call. = FALSE)
    bmi0 <- bw0 / ht^2
    imputed <- rep(FALSE, n)
  } else {
    bmi0 <- as.numeric(bmi)
    if (any(bmi0 <= 0, na.rm = TRUE)) stop("'bmi' must be positive (kg/m^2).", call. = FALSE)
    # Heights: imputation is allowed here; hall_step() refuses kcal changes for
    # imputed rows. EI and fat are absolute and need a real height.
    h <- .hall_bmi_height(ht, female, n, "fraction", ref_height,
                          absolute = !is.na(EI) | !is.na(fat))
    ht <- h$ht; imputed <- h$imputed
    bw0 <- bmi0 * ht^2
  }
  EI0 <- ifelse(is.na(EI), .hall_ei0(bw0, ht, age, female, PAL), EI)
  .hall_check_inputs(bw0, ht, age, PAL, pcarb_base, 0.5, EI0, 0)

  ini <- .hall_init_core(bw0, ht, age, female, PAL, EI, fat, pcarb_base,
                         bw_compat = model == "bw")
  ok <- !is.na(ini$L)
  out <- data.frame(bw0 = bw0, ht = ht, age0 = age, female = female, PAL = PAL,
                    EI = EI, fat = fat, pcarb_base = pcarb_base, EI0 = EI0,
                    bmi0 = bmi0, height_imputed = imputed, bw_compat = model == "bw",
                    L = ini$L, AT = ifelse(ok, 0, NA_real_), ECF = ini$ECF, G = ini$G,
                    dEI = ifelse(ok, 0, NA_real_), h = NA_real_,
                    BW = ifelse(ok, bw0, NA_real_), BMI = ifelse(ok, bmi0, NA_real_),
                    intervals = 0L)
  class(out) <- c("hall_state", "data.frame")
  out
}


#' Advance a carried model state by one interval
#'
#' Advances every individual of a [hall_state()] by one interval (a year by
#' default) with a constant change in energy intake from **baseline**, and
#' returns the new state. Calling it \eqn{Y} times with the columns of an
#' \eqn{n \times Y}{n x Y} `EIchange` matrix gives exactly (bit for bit) the
#' body weights of one call of [hall_bw_annual()] with that matrix, so a
#' microsimulation can carry the state from year to year instead of
#' re-solving each individual's history.
#'
#' @param state A `hall_state` from [hall_state()] or a previous `hall_step()`.
#' @param EIchange Change in energy intake from baseline during this interval:
#'   kcal/day (`EIchange_unit = "kcal"`, the default) or a fraction of each
#'   individual's baseline intake `EI0`. A scalar or a length-\eqn{n} vector.
#'   With `appetite_k > 0` it is the intervention push, as in
#'   [hall_bw_annual()].
#' @param EIchange_unit `"kcal"` (default) or `"fraction"`. Individuals whose
#'   height was imputed by [hall_state()] can only take fractions.
#' @param days Length of the interval in days (default 365). Intervals of
#'   different lengths can follow each other.
#' @param NAchange,pcarb,appetite_k Change in sodium (mg/day), carbohydrate
#'   share and appetite feedback gain during this interval, as in
#'   [hall_bw_annual()]: scalars or length-\eqn{n} vectors (`appetite_k` may
#'   be a preset name). `pcarb` defaults to each individual's `pcarb_base`.
#' @param atol,threads As in [hall_bw_annual()].
#'
#' @return The advanced `hall_state` (same rows, same order): the dynamic
#'   state, `BW` (body weight at the end of the interval, kg), `BMI`
#'   (`bmi0 * BW / bw0`, as in [hall_bmi_annual()]) and `intervals` updated.
#'   An individual with a missing input in this interval, or whose state is
#'   already missing, gets a missing state, which stays missing in later
#'   intervals (as [hall_bw_annual()] stops an individual at the first
#'   missing input).
#'
#' @details A state that has not moved from baseline and is stepped with no
#'   change in intake, sodium or carbohydrate share stays exactly at baseline
#'   (the solver's exact shortcut). Its `BW` is then the model's own baseline
#'   weight, which can differ from `bw0` in the last binary digit; to keep an
#'   untouched individual exactly at `bw0`, do not step them.
#'
#' @seealso [hall_state()]; [hall_bw_annual()] for the model.
#' @examples
#' s <- hall_state(age = c(40, 62), sex = c("male", "female"), bmi = c(29, 33))
#' for (y in 1:3) s <- hall_step(s, -0.05, EIchange_unit = "fraction")
#' s$BMI                                            # after three years of a 5% cut
#' hall_bmi_annual(c(29, 33), c(40, 62), c("male", "female"),
#'                 matrix(-0.05, 2, 3))[, 3]        # the same
#' @export
hall_step <- function(state, EIchange, EIchange_unit = c("kcal", "fraction"),
                      days = 365, NAchange = 0, pcarb = NULL, appetite_k = 0,
                      atol = 1e-4, threads = getOption("hallfast.threads", 1L)) {
  EIchange_unit <- match.arg(EIchange_unit)
  if (!inherits(state, "hall_state") || !all(.HALL_STATE_COLS %in% names(state)))
    stop("'state' must be a hall_state from hall_state() or hall_step().", call. = FALSE)
  .hall_check_controls(days, atol, threads)
  n <- nrow(state)
  if (n == 0L) return(state)
  compat <- unique(state$bw_compat)
  if (length(compat) != 1L)
    stop("All rows of 'state' must use the same 'model'.", call. = FALSE)

  dEI <- .hall_rep_n(EIchange, n, "EIchange")
  if (EIchange_unit == "kcal") {
    if (any(state$height_imputed & !is.na(dEI)))
      stop("Height is required when 'EIchange' is in kcal/day: at a given BMI the ",
           "change in BMI from a fixed kcal change scales roughly with 1/ht^2. ",
           "Some rows of 'state' have an imputed height; give 'EIchange' as a ",
           "fraction of baseline intake (EIchange_unit = \"fraction\").", call. = FALSE)
  } else {
    dEI <- .hall_fraction_to_kcal(dEI, state$EI0)
  }
  dNa <- .hall_rep_n(NAchange, n, "NAchange")
  pc  <- if (is.null(pcarb)) as.numeric(state$pcarb_base) else .hall_rep_n(pcarb, n, "pcarb")
  kk  <- .hall_rep_n(.hall_resolve_k(appetite_k), n, "appetite_k")
  .hall_check_inputs(state$bw0, state$ht, state$age0, state$PAL, state$pcarb_base, pc,
                     state$EI0, dEI, kk)

  res <- .hall_step_core(state$bw0, state$ht, state$age0, as.integer(state$female),
                         state$PAL, state$EI, state$fat, state$pcarb_base,
                         state$L, state$AT, state$ECF, state$G, state$dEI, state$h,
                         dEI, dNa, pc, kk, days = days, atol = atol,
                         bw_compat = isTRUE(compat), nthreads = as.integer(threads))
  state$L <- res$L; state$AT <- res$AT; state$ECF <- res$ECF; state$G <- res$G
  state$dEI <- res$dEI; state$h <- res$h; state$BW <- res$BW
  state$BMI <- state$bmi0 * res$BW / state$bw0
  state$intervals <- state$intervals + 1L
  state
}
