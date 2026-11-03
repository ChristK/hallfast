# BMI interfaces: run the model from BMI, with height optional when the change
# in intake is proportional to baseline intake.
# Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos. GPL-3.

#' Reference adult heights used when height is unknown
#'
#' Sex-specific mean adult heights (m), used by [hall_bmi_annual()] and
#' [hall_bmi_long()] when `ht` is `NULL` or `NA`. The values are the measured
#' means for adults aged 16 and over in the Health Survey for England 2019:
#' 1.755 m for men and 1.617 m for women.
#'
#' The reference height is used only when the change in intake is a fraction
#' of baseline intake (`EIchange_unit = "fraction"`). In that case the
#' proportional change in BMI depends very little on height, so the choice of
#' reference matters little; see the vignette section "BMI and proportional
#' changes in intake". Pass your own values through the `ref_height`
#' argument, for example the mean height of the population being simulated.
#'
#' @format A named numeric vector with elements `female` and `male` (m).
#' @source NHS England (formerly NHS Digital). *Health Survey for England
#'   2019: Adult and child overweight and obesity* (and *Health Survey for
#'   England 2021, Part 4: Trends*), mean measured height of adults.
#' @examples
#' hall_ref_height
#' @export
hall_ref_height <- c(female = 1.617, male = 1.755)


#' BMI under the Hall adult model, by year (wide interface)
#'
#' BMI counterpart of [hall_bw_annual()]. Takes baseline BMI instead of body
#' weight and returns BMI at the end of each interval. Height is optional
#' when the change in intake is given as a fraction of baseline intake
#' (`EIchange_unit = "fraction"`, the default here).
#'
#' @param bmi Baseline BMI (kg/m\eqn{^2}); numeric vector of length \eqn{n}.
#' @param age,sex As in [hall_bw_annual()].
#' @param EIchange Change in energy intake from baseline: a fraction of
#'   baseline energy intake (default, for example `-0.05` for a sustained 5%
#'   cut) or kcal/day with `EIchange_unit = "kcal"`. An
#'   \eqn{n \times Y}{n x Y} matrix or a length-\eqn{n} vector, as in
#'   [hall_bw_annual()].
#' @param ht Height (m): `NULL` (default), a scalar or a length-\eqn{n}
#'   vector, which may contain `NA`. Missing heights are replaced by the
#'   sex-specific `ref_height`; this is allowed only with
#'   `EIchange_unit = "fraction"`.
#' @param EIchange_unit `"fraction"` (default) or `"kcal"`; see `EIchange`.
#'   Note that the default differs from [hall_bw_annual()], whose default is
#'   `"kcal"`.
#' @param EI,fat Optional baseline energy intake (kcal/day) and fat mass (kg),
#'   as in [hall_bw_annual()]. Being absolute quantities, they can be given
#'   only for individuals whose height is supplied.
#' @param ref_height Named numeric vector with elements `female` and `male`:
#'   heights (m) used for individuals without a height. Defaults to
#'   [hall_ref_height].
#' @param output `"bmi"` (default) returns the matrix of BMI. `"all"` returns
#'   the full list from [hall_bw_annual()] (body weight and compartments in kg
#'   for the height used), with `BMI` computed as below, plus `height` (the
#'   height used, m) and `height_imputed` (logical).
#' @param ... Further arguments passed to [hall_bw_annual()], such as `PAL`,
#'   `NAchange`, `pcarb_base`, `pcarb`, `days`, `appetite_k`, `atol`, `model`
#'   and `threads`.
#'
#' @details
#' **Calculation.** Baseline weight is set to `bmi * ht^2` and the model is
#' solved with [hall_bw_annual()]. Because height is fixed in the model, BMI
#' at the end of interval \eqn{y} is
#' \deqn{BMI_y = BMI_0 \times BW_y / BW_0,}{BMI_y = BMI_0 * BW_y / BW_0,}
#' which equals \eqn{BW_y / ht^2} exactly when the true height is supplied.
#'
#' **Why height can be omitted for proportional changes.** At a given BMI,
#' the change in body weight in kg produced by a fixed change in kcal/day
#' hardly depends on height, so the change in BMI scales roughly with
#' \eqn{1/ht^2}: in a woman with BMI 27, a sustained cut of 100 kcal/day
#' lowers BMI after 30 years by 1.92 at 1.50 m but by 1.54 at 1.74 m. A fixed
#' kcal change is a larger share of a smaller person's energy requirement.
#' When the change is instead a fixed *fraction* of baseline intake, which
#' scales with body size, the proportional change in BMI is nearly the same
#' at any height (6.6% to 6.9% for a 5% cut in the same example). The
#' sex-specific reference height then gives almost the same answer as the
#' true height. The vignette quantifies this in a realistic population.
#' Height is therefore required with `EIchange_unit = "kcal"`.
#'
#' **Caveats.** The near-independence from height is an empirical property
#' of the model, not an identity; the vignette quantifies the residual error
#' (in a realistic population, a relative SD of about 1\% in the change in
#' BMI, with no material bias in the mean). Other absolute inputs keep their
#' units: `NAchange` (mg/day) and the appetite feedback gain `appetite_k`
#' (kcal/day per kg of weight change). Feedback did not increase the
#' dependence on height in validation, but it is not scaled to body size.
#' The body-weight and compartment outputs of `output = "all"` refer to the
#' height used, not to the (unknown) true height, when height is imputed.
#'
#' @return With `output = "bmi"`, an \eqn{n \times Y}{n x Y} matrix of BMI
#'   (kg/m\eqn{^2}) at the end of each interval. With `output = "all"`, a list
#'   as described for `output`.
#' @seealso [hall_bw_annual()] for the model, method and other arguments;
#'   [hall_bmi_long()] for long-format data; [hall_ref_height].
#' @examples
#' bmi <- c(24, 31, 36, 28); age <- c(30, 45, 60, 50)
#' sex <- c("female", "male", "female", "male")
#' # Sustained 5% cut in intake; no heights needed
#' B <- hall_bmi_annual(bmi, age, sex, matrix(-0.05, 4, 5))
#' round(B - bmi, 2)                     # change in BMI
#' round(100 * (B / bmi - 1), 2)         # proportional change in BMI (%)
#'
#' # With known heights the result is exact (equal to hall_bw_annual() / ht^2)
#' ht <- c(1.65, 1.75, 1.70, 1.80)
#' hall_bmi_annual(bmi, age, sex, -0.05, ht = ht)
#' hall_bw_annual(bmi * ht^2, ht, age, sex, -0.05, EIchange_unit = "fraction") / ht^2
#'
#' # A change in kcal/day needs height
#' hall_bmi_annual(bmi, age, sex, -100, ht = ht, EIchange_unit = "kcal")
#' @export
hall_bmi_annual <- function(bmi, age, sex, EIchange, ht = NULL,
                            EIchange_unit = c("fraction", "kcal"),
                            EI = NA_real_, fat = NA_real_,
                            ref_height = hall_ref_height,
                            output = c("bmi", "all"), ...) {
  EIchange_unit <- match.arg(EIchange_unit); output <- match.arg(output)
  n <- length(bmi)
  if (n == 0L) stop("'bmi' has length 0.", call. = FALSE)
  .hall_check_dots(list(...))
  if (length(sex) != n)
    stop("'bmi', 'age' and 'sex' must have the same length.", call. = FALSE)
  if (any(bmi <= 0, na.rm = TRUE)) stop("'bmi' must be positive (kg/m^2).", call. = FALSE)
  EI <- .hall_rep_n(EI, n, "EI"); fat <- .hall_rep_n(fat, n, "fat")
  h <- .hall_bmi_height(ht, .hall_code_sex(sex), n, EIchange_unit, ref_height,
                        absolute = !is.na(EI) | !is.na(fat))
  bw0 <- as.numeric(bmi) * h$ht^2

  res <- hall_bw_annual(bw0, h$ht, age, sex, EIchange, EI = EI, fat = fat,
                        EIchange_unit = EIchange_unit,
                        output = if (output == "bmi") "bw" else "all", ...)
  BW  <- if (output == "bmi") res else res$BW
  BMI <- as.numeric(bmi) * BW / bw0       # BMI_0 * BW_y / BW_0 (column-wise recycling)
  dimnames(BMI) <- list(names(bmi), colnames(BW))
  if (output == "bmi") return(BMI)
  res$BMI <- BMI
  res$baseline$BMI <- as.numeric(bmi)
  res$height <- h$ht
  res$height_imputed <- h$imputed
  res
}


#' BMI under the Hall adult model, by year (long interface)
#'
#' Long-format version of [hall_bmi_annual()], for microsimulation tables with
#' one row per person-year. As in [hall_bw_long()], rows of each individual
#' must be contiguous and in chronological order, and baseline
#' characteristics (`bmi`, `ht`, `age`, `sex`, and `PAL`, `EI`, `fat`,
#' `pcarb_base` passed through `...`) are taken from the **first** row of
#' each individual. `EIchange` is read row by row.
#'
#' @param id Individual identifier; rows of an individual must be contiguous.
#' @param ht Height (m): `NULL` (default), a scalar or one value per row,
#'   which may contain `NA`; only the first row of each individual is used.
#'   Missing heights are replaced by the sex-specific `ref_height`, which is
#'   allowed only with `EIchange_unit = "fraction"`.
#' @inheritParams hall_bmi_annual
#' @param ... Further arguments passed to [hall_bw_long()].
#' @return A numeric vector with one value per row: BMI at the **end** of that
#'   row's interval. With `output = "all"`, the list returned by
#'   [hall_bw_long()] with `BMI` computed as in [hall_bmi_annual()], plus
#'   `height` and `height_imputed` (one value per row).
#' @seealso [hall_bmi_annual()] for details.
#' @examples
#' d <- data.frame(pid = rep(1:2, each = 3), year = rep(2025:2027, 2),
#'                 bmi = rep(c(29, 33), each = 3),
#'                 age = rep(c(40, 62), each = 3), sex = rep(c("male", "female"), each = 3),
#'                 dEI = c(-0.02, -0.04, -0.04, 0, -0.03, -0.03))  # fraction of intake
#' d$bmi_end <- with(d, hall_bmi_long(pid, bmi, age, sex, dEI))
#' d
#'
#' # With data.table:
#' # dt[, bmi_end := hall_bmi_long(pid, bmi, age, sex, dEI)]
#' @export
hall_bmi_long <- function(id, bmi, age, sex, EIchange, ht = NULL,
                          EIchange_unit = c("fraction", "kcal"),
                          EI = NA_real_, fat = NA_real_,
                          ref_height = hall_ref_height,
                          output = c("bmi", "all"), ...) {
  EIchange_unit <- match.arg(EIchange_unit); output <- match.arg(output)
  N <- length(id)
  if (N == 0L) stop("'id' has length 0.", call. = FALSE)
  .hall_check_dots(list(...))
  runs <- .hall_id_runs(id)
  bmiN <- .hall_rep_n(bmi, N, "bmi")
  if (any(bmiN <= 0, na.rm = TRUE)) stop("'bmi' must be positive (kg/m^2).", call. = FALSE)
  sexN <- if (length(sex) == 1L) rep(sex, N) else sex
  if (length(sexN) != N) stop(sprintf("'sex' must have length 1 or %d.", N), call. = FALSE)
  # Height of each row; only first rows are used by the solver, but filling
  # every row keeps bw0 = bmi * ht^2 consistent row by row.
  htN <- if (is.null(ht)) NULL else .hall_rep_n(ht, N, "ht")
  EI <- .hall_rep_n(EI, N, "EI"); fat <- .hall_rep_n(fat, N, "fat")
  h <- .hall_bmi_height(htN, .hall_code_sex(sexN), N, EIchange_unit, ref_height,
                        absolute = !is.na(EI) | !is.na(fat), use = runs$first)
  bw0 <- bmiN * h$ht^2

  res <- hall_bw_long(id, bw0, h$ht, age, sexN, EIchange, EI = EI, fat = fat,
                      EIchange_unit = EIchange_unit,
                      output = if (output == "bmi") "bw" else "all", ...)
  BW <- if (output == "bmi") res else res$BW
  # Each person's baseline (first-row) BMI and weight, repeated over their rows
  BMI <- rep.int(bmiN[runs$first], runs$len) * BW / rep.int(bw0[runs$first], runs$len)
  if (output == "bmi") return(BMI)
  res$BMI <- BMI
  res$height <- rep.int(h$ht[runs$first], runs$len)
  res$height_imputed <- rep.int(h$imputed[runs$first], runs$len)
  res
}


# Heights for the BMI interfaces. Returns list(ht, imputed). Missing heights
# (ht NULL or NA) are filled with the sex-specific reference height, which is
# allowed only for proportional changes in intake and for individuals without
# absolute baseline inputs ('absolute': EI or fat supplied). 'use' restricts
# the checks to the rows the solver reads (first row per id).
#' @noRd
.hall_bmi_height <- function(ht, female, n, unit, ref_height, absolute = FALSE,
                             use = seq_len(n)) {
  ht <- if (is.null(ht)) rep(NA_real_, n) else .hall_rep_n(ht, n, "ht")
  imputed <- is.na(ht)
  if (any(imputed[use])) {
    if (unit == "kcal")
      stop("Height is required when 'EIchange' is in kcal/day: at a given BMI the ",
           "change in BMI from a fixed kcal change scales roughly with 1/ht^2. ",
           "Supply 'ht' for everyone, or give 'EIchange' as a fraction of baseline ",
           "intake (EIchange_unit = \"fraction\").", call. = FALSE)
    if (any((imputed & rep_len(absolute, n))[use]))
      stop("'EI' and 'fat' are absolute baseline values and cannot be combined with ",
           "imputed heights. Supply 'ht' for these individuals, or omit 'EI' and 'fat'.",
           call. = FALSE)
    if (!is.numeric(ref_height) || !all(c("female", "male") %in% names(ref_height)) ||
        any(!is.finite(ref_height[c("female", "male")]) | ref_height[c("female", "male")] <= 0))
      stop("'ref_height' must be a numeric vector with positive elements ",
           "'female' and 'male' (m).", call. = FALSE)
    ht[imputed] <- ifelse(female[imputed] == 1L, ref_height[["female"]], ref_height[["male"]])
  }
  list(ht = ht, imputed = imputed)
}

# Arguments of hall_bw_annual()/hall_bw_long() that the BMI interfaces set
# themselves and must not be passed through '...'.
#' @noRd
.hall_check_dots <- function(dots) {
  bad <- intersect(names(dots), c("bw", "output"))
  if (length(bad))
    stop(sprintf("Argument(s) '%s' cannot be passed through '...'; ",
                 paste(bad, collapse = "', '")),
         "the BMI interfaces set them.", call. = FALSE)
  invisible(TRUE)
}
