# Published appetite feedback gains.
# Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos. GPL-3.

#' Published appetite feedback gains (k)
#'
#' Fitted values of the proportional appetite feedback gain \eqn{k} (kcal/day
#' per kg of weight change) in
#' \deqn{\Delta I(t) = -P(t) - k \, \Delta W(t),}{dI(t) = -P(t) - k * dW(t),}
#' the intake-control extension of the Hall et al. (2011) model (Polidori et
#' al. 2016; Hall, Sanghvi and Gobel 2017; Hall 2024). Any entry of
#' `preset` can be passed as `appetite_k` to [hall_bw_annual()] or
#' [hall_bw_long()].
#'
#' Each \eqn{k} is the best fit to the *mean* body-weight trajectory of one
#' trial (Hall 2024, Table 1). The `"baseline"` value comes from covert energy
#' loss under SGLT2 inhibition (Polidori et al. 2016) and is used by Hall
#' (2024) as the strength of the feedback circuit in the absence of an
#' intervention on appetite. The Results text and Table 1 of Hall (2024) give
#' \eqn{k = 82} for CALERIE; the Discussion quotes 83. The semaglutide fit with
#' energy partitioning altered about 3-fold (\eqn{k = 31}) is not included,
#' because this package keeps the Forbes partitioning of the original model.
#'
#' `P0_kcal` and `Pmax_kcal` are the fitted intervention magnitudes (kcal/day;
#' intake cuts are positive) at the start and after dose escalation. They are
#' given for reference only: they were fitted to the trial populations, with
#' baseline PAL assumed to be 1.65, and are not transferable.
#'
#' @format A data frame with 6 rows and 6 columns:
#' \describe{
#'   \item{preset}{Name to pass as `appetite_k`.}
#'   \item{k}{Feedback gain, kcal/day per kg.}
#'   \item{P0_kcal}{Fitted initial intervention magnitude P (kcal/day), or `NA`.}
#'   \item{Pmax_kcal}{Fitted maximum intervention magnitude P (kcal/day), or `NA`.}
#'   \item{setting}{Trial and intervention.}
#'   \item{source}{Publication.}
#' }
#' @references
#' Polidori D, Sanghvi A, Seeley RJ, Hall KD (2016). How strongly does appetite
#' counter weight loss? Quantification of the feedback control of human energy
#' intake. *Obesity* 24(11):2289-2295. \doi{10.1002/oby.21653}
#'
#' Hall KD, Sanghvi A, Gobel B (2017). Proportional feedback control of energy
#' intake during obesity pharmacotherapy. *Obesity* 25(12):2088-2091.
#' \doi{10.1002/oby.21978}
#'
#' Hall KD (2024). Physiology of the weight-loss plateau in response to diet
#' restriction, GLP-1 receptor agonism, and bariatric surgery. *Obesity*
#' 32(6):1163-1168. \doi{10.1002/oby.24027}
#' @seealso [hall_bw_annual()] (argument `appetite_k`).
#' @examples
#' hall_appetite_presets[, c("preset", "k")]
#' @export
hall_appetite_presets <- data.frame(
  preset = c("baseline", "calorie_restriction", "diet_restriction",
             "tirzepatide", "semaglutide", "rygb"),
  k = c(95, 82, 101, 48, 49, 58),
  P0_kcal = c(NA, 830, 1200, 830, 610, 3600),
  Pmax_kcal = c(NA, 830, 1200, 1560, 1300, 3600),
  setting = c(
    "No intervention on appetite: SGLT2 inhibition (canagliflozin 300 mg/d, 52 wk, n = 153, type 2 diabetes)",
    "Intensive calorie restriction, CALERIE phase 2 (24 months)",
    "Diet restriction (low-fat and low-carbohydrate arms combined), DIETFITS (12 months)",
    "Tirzepatide 10 mg weekly (SURMOUNT-1), 12-week dose escalation",
    "Semaglutide 2.4 mg weekly (STEP 1), 16-week dose escalation, standard Forbes partitioning",
    "Roux-en-Y gastric bypass surgery"),
  source = c("Polidori et al. Obesity 2016;24:2289-95, doi:10.1002/oby.21653",
             rep("Hall. Obesity 2024;32:1163-8, Table 1, doi:10.1002/oby.24027", 5)),
  stringsAsFactors = FALSE)

# Resolve appetite_k: numeric (>= 0) or preset name(s) -> numeric vector
# (matrix dimensions, if any, are dropped; column-major order is kept).
#' @noRd
.hall_resolve_k <- function(k) {
  if (is.character(k) || is.factor(k)) {
    k <- as.character(k)
    m <- match(tolower(k), hall_appetite_presets$preset)
    if (anyNA(m))
      stop("Unknown appetite_k preset(s): ", paste(unique(k[is.na(m)]), collapse = ", "),
           ". Use one of: ", paste(hall_appetite_presets$preset, collapse = ", "), ".",
           call. = FALSE)
    return(hall_appetite_presets$k[m])
  }
  k <- as.numeric(k)
  if (any(k < 0, na.rm = TRUE)) stop("'appetite_k' must be >= 0 (kcal/day per kg).", call. = FALSE)
  k
}
