# Main user-facing solvers.
# Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos. GPL-3.

#' Body weight under the Hall adult model, by year (wide interface)
#'
#' Solves the dynamic adult body-weight model of Hall et al. (2011) for \eqn{n}
#' individuals over \eqn{Y} intervals (years by default). The change in energy
#' intake is held constant within each interval. The function returns body
#' weight at the end of each interval. It is the fast counterpart of running a
#' daily solver such as `bw::adult_weight()` for \eqn{365 Y} days and reading
#' days 365, 730, and so on.
#'
#' @param bw Baseline body weight (kg); numeric vector of length \eqn{n}.
#' @param ht Height (m); length \eqn{n}. Height is fixed, so
#'   \eqn{BMI = BW / ht^2}.
#' @param age Baseline age (years); length \eqn{n}. Age enters only the
#'   baseline equations (resting metabolic rate, body fat and extracellular
#'   fluid), as in the appendix.
#' @param sex Sex; length \eqn{n}. Either character or factor (`"male"`,
#'   `"female"`; also `"men"`/`"women"`, `"m"`/`"f"`; case-insensitive),
#'   logical (`TRUE` = female), or numeric 0/1 (1 = female).
#' @param EIchange Change in energy intake from **baseline** (kcal/day). Give an
#'   \eqn{n \times Y}{n x Y} matrix whose column \eqn{y} applies throughout
#'   interval \eqn{y}, or a length-\eqn{n} vector for a single interval. Changes
#'   are relative to baseline, not to the previous year. With
#'   `appetite_k > 0` this is the intervention "push" (\eqn{-P}), not the
#'   realised change; see *Appetite feedback*.
#' @param NAchange Change in sodium intake from baseline (mg/day): a scalar, a
#'   length-\eqn{n} vector or an \eqn{n \times Y}{n x Y} matrix. Default 0.
#' @param EI Optional baseline energy intake (kcal/day), scalar or length
#'   \eqn{n}. `NA` (default) uses the Mifflin-St Jeor resting metabolic rate
#'   times `PAL`.
#' @param fat Optional baseline fat mass (kg), scalar or length \eqn{n}. `NA`
#'   (default) uses the Jackson et al. equations (appendix eq 4).
#' @param PAL Physical activity level, scalar or length \eqn{n}. Default 1.6;
#'   see *Physical activity level*.
#' @param pcarb_base Baseline share of energy from carbohydrate, in (0, 1].
#'   Default 0.5.
#' @param pcarb Share of energy from carbohydrate after the change: a scalar,
#'   a length-\eqn{n} vector or an \eqn{n \times Y}{n x Y} matrix. Defaults to
#'   `pcarb_base`.
#' @param days Length of each interval in days. Default 365; use for example
#'   `365/4` for quarterly inputs.
#' @param appetite_k Appetite feedback gain \eqn{k} (kcal/day per kg of weight
#'   change; Hall 2024). The default 0 means no feedback (the Lancet 2011
#'   model). It may be a number, a preset name from
#'   [hall_appetite_presets] (for example `"baseline"`, which is 95), a
#'   length-\eqn{n} vector (numeric or preset names), or an
#'   \eqn{n \times Y}{n x Y} numeric matrix so that \eqn{k} can change over
#'   time, for example when a drug is started.
#' @param atol Absolute local error tolerance (kg of lean mass) of the
#'   adaptive integrator. The default `1e-4` gives errors below 0.001 kg per
#'   year in validation; smaller values cost more evaluations.
#' @param model `"hall2011"` (default) follows the web appendix of Hall et al.
#'   (2011): the glycogen energy term is \eqn{\rho_G\,dG/dt} and \eqn{K}
#'   satisfies \eqn{EE(0) = EI(0)}. `"bw"` reproduces the R package 'bw'
#'   v1.0.0, which adds \eqn{dG/dt} in kg/day and always sets \eqn{K} from
#'   RMR \eqn{\times}{x} PAL. Its default PAL is 1.5, so also set `PAL = 1.5`
#'   to reproduce `bw` results.
#' @param output `"bw"` (default) returns the matrix of body weight. `"all"`
#'   returns a list that also contains BMI, the body compartments, the
#'   realised intake change and the baseline values.
#' @param threads Number of OpenMP threads. Defaults to
#'   `getOption("hallfast.threads", 1L)`. It is ignored if the package was
#'   built without OpenMP (see [hall_has_openmp()]).
#'
#' @details
#' **Model.** The model is the four-compartment system of the appendix of
#' Hall et al. (2011): adaptive thermogenesis (AT; eq 7), extracellular fluid
#' (ECF; eq 2), glycogen with 2.7 g of bound water per g (G; eq 1), and lean
#' and fat mass partitioned by \eqn{p = C/(C+F)} with
#' \eqn{C = 10.4\,\mathrm{kg} \times \rho_L/\rho_F}{C = 10.4 kg x rho_L/rho_F}
#' (eq 3). Energy expenditure follows eq 5, the thermic effect of feeding
#' eq 6, and physical activity eq 8. Baseline RMR comes from the Mifflin-St
#' Jeor equations, baseline fat mass from the Jackson et al. equations, and
#' baseline ECF from the sex-specific regressions of Silva et al. (as in 'bw';
#' the appendix does not specify these).
#'
#' **Method.** Within an interval the inputs are constant. AT, ECF and
#' glycogen do not depend on fat or lean mass and have exact closed-form
#' solutions (exponential relaxation, and a Riccati equation for glycogen).
#' Forbes partitioning gives \eqn{F = F_0 \exp((L-L_0)/10.4)}, so fat and lean
#' mass reduce to one scalar ODE. The effect of the 1-2-day ECF and glycogen
#' transients is added as an exact time-integral. Most of the 14-day AT
#' transient is removed by an exact change of variables. The remaining
#' equation is solved with an adaptive Dormand-Prince 5(4) method. This takes
#' about 10-30 right-hand-side evaluations per person-year, against more than
#' 5,000 elementary evaluations for daily fourth-order Runge-Kutta. See
#' `vignette("hallfast")`.
#'
#' **Physical activity level.** The appendix assumed a sedentary PAL of 1.5
#' when estimating the activity parameter (eq 8). The default here is 1.6,
#' because SACN (2011) estimated a median adult PAL of 1.63 from doubly
#' labelled water data (25th and 75th centiles 1.49 and 1.78), and Hall
#' (2024) assumed 1.65, citing SACN. Supply individual or subgroup values
#' where known.
#'
#' **Appetite feedback.** With `appetite_k = k > 0` the change in intake
#' follows Hall (2024),
#' \eqn{\Delta I(t) = -P(t) - k\,\Delta W(t)}{dI(t) = -P(t) - k * dW(t)}.
#' Here `EIchange` plays the role of \eqn{-P}, and the realised change in
#' intake is `EIchange - k * (BW(t) - bw)`, truncated so that total intake
#' cannot become negative. The thermic effect of feeding, AT, carbohydrate
#' intake, glycogen and ECF all respond to the realised intake. Use feedback
#' only when `EIchange` is an intervention's persistent push (for example a
#' reformulation or a drug effect). If `EIchange` is already a net, observed
#' change in intake, compensation is already included and feedback would
#' count it twice. The published \eqn{k} values are fits to mean trajectories
#' in weight-loss trials, so applying them to weight gain or to small
#' population-wide shifts is an extrapolation. See [hall_appetite_presets]
#' and the vignette for details and limitations. With feedback, AT is
#' integrated as a state and ECF and glycogen are treated as quasi-steady at
#' the current intake plus their exact step transients. This path is about
#' 6 times slower.
#'
#' **Accuracy.** Against an independent high-precision solution of the full
#' four-compartment system ([hall_bw_reference()]), the maximum error was
#' 8.5e-5 kg over 10 years in a realistic population, and below 7e-4 kg in a
#' stress test (BMI 16-55, PAL up to 2.4, intake changes from -800 to
#' +500 kcal/day). With feedback, the maximum error was 8.3e-5 kg in the
#' realistic population and 1.8e-3 kg in the stress test. The validation
#' scripts are in the `validation/` folder of the source repository.
#'
#' @return With `output = "bw"`, an \eqn{n \times Y}{n x Y} matrix of body
#'   weight (kg) at the end of each interval. With `output = "all"`, a list
#'   with \eqn{n \times Y}{n x Y} matrices `BW`, `BMI`, `fat_mass`,
#'   `lean_mass`, `ECF`, `glycogen`, `adaptive_thermogenesis` (kcal/day),
#'   `realised_EIchange` (kcal/day, at the end of each interval) and
#'   `rhs_evaluations` (a diagnostic), plus a data frame `baseline` (BW, BMI,
#'   EI, fat_mass, lean_mass). Values are `NA` for individuals with missing
#'   inputs, and from the first interval with a missing input onwards.
#'
#' @references
#' Hall KD, Sacks G, Chandramohan D, Chow CC, Wang YC, Gortmaker SL, Swinburn
#' BA (2011). Quantification of the effect of energy imbalance on bodyweight.
#' *Lancet* 378(9793):826-837. \doi{10.1016/S0140-6736(11)60812-X}
#'
#' Hall KD (2024). Physiology of the weight-loss plateau in response to diet
#' restriction, GLP-1 receptor agonism, and bariatric surgery. *Obesity*
#' 32(6):1163-1168. \doi{10.1002/oby.24027}
#'
#' Scientific Advisory Committee on Nutrition (2011). *Dietary Reference Values
#' for Energy.* London: The Stationery Office.
#'
#' @seealso [hall_bw_long()] for long-format data, [hall_bw_reference()] for
#'   the reference solver, [hall_appetite_presets] for published values of
#'   \eqn{k}.
#' @examples
#' # Four adults; intake 50 kcal/day lower in year 1, then 100 kcal/day lower
#' bw  <- c(70, 85, 100, 120); ht <- c(1.65, 1.75, 1.70, 1.80)
#' age <- c(30, 45, 60, 50);   sex <- c("female", "male", "female", "male")
#' dEI <- matrix(c(-50, -100, -100, -100, -100), 4, 5, byrow = TRUE)
#'
#' BW <- hall_bw_annual(bw, ht, age, sex, dEI)
#' round(BW - bw, 2)                       # change in body weight (kg)
#' round((BW - bw) / ht^2, 3)              # change in BMI
#'
#' # The same push with normal appetite feedback (k = 95 kcal/day per kg)
#' round(hall_bw_annual(bw, ht, age, sex, dEI, appetite_k = "baseline") - bw, 2)
#'
#' # Full output
#' res <- hall_bw_annual(bw, ht, age, sex, dEI, output = "all")
#' names(res)
#'
#' # Reproduce the 'bw' package (v1.0.0)
#' hall_bw_annual(80, 1.8, 40, "female", -100, PAL = 1.5, model = "bw")
#' @export
hall_bw_annual <- function(bw, ht, age, sex, EIchange, NAchange = 0,
                           EI = NA_real_, fat = NA_real_, PAL = 1.6,
                           pcarb_base = 0.5, pcarb = pcarb_base, days = 365,
                           appetite_k = 0, atol = 1e-4,
                           model = c("hall2011", "bw"),
                           output = c("bw", "all"),
                           threads = getOption("hallfast.threads", 1L)) {
  model <- match.arg(model); output <- match.arg(output)
  .hall_check_controls(days, atol, threads)
  n <- length(bw)
  if (n == 0L) stop("'bw' has length 0.", call. = FALSE)
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
  for (nm in c("NAchange", "pcarb", "appetite_k")) {
    v <- get(nm)
    if (is.matrix(v) && !all(dim(v) == c(n, Y)))
      stop(sprintf("'%s' must be %d x %d.", nm, n, Y), call. = FALSE)
  }
  dNa <- .hall_kind(NAchange, n, n * Y, "NAchange")
  pc  <- .hall_kind(pcarb, n, n * Y, "pcarb")
  kk  <- .hall_kind(.hall_resolve_k(appetite_k), n, n * Y, "appetite_k")
  EI0 <- ifelse(is.na(EI), .hall_ei0(bw, ht, age, female, PAL), EI)
  .hall_check_inputs(bw, ht, age, PAL, pcarb_base, pc$v, EI0, dEI, kk$v)

  res <- .hall_core(as.numeric(bw), as.numeric(ht), as.numeric(age), female, PAL,
                    EI, fat, pcarb_base, dEI, dNa$v, dNa$kind, pc$v, pc$kind,
                    kk$v, kk$kind,
                    start = seq.int(0L, n - 1L), len = rep.int(Y, n), stride = n,
                    days = days, atol = atol, bw_compat = model == "bw",
                    full_output = output == "all", nthreads = as.integer(threads))
  BW <- matrix(res$BW, n, Y, dimnames = list(names(bw), paste0("year", seq_len(Y))))
  if (output == "bw") return(BW)

  shape <- function(v) matrix(v, n, Y, dimnames = dimnames(BW))
  list(BW = BW, BMI = BW / ht^2,
       fat_mass = shape(res$FM), lean_mass = shape(res$LM),
       ECF = shape(res$ECF), glycogen = shape(res$GLY),
       adaptive_thermogenesis = shape(res$AT),
       realised_EIchange = shape(res$EIR), rhs_evaluations = shape(res$nevals),
       baseline = data.frame(BW = as.numeric(bw), BMI = as.numeric(bw) / ht^2,
                             EI = res$EI0, fat_mass = res$F0, lean_mass = res$L0))
}


#' Body weight under the Hall adult model, by year (long interface)
#'
#' Long-format version of [hall_bw_annual()], for microsimulation tables with
#' one row per person-year (for example a `data.table` keyed by person and
#' year). Rows of each individual must be contiguous and in chronological
#' order. Baseline characteristics (`bw`, `ht`, `age`, `sex`, `PAL`, `EI`,
#' `fat`, `pcarb_base`) are taken from the **first** row of each individual;
#' values on later rows are ignored. `EIchange`, `NAchange`, `pcarb` and
#' `appetite_k` are read row by row, one interval per row (a single value is
#' recycled).
#'
#' @param id Individual identifier; rows of an individual must be contiguous.
#' @inheritParams hall_bw_annual
#' @return A numeric vector with one value per row: body weight (kg) at the
#'   **end** of that row's interval. With `output = "all"`, a list of such
#'   vectors (`BW`, `BMI`, `fat_mass`, `lean_mass`, `ECF`, `glycogen`,
#'   `adaptive_thermogenesis`, `realised_EIchange`, `rhs_evaluations`).
#' @seealso [hall_bw_annual()] for details of the model, method and arguments.
#' @examples
#' d <- data.frame(pid = rep(1:2, each = 3), year = rep(2025:2027, 2),
#'                 bw = rep(c(80, 95), each = 3), ht = rep(c(1.75, 1.68), each = 3),
#'                 age = rep(c(40, 62), each = 3), sex = rep(c("male", "female"), each = 3),
#'                 dEI = c(-50, -100, -100, 0, -60, -60))
#' d$bw_end <- with(d, hall_bw_long(pid, bw, ht, age, sex, dEI))
#' d$dBMI   <- with(d, (bw_end - bw) / ht^2)
#' d
#'
#' # With data.table:
#' # dt[, bw_end := hall_bw_long(pid, bw, ht, age, sex, dEI)]
#' @export
hall_bw_long <- function(id, bw, ht, age, sex, EIchange, NAchange = 0,
                         EI = NA_real_, fat = NA_real_, PAL = 1.6,
                         pcarb_base = 0.5, pcarb = pcarb_base, days = 365,
                         appetite_k = 0, atol = 1e-4,
                         model = c("hall2011", "bw"),
                         output = c("bw", "all"),
                         threads = getOption("hallfast.threads", 1L)) {
  model <- match.arg(model); output <- match.arg(output)
  .hall_check_controls(days, atol, threads)
  N <- length(id)
  if (N == 0L) stop("'id' has length 0.", call. = FALSE)
  r <- rle(if (is.factor(id)) as.integer(id) else id)
  if (anyDuplicated(r$values))
    stop("Rows of each id must be contiguous (sort by id, then by time).", call. = FALSE)
  len   <- as.integer(r$lengths)
  first <- c(1L, cumsum(len)[-length(len)] + 1L)          # 1-based first row per id
  rowN  <- function(x, name) .hall_rep_n(x, N, name)
  sub   <- function(x, name) rowN(x, name)[first]

  bw0  <- sub(bw, "bw"); ht0 <- sub(ht, "ht"); age0 <- sub(age, "age")
  female <- .hall_code_sex(if (length(sex) == 1L) rep(sex, N) else sex)[first]
  PAL0 <- sub(PAL, "PAL"); EIb <- sub(EI, "EI"); fat0 <- sub(fat, "fat")
  pb0  <- sub(pcarb_base, "pcarb_base")
  dEI  <- rowN(EIchange, "EIchange")
  per_row <- function(x, name)
    if (length(x) == 1L) list(v = as.numeric(x), kind = 0L) else list(v = rowN(x, name), kind = 2L)
  dNa <- per_row(NAchange, "NAchange")
  pc  <- per_row(pcarb, "pcarb")
  kk  <- per_row(.hall_resolve_k(appetite_k), "appetite_k")
  EI0 <- ifelse(is.na(EIb), .hall_ei0(bw0, ht0, age0, female, PAL0), EIb)
  .hall_check_inputs(bw0, ht0, age0, PAL0, pb0, pc$v, rep.int(EI0, len), dEI, kk$v)

  res <- .hall_core(bw0, ht0, age0, female, PAL0, EIb, fat0, pb0, dEI,
                    dNa$v, dNa$kind, pc$v, pc$kind, kk$v, kk$kind,
                    start = first - 1L, len = len, stride = 1L,
                    days = days, atol = atol, bw_compat = model == "bw",
                    full_output = output == "all", nthreads = as.integer(threads))
  if (output == "bw") return(res$BW)
  list(BW = res$BW, BMI = res$BW / rep.int(ht0, len)^2, fat_mass = res$FM,
       lean_mass = res$LM, ECF = res$ECF, glycogen = res$GLY,
       adaptive_thermogenesis = res$AT, realised_EIchange = res$EIR,
       rhs_evaluations = res$nevals)
}
