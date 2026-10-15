# Internal helpers for argument handling and checks.
# Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos. GPL-3.

# Code sex as 1 = female, 0 = male. Accepts character/factor
# ("female"/"male", "women"/"men", "f"/"m"; case-insensitive), logical
# (TRUE = female) or numeric 0/1 (1 = female, as in bw's C++ code).
#' @noRd
.hall_code_sex <- function(sex) {
  if (is.factor(sex)) sex <- as.character(sex)
  if (is.character(sex)) {
    s <- tolower(trimws(sex))
    out <- rep(NA_integer_, length(s))
    out[s %in% c("female", "women", "woman", "f", "w")] <- 1L
    out[s %in% c("male", "men", "man", "m")] <- 0L
    if (anyNA(out[!is.na(s)]))
      stop("Unrecognised values in 'sex'. Use 'male'/'female' or 0/1 (1 = female).", call. = FALSE)
    return(out)
  }
  if (is.logical(sex)) return(as.integer(sex))
  if (is.numeric(sex)) {
    if (any(!sex %in% c(0, 1, NA)))
      stop("Numeric 'sex' must be 0 (male) or 1 (female).", call. = FALSE)
    return(as.integer(sex))
  }
  stop("'sex' must be character, factor, logical or numeric 0/1.", call. = FALSE)
}

# Recycle a scalar or length-n vector to length n (numeric).
#' @noRd
.hall_rep_n <- function(x, n, name) {
  if (length(x) == 1L) return(rep(as.numeric(x), n))
  if (length(x) != n) stop(sprintf("'%s' must have length 1 or %d.", name, n), call. = FALSE)
  as.numeric(x)
}

# Expand scalar / length-n vector / n x Y matrix to an n x Y numeric matrix.
#' @noRd
.hall_as_matrix <- function(x, n, Y, name) {
  if (is.matrix(x)) {
    if (nrow(x) != n || ncol(x) != Y)
      stop(sprintf("'%s' must be %d x %d.", name, n, Y), call. = FALSE)
    storage.mode(x) <- "double"
    return(x)
  }
  if (length(x) == 1L || length(x) == n) return(matrix(as.numeric(x), n, Y))
  stop(sprintf("'%s' must be a scalar, a length-%d vector or a %d x %d matrix.", name, n, n, Y),
       call. = FALSE)
}

# Classify an input that may be constant (kind 0), per individual (kind 1,
# length n) or per person-year (kind 2, length nrow or an n x Y matrix) without
# expanding it to a full matrix (saves memory for very large populations).
#' @noRd
.hall_kind <- function(x, n, nrow, name) {
  x <- as.numeric(x)   # drops matrix dim; column-major order is preserved
  if (length(x) == 1L)   return(list(v = x, kind = 0L))
  if (length(x) == nrow) return(list(v = x, kind = 2L))
  if (length(x) == n)    return(list(v = x, kind = 1L))
  stop(sprintf("'%s' must have length 1, %d (per individual) or %d (per person-year).",
               name, n, nrow), call. = FALSE)
}

# Input checks shared by all interfaces.
#' @noRd
.hall_check_inputs <- function(bw, ht, age, PAL, pcarb_base, pcarb, EI0, dEI, k = 0) {
  if (any(bw <= 0, na.rm = TRUE) || any(ht <= 0, na.rm = TRUE))
    stop("'bw' and 'ht' must be positive (kg, m).", call. = FALSE)
  if (any(ht > 3, na.rm = TRUE))
    warning("Some heights exceed 3: 'ht' must be in metres.", call. = FALSE)
  if (any(age < 18, na.rm = TRUE))
    warning("Some ages are < 18: the adult model is not intended for children.", call. = FALSE)
  if (any(PAL <= 0, na.rm = TRUE)) stop("'PAL' must be positive.", call. = FALSE)
  if (any(pcarb_base <= 0 | pcarb_base > 1, na.rm = TRUE))
    stop("'pcarb_base' must be in (0, 1].", call. = FALSE)
  if (any(pcarb < 0 | pcarb > 1, na.rm = TRUE)) stop("'pcarb' must be in [0, 1].", call. = FALSE)
  # Without feedback a negative total intake is an input error. With feedback,
  # EIchange is a push that may exceed baseline intake (e.g. bariatric surgery);
  # the realised intake is then truncated at zero inside the solver.
  if (any(EI0 + dEI < 0 & k == 0, na.rm = TRUE))
    stop("Baseline EI + EIchange is negative for some person-years (appetite_k = 0).",
         call. = FALSE)
  invisible(TRUE)
}

# Mifflin-St Jeor RMR x PAL (the default baseline energy intake).
#' @noRd
.hall_ei0 <- function(bw, ht, age, female, PAL)
  (9.99 * bw + 625 * ht - 4.92 * age + ifelse(female == 1L, -161, 5)) * PAL

# Validate scalar control arguments.
#' @noRd
.hall_check_controls <- function(days, atol, threads) {
  if (!is.numeric(days) || length(days) != 1L || !is.finite(days) || days <= 0)
    stop("'days' must be a single positive number.", call. = FALSE)
  if (!is.numeric(atol) || length(atol) != 1L || !is.finite(atol) || atol <= 0)
    stop("'atol' must be a single positive number.", call. = FALSE)
  if (!is.numeric(threads) || length(threads) != 1L || threads < 1)
    stop("'threads' must be a single integer >= 1.", call. = FALSE)
  invisible(TRUE)
}
