# =============================================================================
#  validate_vs_bw.R
#
#  Validates hall_fast (annual-resolution solver) against
#    (1) an independent high-precision reference: the full 4-compartment ODE
#        system of the Hall et al. (2011) web appendix, integrated with RK4 at
#        dt = 1/20 day (hall_bw_reference()). This is the gold standard.
#    (2) bw::adult_weight() (dt = 1 day, and dt = 0.1 day on a subsample), for
#        which hall_fast(model = "bw") is the like-for-like comparator.
#  It also quantifies how much bw's deviations from the appendix matter,
#  validates the optional appetite-feedback extension (Hall 2024) against the
#  same reference ODE with feedback, and benchmarks speed. Results are written
#  to validation/results/.
#
#  Usage: from the root of the source repository,
#    Rscript validation/validate_vs_bw.R
#  Requirements: the installed 'hallfast' package, 'bw' (CRAN, used only as a
#  comparator) and 'data.table'. Run time: about 15 minutes, dominated by the
#  bw runs. Results are written to validation/results/ and copied to
#  inst/extdata/validation/ (used by the vignette).
# =============================================================================

suppressPackageStartupMessages({ library(hallfast); library(bw); library(data.table) })

# ---- locate files ----------------------------------------------------------
args    <- commandArgs(trailingOnly = FALSE)
this    <- sub("^--file=", "", args[grep("^--file=", args)])
HF_DIR  <- if (length(this)) normalizePath(file.path(dirname(this), "..")) else getwd()
OUT_DIR <- file.path(HF_DIR, "validation", "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

has_omp <- hall_has_openmp()
cat("hallfast", as.character(packageVersion("hallfast")), "- OpenMP available:", has_omp, "\n")

# ---- helpers ---------------------------------------------------------------
# Run bw::adult_weight with annual step inputs; return end-of-year BW and tissue
# (F + L). bw needs one more column than days, otherwise it stops a day early.
run_bw <- function(p, dEI, dNa, dt = 1, EI = NA, fat = rep(NA, nrow(p))) {
  Y <- ncol(dEI); spd <- round(1 / dt)
  idx <- rep(seq_len(Y), each = 365 * spd)
  EId <- cbind(dEI[, idx, drop = FALSE], dEI[, Y])
  NAd <- cbind(dNa[, idx, drop = FALSE], dNa[, Y])
  r <- suppressWarnings(adult_weight(p$bw, p$ht, p$age, ifelse(p$female == 1, "female", "male"),
                                     EId, NAd, EI = EI, fat = fat, PAL = p$PAL,
                                     pcarb_base = p$pcarb_base, pcarb = p$pcarb,
                                     days = 365 * Y, dt = dt))
  cols <- 365 * spd * seq_len(Y) + 1
  list(BW = r$Body_Weight[, cols, drop = FALSE],
       tissue = r$Fat_Mass[, cols, drop = FALSE] + r$Lean_Mass[, cols, drop = FALSE])
}

# Independent reference: full 4-compartment system, RK4 at dt = 0.05 day
ref_ode <- function(p, dEI, dNa, EI = rep(NA_real_, nrow(p)), fat = rep(NA_real_, nrow(p)),
                    bw_compat = FALSE, tissue_only = FALSE, k = 0) {
  r <- hall_bw_reference(p$bw, p$ht, p$age, p$female, dEI, NAchange = dNa, EI = EI, fat = fat,
                         PAL = p$PAL, pcarb_base = p$pcarb_base, pcarb = p$pcarb,
                         appetite_k = k, model = if (bw_compat) "bw" else "hall2011",
                         steps_per_day = 20, output = "all")
  unname(if (tissue_only) r$tissue else r$BW)
}

fast <- function(p, dEI, dNa = 0, EI = NA_real_, fat = NA_real_, output = "bw", ...) {
  r <- hall_bw_annual(p$bw, p$ht, p$age, p$female, dEI, NAchange = dNa, EI = EI, fat = fat,
                      PAL = p$PAL, pcarb_base = p$pcarb_base, pcarb = p$pcarb, output = output, ...)
  if (output == "bw") unname(r) else lapply(r, function(z) if (is.matrix(z)) unname(z) else z)
}

err_row <- function(section, population, comparison, est, ref, ht, base = NULL) {
  e  <- abs(est - ref)
  out <- data.table(section, population, comparison,
                    max_abs_kg = max(e), p99_abs_kg = unname(quantile(e, 0.99)),
                    mean_abs_kg = mean(e), max_abs_BMI = max(e / ht^2))
  if (!is.null(base)) {                 # error relative to the simulated change
    chg <- abs(ref - base); k <- chg > 0.5
    out[, max_rel_err_of_change_pct := 100 * max(e[k] / chg[k])]
  }
  out
}

# Linearised closed form of the appendix (eqs 10-15), for context only.
# Piecewise-constant dEI: x_y = x_{y-1} e^{-T/tau} + (dEI_y / eps)(1 - e^{-T/tau}).
hall_linear <- function(p, dEI, F0) {
  roF <- 9440.727; roL <- 1816.444; gF <- 3.107075; gL <- 21.98853
  eF <- 179.2543; eL <- 229.4455; beta <- 0.1 + 0.14
  rmr <- 9.99 * p$bw + 625 * p$ht - 4.92 * p$age + ifelse(p$female == 1, -161, 5)
  delta <- ((1 - 0.1) * p$PAL - 1) * rmr / p$bw
  a   <- 10.4 / F0
  rho <- (eF + roF + a * eL + a * roL) / ((1 - beta) * (1 + a))     # eq 12
  eps <- ((gF + a * gL) / (1 + a) + delta) / (1 - beta)              # eq 13
  tau <- rho / eps; x <- 0; out <- dEI
  for (y in seq_len(ncol(dEI))) {
    x <- x * exp(-365 / tau) + dEI[, y] / eps * (1 - exp(-365 / tau)); out[, y] <- p$bw + x
  }
  out
}

make_pop <- function(n, type, seed) {
  set.seed(seed)
  female <- rbinom(n, 1, 0.5)
  ht <- ifelse(female == 1, rnorm(n, 1.62, 0.065), rnorm(n, 1.76, 0.07))
  if (type == "realistic") {
    age <- runif(n, 20, 90)
    bmi <- pmin(pmax(rlnorm(n, log(27), 0.18), 16), 60)
    PAL <- runif(n, 1.4, 1.9); pcb <- rep(0.5, n); pc <- pcb
  } else {
    age <- runif(n, 18, 95); bmi <- runif(n, 16, 55); PAL <- runif(n, 1.4, 2.4)
    pcb <- runif(n, 0.4, 0.6); pc <- runif(n, 0.3, 0.7)
  }
  data.table(female, age, ht, bw = bmi * ht^2, PAL, pcarb_base = pcb, pcarb = pc)
}

results <- list()
add <- function(...) results[[length(results) + 1L]] <<- err_row(...)

# =============================================================================
# Realistic population, 10 years, policy-like intake changes (random walk)
# =============================================================================
n <- 3000; Y <- 10; popR <- "realistic (n=3000, 10 y)"
p <- make_pop(n, "realistic", 2026)
dEI <- matrix(0, n, Y); dEI[, 1] <- rnorm(n, -50, 50)
for (y in 2:Y) dEI[, y] <- dEI[, y - 1] + rnorm(n, 0, 30)
dEI <- pmin(pmax(dEI, -400), 300)
dNa <- matrix(0, n, Y)
base <- matrix(p$bw, n, Y)

f_h  <- fast(p, dEI, output = "all")                         # appendix model
f_b  <- fast(p, dEI, model = "bw", output = "all")           # bw-compatible model
r_h  <- ref_ode(p, dEI, dNa)
r_ht <- ref_ode(p, dEI, dNa, tissue_only = TRUE)
r_b  <- ref_ode(p, dEI, dNa, bw_compat = TRUE)
t_bw <- system.time(b1 <- run_bw(p, dEI, dNa))[["elapsed"]]

# A. Accuracy against the gold standard (appendix)
add("A. appendix (gold standard)", popR, "hall_fast vs reference ODE", f_h$BW, r_h, p$ht, base)
add("A. appendix (gold standard)", popR, "fat+lean: hall_fast vs reference ODE",
    f_h$fat_mass + f_h$lean_mass, r_ht, p$ht)
add("A. appendix (gold standard)", popR, "linearised appendix eq 14 vs reference ODE (context)",
    hall_linear(p, dEI, f_h$baseline$fat_mass), r_h, p$ht, base)
add("A. appendix (gold standard)", popR, "bw (dt = 1) vs reference ODE (appendix)", b1$BW, r_h, p$ht, base)

# B. Like-for-like reproduction of bw (model = "bw")
add("B. bw compatibility", popR, "hall_fast(model='bw') vs reference ODE (bw variant)", f_b$BW, r_b, p$ht, base)
add("B. bw compatibility", popR, "bw (dt = 1) vs reference ODE (bw variant)", b1$BW, r_b, p$ht, base)
add("B. bw compatibility", popR, "hall_fast(model='bw') vs bw (dt = 1)", f_b$BW, b1$BW, p$ht, base)
add("B. bw compatibility", popR, "fat+lean: hall_fast(model='bw') vs bw (dt = 1)",
    f_b$fat_mass + f_b$lean_mass, b1$tissue, p$ht)

# bw's year-boundary look-ahead: error explained by next year's change in intake
jump <- dEI[, -1] - dEI[, -Y]
lk <- lm(as.vector((b1$BW - r_b)[, -Y]) ~ as.vector(jump))
lookahead <- data.table(slope_kg_per_1000kcal_jump = 1000 * unname(coef(lk)[2]),
                        r_squared = summary(lk)$r.squared)
fwrite(lookahead, file.path(OUT_DIR, "bw_lookahead_artefact.csv"))

# bw convergence with dt = 0.1 day on a subsample
ss <- 1:150; popS <- "realistic subsample (n=150, 10 y)"
b01 <- run_bw(p[ss], dEI[ss, ], dNa[ss, ], dt = 0.1)
add("B. bw compatibility", popS, "bw (dt = 0.1) vs reference ODE (bw variant)", b01$BW, r_b[ss, ], p$ht[ss])
add("B. bw compatibility", popS, "hall_fast(model='bw') vs bw (dt = 0.1)", f_b$BW[ss, ], b01$BW, p$ht[ss])

# Single-year bw runs (no intake change inside the run)
b1y <- run_bw(p, dEI[, 1, drop = FALSE], dNa[, 1, drop = FALSE])
add("B. bw compatibility", "realistic (n=3000, 1 y)", "hall_fast(model='bw') vs bw (dt = 1), single year",
    f_b$BW[, 1], b1y$BW[, 1], p$ht)

# C. How much do bw's deviations from the appendix matter?
add("C. size of bw deviations", popR, "glycogen term: appendix vs bw variant (both exact solutions)", r_h, r_b, p$ht)
# K initialisation: baseline EI supplied 5% above RMR x PAL and NO change in intake.
# Appendix: K satisfies EE(0) = EI(0), so weight stays at baseline.
# bw: K from RMR x PAL, so the 5% surplus acts as a permanent energy imbalance.
EI5 <- 1.05 * (9.99 * p$bw + 625 * p$ht - 4.92 * p$age + ifelse(p$female == 1, -161, 5)) * p$PAL
zero <- matrix(0, n, Y)
drift_h <- fast(p, zero, EI = EI5)
drift_b <- fast(p, zero, EI = EI5, model = "bw")
drift <- data.table(model = c("appendix (hall2011)", "bw variant"),
                    mean_EI_surplus_kcal = mean(EI5 / 1.05 * 0.05),
                    mean_drift_year1_kg  = c(mean(drift_h[, 1] - p$bw), mean(drift_b[, 1] - p$bw)),
                    mean_drift_year10_kg = c(mean(drift_h[, Y] - p$bw), mean(drift_b[, Y] - p$bw)),
                    max_abs_drift_year10_kg = c(max(abs(drift_h[, Y] - p$bw)), max(abs(drift_b[, Y] - p$bw))))
fwrite(drift, file.path(OUT_DIR, "bw_K_initialisation_drift.csv"))

# =============================================================================
# Stress test: extreme BMI/PAL, large random intake swings, sodium and
# carbohydrate changes, and all four initialisation paths (EI and/or fat given)
# =============================================================================
n2 <- 2000; Y2 <- 5
q <- make_pop(n2, "stress", 7)
dEI2 <- matrix(runif(n2 * Y2, -800, 500), n2, Y2)
dNa2 <- matrix(runif(n2 * Y2, -1500, 1500), n2, Y2)
paths <- list(none = 1:500, EI = 501:1000, fat = 1001:1500, both = 1501:2000)
EIin  <- rep(NA_real_, n2); fatin <- rep(NA_real_, n2)
EIin[c(paths$EI, paths$both)]   <- runif(1000, 1600, 3200)
fatin[c(paths$fat, paths$both)] <- q$bw[c(paths$fat, paths$both)] * runif(1000, 0.12, 0.45)
f2h <- fast(q, dEI2, dNa2, EI = EIin, fat = fatin, output = "all")
f2b <- fast(q, dEI2, dNa2, EI = EIin, fat = fatin, model = "bw", output = "all")
r2h <- ref_ode(q, dEI2, dNa2, EIin, fatin)
r2b <- ref_ode(q, dEI2, dNa2, EIin, fatin, bw_compat = TRUE)
stiff <- list()
for (nm in names(paths)) {
  ii <- paths[[nm]]
  bq <- run_bw(q[ii], dEI2[ii, ], dNa2[ii, ],
               EI = if (nm %in% c("EI", "both")) EIin[ii] else NA,
               fat = if (nm %in% c("fat", "both")) fatin[ii] else rep(NA, length(ii)))
  # bw's glycogen equation is stiff for high carbohydrate intakes: its linearised
  # decay rate at baseline is 2 k_G G0 / rho_G = 4 CI_b / rho_G per day (G0 = 0.5 kg),
  # and classical RK4 with dt = 1 is unstable when this exceeds ~2.785.
  EI0q <- ifelse(is.na(EIin[ii]), (9.99 * q$bw[ii] + 625 * q$ht[ii] - 4.92 * q$age[ii] +
                                   ifelse(q$female[ii] == 1, -161, 5)) * q$PAL[ii], EIin[ii])
  lam  <- 4 * q$pcarb_base[ii] * EI0q / 4206.501
  stiff[[nm]] <- data.table(rate = lam, err = apply(abs(f2b$BW[ii, ] - bq$BW), 1, max))
  lab <- sprintf("stress, baseline given: %s (n=500, 5 y)", nm)
  add("A. appendix (gold standard)", lab, "hall_fast vs reference ODE", f2h$BW[ii, ], r2h[ii, ], q$ht[ii])
  add("B. bw compatibility", lab, "hall_fast(model='bw') vs reference ODE (bw variant)", f2b$BW[ii, ], r2b[ii, ], q$ht[ii])
  add("B. bw compatibility", lab, "fat+lean: hall_fast(model='bw') vs bw (dt = 1)",
      (f2b$fat_mass + f2b$lean_mass)[ii, ], bq$tissue, q$ht[ii])
  add("B. bw compatibility", lab, "hall_fast(model='bw') vs bw (dt = 1)", f2b$BW[ii, ], bq$BW, q$ht[ii])
  add("C. size of bw deviations", lab, "appendix vs bw variant (both exact solutions)", r2h[ii, ], r2b[ii, ], q$ht[ii])
}
stiff <- rbindlist(stiff)
stiff[, band := cut(rate, c(0, 1, 2, 2.785, Inf), right = FALSE,
                    labels = c("<1 /day", "1-2 /day", "2-2.785 /day", ">=2.785 /day (RK4 dt=1 unstable)"))]
stiff_tab <- stiff[, .(n = .N, max_abs_err_kg_hall_fast_bw_vs_bw = max(err),
                       mean_abs_err_kg = mean(err)), keyby = band]
fwrite(stiff_tab, file.path(OUT_DIR, "bw_glycogen_stiffness.csv"))

evals <- data.table(population = c("realistic", "stress"),
                    mean_rhs_evals_per_person_year = c(mean(f_h$rhs_evaluations), mean(f2h$rhs_evaluations)),
                    max_rhs_evals_per_person_year  = c(max(f_h$rhs_evaluations), max(f2h$rhs_evaluations)),
                    range_BW_kg = c(paste(round(range(r_h), 1), collapse = " to "),
                                    paste(round(range(r2h), 1), collapse = " to ")))
fwrite(evals, file.path(OUT_DIR, "integrator_cost.csv"))

# Tolerance sensitivity on the stress population (appendix model)
tol <- rbindlist(lapply(c(1e-3, 1e-4, 1e-5, 1e-6), function(a) {
  f <- fast(q, dEI2, dNa2, EI = EIin, fat = fatin, atol = a, output = "all")
  data.table(atol = a, max_abs_err_kg = max(abs(f$BW - r2h)), mean_abs_err_kg = mean(abs(f$BW - r2h)),
             mean_rhs_evals = mean(f$rhs_evaluations))
}))
fwrite(tol, file.path(OUT_DIR, "tolerance_sensitivity.csv"))

acc <- rbindlist(results, fill = TRUE)
fwrite(acc, file.path(OUT_DIR, "accuracy.csv"))


# =============================================================================
# D. Appetite feedback (Hall 2024): hall_fast(appetite_k = k) vs the reference
#    ODE in which intake = EIchange - k (BW - BW0) is evaluated from the full
#    4-compartment state at every RK4 stage
# =============================================================================
fb <- list()
addfb <- function(...) fb[[length(fb) + 1L]] <<- err_row(...)
for (kn in c("baseline", "calorie_restriction", "tirzepatide", "rygb")) {
  kval <- hall_appetite_presets$k[hall_appetite_presets$preset == kn]
  fk <- fast(p, dEI, appetite_k = kn, output = "all")
  rk <- ref_ode(p, dEI, dNa, k = kval)
  addfb("D. appetite feedback", popR, sprintf("hall_fast vs reference ODE, k = %g (%s)", kval, kn), fk$BW, rk, p$ht, base)
}
# Time-varying k: baseline appetite control (95) switching to tirzepatide (48) in year 4
kmat <- cbind(matrix(95, n, 3), matrix(48, n, Y - 3))
fkv <- fast(p, dEI, appetite_k = kmat)
rkv <- ref_ode(p, dEI, dNa, k = kmat)
addfb("D. appetite feedback", popR, "hall_fast vs reference ODE, k = 95 then 48 from year 4", fkv, rkv, p$ht, base)
# Consistency of the two code paths: k -> 0 through the feedback path vs k = 0
fk0 <- fast(p, dEI, appetite_k = 1e-9)
addfb("D. appetite feedback", popR, "feedback path with k = 1e-9 vs reference ODE with k = 0", fk0, r_h, p$ht, base)
# Stress: random k, pushes up to 1,500 kcal/d (some beyond baseline intake), sodium swings
kq  <- runif(n2, 30, 110)
dEIq <- matrix(runif(n2 * Y2, -1500, 500), n2, Y2)
fq <- fast(q, dEIq, dNa2, appetite_k = kq, output = "all")
rq <- ref_ode(q, dEIq, dNa2, k = kq)
addfb("D. appetite feedback", "stress (n=2000, 5 y), k ~ U(30, 110), push -1500 to +500 kcal/d",
      "hall_fast vs reference ODE", fq$BW, rq, q$ht)
okq <- apply(dEIq, 1, min) > -0.8 * fq$baseline$EI
addfb("D. appetite feedback", sprintf("stress subset, push never > 80%% of baseline intake (n=%d)", sum(okq)),
      "hall_fast vs reference ODE", fq$BW[okq, ], rq[okq, ], q$ht[okq])
fbacc <- rbindlist(fb, fill = TRUE)
fwrite(fbacc, file.path(OUT_DIR, "appetite_feedback_accuracy.csv"))
fbcost <- data.table(population = c("realistic, k = 95", "stress"),
                     mean_rhs_evals_per_person_year = c(mean(fast(p, dEI, appetite_k = 95, output = "all")$rhs_evaluations),
                                                         mean(fq$rhs_evaluations)))
fwrite(fbcost, file.path(OUT_DIR, "appetite_feedback_cost.csv"))

# =============================================================================
# Speed
# =============================================================================
med3 <- function(f) median(replicate(3, system.time(f())[["elapsed"]]))
bench <- list()
# (a) one year, n = 10,000
nb <- 10000; pb <- make_pop(nb, "realistic", 11); dEb <- matrix(rnorm(nb, -50, 50), nb, 1)
tb <- system.time(run_bw(pb, dEb, dEb * 0))[["elapsed"]]
tf <- med3(function() fast(pb, dEb))
bench[[1]] <- data.table(task = "n = 10,000, 1 year", person_years = 1e4, bw_s = tb, hall_fast_s = tf)
# (b) 10 years, n = 3,000
tf <- med3(function() fast(p, dEI))
bench[[2]] <- data.table(task = "n = 3,000, 10 years", person_years = 3e4, bw_s = t_bw, hall_fast_s = tf)
# (c) 1 million people x 10 years, hall_fast only (bw would need ~30 GB of daily matrices)
nm <- 1e6; pm <- make_pop(nm, "realistic", 12)
dEm <- matrix(rnorm(nm * 10, -50, 50), nm, 10)
tf1 <- med3(function() fast(pm, dEm, threads = 1))
bench[[3]] <- data.table(task = "n = 1,000,000, 10 years, 1 thread", person_years = 1e7, bw_s = NA_real_, hall_fast_s = tf1)
if (has_omp && parallel::detectCores() > 1) {
  nth <- parallel::detectCores()
  tfn <- med3(function() fast(pm, dEm, threads = nth))
  bench[[4]] <- data.table(task = sprintf("n = 1,000,000, 10 years, %d threads", nth),
                           person_years = 1e7, bw_s = NA_real_, hall_fast_s = tfn)
}
# (c2) appetite feedback (k = 95), 1 million people x 10 years
tfk <- med3(function() fast(pm, dEm, appetite_k = 95, threads = 1))
bench[[length(bench) + 1L]] <- data.table(task = "n = 1,000,000, 10 years, appetite feedback k = 95, 1 thread",
                                          person_years = 1e7, bw_s = NA_real_, hall_fast_s = tfk)
# (c3) sustained intake change from year 4 (typical policy scenario), no feedback
dEs <- cbind(matrix(0, nm, 3), matrix(rnorm(nm, -50, 50), nm, 7))
tfs <- med3(function() fast(pm, dEs, threads = 1))
bench[[length(bench) + 1L]] <- data.table(task = "n = 1,000,000, 10 years, sustained change from year 4, 1 thread",
                                          person_years = 1e7, bw_s = NA_real_, hall_fast_s = tfs)
# (d) long format via data.table, 1 million people x 10 years
DT <- data.table(pid = rep(seq_len(nm), each = 10), year = rep(1:10, nm),
                 bw = rep(pm$bw, each = 10), ht = rep(pm$ht, each = 10), age = rep(pm$age, each = 10),
                 sex = rep(pm$female, each = 10), PAL = rep(pm$PAL, each = 10), dEI = as.vector(t(dEm)))
tl <- med3(function() DT[, bw_hall := hall_bw_long(pid, bw, ht, age, sex, dEI, PAL = PAL)])
bench[[length(bench) + 1L]] <- data.table(task = "n = 1,000,000, 10 years, long format in data.table, 1 thread",
                                          person_years = 1e7, bw_s = NA_real_, hall_fast_s = tl)
stopifnot(isTRUE(all.equal(DT$bw_hall, as.vector(t(fast(pm, dEm))))))   # long == wide

bench <- rbindlist(bench)
bench[, `:=`(bw_us_per_person_year = 1e6 * bw_s / person_years,
             hall_fast_us_per_person_year = 1e6 * hall_fast_s / person_years,
             speed_up = bw_s / hall_fast_s)]
fwrite(bench, file.path(OUT_DIR, "benchmark.csv"))

# Copy the results used by the vignette into the package
EXT <- file.path(HF_DIR, "inst", "extdata", "validation")
dir.create(EXT, showWarnings = FALSE, recursive = TRUE)
invisible(file.copy(list.files(OUT_DIR, pattern = "\\.csv$", full.names = TRUE), EXT, overwrite = TRUE))

sink(file.path(OUT_DIR, "session_info.txt")); print(sessionInfo())
cat("\nCPU cores:", parallel::detectCores(), "\nOpenMP:", has_omp, "\n"); sink()

options(width = 250)
cat("\n==== ACCURACY ====\n");  print(acc, digits = 3)
cat("\n==== APPETITE FEEDBACK ====\n"); print(fbacc, digits = 3); print(fbcost)
cat("\n==== bw LOOK-AHEAD ====\n"); print(lookahead)
cat("\n==== bw K INITIALISATION (EI supplied 5% above RMR x PAL, no intake change) ====\n"); print(drift)
cat("\n==== bw GLYCOGEN STIFFNESS (stress population) ====\n"); print(stiff_tab)
cat("\n==== INTEGRATOR COST ====\n"); print(evals)
cat("\n==== TOLERANCE ====\n"); print(tol)
cat("\n==== SPEED ====\n"); print(bench, digits = 3)
