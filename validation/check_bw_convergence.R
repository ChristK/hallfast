# =============================================================================
#  check_bw_convergence.R
#
#  Takes the 6 individuals of the stress population (first initialisation path,
#  n = 500) where bw (dt = 1) differs most from hall_fast(model = "bw"), and
#  re-runs bw with dt = 0.1 and dt = 0.01. If bw converges to hall_fast as
#  dt -> 0, the discrepancy is bw's daily-step discretisation error, not an
#  error in hall_fast. Population and inputs are generated exactly as in
#  validate_vs_bw.R (same seed).
#
#  Usage (from the root of the source repository):
#    Rscript validation/check_bw_convergence.R   (~2 min)
# =============================================================================
suppressPackageStartupMessages({ library(hallfast); library(bw); library(data.table) })
args   <- commandArgs(trailingOnly = FALSE)
this   <- sub("^--file=", "", args[grep("^--file=", args)])
HF_DIR <- if (length(this)) normalizePath(file.path(dirname(this), "..")) else getwd()

# Stress population, identical to validate_vs_bw.R
set.seed(7); n2 <- 2000; Y2 <- 5
female <- rbinom(n2, 1, 0.5)
ht <- ifelse(female == 1, rnorm(n2, 1.62, 0.065), rnorm(n2, 1.76, 0.07))
age <- runif(n2, 18, 95); bmi <- runif(n2, 16, 55); PAL <- runif(n2, 1.4, 2.4)
pcb <- runif(n2, 0.4, 0.6); pc <- runif(n2, 0.3, 0.7)
q <- data.table(female, age, ht, bw = bmi * ht^2, PAL, pcarb_base = pcb, pcarb = pc)
dEI2 <- matrix(runif(n2 * Y2, -800, 500), n2, Y2)
dNa2 <- matrix(runif(n2 * Y2, -1500, 1500), n2, Y2)

ii <- 1:500
fb <- hall_bw_annual(q$bw[ii], q$ht[ii], q$age[ii], q$female[ii], dEI2[ii, ], NAchange = dNa2[ii, ],
                     PAL = q$PAL[ii], pcarb_base = q$pcarb_base[ii], pcarb = q$pcarb[ii],
                     model = "bw", output = "all")
run_bw <- function(k, dt) {
  spd <- round(1 / dt); idx <- rep(1:Y2, each = 365 * spd)
  EId <- cbind(dEI2[k, idx, drop = FALSE], dEI2[k, Y2]); NAd <- cbind(dNa2[k, idx, drop = FALSE], dNa2[k, Y2])
  r <- suppressWarnings(adult_weight(q$bw[k], q$ht[k], q$age[k], ifelse(q$female[k] == 1, "female", "male"),
                                     EId, NAd, PAL = q$PAL[k], pcarb_base = q$pcarb_base[k],
                                     pcarb = q$pcarb[k], days = 365 * Y2, dt = dt))
  r$Body_Weight[, 365 * spd * (1:Y2) + 1, drop = FALSE]
}
b1 <- run_bw(ii, 1)
e  <- apply(abs(b1 - fb$BW), 1, max)
worst <- order(-e)[1:6]
EI0 <- (9.99 * q$bw + 625 * q$ht - 4.92 * q$age + ifelse(q$female == 1, -161, 5)) * q$PAL
cat("Worst 6 individuals (bw dt = 1 vs hall_fast(model='bw')):\n")
print(data.table(id = worst, BMI0 = round(bmi[worst], 1), PAL = round(PAL[worst], 2),
                 EI0_kcal = round(EI0[worst]),
                 glycogen_rate_per_day = round(4 * pcb[worst] * EI0[worst] / 4206.501, 2),
                 max_err_kg = round(e[worst], 3)))
conv <- rbindlist(lapply(c(1, 0.1, 0.01), function(dt)
  data.table(bw_dt_days = dt, max_abs_diff_kg = max(abs(run_bw(worst, dt) - fb$BW[worst, ])))))
cat("\nConvergence of bw towards hall_fast(model='bw') on these individuals:\n"); print(conv)
fwrite(conv, file.path(HF_DIR, "validation", "results", "bw_convergence_worst_cases.csv"))
fwrite(conv, file.path(HF_DIR, "inst", "extdata", "validation", "bw_convergence_worst_cases.csv"))
