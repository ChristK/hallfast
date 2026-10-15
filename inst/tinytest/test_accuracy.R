# Accuracy of the fast solver against the independent reference solver
# (full four-compartment system, fixed-step RK4), which shares no algebra
# with the fast solver beyond the model equations.
source("helpers.R")

# Realistic population, intake changing every year
p <- make_pop(60, "realistic", seed = 11)
d <- random_walk(60, 4, seed = 12)
fast <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL)
ref  <- hall_bw_reference(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, steps_per_day = 10)
expect_true(maxabs(fast, ref) < 5e-4, info = "realistic, appendix model")

# Compatibility model against the reference with bw's equations
fastb <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, model = "bw")
refb  <- hall_bw_reference(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, model = "bw",
                           steps_per_day = 10)
expect_true(maxabs(fastb, refb) < 5e-4, info = "realistic, bw model")

# Stress: extreme BMI and PAL, large intake swings, sodium and carbohydrate
# changes, supplied baseline EI and fat mass
q <- make_pop(40, "stress", seed = 21)
set.seed(22)
dq  <- matrix(runif(40 * 3, -800, 500), 40, 3)
nq  <- matrix(runif(40 * 3, -1500, 1500), 40, 3)
EIq <- ifelse(seq_len(40) %% 2 == 0, runif(40, 1600, 3200), NA)
fq  <- ifelse(seq_len(40) %% 3 == 0, q$bw * runif(40, 0.15, 0.4), NA)
fs <- hall_bw_annual(q$bw, q$ht, q$age, q$female, dq, NAchange = nq, EI = EIq, fat = fq,
                     PAL = q$PAL, pcarb_base = q$pcarb_base, pcarb = q$pcarb, output = "all")
rs <- hall_bw_reference(q$bw, q$ht, q$age, q$female, dq, NAchange = nq, EI = EIq, fat = fq,
                        PAL = q$PAL, pcarb_base = q$pcarb_base, pcarb = q$pcarb,
                        steps_per_day = 10, output = "all")
expect_true(maxabs(fs$BW, rs$BW) < 2e-3, info = "stress, body weight")
expect_true(maxabs(fs$fat_mass + fs$lean_mass, rs$tissue) < 2e-3, info = "stress, fat + lean")

# Internal consistency of the full output
expect_equal(fs$BW, fs$fat_mass + fs$lean_mass + fs$ECF + 3.7 * fs$glycogen, tolerance = 1e-12)
expect_equal(fs$BMI, fs$BW / q$ht^2, tolerance = 1e-12)
expect_equal(fs$realised_EIchange, dq, tolerance = 1e-12, check.attributes = FALSE)

# Tighter tolerance reduces (or does not increase) the error
f5 <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, atol = 1e-6)
expect_true(maxabs(f5, ref) <= maxabs(fast, ref) + 1e-6)
