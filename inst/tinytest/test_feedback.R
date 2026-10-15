# Appetite feedback (Hall 2024): realised intake = EIchange - k (BW - BW0)
source("helpers.R")
p <- make_pop(40, "realistic", seed = 31)
d <- random_walk(40, 3, seed = 32)

# Against the reference solver, which evaluates intake from the full state
for (pr in c("baseline", "tirzepatide")) {
  f <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, appetite_k = pr)
  r <- hall_bw_reference(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, appetite_k = pr,
                         steps_per_day = 10)
  expect_true(maxabs(f, r) < 5e-4, info = paste("feedback preset", pr))
}

# Time-varying k (normal appetite control, then a drug from year 2)
km <- cbind(matrix(95, 40, 1), matrix(48, 40, 2))
f <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, appetite_k = km)
r <- hall_bw_reference(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, appetite_k = km,
                       steps_per_day = 10)
expect_true(maxabs(f, r) < 5e-4, info = "time-varying k")

# The feedback code path with k -> 0 agrees with the no-feedback path
f0 <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL)
fe <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, appetite_k = 1e-9)
expect_true(maxabs(f0, fe) < 3e-4)

# Presets resolve to their published values, per person too
expect_identical(hall_bw_annual(80, 1.75, 40, "male", -100, appetite_k = "baseline"),
                 hall_bw_annual(80, 1.75, 40, "male", -100, appetite_k = 95))
expect_identical(
  hall_bw_annual(c(80, 90), c(1.75, 1.7), c(40, 50), c("male", "female"), c(-100, -100),
                 appetite_k = c("baseline", "rygb")),
  hall_bw_annual(c(80, 90), c(1.75, 1.7), c(40, 50), c("male", "female"), c(-100, -100),
                 appetite_k = c(95, 58)))
expect_equal(hall_appetite_presets$k[hall_appetite_presets$preset == "calorie_restriction"], 82)

# Definition of the realised intake and a weight plateau
a <- hall_bw_annual(80, 1.75, 40, "male", matrix(-100, 1, 8), appetite_k = 95, output = "all")
expect_equal(a$realised_EIchange[1, 8], -100 - 95 * (a$BW[1, 8] - 80), tolerance = 1e-6, scale = 1)
expect_true(abs(a$BW[1, 8] - a$BW[1, 7]) < 1e-3)               # plateau within years
noFb <- hall_bw_annual(80, 1.75, 40, "male", matrix(-100, 1, 8))
expect_true(abs(a$BW[1, 8] - 80) < abs(noFb[1, 8] - 80) / 3)   # strongly attenuated

# Intake is truncated at zero for pushes larger than baseline intake
b <- hall_bw_annual(90, 1.75, 45, "male", matrix(-5000, 1, 2), appetite_k = "rygb", output = "all")
expect_true(all(b$baseline$EI + b$realised_EIchange >= -1e-9))
expect_true(all(is.finite(b$BW)))
