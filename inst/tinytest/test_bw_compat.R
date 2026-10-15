# model = "bw" reproduces the R package 'bw' (v1.0.0, adult_weight(), dt = 1 day).
# Expected values are bw's body weight at day 365 for single-year runs
# (computed with bw 1.0.0; bw's own daily-step error is about 1e-4 kg).
bwv <- function(...) unname(hall_bw_annual(..., model = "bw")[1, 1])

expect_equal(bwv(80, 1.80, 40, "female", -100, PAL = 1.5), 77.284985, tolerance = 1e-3, scale = 1)
expect_equal(bwv(100, 1.75, 55, "male", 200, PAL = 1.7), 104.970861, tolerance = 1e-3, scale = 1)
expect_equal(bwv(65, 1.60, 30, "female", -300, NAchange = -500, PAL = 1.5, pcarb = 0.4),
             56.273458, tolerance = 1e-3, scale = 1)
expect_equal(bwv(90, 1.80, 45, "male", -150, fat = 25, PAL = 1.5), 85.946972,
             tolerance = 1e-3, scale = 1)
# bw sets K from RMR x PAL even when baseline EI is supplied, so a supplied EI
# that differs from RMR x PAL causes drift with no change in intake:
expect_equal(bwv(75, 1.65, 50, "female", 0, EI = 2300, PAL = 1.5), 83.218941,
             tolerance = 1e-3, scale = 1)
# ... whereas the appendix model (K from EE(0) = EI(0)) does not drift
expect_identical(unname(hall_bw_annual(75, 1.65, 50, "female", 0, EI = 2300, PAL = 1.5)[1, 1]), 75)
