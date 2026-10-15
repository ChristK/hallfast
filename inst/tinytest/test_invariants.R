# Model invariants and sanity checks
# No change in intake: weight stays exactly at baseline (appendix model)
z <- hall_bw_annual(c(60, 80, 120), c(1.6, 1.75, 1.8), c(25, 45, 70),
                    c("female", "male", "female"), matrix(0, 3, 5))
expect_equal(as.vector(z), rep(c(60, 80, 120), 5), tolerance = 1e-10, scale = 1)

# Default PAL is 1.6
expect_identical(hall_bw_annual(80, 1.75, 40, "male", -100),
                 hall_bw_annual(80, 1.75, 40, "male", -100, PAL = 1.6))

# Monotone in the size of the change
r <- hall_bw_annual(rep(85, 4), rep(1.75, 4), rep(45, 4), rep("male", 4),
                    matrix(c(-300, -100, 100, 300), 4, 3))
expect_true(all(diff(r[, 3]) > 0))

# Hall's rule of thumb: roughly half of the long-run change occurs in year 1
s <- hall_bw_annual(80, 1.75, 40, "male", matrix(-100, 1, 15))
frac <- (s[1, 1] - 80) / (s[1, 15] - 80)
expect_true(frac > 0.4 && frac < 0.75)

# Higher PAL gives a smaller long-run change (larger energy cost per kg)
lo <- hall_bw_annual(80, 1.75, 40, "male", matrix(-100, 1, 10), PAL = 1.5)
hi <- hall_bw_annual(80, 1.75, 40, "male", matrix(-100, 1, 10), PAL = 1.8)
expect_true(abs(hi[1, 10] - 80) < abs(lo[1, 10] - 80))

# Starting a change after some unchanged years equals starting it at year 1
late  <- hall_bw_annual(80, 1.75, 40, "male", matrix(c(0, 0, -100, -100), 1, 4))
early <- hall_bw_annual(80, 1.75, 40, "male", matrix(c(-100, -100), 1, 2))
expect_equal(late[1, 3:4], early[1, 1:2], tolerance = 1e-12, check.attributes = FALSE)
