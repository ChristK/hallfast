# Long vs wide interfaces, interval length and threads
source("helpers.R")
p <- make_pop(30, "realistic", seed = 41)
d <- random_walk(30, 4, seed = 42)
na <- matrix(seq(-300, 300, length.out = 30 * 4), 30, 4)

w <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, NAchange = na, PAL = p$PAL,
                    appetite_k = "baseline")
long <- data.frame(id = rep(seq_len(30), each = 4),
                   bw = rep(p$bw, each = 4), ht = rep(p$ht, each = 4),
                   age = rep(p$age, each = 4), sex = rep(p$female, each = 4),
                   PAL = rep(p$PAL, each = 4), dEI = as.vector(t(d)), dNa = as.vector(t(na)))
l <- with(long, hall_bw_long(id, bw, ht, age, sex, dEI, NAchange = dNa, PAL = PAL,
                             appetite_k = "baseline"))
expect_equal(l, as.vector(t(w)), tolerance = 1e-12)

# Rows of an id must be contiguous
expect_error(hall_bw_long(c(1, 2, 1), rep(70, 3), rep(1.7, 3), rep(40, 3), rep("male", 3),
                          rep(-100, 3)), "contiguous")

# A constant change gives the same end-of-year weights with quarterly intervals
a <- hall_bw_annual(80, 1.75, 40, "male", matrix(-150, 1, 2))
q <- hall_bw_annual(80, 1.75, 40, "male", matrix(-150, 1, 8), days = 365 / 4)
expect_true(maxabs(a[1, ], q[1, c(4, 8)]) < 2e-4)

# Results do not depend on the number of threads
if (hall_has_openmp()) {
  w1 <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, threads = 1)
  w2 <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, threads = 2)
  expect_identical(w1, w2)
}
