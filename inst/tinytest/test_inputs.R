# Argument handling, recycling, errors and missing values
x <- hall_bw_annual(c(70, 90), c(1.65, 1.8), c(30, 50), c("female", "male"), c(-100, -100))
for (s in list(c("Female", "MALE"), c("f", "m"), c("women", "men"), c(TRUE, FALSE), c(1, 0),
               factor(c("female", "male"))))
  expect_identical(hall_bw_annual(c(70, 90), c(1.65, 1.8), c(30, 50), s, c(-100, -100)), x)

expect_error(hall_bw_annual(70, 1.65, 30, "x", -100), "sex")
expect_error(hall_bw_annual(70, 1.65, 30, 2, -100), "sex")
expect_error(hall_bw_annual(c(70, 80), 1.65, c(30, 40), c("male", "male"), -100), "same length")
expect_error(hall_bw_annual(70, 1.65, 30, "male", matrix(-100, 2, 1)), "one row per")
expect_error(hall_bw_annual(-70, 1.65, 30, "male", -100), "positive")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, PAL = 0), "PAL")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, pcarb_base = 0), "pcarb_base")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, pcarb = 1.5), "pcarb")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, appetite_k = -1), "appetite_k")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, appetite_k = "nope"), "preset")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -5000), "negative")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, days = 0), "days")
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, atol = -1), "atol")
expect_error(hall_bw_annual(70, 1.65, 30, "male", matrix(-100, 1, 2),
                            NAchange = matrix(0, 1, 3)), "NAchange")
expect_warning(hall_bw_annual(70, 1.65, 16, "male", -100), "children")
expect_warning(hall_bw_annual(70, 165, 30, "male", -100), "metres")

# Missing values: a missing baseline gives NA for that person only; a missing
# change gives NA from that year onwards
m <- hall_bw_annual(c(70, NA), c(1.65, 1.8), c(30, 50), c("female", "male"),
                    matrix(-100, 2, 3))
expect_true(all(is.na(m[2, ])) && all(!is.na(m[1, ])))
m2 <- hall_bw_annual(70, 1.65, 30, "female", matrix(c(-100, NA, -100), 1, 3))
expect_true(!is.na(m2[1, 1]) && all(is.na(m2[1, 2:3])))

# Names and dimensions
nm <- hall_bw_annual(c(a = 70, b = 90), c(1.65, 1.8), c(30, 50), c("female", "male"),
                     matrix(-100, 2, 2))
expect_identical(dimnames(nm), list(c("a", "b"), c("year1", "year2")))
