# Proportional changes in intake (EIchange_unit = "fraction") and the BMI
# interfaces hall_bmi_annual() / hall_bmi_long().
source("helpers.R")
p  <- make_pop(40, "realistic", seed = 51)
Y  <- 4
fr <- matrix(seq(-0.08, 0.04, length.out = 40 * Y), 40, Y)   # fractions of baseline EI
EI0 <- (9.99 * p$bw + 625 * p$ht - 4.92 * p$age + ifelse(p$female == 1, -161, 5)) * p$PAL

# ---- Fraction of baseline intake = the same change in kcal/day -------------
w_fr <- hall_bw_annual(p$bw, p$ht, p$age, p$female, fr, PAL = p$PAL,
                       EIchange_unit = "fraction")
w_kc <- hall_bw_annual(p$bw, p$ht, p$age, p$female, fr * EI0, PAL = p$PAL)
expect_equal(w_fr, w_kc, tolerance = 1e-12)

# ... also with a supplied baseline EI (the fraction is then of that EI)
EIs <- EI0 * 1.03
expect_equal(hall_bw_annual(p$bw, p$ht, p$age, p$female, fr, EI = EIs, PAL = p$PAL,
                            EIchange_unit = "fraction"),
             hall_bw_annual(p$bw, p$ht, p$age, p$female, fr * EIs, EI = EIs, PAL = p$PAL),
             tolerance = 1e-12)

# ... with appetite feedback (the push is then a fraction of baseline intake)
expect_equal(hall_bw_annual(p$bw, p$ht, p$age, p$female, fr, PAL = p$PAL,
                            appetite_k = "baseline", EIchange_unit = "fraction"),
             hall_bw_annual(p$bw, p$ht, p$age, p$female, fr * EI0, PAL = p$PAL,
                            appetite_k = "baseline"),
             tolerance = 1e-12)

# ... in the long interface (each person's EI0 applies to all of their rows)
id <- rep(seq_len(40), each = Y)
l_fr <- hall_bw_long(id, rep(p$bw, each = Y), rep(p$ht, each = Y), rep(p$age, each = Y),
                     rep(p$female, each = Y), as.vector(t(fr)), PAL = rep(p$PAL, each = Y),
                     EIchange_unit = "fraction")
expect_equal(l_fr, as.vector(t(w_fr)), tolerance = 1e-12)

# ... and in the reference solver, which agrees with the fast solver
r_fr <- hall_bw_reference(p$bw[1:5], p$ht[1:5], p$age[1:5], p$female[1:5], fr[1:5, 1:2],
                          PAL = p$PAL[1:5], EIchange_unit = "fraction", steps_per_day = 4)
expect_true(maxabs(r_fr, w_fr[1:5, 1:2]) < 1e-3)

# Fractions outside [-1, 1] are rejected as probable kcal/day inputs
expect_error(hall_bw_annual(70, 1.65, 30, "male", -100, EIchange_unit = "fraction"),
             "fraction")
expect_error(hall_bw_long(1, 70, 1.65, 30, "male", 2, EIchange_unit = "fraction"),
             "fraction")

# ---- BMI interface with known height = body-weight interface / ht^2 --------
bmi <- p$bw / p$ht^2
b_fr <- hall_bmi_annual(bmi, p$age, p$female, fr, ht = p$ht, PAL = p$PAL)
expect_equal(unname(b_fr), unname(w_fr / p$ht^2), tolerance = 1e-12)
dk <- random_walk(40, Y, seed = 52)
b_kc <- hall_bmi_annual(bmi, p$age, p$female, dk, ht = p$ht, EIchange_unit = "kcal",
                        PAL = p$PAL)
expect_equal(unname(b_kc),
             unname(hall_bw_annual(p$bw, p$ht, p$age, p$female, dk, PAL = p$PAL) / p$ht^2),
             tolerance = 1e-12)

# No change in intake leaves BMI exactly at baseline, with or without height
expect_equal(unname(hall_bmi_annual(bmi, p$age, p$female, matrix(0, 40, 2))),
             matrix(bmi, 40, 2), tolerance = 1e-12)

# Proportional change in BMI = proportional change in weight
expect_equal(unname(b_fr / bmi), unname(w_fr / p$bw), tolerance = 1e-12)

# ---- Missing height: reference height ---------------------------------------
ref <- ifelse(p$female == 1, hall_ref_height[["female"]], hall_ref_height[["male"]])
b_ref <- hall_bmi_annual(bmi, p$age, p$female, fr, PAL = p$PAL)
expect_equal(b_ref, hall_bmi_annual(bmi, p$age, p$female, fr, ht = ref, PAL = p$PAL),
             tolerance = 1e-12)
# NA heights are filled individually; supplied heights are kept
ht_na <- p$ht; ht_na[c(2, 5)] <- NA
b_na <- hall_bmi_annual(bmi, p$age, p$female, fr, ht = ht_na, PAL = p$PAL)
expect_equal(b_na[-c(2, 5), ], b_fr[-c(2, 5), ], tolerance = 1e-12)
expect_equal(b_na[c(2, 5), ], b_ref[c(2, 5), ], tolerance = 1e-12)
# User-supplied reference heights
expect_equal(hall_bmi_annual(25, 50, "male", -0.05, ref_height = c(female = 1.6, male = 1.8)),
             hall_bmi_annual(25, 50, "male", -0.05, ht = 1.8))

# The point of the proportional formulation: at a fixed BMI, the proportional
# change in BMI hardly depends on height, unlike a fixed change in kcal/day.
hs <- c(1.50, 1.62, 1.74)
rel <- function(unit, d) {
  B <- hall_bmi_annual(rep(27, 3), rep(45, 3), rep("female", 3), matrix(d, 3, 30),
                       ht = hs, EIchange_unit = unit)
  B[, 30] / 27 - 1
}
r_fr <- rel("fraction", -0.05); r_kc <- rel("kcal", -100)
expect_true(diff(range(r_fr)) / abs(mean(r_fr)) < 0.06)   # within ~5% across 24 cm
expect_true(diff(range(r_kc)) / abs(mean(r_kc)) > 0.15)   # ~22% for a fixed kcal change

# In a realistic population, the reference height gives an unbiased mean
# change in BMI and small individual errors with proportional changes.
q <- make_pop(400, "realistic", seed = 53)
qb <- q$bw / q$ht^2
qr <- ifelse(q$female == 1, hall_ref_height[["female"]], hall_ref_height[["male"]])
t5 <- hall_bmi_annual(qb, q$age, q$female, matrix(-0.05, 400, 5), ht = q$ht, PAL = q$PAL)
i5 <- hall_bmi_annual(qb, q$age, q$female, matrix(-0.05, 400, 5), ht = qr, PAL = q$PAL)
expect_true(abs(mean(i5[, 5] - qb) / mean(t5[, 5] - qb) - 1) < 0.01)  # mean bias < 1%
expect_true(sd((i5[, 5] - qb) / (t5[, 5] - qb)) < 0.03)               # individual SD < 3%

# ---- Errors -----------------------------------------------------------------
expect_error(hall_bmi_annual(27, 45, "female", -100, EIchange_unit = "kcal"), "Height is required")
expect_error(hall_bmi_annual(27, 45, "female", -100, ht = NA, EIchange_unit = "kcal"),
             "Height is required")
expect_error(hall_bmi_annual(27, 45, "female", -0.05, EI = 2000), "imputed")
expect_error(hall_bmi_annual(27, 45, "female", -0.05, fat = 25), "imputed")
expect_silent(hall_bmi_annual(27, 45, "female", -0.05, ht = 1.6, EI = 2000))
expect_error(hall_bmi_annual(27, 45, "female", -0.05, bw = 70), "cannot be passed")
expect_error(hall_bmi_annual(27, 45, "female", -0.05, output = "bw"), "should be one of")
expect_error(hall_bmi_annual(-27, 45, "female", -0.05), "positive")
expect_error(hall_bmi_annual(c(27, 30), 45, "female", -0.05), "same length")
expect_error(hall_bmi_annual(27, 45, "female", -0.05, ref_height = c(1.6, 1.7)), "ref_height")

# ---- Full output ------------------------------------------------------------
a <- hall_bmi_annual(bmi[1:3], p$age[1:3], p$female[1:3], fr[1:3, ], ht = c(p$ht[1], NA, NA),
                     PAL = p$PAL[1:3], output = "all")
expect_equal(a$height_imputed, c(FALSE, TRUE, TRUE))
expect_equal(a$height[1], p$ht[1])
expect_equal(a$baseline$BMI, bmi[1:3])
expect_equal(a$BMI, a$BW / a$height^2, tolerance = 1e-12)

# ---- Long interface = wide interface ----------------------------------------
lb <- hall_bmi_long(id, rep(bmi, each = Y), rep(p$age, each = Y), rep(p$female, each = Y),
                    as.vector(t(fr)), PAL = rep(p$PAL, each = Y))
expect_equal(lb, as.vector(t(b_ref)), tolerance = 1e-12)
lh <- hall_bmi_long(id, rep(bmi, each = Y), rep(p$age, each = Y), rep(p$female, each = Y),
                    as.vector(t(fr)), ht = rep(ht_na, each = Y), PAL = rep(p$PAL, each = Y),
                    output = "all")
expect_equal(lh$BMI, as.vector(t(b_na)), tolerance = 1e-12)
expect_equal(lh$height_imputed, rep(is.na(ht_na), each = Y))
expect_error(hall_bmi_long(c(1, 2, 1), rep(27, 3), rep(45, 3), "male", rep(-0.05, 3)),
             "contiguous")
expect_error(hall_bmi_long(1:2, c(27, 30), c(45, 50), c("male", "female"), c(-100, -100),
                           EIchange_unit = "kcal"), "Height is required")
