# One interval at a time: hall_state() + hall_step() (0.3.0).
# The defining property: Y calls of hall_step() give exactly (identical(), not
# a tolerance) the weights of one hall_bw_annual() call over the Y years.
source("helpers.R")
p <- make_pop(40, "stress", seed = 61)
Y <- 5
d <- random_walk(40, Y, seed = 62)

step_all <- function(s, D, ...) {
  out <- matrix(NA_real_, nrow(s), ncol(D))
  for (y in seq_len(ncol(D))) { s <- hall_step(s, D[, y], ...); out[, y] <- s$BW }
  list(state = s, BW = out)
}
strip <- function(m) { dimnames(m) <- NULL; m }

# ---- 1. Bit-identical to hall_bw_annual() ------------------------------------
s0 <- hall_state(p$age, p$female, bw = p$bw, ht = p$ht, PAL = p$PAL,
                 pcarb_base = p$pcarb_base)
w  <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL,
                     pcarb_base = p$pcarb_base, pcarb = p$pcarb)
st <- step_all(s0, d, pcarb = p$pcarb)
expect_identical(st$BW, strip(w))
expect_identical(st$state$intervals, rep(as.integer(Y), 40))

# ... with appetite feedback (a different solver path)
wf <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = p$PAL, appetite_k = "baseline")
sf <- step_all(hall_state(p$age, p$female, bw = p$bw, ht = p$ht, PAL = p$PAL), d,
               appetite_k = "baseline")
expect_identical(sf$BW, strip(wf))

# ... with sodium, a known baseline intake and fat mass, and the 'bw' model
na <- matrix(seq(-300, 300, length.out = 40 * Y), 40, Y)
EIs <- 2300 + 10 * seq_len(40); fats <- p$bw * 0.3
wx <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, NAchange = na, EI = EIs, fat = fats,
                     PAL = p$PAL)
sx <- hall_state(p$age, p$female, bw = p$bw, ht = p$ht, PAL = p$PAL, EI = EIs, fat = fats)
bwx <- matrix(NA_real_, 40, Y)
for (y in seq_len(Y)) { sx <- hall_step(sx, d[, y], NAchange = na[, y]); bwx[, y] <- sx$BW }
expect_identical(bwx, strip(wx))
wb <- hall_bw_annual(p$bw, p$ht, p$age, p$female, d, PAL = 1.5, model = "bw")
sb <- step_all(hall_state(p$age, p$female, bw = p$bw, ht = p$ht, PAL = 1.5, model = "bw"), d)
expect_identical(sb$BW, strip(wb))

# ... fractions of baseline intake, from BMI without heights (= hall_bmi_annual)
fr <- matrix(seq(-0.08, 0.04, length.out = 40 * Y), 40, Y)
bmi <- p$bw / p$ht^2
wm <- hall_bmi_annual(bmi, p$age, p$female, fr)
sm <- hall_state(p$age, p$female, bmi = bmi)
bm <- matrix(NA_real_, 40, Y)
for (y in seq_len(Y)) { sm <- hall_step(sm, fr[, y], EIchange_unit = "fraction"); bm[, y] <- sm$BMI }
expect_identical(bm, strip(wm))
expect_true(all(sm$height_imputed))

# ... and with real heights, from BMI = from body weight
sh <- hall_state(p$age, p$female, bmi = bmi, ht = p$ht)
expect_false(any(sh$height_imputed))
expect_equal(sh$bw0, p$bw, tolerance = 1e-14)

# ---- 2. Rows are independent: subset, reorder, bind --------------------------
s2 <- s0
for (y in 1:2) s2 <- hall_step(s2, d[, y], pcarb = p$pcarb)
i <- c(7, 3, 21)
a <- hall_step(s2[i, ], d[i, 3], pcarb = p$pcarb[i])
b <- hall_step(s2, d[, 3], pcarb = p$pcarb)[i, ]
expect_identical(a$BW, b$BW)
expect_identical(a$L, b$L)
expect_true(inherits(s2[i, ], "hall_state"))
# a newcomer bound to a running state starts from its own baseline
nb <- hall_state(50, "male", bw = 90, ht = 1.8)
both <- rbind(s2[1:3, ], nb)
expect_true(inherits(both, "hall_state"))
nxt <- hall_step(both, -100)
expect_identical(nxt$BW[4], unname(hall_bw_annual(90, 1.8, 50, "male", -100)[1, 1]))

# ---- 3. No change at baseline stays exactly at baseline ----------------------
z <- hall_step(s0, 0)
expect_identical(z$L, s0$L)
expect_identical(z$ECF, s0$ECF)
expect_identical(z$AT, s0$AT)
# ... and the realised change is exactly 0
expect_identical(z$dEI, rep(0, 40))

# ---- 4. Missing inputs: NA from that interval on, others untouched -----------
dn <- d[, 1]; dn[5] <- NA
s5 <- hall_step(s0, dn, pcarb = p$pcarb)
expect_true(is.na(s5$BW[5]) && is.na(s5$L[5]))
expect_identical(s5$BW[-5], strip(w)[-5, 1])
s6 <- hall_step(s5, d[, 2], pcarb = p$pcarb)
expect_true(is.na(s6$BW[5]))
expect_identical(s6$BW[-5], strip(w)[-5, 2])

# ---- 5. Units and checks -------------------------------------------------------
expect_error(hall_step(hall_state(40, "male", bmi = 30), -100), "Height is required")
expect_silent(hall_step(hall_state(40, "male", bmi = 30), -0.05, EIchange_unit = "fraction"))
expect_error(hall_step(hall_state(40, "male", bmi = 30), -5, EIchange_unit = "fraction"),
             "fraction")
expect_error(hall_state(40, "male"), "exactly one")
expect_error(hall_state(40, "male", bw = 80), "'ht' is required")
expect_error(hall_step(data.frame(a = 1), -100), "hall_state")
expect_error(hall_state(40, "male", bmi = 30, EI = 2500), "imputed heights")

# ---- 6. Interval length: two half-years ~ one year ----------------------------
one <- hall_step(s0, -150, days = 365)
two <- hall_step(hall_step(s0, -150, days = 182.5), -150, days = 182.5)
expect_true(maxabs(one$BW, two$BW) < 2e-4)

# ---- 7. Threads do not change results -----------------------------------------
if (hall_has_openmp()) {
  expect_identical(hall_step(s0, d[, 1], threads = 1)$BW, hall_step(s0, d[, 1], threads = 2)$BW)
}

# ---- 8. Controls: the identity in 1 must FAIL when the carried state is broken -
# Reset one carried compartment to its baseline value before every interval:
# each must break the identity. (The step-size hint h is carried as the core
# does, but resetting it was measured to change nothing here, so it is no
# control.)
for (cmp in c("L", "AT", "ECF", "G")) {
  bad <- s0; bwbad <- matrix(NA_real_, 40, Y)
  for (y in seq_len(Y)) {
    bad[[cmp]] <- s0[[cmp]]
    bad <- hall_step(bad, d[, y], pcarb = p$pcarb); bwbad[, y] <- bad$BW
  }
  expect_false(identical(bwbad, strip(w)), info = paste("control: reset", cmp))
}
# (`dEI`, the realised change, is an OUTPUT of each interval: the solver never
# reads it, so resetting it changes nothing and is no control either.)
# lean mass not carried is a large error, not a rounding one
bad <- s0; bwbad <- matrix(NA_real_, 40, Y)
for (y in seq_len(Y)) { bad$L <- s0$L; bad <- hall_step(bad, d[, y], pcarb = p$pcarb); bwbad[, y] <- bad$BW }
expect_true(maxabs(bwbad, strip(w)) > 0.5)
