# Shared helpers for the tinytest suite (sourced by each test file).
maxabs <- function(a, b) max(abs(a - b))

# Small reproducible adult populations
make_pop <- function(n, type = c("realistic", "stress"), seed = 1L) {
  type <- match.arg(type)
  set.seed(seed)
  female <- rbinom(n, 1, 0.5)
  ht <- ifelse(female == 1, rnorm(n, 1.62, 0.065), rnorm(n, 1.76, 0.07))
  if (type == "realistic") {
    age <- runif(n, 20, 90)
    bmi <- pmin(pmax(exp(rnorm(n, log(27), 0.18)), 16), 60)
    PAL <- runif(n, 1.4, 1.9); pcb <- rep(0.5, n); pc <- pcb
  } else {
    age <- runif(n, 18, 95); bmi <- runif(n, 16, 55); PAL <- runif(n, 1.4, 2.4)
    pcb <- runif(n, 0.4, 0.6); pc <- runif(n, 0.3, 0.7)
  }
  data.frame(female = female, age = age, ht = ht, bw = bmi * ht^2, PAL = PAL,
             pcarb_base = pcb, pcarb = pc)
}
random_walk <- function(n, Y, seed = 2L) {
  set.seed(seed)
  d <- matrix(0, n, Y); d[, 1] <- rnorm(n, -50, 50)
  if (Y > 1) for (y in 2:Y) d[, y] <- d[, y - 1] + rnorm(n, 0, 30)
  pmin(pmax(d, -400), 300)
}
