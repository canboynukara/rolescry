test_that(".pair_admitted reproduces the pair detector's admission", {
  set.seed(11L); n <- 200L; t <- rnorm(n); a <- t + rnorm(n, 0, .5); b <- t + 1 + rnorm(n, 0, .5)
  p <- rolescry:::.pair_admitted(a, b)
  expect_true(p$admitted); expect_equal(p$route, "closed-form"); expect_gt(p$bits, 0)
  q <- rolescry:::.pair_admitted(rnorm(n), rnorm(n))
  expect_false(q$admitted); expect_equal(q$route, "none")
})
test_that(".survival_assoc reproduces the survival association routes", {
  set.seed(12L); n <- 300L; x <- rnorm(n); t <- rweibull(n, 1.3, exp(5 - .5 * x)); c <- runif(n, 0, quantile(t, .9))
  time <- pmin(t, c); ev <- as.integer(t <= c)
  s <- rolescry:::.survival_assoc(time, ev)
  expect_true(s$admitted); expect_true(s$monotone || s$nonmono)
  s0 <- rolescry:::.survival_assoc(runif(n, 1, 100), rbinom(n, 1, .5))
  expect_false(s0$admitted)
  expect_false(rolescry:::.survival_assoc(time[1:8], ev[1:8])$admitted)
})
test_that(".survival_assoc n = 10 boundary: default min_n skips, min_n = 9 evaluates", {
  set.seed(13L); tt <- rexp(10) + .1; ev <- rep(0:1, 5)
  expect_false(rolescry:::.survival_assoc(tt, ev)$admitted)
  expect_true(is.na(rolescry:::.survival_assoc(tt, ev)$bits))
  expect_true(is.finite(rolescry:::.survival_assoc(tt, ev, min_n = 9L)$bits))
})
test_that(".id_like and .scale_item_type follow the detector rules", {
  set.seed(14L)
  expect_true(rolescry:::.id_like(1:50)); expect_false(rolescry:::.id_like(rep(1:5, 10)))
  expect_false(rolescry:::.id_like(c(1:20, rep(NA, 10))))      # > 20 % NA
  expect_equal(rolescry:::.scale_item_type(sample(1:5, 100, TRUE)), "ordinal")
  expect_equal(rolescry:::.scale_item_type(rbinom(100, 1, .5)), "binary")
  expect_true(is.na(rolescry:::.scale_item_type(rnorm(100))))
  expect_true(is.na(rolescry:::.scale_item_type(sample(1:5, 15, TRUE))))   # < 20 observed
})
