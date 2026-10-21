test_that("one-factor bits agree in sign with the pair criterion at k = 2", {
  for (sd_e in c(.7, 1.5)) {   # r ~ 0.67 and r ~ 0.3
    set.seed(61L); n <- 300L; t <- rnorm(n); a <- t + rnorm(n, 0, sd_e); b <- t + rnorm(n, 0, sd_e)
    f <- rolescry:::.mdl_factor_bits(cbind(a, b)); p <- rolescry:::.mdl_pair_bits(a, b)
    expect_true(f$bits > 0); expect_true(p$bits - 0.5 * log2(n) > 0)
    expect_equal(sign(f$bits), sign(p$bits - rolescry:::.mdl_param_bits(1, n)))
  }
  set.seed(61L); n <- 300L; t <- rnorm(n); a <- t + rnorm(n, 0, .7); b <- t + rnorm(n, 0, .7)
  expect_equal(rolescry:::.mdl_factor_bits(cbind(a, b))$same_quantity, rolescry:::.mdl_pair_bits(a, b)$agree)
})
test_that("a 3-item battery with inter-item r ~ 0.3 is recalled at n = 300", {
  set.seed(65L); n <- 300L; t <- rnorm(n)
  M <- sapply(1:3, function(i) t + rnorm(n, 0, sqrt(1 / .3 - 1) * 1))   # r = 1/(1+s^2) = .3
  expect_true(rolescry:::.mdl_factor_bits(M)$bits > 0)
})
test_that("independent columns do not compress; a one-factor block does", {
  set.seed(62L); n <- 300L
  ind <- replicate(4, rnorm(n)); lat <- rnorm(n); blk <- sapply(1:4, function(i) lat + rnorm(n, 0, .8))
  expect_true(rolescry:::.mdl_factor_bits(ind)$bits <= 0)
  expect_true(rolescry:::.mdl_factor_bits(blk)$bits > 0)
  expect_true(rolescry:::.mdl_factor_bits(blk)$same_quantity)
})
test_that("a two-factor block compresses less under one factor than a one-factor block of equal size", {
  set.seed(63L); n <- 300L; l1 <- rnorm(n); l2 <- rnorm(n)
  two <- cbind(l1 + rnorm(n,0,.5), l1 + rnorm(n,0,.5), l2 + rnorm(n,0,.5), l2 + rnorm(n,0,.5))
  one <- sapply(1:4, function(i) l1 + rnorm(n, 0, .5))
  expect_lt(rolescry:::.mdl_factor_bits(two)$bits, rolescry:::.mdl_factor_bits(one)$bits)
})
test_that("collinear columns give a finite result", {
  a <- rnorm(60); M <- cbind(a, 2 * a, a + 1)
  f <- rolescry:::.mdl_factor_bits(M)
  expect_true(is.finite(f$bits) || identical(f$bits, -Inf))
})
test_that("wide-scale blocks do not error", {
  set.seed(64L); t <- rnorm(80); e <- rnorm(80)
  M <- cbind(1e-5 * (t + e), 1e5 * (t + e), 1e-5 * e)
  expect_error(f <- rolescry:::.mdl_factor_bits(M), NA)
  expect_true(is.finite(f$bits) || identical(f$bits, -Inf))
})
