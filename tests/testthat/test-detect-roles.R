test_that("detect_roles returns a structured role_detection object", {
  set.seed(3L)
  d <- data.frame(a = stats::rnorm(50), b = stats::rbinom(50, 1, 0.5))
  r <- detect_roles(d)
  expect_s3_class(r, "role_detection")
  expect_named(r, c("var_info", "roles", "value_types", "potential_pairs", "n_obs", "n_var"))
  expect_equal(r$n_obs, 50L)
  expect_equal(r$n_var, 2L)
  expect_s3_class(r$var_info, "data.frame")
  expect_identical(names(r$var_info), c("column", "type"))
})

test_that("summary and print methods work", {
  set.seed(3L)
  d <- data.frame(a = stats::rnorm(50), b = stats::rbinom(50, 1, 0.5))
  r <- detect_roles(d)
  s <- summary(r)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("role", "found", "columns", "pct") %in% names(s)))
  expect_output(print(r), "role_detection")
})

test_that("detect_roles validates its inputs", {
  expect_error(detect_roles(1:10), "data.frame")
  expect_error(detect_roles(data.frame()), "no columns")
  expect_error(detect_roles(data.frame(a = 1), name_bonus = "nope"), "name_bonus")
})
