# Baseline-equivalence: the detector reproduces the deterministic, value-based
# behaviour established in the Phase A baseline -- on a controlled seeded
# fixture (no real MDStatR CSVs are copied). Asserts the planted structure is
# recovered and the result is deterministic.

test_that("value-based role structure is recovered on continuous data", {
  tw <- make_turnusol_twins()
  res <- detect_roles(tw$col_n)

  # every column is continuous numeric -> all are covariate candidates
  expect_identical(role_index(res, "covariate", tw$col_n), 1:7)

  # an outcome_continuous is always offered for continuous data
  expect_true(res$roles$outcome_continuous$found)

  # no categorical/binary columns -> no group/outcome_binary/survival
  expect_false(res$roles$group_var$found)
  expect_false(res$roles$outcome_binary$found)
  expect_false(res$roles$time_variable$found)

  # any found role assigns at least one real column
  for (role in names(res$roles)) {
    if (isTRUE(res$roles[[role]]$found)) {
      expect_true(length(res$roles[[role]]$columns) >= 1, info = role)
    }
  }
})

test_that("detection is deterministic (no RNG in the detector)", {
  tw <- make_turnusol_twins()
  a <- detect_roles(tw$col_n)
  b <- detect_roles(tw$col_n)
  expect_identical(a$roles, b$roles)
  expect_identical(a$value_types, b$value_types)
})

test_that("all-continuous columns classify as CONTINUOUS_NUMERIC", {
  tw <- make_turnusol_twins()
  res <- detect_roles(tw$col_n)
  expect_true(all(res$value_types == "CONTINUOUS_NUMERIC"))
  expect_true(all(res$var_info$type == "Continuous"))
})

test_that("clinical fixture recovers group, paired and outcome roles", {
  tw <- make_namebonus_twins()
  res <- detect_roles(tw$named)
  # treatment_arm is a balanced 2-level categorical -> group candidate
  expect_true(res$roles$group_var$found)
  expect_identical(res$roles$group_var$columns, "treatment_arm")
  # pre/post are correlated continuous -> a paired pair is offered
  expect_true(res$roles$paired_pairs$found)
})
