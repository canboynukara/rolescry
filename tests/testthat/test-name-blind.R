# Name-blindness contract (Ax1) and the capped name_bonus tie-breaker.

test_that("classify_value_type classifies name-traps by DATA, not name", {
  set.seed(2026L)
  # name suggests one type, data is another
  age_trap <- rep(c("Yes", "No"), 60)                 # name=>continuous, DATA=categorical
  code_trap <- round(stats::rnorm(120, 50, 12), 2)    # name=>categorical, DATA=continuous
  expect_match(classify_value_type(age_trap)$type, "BINARY|CATEGORICAL")
  expect_match(classify_value_type(code_trap)$type, "CONTINUOUS")
})

test_that("score_gap_ok enforces math-dominant, name-capped", {
  expect_true(score_gap_ok(0.95, 0.05))   # legit
  expect_false(score_gap_ok(0.80, 0.50))  # name-dominant -> rejected
  expect_false(score_gap_ok(0.89, 0.05))  # math below 0.90 -> rejected
})

test_that("name_bonus tie-breaks outcome selection (Phase A caveat exercised)", {
  tw <- make_namebonus_twins()
  nb <- rolescry_default_name_bonus()

  # Pure (NULL): the first 2-level column wins positionally.
  pure <- detect_roles(tw$named)
  expect_true(pure$roles$outcome_binary$found)

  # With the dictionary, "death" is selected over the positional default.
  hinted <- detect_roles(tw$named, name_bonus = nb)
  expect_identical(hinted$roles$outcome_binary$columns, "death")
  expect_false(identical(pure$roles$outcome_binary$columns,
                         hinted$roles$outcome_binary$columns))
})

test_that("name_bonus is inert on col_N columns (no names to match)", {
  tw <- make_namebonus_twins()
  nb <- rolescry_default_name_bonus()
  pure <- detect_roles(tw$col_n)
  hinted <- detect_roles(tw$col_n, name_bonus = nb)
  expect_identical(pure$roles$outcome_binary$columns, hinted$roles$outcome_binary$columns)
  expect_identical(pure$roles$group_var$columns, hinted$roles$group_var$columns)
})
