# Name-blindness contract (Ax1) and the capped name_bonus tie-breaker.

test_that("classify_value_type classifies name-traps by DATA, not name", {
  set.seed(2026L)
  # name suggests one type, data is another
  age_trap <- rep(c("Yes", "No"), 60)                 # name=>continuous, DATA=categorical
  code_trap <- round(stats::rnorm(120, 50, 12), 2)    # name=>categorical, DATA=continuous
  expect_match(classify_value_type(age_trap)$type, "BINARY|CATEGORICAL")
  expect_match(classify_value_type(code_trap)$type, "CONTINUOUS")
})

test_that("outcome_binary is detected name-blind; name_bonus agrees with a decisive math pick", {
  tw <- make_namebonus_twins()
  nb <- rolescry_default_name_bonus()

  # The signature alone selects `death` as the binary outcome; because the math
  # pick is decisive, the capped name_bonus agrees rather than overrides it.
  pure   <- detect_roles(tw$named)
  hinted <- detect_roles(tw$named, name_bonus = nb)
  expect_true(pure$roles$outcome_binary$found)
  expect_identical(pure$roles$outcome_binary$columns, "death")
  expect_identical(hinted$roles$outcome_binary$columns, "death")

  # Turnusol: the name-stripped twin selects the same column by position.
  blind <- detect_roles(tw$col_n)
  expect_identical(match(pure$roles$outcome_binary$columns,  names(tw$named)),
                   match(blind$roles$outcome_binary$columns, names(tw$col_n)))
})

test_that("name_bonus is a capped tie-breaker when the math margin is within the cap", {
  # Two plausible groupings: perfectly balanced `site` is the mathematical pick;
  # a keyword dictionary nudges the choice to the intended `treated` arm
  # (a <= 10% score bump), exercising the capped name channel.
  set.seed(4L)
  clin <- data.frame(
    site    = rep(c("north", "south"), length.out = 160L),
    treated = sample(c("no", "yes"), 160L, replace = TRUE, prob = c(0.62, 0.38))
  )
  nb <- rolescry_default_name_bonus()
  expect_identical(detect_roles(clin)$roles$group_var$columns, "site")
  expect_identical(detect_roles(clin, name_bonus = nb)$roles$group_var$columns, "treated")
})

test_that("name_bonus is inert on col_N columns (no names to match)", {
  tw <- make_namebonus_twins()
  nb <- rolescry_default_name_bonus()
  pure <- detect_roles(tw$col_n)
  hinted <- detect_roles(tw$col_n, name_bonus = nb)
  expect_identical(pure$roles$outcome_binary$columns, hinted$roles$outcome_binary$columns)
  expect_identical(pure$roles$group_var$columns, hinted$roles$group_var$columns)
})
