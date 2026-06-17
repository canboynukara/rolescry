# KEYSTONE: the name-blindness (turnusol) invariant. A col_N dataset and its
# named twin carry byte-identical data, so every role must be assigned to the
# same columns BY POSITION regardless of the names.

test_that("col_N twin and named twin yield identical roles by index", {
  tw <- make_turnusol_twins()
  rc <- detect_roles(tw$col_n)
  rn <- detect_roles(tw$named)

  expect_identical(names(rc$roles), names(rn$roles))
  for (role in names(rc$roles)) {
    expect_identical(
      role_index(rc, role, tw$col_n),
      role_index(rn, role, tw$named),
      info = paste("role:", role)
    )
    # the found flag must also match
    expect_identical(rc$roles[[role]]$found, rn$roles[[role]]$found, info = role)
  }
})

test_that("the twin data is genuinely identical (only names differ)", {
  tw <- make_turnusol_twins()
  expect_equal(unname(as.matrix(tw$col_n)), unname(as.matrix(tw$named)))
  expect_false(identical(names(tw$col_n), names(tw$named)))
})

test_that("name_bonus is inert on col_N columns (turnusol preserved without names)", {
  # On col_N columns no keyword can match, so supplying a name dictionary must
  # leave every role assignment identical to the pure-mathematical result.
  # (On the NAMED twin a bonus may intentionally change selection -- that is the
  # whole point of the bonus -- so named-vs-col_N equality holds only for
  # name_bonus = NULL, which the first test already checks.)
  tw <- make_turnusol_twins()
  nb <- rolescry_default_name_bonus()
  pure <- detect_roles(tw$col_n)
  hint <- detect_roles(tw$col_n, name_bonus = nb)
  for (role in names(pure$roles)) {
    expect_identical(
      role_index(pure, role, tw$col_n),
      role_index(hint, role, tw$col_n),
      info = role
    )
  }
})

test_that("classify_value_type is name-blind by construction", {
  tw <- make_turnusol_twins()
  tc <- vapply(tw$col_n, function(x) classify_value_type(x)$type, character(1))
  tn <- vapply(tw$named, function(x) classify_value_type(x)$type, character(1))
  expect_equal(unname(tc), unname(tn))
})
