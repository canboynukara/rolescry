test_that("compute_nmi is 0 for independence and ~1 for a deterministic map", {
  set.seed(1L)
  g <- sample(c("A", "B", "C"), 400, replace = TRUE)
  expect_lt(compute_nmi(g, sample(g)), 0.05)            # independent
  expect_gt(compute_nmi(g, paste0("y", g)), 0.99)       # deterministic 1:1
})

test_that("compute_nmi accepts a contingency table or two vectors", {
  x <- c("a", "a", "b", "b"); y <- c(1, 1, 2, 2)
  expect_equal(compute_nmi(x, y), compute_nmi(table(x, y)))
})

test_that("compute_nmi is bounded and degenerate-safe", {
  v <- compute_nmi(sample(letters[1:3], 100, TRUE), sample(1:3, 100, TRUE))
  expect_gte(v, 0); expect_lte(v, 1)
  expect_equal(compute_nmi(matrix(5, 1, 3)), 0)         # < 2 rows
  expect_equal(compute_nmi(rep("a", 10), rep("b", 10)), 0)  # constant marginals
})
