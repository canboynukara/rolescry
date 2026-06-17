test_that("classify_value_type covers the core types", {
  expect_identical(classify_value_type(character(0))$type, "EMPTY")
  expect_identical(classify_value_type(rep("x", 10))$type, "CONSTANT")
  expect_identical(classify_value_type(c("a", "b", "a", "b"))$type, "BINARY_CATEGORICAL")
  expect_identical(classify_value_type(c(0L, 1L, 0L, 1L))$type, "BINARY_NUMERIC_ENCODED")
  expect_identical(classify_value_type(rep(1:4, 30))$type, "ORDINAL_NUMERIC")
  expect_identical(classify_value_type(stats::rnorm(200))$type, "CONTINUOUS_NUMERIC")
})

test_that("non-finite numerics do not crash the integer check", {
  expect_identical(classify_value_type(c(1, 2, Inf, NaN, 3))$type, "CONTINUOUS_NUMERIC")
})
