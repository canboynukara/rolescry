test_that("detect_header finds row 1 and returns repaired names", {
  raw <- data.frame(
    X1 = c("patient_age", "34", "51", "29"),
    X2 = c("sex", "M", "F", "M"),
    X3 = c("lab_score", "8.1", "7.4", "9.0"),
    stringsAsFactors = FALSE
  )
  h <- detect_header(raw)
  expect_equal(h$header_row, 1L)
  expect_equal(h$names, c("patient_age", "sex", "lab_score"))
})

test_that("detect_header skips a sparse title row", {
  raw <- data.frame(
    X1 = c("Study Export", "age", "30", "40", "50"),
    X2 = c(NA, "weight", "70", "80", "90"),
    X3 = c(NA, "height", "170", "180", "175"),
    stringsAsFactors = FALSE
  )
  h <- detect_header(raw)
  expect_equal(h$header_row, 2L)
  expect_equal(h$names, c("age", "weight", "height"))
})

test_that("empty header cells fall back to col_N", {
  raw <- data.frame(
    X1 = c("id", "1", "2"),
    X2 = c("name", "a", "b"),
    X3 = c("", "x", "y"),
    X4 = c("grp", "P", "Q"),
    stringsAsFactors = FALSE
  )
  h <- detect_header(raw)
  expect_equal(h$header_row, 1L)
  expect_true("col_3" %in% h$names)
  expect_equal(anyDuplicated(h$names), 0L)
})
