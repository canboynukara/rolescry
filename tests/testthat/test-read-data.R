test_that("read_data reads a CSV with header detection and type conversion", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("age,sex,score", "34,M,8.1", "51,F,7.4"), tmp)
  d <- read_data(tmp)
  expect_equal(names(d), c("age", "sex", "score"))
  expect_equal(nrow(d), 2L)
  expect_type(d$age, "integer")
  expect_type(d$score, "double")
})

test_that("read_data honors an explicit header row", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("export,note", "age,sex", "34,M", "51,F"), tmp)
  d <- read_data(tmp, header = 2)
  expect_equal(names(d), c("age", "sex"))
  expect_equal(nrow(d), 2L)
})

test_that("read_data errors on a missing file", {
  expect_error(read_data(file.path(tempdir(), "no_such_file_zzz.csv")), "not found")
})
