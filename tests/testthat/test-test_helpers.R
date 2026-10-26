test_that("published comparisons enforce each rounding interval", {
  published <- c(-2.1898, 1.8602, 0.0317, 0.3478, 0.9513)
  actual <- published
  actual[1L] <- actual[1L] + 1e-4
  result <- tryCatch(
    expect_published(actual, published, 4L),
    expectation_failure = identity
  )
  expect_s3_class(result, "expectation_failure")
})

test_that("published comparisons validate their arguments", {
  expect_error(expect_published(NA_real_, 1, 4L), "'object' must")
  expect_error(expect_published(1, NA_real_, 4L), "'published' must")
  expect_error(expect_published(1:2, 1, 4L), "must have matching shapes")
  expect_error(expect_published(1, 1, -1L), "'digits' must")
  expect_error(expect_published(1, 1, 1000L), "positive rounding unit")
})

test_that("subprocess evaluation uses a distinct R process", {
  result <- in_subprocess("cat(Sys.getpid())")
  expect_identical(result$status, 0L)
  expect_false(identical(result$output, as.character(Sys.getpid())))
  expect_error(in_subprocess(1), "'code' must be a character string")
})
