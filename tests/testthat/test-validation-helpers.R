testthat::local_edition(3)

test_that("check_geer_object rejects non-geer inputs", {
  expect_error(
    check_geer_object(list()),
    "'object' must be of 'geer' class"
  )
})

test_that("validate_step_thresholds returns validated thresholds", {
  out <- validate_step_thresholds(0.10, 0.20)

  expect_type(out, "list")
  expect_identical(out$p_enter, 0.10)
  expect_identical(out$p_remove, 0.20)
})

test_that("validate_step_thresholds rejects invalid p_enter values", {
  expect_error(
    validate_step_thresholds(0, 0.20),
    "'p_enter' must be strictly between 0 and 1"
  )
  expect_error(
    validate_step_thresholds(c(0.10, 0.20), 0.20),
    "'p_enter' must be a single finite numeric value"
  )
  expect_error(
    validate_step_thresholds("0.10", 0.20),
    "'p_enter' must be a single finite numeric value"
  )
})

test_that("validate_step_thresholds rejects invalid p_remove values", {
  expect_error(
    validate_step_thresholds(0.10, 1.10),
    "'p_remove' must be strictly between 0 and 1"
  )
  expect_error(
    validate_step_thresholds(0.10, c(0.20, 0.30)),
    "'p_remove' must be a single finite numeric value"
  )
  expect_error(
    validate_step_thresholds(0.10, "0.20"),
    "'p_remove' must be a single finite numeric value"
  )
})

test_that("validate_step_count returns a non-negative integer count", {
  expect_identical(validate_step_count(10), 10L)
  expect_identical(validate_step_count(0), 0L)
})

test_that("validate_step_count rejects invalid step values", {
  expect_error(
    validate_step_count(-1),
    "'steps' must be a single non-negative integer"
  )
  expect_error(
    validate_step_count(1.5),
    "'steps' must be a single non-negative integer"
  )
  expect_error(
    validate_step_count(c(1, 2)),
    "'steps' must be a single finite numeric value"
  )
})

test_that("normalize_test_options drops pmethod for non-working tests", {
  out <- normalize_test_options("wald", "robust", "rao-scott")

  expect_type(out, "list")
  expect_identical(out$test, "wald")
  expect_identical(out$cov_type, "robust")
  expect_null(out$pmethod)
})

test_that("normalize_test_options keeps pmethod for working tests", {
  out <- normalize_test_options("working-score", "robust", "satterthwaite")

  expect_type(out, "list")
  expect_identical(out$test, "working-score")
  expect_identical(out$cov_type, "robust")
  expect_identical(out$pmethod, "satterthwaite")
})
