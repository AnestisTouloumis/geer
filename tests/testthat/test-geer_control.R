testthat::local_edition(3)

test_that("geer_control returns the expected defaults", {
  ctrl <- geer_control()

  expect_type(ctrl, "list")
  expect_named(
    ctrl,
    c("tolerance", "maxiter", "or_adding", "step_maxit", "step_multi", "jeffreys_power")
  )

  expect_equal(ctrl$tolerance, 1e-6)
  expect_equal(ctrl$maxiter, 500L)
  expect_equal(ctrl$or_adding, 0.5)
  expect_equal(ctrl$step_maxit, 10L)
  expect_equal(ctrl$step_multi, 1L)
  expect_equal(ctrl$jeffreys_power, 0.5)

  expect_true(is.numeric(ctrl$tolerance) && length(ctrl$tolerance) == 1L)
  expect_true(is.numeric(ctrl$or_adding) && length(ctrl$or_adding) == 1L)
  expect_true(is.numeric(ctrl$jeffreys_power) && length(ctrl$jeffreys_power) == 1L)
  expect_true(is.integer(ctrl$maxiter) && length(ctrl$maxiter) == 1L)
  expect_true(is.integer(ctrl$step_maxit) && length(ctrl$step_maxit) == 1L)
  expect_true(is.integer(ctrl$step_multi) && length(ctrl$step_multi) == 1L)
})

test_that("geer_control accepts valid inputs and coerces integer fields", {
  ctrl <- geer_control(
    tolerance = 1e-4,
    maxiter = 123,
    or_adding = 0.25,
    step_maxit = 7,
    step_multi = 2,
    jeffreys_power = 0.75
  )

  expect_equal(ctrl$tolerance, 1e-4)
  expect_equal(ctrl$or_adding, 0.25)
  expect_equal(ctrl$jeffreys_power, 0.75)

  expect_true(is.integer(ctrl$maxiter))
  expect_true(is.integer(ctrl$step_maxit))
  expect_true(is.integer(ctrl$step_multi))

  expect_equal(ctrl$maxiter, 123L)
  expect_equal(ctrl$step_maxit, 7L)
  expect_equal(ctrl$step_multi, 2L)
})

test_that("geer_control validates positive numeric controls", {
  expect_error(
    geer_control(tolerance = 0),
    "'tolerance' must be a positive number"
  )
  expect_error(
    geer_control(or_adding = -0.5),
    "'or_adding' must be a positive number"
  )
  expect_error(
    geer_control(jeffreys_power = NA_real_),
    "'jeffreys_power' must be a positive number"
  )
})

test_that("geer_control validates positive integer controls", {
  expect_error(
    geer_control(maxiter = 1.5),
    "'maxiter' must be a positive integer"
  )
  expect_error(
    geer_control(step_maxit = 0),
    "'step_maxit' must be a positive integer"
  )
  expect_error(
    geer_control(step_multi = NA_integer_),
    "'step_multi' must be a positive integer"
  )
})
