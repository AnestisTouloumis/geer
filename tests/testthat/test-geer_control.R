testthat::local_edition(3)

test_that("normalize_geer_control fills defaults for partial control lists", {
  ctrl <- normalize_geer_control(list(maxiter = 25L, tolerance = 1e-5))
  expect_type(ctrl, "list")
  expect_named(
    ctrl,
    c(
      "tolerance", "maxiter", "or_adding",
      "step_maxiter", "step_multiplier", "jeffreys_power"
    )
  )
  expect_equal(ctrl$tolerance, 1e-5)
  expect_equal(ctrl$maxiter, 25L)
  expect_equal(ctrl$or_adding, 0.5)
  expect_equal(ctrl$step_maxiter, 10L)
  expect_equal(ctrl$step_multiplier, 1L)
  expect_equal(ctrl$jeffreys_power, 0.5)
})


test_that("normalize_geer_control validates malformed partial control lists", {
  expect_error(
    normalize_geer_control(list(maxiter = "bad", tolerance = -1)),
    "positive"
  )
})


test_that("geer_control returns the expected defaults", {
  ctrl <- geer_control()
  expect_type(ctrl, "list")
  expect_named(
    ctrl,
    c("tolerance", "maxiter", "or_adding", "step_maxiter", "step_multiplier", "jeffreys_power")
  )
  expect_equal(ctrl$tolerance, 1e-6)
  expect_equal(ctrl$maxiter, 500L)
  expect_equal(ctrl$or_adding, 0.5)
  expect_equal(ctrl$step_maxiter, 10L)
  expect_equal(ctrl$step_multiplier, 1L)
  expect_equal(ctrl$jeffreys_power, 0.5)
  expect_true(is.numeric(ctrl$tolerance) && length(ctrl$tolerance) == 1L)
  expect_true(is.numeric(ctrl$or_adding) && length(ctrl$or_adding) == 1L)
  expect_true(is.numeric(ctrl$jeffreys_power) && length(ctrl$jeffreys_power) == 1L)
  expect_true(is.integer(ctrl$maxiter) && length(ctrl$maxiter) == 1L)
  expect_true(is.integer(ctrl$step_maxiter) && length(ctrl$step_maxiter) == 1L)
  expect_true(is.integer(ctrl$step_multiplier) && length(ctrl$step_multiplier) == 1L)
})


test_that("geer_control accepts valid inputs and coerces integer fields", {
  ctrl <- geer_control(
    tolerance = 1e-4,
    maxiter = 123,
    or_adding = 0.25,
    step_maxiter = 7,
    step_multiplier = 2,
    jeffreys_power = 0.75
  )
  expect_equal(ctrl$tolerance, 1e-4)
  expect_equal(ctrl$or_adding, 0.25)
  expect_equal(ctrl$jeffreys_power, 0.75)
  expect_true(is.integer(ctrl$maxiter))
  expect_true(is.integer(ctrl$step_maxiter))
  expect_true(is.integer(ctrl$step_multiplier))
  expect_equal(ctrl$maxiter, 123L)
  expect_equal(ctrl$step_maxiter, 7L)
  expect_equal(ctrl$step_multiplier, 2L)
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
    geer_control(step_maxiter = 0),
    "'step_maxiter' must be a positive integer"
  )
  expect_error(
    geer_control(step_multiplier = NA_integer_),
    "'step_multiplier' must be a positive integer"
  )
})
