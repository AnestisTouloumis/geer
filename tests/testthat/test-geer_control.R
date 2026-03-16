test_that("geer_control returns a list with expected names and types", {
  ctrl <- geer_control()
  expect_type(ctrl, "list")
  expect_named(
    ctrl,
    c("tolerance", "maxiter", "or_adding", "step_maxit", "step_multi", "jeffreys_power")
  )
  expect_true(is.numeric(ctrl$tolerance) && length(ctrl$tolerance) == 1L)
  expect_true(is.numeric(ctrl$or_adding) && length(ctrl$or_adding) == 1L)
  expect_true(is.numeric(ctrl$jeffreys_power) && length(ctrl$jeffreys_power) == 1L)
  expect_true(is.integer(ctrl$maxiter) && length(ctrl$maxiter) == 1L)
  expect_true(is.integer(ctrl$step_maxit) && length(ctrl$step_maxit) == 1L)
  expect_true(is.integer(ctrl$step_multi) && length(ctrl$step_multi) == 1L)
})


test_that("geer_control accepts valid user inputs and coerces integer fields", {
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
  expect_equal(ctrl$maxiter, as.integer(123))
  expect_equal(ctrl$step_maxit, as.integer(7))
  expect_equal(ctrl$step_multi, as.integer(2))
})


test_that("geer_control errors on non-positive or invalid scalars", {
  expect_error(geer_control(tolerance = 0), "tolerance")
  expect_error(geer_control(tolerance = -1e-6), "tolerance")
  expect_error(geer_control(tolerance = NA_real_), "tolerance")
  expect_error(geer_control(or_adding = 0), "or_adding")
  expect_error(geer_control(or_adding = -0.5), "or_adding")
  expect_error(geer_control(or_adding = NA_real_), "or_adding")
  expect_error(geer_control(jeffreys_power = 0), "jeffreys_power")
  expect_error(geer_control(jeffreys_power = -0.1), "jeffreys_power")
  expect_error(geer_control(jeffreys_power = NA_real_), "jeffreys_power")
})


test_that("geer_control errors on invalid integer controls", {
  expect_error(geer_control(maxiter = 0), "maxiter")
  expect_error(geer_control(maxiter = -1), "maxiter")
  expect_error(geer_control(maxiter = 1.5), "maxiter")
  expect_error(geer_control(maxiter = NA_integer_), "maxiter")
  expect_error(geer_control(step_maxit = 0), "step_maxit")
  expect_error(geer_control(step_maxit = -2), "step_maxit")
  expect_error(geer_control(step_maxit = 2.2), "step_maxit")
  expect_error(geer_control(step_maxit = NA_integer_), "step_maxit")
  expect_error(geer_control(step_multi = 0), "step_multi")
  expect_error(geer_control(step_multi = -3), "step_multi")
  expect_error(geer_control(step_multi = 1.1), "step_multi")
  expect_error(geer_control(step_multi = NA_integer_), "step_multi")
})
