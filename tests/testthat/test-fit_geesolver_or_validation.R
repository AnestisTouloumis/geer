testthat::local_edition(3)


## Direct calls to the internal odds-ratio solver: the R interface already
## rejects responses outside [0, 1], so these tests guard the C++ routine
## against being fed such values by other callers.
or_solver_x <- c(
  -1.2, 0.3, 0.8, -0.5, 1.1, -0.2, 0.9, -1.4, 0.1, 0.4, -0.7, 1.3,
  -0.9, 0.6, -0.1, 1.5, -0.3, 0.2, -1.1, 0.7, -0.6, 0.0, 1.0, -0.8
)
or_solver_y <- c(
  0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1,
  0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0
)


call_or_solver <- function(y) {
  geer:::fit_geesolver_or(
    y, cbind(1, or_solver_x), rep(1:8, each = 3), rep(1:3, times = 8),
    rep(1, 24), "logit", c(0, 0), rep(0, 24), 100L, 1e-6, 10L, 1, 0.5,
    "gee", rep(1.5, 3)
  )
}


test_that("the odds-ratio solver fits binary and proportion responses", {
  binary <- call_or_solver(or_solver_y)
  expect_identical(binary$failure, "")
  expect_true(all(is.finite(binary$beta_hat)))
  iterations <- sum(binary$criterion != 0)
  expect_lte(binary$criterion[iterations], 1e-6)

  proportions <- or_solver_y
  proportions[c(3, 8)] <- c(0.3, 0.75)
  fractional <- call_or_solver(proportions)
  expect_identical(fractional$failure, "")
  expect_true(all(is.finite(fractional$beta_hat)))
})


test_that("the odds-ratio solver rejects responses outside [0, 1]", {
  bad_above <- or_solver_y
  bad_above[4] <- 2
  expect_error(
    call_or_solver(bad_above),
    "fit_geesolver_or: the response must be finite and lie in \\[0, 1\\]; observation 4 has value 2"
  )

  bad_below <- or_solver_y
  bad_below[11] <- -1
  expect_error(call_or_solver(bad_below), "observation 11 has value -1")

  for (bad in c(NaN, Inf, -Inf)) {
    bad_value <- or_solver_y
    bad_value[6] <- bad
    expect_error(call_or_solver(bad_value), "observation 6 has value")
  }
})
