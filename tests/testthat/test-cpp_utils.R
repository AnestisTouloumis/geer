skip_on_cran()

test_that("asymmetric solves are invariant to matrix scale", {
  compute <- get_cached_cpp_utils()
  for (scale in c(1e-200, 1e-12, 1, 1e12, 1e200)) {
    x <- scale * test_data$cpp_asymmetric
    b <- rep(scale, 2L)
    observed <- as.numeric(compute(x, b, 0L))
    expect_equal(observed, c(10 / 19, 18 / 19), tolerance = 1e-12)
    residual <- max(abs(x %*% observed - b)) / scale
    expect_lt(residual, 1e-12)
    expect_error(compute(x, b, 2L))
  }
})

test_that("symmetrization preserves zero and extreme symmetric matrices", {
  compute <- get_cached_cpp_utils()
  for (scale in c(0, 1e-320, 1e-200, 1, 1e200, 1e308)) {
    x <- diag(scale, 2L)
    expect_identical(compute(x, numeric(2L), 2L), x)
  }
  expect_error(
    compute(test_data$cpp_asymmetric, numeric(2L), 2L),
    "non-trivial asymmetry"
  )
})

test_that("symmetry checks retain undefined covariance sentinels", {
  compute <- get_cached_cpp_utils()
  x <- matrix(NA_real_, 2L, 2L)
  expect_identical(compute(x, numeric(2L), 2L), x)
  x[1L, 1L] <- 1
  expect_error(compute(x, numeric(2L), 2L))
  expect_error(
    compute(x, numeric(2L), 0L), "the matrix must be finite"
  )
})


test_that("generic indefinite matrices retain the documented LU warning", {
  compute <- get_cached_cpp_utils()
  expect_warning(
    out <- compute(matrix(c(1, 2, 2, 1), 2L), rep(1, 2L), 0L),
    "a symmetric matrix that is not positive definite was solved by LU",
    fixed = TRUE
  )
  expect_equal(as.numeric(out), rep(1 / 3, 2L), tolerance = 1e-14)
})

test_that("cluster wrappers distinguish input errors from numerical failures", {
  compute <- get_cached_cpp_utils()
  input <- tryCatch(compute(diag(2L), rep(1, 2L), 3L), error = identity)
  expect_identical(
    conditionMessage(input), "invalid repeated index; cluster 2 (id = 42)"
  )
  expect_null(conditionCall(input))
  numerical <- tryCatch(compute(diag(2L), rep(1, 2L), 4L), error = identity)
  expect_identical(
    conditionMessage(numerical),
    "compute_cpp_utils failed: cluster 2 (id = 42): the matrix is singular"
  )
})
