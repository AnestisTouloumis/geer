test_that("MCAR coefficient tables reject invalid variances with context", {
  fit <- fit_geewa_bin_exch
  covariance <- diag(length(coef(fit)))
  for (value in c(-1, NA_real_, Inf)) {
    covariance[1L, 1L] <- value
    expect_error(
      build_mcar_coefficient_table(fit, covariance, character(0)),
      "^mcar_logistic_test failed:"
    )
  }
})

test_that("MCAR coefficient tables retain zero-variance behavior", {
  fit <- fit_geewa_bin_exch
  covariance <- diag(length(coef(fit)))
  covariance[1L, 1L] <- 0
  result <- build_mcar_coefficient_table(fit, covariance, character(0))
  expect_identical(result$std.error[1L], 0)
  expect_true(is.infinite(result$statistic[1L]))
  expect_identical(result$p.value[1L], 0)
})
