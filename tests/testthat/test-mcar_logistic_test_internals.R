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


test_that("MCAR reconstruction reports an empty original response", {
  local_mocked_bindings(
    normalize_mcar_response = function(response, object) numeric(0),
    .package = "geer"
  )
  expect_error(
    build_mcar_transition_frame(fit_geewa_bin_exch, ~1),
    "mcar_logistic_test failed: the reconstructed original response is empty",
    fixed = TRUE
  )
})


test_that("MCAR reconstruction checks the risk-set response length", {
  calls <- 0L
  local_mocked_bindings(
    model.response = function(data, ...) {
      calls <<- calls + 1L
      rep(0, nrow(data) + as.integer(calls > 1L))
    },
    .package = "stats"
  )
  expect_error(
    build_mcar_transition_frame(fit_geewa_bin_exch, ~1),
    paste0(
      "mcar_logistic_test failed: the reconstructed response must have ",
      "one value per risk-set row"
    ),
    fixed = TRUE
  )
})
