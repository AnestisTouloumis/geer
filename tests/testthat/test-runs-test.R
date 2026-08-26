testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch


count_runs <- function(signs) {
  signs <- signs[signs != 0]
  1L + sum(signs[-1L] != signs[-length(signs)])
}


test_that("runs-test calculation matches the Chang sign-sequence example", {
  signs <- c(1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1)
  out <- geer:::compute_runs_statistics(signs)

  expected_runs <- 71 / 11
  variance_runs <- 294 / 121
  expected_z <- (6 - expected_runs) / sqrt(variance_runs)

  expect_identical(out$runs, 6L)
  expect_identical(out$positive, 5L)
  expect_identical(out$negative, 6L)
  expect_identical(out$zero, 0L)
  expect_equal(out$expected_runs, expected_runs, tolerance = 1e-12)
  expect_equal(out$variance_runs, variance_runs, tolerance = 1e-12)
  expect_equal(out$statistic, expected_z, tolerance = 1e-12)
  expect_equal(
    out$p_value,
    2 * stats::pnorm(-abs(expected_z)),
    tolerance = 1e-12
  )
})


test_that("runs_test uses natural cluster/repeated ordering by default", {
  out <- runs_test(count_fit, type = "pearson")
  residual_values <- residuals(count_fit, type = "pearson")
  ord <- order(
    count_fit$id,
    count_fit$repeated,
    seq_along(residual_values)
  )
  signs <- sign(residual_values[ord])

  expect_s3_class(out, "htest")
  expect_identical(out$alternative, "two.sided")
  expect_identical(out$residual_type, "pearson")
  expect_identical(out$order_by, "natural cluster/repeated order")
  expect_identical(out$runs, as.integer(count_runs(signs)))
  expect_identical(out$positive, as.integer(sum(signs > 0)))
  expect_identical(out$negative, as.integer(sum(signs < 0)))
  expect_identical(out$zero, as.integer(sum(signs == 0)))
  expect_equal(
    unname(out$statistic),
    (out$runs - out$expected_runs) / sqrt(out$variance_runs),
    tolerance = 1e-12
  )
  expect_equal(
    out$p.value,
    2 * stats::pnorm(-abs(unname(out$statistic))),
    tolerance = 1e-12
  )
})


test_that("runs_test can order residuals by fitted values", {
  out <- runs_test(count_fit, type = "deviance", order_by = "fitted")
  residual_values <- residuals(count_fit, type = "deviance")
  natural <- order(count_fit$id, count_fit$repeated, seq_along(residual_values))
  ord <- natural[order(count_fit$fitted.values[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "fitted values")
  expect_identical(out$runs, as.integer(count_runs(signs)))
  expect_identical(out$positive, as.integer(sum(signs > 0)))
  expect_identical(out$negative, as.integer(sum(signs < 0)))
})


test_that("runs_test can order residuals by a model-matrix covariate", {
  out <- runs_test(count_fit, order_by = "lnage")
  residual_values <- residuals(count_fit, type = "pearson")
  natural <- order(count_fit$id, count_fit$repeated, seq_along(residual_values))
  key <- count_fit$x[, "lnage"]
  ord <- natural[order(key[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "model-matrix column 'lnage'")
  expect_identical(out$runs, as.integer(count_runs(signs)))
})


test_that("runs_test accepts a supplied ordering vector", {
  key <- rep(c(2, 1, 3), length.out = count_fit$obs_no)
  out <- runs_test(count_fit, type = "working", order_by = key)
  residual_values <- residuals(count_fit, type = "working")
  natural <- order(count_fit$id, count_fit$repeated, seq_along(residual_values))
  ord <- natural[order(key[natural], seq_along(natural))]
  signs <- sign(residual_values[ord])

  expect_identical(out$order_by, "supplied ordering vector")
  expect_identical(out$runs, as.integer(count_runs(signs)))
})


test_that("zero residuals are omitted from the sign sequence", {
  fit <- count_fit
  fit$residuals <- rep(c(1, 0, -1, 0), length.out = fit$obs_no)
  out <- runs_test(fit, type = "working")

  expect_identical(out$zero, as.integer(sum(fit$residuals == 0)))
  expect_identical(out$nonzero, fit$obs_no - out$zero)
  expect_identical(out$positive + out$negative, out$nonzero)
})


test_that("runs_test warns when the normal approximation is based on small sign counts", {
  fit <- count_fit
  fit$residuals <- c(
    rep(1, 10),
    rep(-1, 10),
    rep(0, fit$obs_no - 20)
  )

  expect_warning(
    out <- runs_test(fit, type = "working"),
    "normal approximation may be unreliable"
  )
  expect_identical(out$positive, 10L)
  expect_identical(out$negative, 10L)
})


test_that("runs_test validates ordering and residual sign requirements", {
  expect_error(
    runs_test(count_fit, order_by = "not-a-variable"),
    "unknown 'order_by'"
  )
  expect_error(
    runs_test(count_fit, order_by = seq_len(count_fit$obs_no - 1L)),
    "one value per fitted observation"
  )

  bad_order <- seq_len(count_fit$obs_no)
  bad_order[1] <- NA_real_
  expect_error(
    runs_test(count_fit, order_by = bad_order),
    "finite and non-missing"
  )

  fit <- count_fit
  fit$residuals <- rep(1, fit$obs_no)
  expect_error(
    runs_test(fit, type = "working"),
    "at least one positive and one negative residual"
  )
})
