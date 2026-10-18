respiratory <- test_data$respiratory

test_that("summary returns a coherent summary object", {
  out <- summary(fit_geewa_pois_exch)
  expect_s3_class(out, "summary.geer")
  expect_true(is.matrix(out$coefficients) || is.data.frame(out$coefficients))
  expect_equal(
    colnames(out$coefficients),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  )
  expect_equal(nrow(out$coefficients), length(coef(fit_geewa_pois_exch)))
  expect_equal(out$family$family, "poisson")
  expect_equal(out$family$link, "log")
})

test_that("summary respects the requested covariance type", {
  for (cov_type in c("robust", "naive", "bias-corrected", "df-adjusted",
                     "kauermann-carroll")) {
    out <- summary(fit_geewa_pois_exch, cov_type = cov_type)

    expect_s3_class(out, "summary.geer")
    expect_identical(out$cov_type, cov_type)
    expect_equal(nrow(out$coefficients), length(coef(fit_geewa_pois_exch)))
    expect_equal(
      unname(out$coefficients[, "Std. Error"]),
      unname(sqrt(diag(vcov(fit_geewa_pois_exch, cov_type = cov_type)))),
      tolerance = 1e-8
    )
  }
})

test_that("summary handles zero standard errors without crashing", {
  fit <- fit_geewa_pois_exch
  fit$naive_covariance[, ] <- 0
  fit$robust_covariance[, ] <- 0
  fit$bias_corrected_covariance[, ] <- 0
  out <- summary(fit)
  expect_s3_class(out, "summary.geer")
  expect_true(all(is.na(out$coefficients[, "z value"])))
})

test_that("print methods return invisibly", {
  expect_invisible(print(fit_geewa_pois_exch))
  expect_invisible(print(summary(fit_geewa_pois_exch)))
})

test_that("print and summary display a notice for a fit that did not converge", {
  res <- collect_warnings(
    geewa(status ~ treatment + baseline,
      family = binomial(), data = respiratory,
      id = id, repeated = visit, corstr = "ar1", method = "gee",
      control = geer_control(maxiter = 1)
    )
  )
  fit <- res$value
  expect_false(fit$converged)
  expect_true(any(grepl("did not converge", res$warnings, fixed = TRUE)))
  out_print <- paste(capture.output(print(fit)), collapse = "\n")
  out_summary <- paste(capture.output(print(summary(fit))), collapse = "\n")
  expect_match(out_print, "WARNING: the algorithm did not converge", fixed = TRUE)
  expect_match(out_summary, "WARNING: the algorithm did not converge", fixed = TRUE)
})

test_that("print and summary of a converged fit carry no convergence notice", {
  fit <- collect_warnings(
    geewa(status ~ treatment + baseline,
      family = binomial(), data = respiratory,
      id = id, repeated = visit, corstr = "exchangeable"
    )
  )$value
  expect_true(fit$converged)
  expect_false(any(grepl("WARNING", capture.output(print(fit)), fixed = TRUE)))
  expect_false(any(grepl("WARNING", capture.output(print(summary(fit))), fixed = TRUE)))
})

test_that("summary rejects unused or misspelled arguments", {
  expect_error(
    summary(fit_geewa_pois_exch, covtype = "naive"),
    "summary.geer does not use the argument(s): 'covtype'",
    fixed = TRUE
  )
  expect_error(
    summary(fit_geewa_bin_exch, "naive", TRUE),
    "summary.geer does not use the argument(s): <unnamed>",
    fixed = TRUE
  )
})
