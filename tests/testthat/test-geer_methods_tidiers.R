test_that("tidy rejects invalid inputs", {
  expect_error(tidy(fit_geewa_pois_exch, conf.int = NA), "'conf.int' must be")
  expect_error(
    tidy(fit_geewa_pois_exch, conf.int = TRUE, conf.level = 0),
    "'conf.level' must be"
  )
  expect_error(
    tidy(fit_geewa_pois_exch, cov_type = "sandwich"),
    "cov_type|arg",
    ignore.case = TRUE
  )
})

test_that("tidy returns coefficient-level output with expected columns", {
  out <- tidy(fit_geewa_pois_exch)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), length(coef(fit_geewa_pois_exch)))
  expect_named(
    out,
    c("term", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(out$term, names(coef(fit_geewa_pois_exch)))
  expect_equal(out$estimate, unname(coef(fit_geewa_pois_exch)), tolerance = 1e-6)
})

test_that("tidy standard errors use the requested covariance type", {
  out <- tidy(fit_geewa_pois_exch, cov_type = "naive")
  expected_se <- unname(sqrt(diag(vcov(fit_geewa_pois_exch, cov_type = "naive"))))
  expect_equal(out$std.error, expected_se, tolerance = 1e-6)
})

test_that("tidy adds confidence intervals when requested", {
  out <- tidy(fit_geewa_pois_exch, conf.int = TRUE, conf.level = 0.95, cov_type = "robust")
  ci <- confint(fit_geewa_pois_exch, level = 0.95, cov_type = "robust")
  expect_true(all(c("conf.low", "conf.high") %in% names(out)))
  expect_equal(out$conf.low, unname(ci[, 1L]), tolerance = 1e-6)
  expect_equal(out$conf.high, unname(ci[, 2L]), tolerance = 1e-6)
})

test_that("tidy exponentiate transforms estimates and intervals only", {
  raw <- tidy(fit_geewa_pois_exch, conf.int = TRUE)
  exp_out <- tidy(fit_geewa_pois_exch, conf.int = TRUE, exponentiate = TRUE)
  expect_equal(exp_out$estimate, exp(raw$estimate), tolerance = 1e-10)
  expect_equal(exp_out$conf.low, exp(raw$conf.low), tolerance = 1e-10)
  expect_equal(exp_out$conf.high, exp(raw$conf.high), tolerance = 1e-10)
  expect_equal(exp_out$std.error, raw$std.error, tolerance = 1e-10)
  expect_equal(exp_out$statistic, raw$statistic, tolerance = 1e-10)
  expect_equal(exp_out$p.value, raw$p.value, tolerance = 1e-10)
})

test_that("tidy works for geewa_binary fits", {
  out <- tidy(fit_geewa_bin_exch)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), length(coef(fit_geewa_bin_exch)))
  expect_equal(out$estimate, unname(coef(fit_geewa_bin_exch)), tolerance = 1e-6)
})

test_that("glance returns a one-row summary with expected fields", {
  out <- glance(fit_geewa_pois_exch)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_true(all(c(
    "family", "link", "method", "wastr",
    "nobs", "nclusters", "min.cluster.size", "max.cluster.size",
    "npar", "df.residual", "phi",
    "QIC", "QICu", "CIC",
    "converged", "niter"
  ) %in% names(out)))
})

test_that("glance summary fields match the fitted object", {
  out <- glance(fit_geewa_pois_exch)
  expect_equal(out$family, fit_geewa_pois_exch$family$family)
  expect_equal(out$link, fit_geewa_pois_exch$family$link)
  expect_equal(out$method, fit_geewa_pois_exch$method)
  expect_equal(out$wastr, fit_geewa_pois_exch$association_structure)
  expect_equal(out$nobs, fit_geewa_pois_exch$obs_no)
  expect_equal(out$nclusters, fit_geewa_pois_exch$clusters_no)
  expect_equal(out$npar, length(coef(fit_geewa_pois_exch)))
  expect_equal(out$df.residual, fit_geewa_pois_exch$df.residual)
  expect_equal(out$phi, fit_geewa_pois_exch$phi, tolerance = 1e-6)
  expect_true(out$converged)
  expect_gte(out$niter, 1L)
})

test_that("glance QIC-based criteria match geecriteria", {
  out <- glance(fit_geewa_pois_exch)
  crit <- geecriteria(fit_geewa_pois_exch, cov_type = "robust", digits = 15)
  expect_equal(out$QIC, crit$QIC, tolerance = 1e-8)
  expect_equal(out$QICu, crit$QICu, tolerance = 1e-8)
  expect_equal(out$CIC, crit$CIC, tolerance = 1e-8)
})

test_that("glance works for geewa_binary fits", {
  out <- glance(fit_geewa_bin_exch)
  expect_equal(out$phi, 1, tolerance = 1e-10)
  expect_equal(out$npar, nrow(tidy(fit_geewa_bin_exch)))
})

test_that("tidy rejects unused or misspelled arguments", {
  expect_error(
    tidy(fit_geewa_pois_exch, covtype = "naive"),
    "tidy.geer does not use the argument(s): 'covtype'",
    fixed = TRUE
  )
  expect_error(
    tidy(fit_geewa_pois_exch, conf_int = TRUE),
    "tidy.geer does not use the argument(s): 'conf_int'",
    fixed = TRUE
  )
})
