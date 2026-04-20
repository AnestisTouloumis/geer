testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch
binary_fit <- fit_geewa_bin_exch

test_that("tidy.geer rejects invalid inputs", {
  expect_error(tidy.geer(list()), "'object' must be of 'geer' class")
  expect_error(tidy(count_fit, conf.int = NA), "'conf.int' must be")
  expect_error(
    tidy(count_fit, conf.int = TRUE, conf.level = 0),
    "'conf.level' must be"
  )
  expect_error(
    tidy(count_fit, cov_type = "sandwich"),
    "cov_type|arg",
    ignore.case = TRUE
  )
})


test_that("tidy.geer returns coefficient-level output with expected columns", {
  out <- tidy(count_fit)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), length(coef(count_fit)))
  expect_identical(
    names(out),
    c("term", "estimate", "std.error", "statistic", "p.value")
  )
  expect_identical(out$term, names(coef(count_fit)))
  expect_equal(out$estimate, unname(coef(count_fit)), tolerance = 1e-6)
})


test_that("tidy.geer standard errors use the requested covariance type", {
  out <- tidy(count_fit, cov_type = "naive")
  expected_se <- unname(sqrt(diag(vcov(count_fit, cov_type = "naive"))))
  expect_equal(out$std.error, expected_se, tolerance = 1e-6)
})


test_that("tidy.geer adds confidence intervals when requested", {
  out <- tidy(count_fit, conf.int = TRUE, conf.level = 0.95)
  ci <- confint(count_fit, level = 0.95, cov_type = "robust")
  expect_true(all(c("conf.low", "conf.high") %in% names(out)))
  expect_equal(out$conf.low, unname(ci[, 1L]), tolerance = 1e-6)
  expect_equal(out$conf.high, unname(ci[, 2L]), tolerance = 1e-6)
})


test_that("tidy.geer exponentiate transforms estimates and intervals only", {
  raw <- tidy(count_fit, conf.int = TRUE)
  exp_out <- tidy(count_fit, conf.int = TRUE, exponentiate = TRUE)
  expect_equal(exp_out$estimate, exp(raw$estimate), tolerance = 1e-10)
  expect_equal(exp_out$conf.low, exp(raw$conf.low), tolerance = 1e-10)
  expect_equal(exp_out$conf.high, exp(raw$conf.high), tolerance = 1e-10)
  expect_equal(exp_out$std.error, raw$std.error, tolerance = 1e-10)
  expect_equal(exp_out$statistic, raw$statistic, tolerance = 1e-10)
  expect_equal(exp_out$p.value, raw$p.value, tolerance = 1e-10)
})


test_that("tidy.geer works for geewa_binary fits", {
  out <- tidy(binary_fit)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), length(coef(binary_fit)))
  expect_equal(out$estimate, unname(coef(binary_fit)), tolerance = 1e-6)
})


test_that("glance.geer rejects non-geer input", {
  expect_error(glance.geer(list()), "'object' must be of 'geer' class")
})


test_that("glance.geer returns a one-row summary with expected fields", {
  out <- glance(count_fit)
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


test_that("glance.geer summary fields match the fitted object", {
  out <- glance(count_fit)
  expect_equal(out$family, count_fit$family$family)
  expect_equal(out$link, count_fit$family$link)
  expect_equal(out$method, count_fit$method)
  expect_equal(out$wastr, count_fit$association_structure)
  expect_equal(out$nobs, count_fit$obs_no)
  expect_equal(out$nclusters, count_fit$clusters_no)
  expect_equal(out$npar, length(coef(count_fit)))
  expect_equal(out$df.residual, count_fit$df.residual)
  expect_equal(out$phi, count_fit$phi, tolerance = 1e-6)
  expect_true(out$converged)
  expect_true(out$niter >= 1L)
})


test_that("glance.geer QIC-based criteria match geecriteria()", {
  out <- glance(count_fit)
  crit <- geecriteria(count_fit, cov_type = "robust", digits = 15)
  expect_equal(out$QIC, crit$QIC, tolerance = 1e-2)
  expect_equal(out$QICu, crit$QICu, tolerance = 1e-2)
  expect_equal(out$CIC, crit$CIC, tolerance = 1e-2)
})


test_that("glance.geer works for geewa_binary fits", {
  out <- glance(binary_fit)
  expect_equal(out$phi, 1, tolerance = 1e-10)
  expect_equal(out$npar, nrow(tidy(binary_fit)))
})
