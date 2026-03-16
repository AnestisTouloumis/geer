.fit_cc <- local({
  data("epilepsy", package = "geer")
  geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data    = epilepsy,
    id      = id,
    family  = poisson(link = "log"),
    corstr  = "exchangeable",
    method  = "gee"
  )
})

.fit_or <- local({
  data("cerebrovascular", package = "geer")
  geewa_binary(
    formula = ecg ~ treatment + factor(period),
    id      = id,
    data    = cerebrovascular,
    link    = "logit",
    orstr   = "exchangeable",
    method  = "gee"
  )
})

test_that("tidy.geer errors on non-geer input", {
  expect_error(tidy.geer(list()),               "'x' must be a 'geer' object")
  expect_error(tidy.geer(lm(mpg ~ cyl, mtcars)), "'x' must be a 'geer' object")
})

test_that("tidy.geer errors on invalid arguments", {
  expect_error(tidy(.fit_cc, conf.int = NA),           "'conf.int' must be")
  expect_error(tidy(.fit_cc, conf.int = "yes"),        "'conf.int' must be")
  expect_error(tidy(.fit_cc, conf.int = TRUE, conf.level = 0),   "'conf.level' must be")
  expect_error(tidy(.fit_cc, conf.int = TRUE, conf.level = 1.5), "'conf.level' must be")
  expect_error(tidy(.fit_cc, exponentiate = NA),       "'exponentiate' must be")
  expect_error(tidy(.fit_cc, cov_type = "sandwich"))
})

test_that("tidy.geer returns one row per coefficient with mandatory columns", {
  out <- tidy(.fit_cc)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), length(coef(.fit_cc)))
  expect_true(all(c("term", "estimate", "std.error",
                    "statistic", "p.value") %in% names(out)))
  expect_false("conf.low"  %in% names(out))
  expect_false("conf.high" %in% names(out))
})

test_that("tidy.geer term and estimate match coef()", {
  out <- tidy(.fit_cc)
  expect_identical(out$term, names(coef(.fit_cc)))
  expect_equal(out$estimate, unname(coef(.fit_cc)), tolerance = 1e-6)
})

test_that("tidy.geer std.error reflects cov_type correctly", {
  for (ct in c("robust", "bias-corrected", "df-adjusted", "naive")) {
    out <- tidy(.fit_cc, cov_type = ct)
    expected_se <- unname(sqrt(diag(vcov(.fit_cc, cov_type = ct))))
    expect_equal(out$std.error, expected_se, tolerance = 1e-6,
                 label = paste("std.error for cov_type =", ct))
  }
})

test_that("tidy.geer robust and naive SEs differ", {
  se_r <- tidy(.fit_cc, cov_type = "robust")$std.error
  se_n <- tidy(.fit_cc, cov_type = "naive")$std.error
  expect_false(isTRUE(all.equal(se_r, se_n, tolerance = 1e-4)))
})

test_that("tidy.geer conf.int adds correct CI columns", {
  out <- tidy(.fit_cc, conf.int = TRUE, conf.level = 0.95)
  expect_true(all(c("conf.low", "conf.high") %in% names(out)))
  ci <- confint(.fit_cc, level = 0.95, cov_type = "robust")
  expect_equal(out$conf.low,  unname(ci[, 1L]), tolerance = 1e-6)
  expect_equal(out$conf.high, unname(ci[, 2L]), tolerance = 1e-6)
})

test_that("tidy.geer narrower CI for smaller conf.level", {
  w90 <- with(tidy(.fit_cc, conf.int = TRUE, conf.level = 0.90),
              conf.high - conf.low)
  w95 <- with(tidy(.fit_cc, conf.int = TRUE, conf.level = 0.95),
              conf.high - conf.low)
  expect_true(all(w90 < w95))
})

test_that("tidy.geer exponentiate transforms estimate and CI, not SE/statistic", {
  raw <- tidy(.fit_cc, conf.int = TRUE)
  exp <- tidy(.fit_cc, conf.int = TRUE, exponentiate = TRUE)
  expect_equal(exp$estimate,  base::exp(raw$estimate),  tolerance = 1e-10)
  expect_equal(exp$conf.low,  base::exp(raw$conf.low),  tolerance = 1e-10)
  expect_equal(exp$conf.high, base::exp(raw$conf.high), tolerance = 1e-10)
  expect_equal(exp$std.error, raw$std.error,            tolerance = 1e-10)
  expect_equal(exp$statistic, raw$statistic,            tolerance = 1e-10)
})

test_that("tidy.geer works on geewa_binary() fit", {
  out <- tidy(.fit_or)
  expect_equal(nrow(out), length(coef(.fit_or)))
  expect_equal(out$estimate, unname(coef(.fit_or)), tolerance = 1e-6)
})

test_that("glance.geer errors on non-geer input", {
  expect_error(glance.geer(list()), "'x' must be a 'geer' object")
  expect_error(glance.geer(42L),    "'x' must be a 'geer' object")
})

test_that("glance.geer returns one row with all expected columns", {
  out <- glance(.fit_cc)
  expect_equal(nrow(out), 1L)
  expect_true(all(c("family", "link", "method", "corstr",
                    "nobs", "nclusters", "min.cluster.size", "max.cluster.size",
                    "npar", "df.residual", "phi",
                    "QIC", "QICu", "CIC",
                    "converged", "niter") %in% names(out)))
})

test_that("glance.geer scalar fields match object slots", {
  g <- glance(.fit_cc)
  expect_equal(g$family,           .fit_cc$family$family)
  expect_equal(g$link,             .fit_cc$family$link)
  expect_equal(g$method,           .fit_cc$method)
  expect_equal(g$corstr,           .fit_cc$association_structure)
  expect_equal(g$nobs,             .fit_cc$obs_no)
  expect_equal(g$nclusters,        .fit_cc$clusters_no)
  expect_equal(g$min.cluster.size, .fit_cc$min_cluster_size)
  expect_equal(g$max.cluster.size, .fit_cc$max_cluster_size)
  expect_equal(g$npar,             length(coef(.fit_cc)))
  expect_equal(g$df.residual,      .fit_cc$df.residual)
  expect_equal(g$phi,              .fit_cc$phi,    tolerance = 1e-6)
  expect_true(g$converged)
  expect_true(g$niter >= 1L)
})

test_that("glance.geer QIC/QICu/CIC match geecriteria()", {
  g    <- glance(.fit_cc)
  crit <- geecriteria(.fit_cc, cov_type = "robust", digits = 15)
  expect_true(is.finite(g$QIC))
  expect_true(is.finite(g$QICu))
  expect_true(is.finite(g$CIC) && g$CIC > 0)
  expect_equal(g$QIC,  crit$QIC,  tolerance = 1e-2)
  expect_equal(g$QICu, crit$QICu, tolerance = 1e-2)
  expect_equal(g$CIC,  crit$CIC,  tolerance = 1e-2)
})

test_that("glance.geer phi is 1 for geewa_binary() fit", {
  expect_equal(glance(.fit_or)$phi, 1, tolerance = 1e-10)
})

test_that("glance npar equals tidy nrow", {
  expect_equal(glance(.fit_cc)$npar, nrow(tidy(.fit_cc)))
})
