testthat::local_edition(3)


test_that("little_mcar_test reproduces Little's bivariate monotone result", {
  x <- cbind(
    y1 = 1:8,
    y2 = c(2.0, 3.2, 4.1, 5.0, NA, NA, NA, NA)
  )

  out <- little_mcar_test(x)

  expect_s3_class(out, "htest")
  expect_equal(unname(out$statistic), 5.333333333333333, tolerance = 1e-6)
  expect_equal(unname(out$parameter), c(1, 6))
  expect_equal(out$exact.f.statistic, 19.2, tolerance = 1e-6)
  expect_equal(
    out$p.value,
    stats::pf(19.2, df1 = 1, df2 = 6, lower.tail = FALSE),
    tolerance = 1e-8
  )
  expect_equal(
    out$asymptotic.p.value,
    stats::pchisq(5.333333333333333, df = 1, lower.tail = FALSE),
    tolerance = 1e-8
  )
  expect_identical(out$reference, "exact bivariate F")
  expect_identical(out$missing.patterns, 2L)
  expect_identical(out$n, 8L)
  expect_identical(out$variables, 2L)
  expect_true(out$converged)
  expect_gt(out$iterations, 0L)
})


test_that("little_mcar_test can force the asymptotic reference", {
  x <- cbind(
    y1 = 1:8,
    y2 = c(2.0, 3.2, 4.1, 5.0, NA, NA, NA, NA)
  )

  out <- little_mcar_test(x, reference = "asymptotic")

  expect_equal(unname(out$statistic), 5.333333333333333, tolerance = 1e-6)
  expect_identical(unname(out$parameter), 1L)
  expect_equal(
    out$p.value,
    stats::pchisq(5.333333333333333, df = 1, lower.tail = FALSE),
    tolerance = 1e-8
  )
  expect_identical(out$reference, "asymptotic chi-squared")
  expect_equal(out$exact.p.value, stats::pf(19.2, 1, 6, lower.tail = FALSE), tolerance = 1e-8)
})


test_that("little_mcar_test applies Little's covariance correction", {
  x <- cbind(
    y1 = 1:8,
    y2 = c(2.0, 3.2, 4.1, 5.0, NA, NA, NA, NA)
  )
  out <- little_mcar_test(x, reference = "asymptotic")

  expect_equal(out$covariance, nrow(x) / (nrow(x) - 1) * out$ml.covariance, tolerance = 1e-10)
})


test_that("little_mcar_test agrees with the corrected airquality benchmark", {
  out <- little_mcar_test(datasets::airquality)

  expect_s3_class(out, "htest")
  expect_equal(unname(out$statistic), 34.8954, tolerance = 0.05)
  expect_identical(unname(out$parameter), 14L)
  expect_lt(abs(out$p.value - 0.00152), 7e-5)
  expect_identical(out$missing.patterns, 4L)
  expect_identical(out$reference, "asymptotic chi-squared")
})


test_that("little_mcar_test reconstructs repeated responses from a geer fit", {
  id <- rep(seq_len(30), each = 3)
  visit <- rep(seq_len(3), times = 30)
  xcov <- rep(c(0, 1), length.out = length(id))
  y <- 2 + 0.15 * id + 0.4 * visit + 0.25 * xcov + sin(id * visit / 5)
  y[c(8, 17, 29, 44, 62, 77)] <- NA_real_
  dat <- data.frame(id = id, visit = visit, xcov = xcov, y = y)

  fit <- geewa(
    y ~ xcov + factor(visit),
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )

  wide <- matrix(NA_real_, nrow = 30, ncol = 3)
  wide[cbind(dat$id, dat$visit)] <- dat$y

  from_fit <- little_mcar_test(fit)
  from_wide <- little_mcar_test(wide)

  expect_equal(from_fit$statistic, from_wide$statistic, tolerance = 1e-8)
  expect_identical(from_fit$parameter, from_wide$parameter)
  expect_equal(from_fit$p.value, from_wide$p.value, tolerance = 1e-10)
  expect_identical(from_fit$missing.patterns, from_wide$missing.patterns)
  expect_identical(from_fit$n, 30L)
  expect_identical(from_fit$variables, 3L)
})


test_that("little_mcar_test uses within-cluster row order when repeated is omitted", {
  id <- rep(seq_len(24), each = 3)
  visit <- rep(seq_len(3), times = 24)
  y <- 1.5 + 0.1 * id + 0.3 * visit + sin(id / 3 + visit)
  y[c(11, 28, 46, 65)] <- NA_real_
  dat <- data.frame(id = id, visit = visit, y = y)

  fit <- geewa(
    y ~ factor(visit),
    data = dat,
    id = id,
    family = gaussian(),
    corstr = "independence"
  )

  wide <- matrix(NA_real_, nrow = 24, ncol = 3)
  wide[cbind(dat$id, dat$visit)] <- dat$y

  expect_equal(
    little_mcar_test(fit)$statistic,
    little_mcar_test(wide)$statistic,
    tolerance = 1e-8
  )
})


test_that("little_mcar_test accepts original data supplied explicitly", {
  id <- rep(seq_len(20), each = 3)
  visit <- rep(seq_len(3), times = 20)
  y <- 1 + 0.2 * id + cos(id + visit)
  y[c(6, 23, 41, 54)] <- NA_real_
  dat <- data.frame(id = id, visit = visit, y = y)

  fit <- geewa(
    y ~ factor(visit),
    data = dat,
    id = id,
    repeated = visit,
    family = gaussian(),
    corstr = "independence"
  )
  fit$data <- NULL

  expect_error(
    little_mcar_test(fit),
    "original data are not available"
  )

  out <- little_mcar_test(fit, data = dat)
  expect_s3_class(out, "htest")
  expect_identical(out$n, 20L)
  expect_identical(out$variables, 3L)
})


test_that("little_mcar_test handles complete data explicitly", {
  x <- cbind(
    a = seq_len(10),
    b = seq_len(10)^2
  )
  out <- little_mcar_test(x)

  expect_equal(unname(out$statistic), 0)
  expect_identical(unname(out$parameter), 0L)
  expect_equal(out$p.value, 1)
  expect_identical(out$missing.patterns, 1L)
  expect_identical(out$iterations, 0L)
  expect_true(out$converged)
  expect_equal(out$covariance, stats::cov(x), tolerance = 1e-12)
  expect_equal(out$ml.covariance, stats::cov(x) * 9 / 10, tolerance = 1e-12)
})


test_that("little_mcar_test warns for binary variables", {
  x <- cbind(
    binary = rep(c(0, 1), 6),
    y2 = c(1.1, 2.4, 1.8, 3.7, 2.2, NA, 4.6, 3.1, NA, 5.2, 4.1, 6.0),
    y3 = c(3.2, 1.7, 4.4, 2.9, 5.1, 3.8, 6.2, 4.7, 7.0, 5.5, 8.1, 6.6)
  )

  expect_warning(
    little_mcar_test(x),
    "most appropriate for quantitative variables"
  )
})


test_that("little_mcar_test validates inputs", {
  expect_error(
    little_mcar_test(data.frame(a = c(1, NA, 3), b = letters[1:3])),
    "must be numeric"
  )
  expect_error(
    little_mcar_test(matrix(1:5, ncol = 1)),
    "at least two variables"
  )
  expect_error(
    little_mcar_test(cbind(a = c(1, 2, 3), b = c(NA, NA, NA))),
    "at least two observed values"
  )
  expect_error(
    little_mcar_test(cbind(a = c(1, 2, NA, NA), b = c(NA, NA, 3, 4))),
    "at least one jointly observed value"
  )
  expect_error(
    little_mcar_test(matrix(1:12, ncol = 2), maxit = 0),
    "positive integer"
  )
  expect_error(
    little_mcar_test(matrix(1:12, ncol = 2), tol = 0),
    "positive finite"
  )
  expect_error(
    little_mcar_test(matrix(1:12, ncol = 2), data = data.frame()),
    "only used when 'object'"
  )
  expect_error(
    little_mcar_test(matrix(1:12, ncol = 2), reference = "invalid"),
    "'arg' should be one of"
  )
})
