testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch
binary_fit <- fit_geewa_bin_exch

test_that("vcov returns well-formed covariance matrices for supported cov_type values", {
  p <- length(coef(count_fit))

  for (cov_type in c("robust", "naive", "bias-corrected", "df-adjusted")) {
    v <- vcov(count_fit, cov_type = cov_type)

    expect_true(is.matrix(v))
    expect_equal(dim(v), c(p, p))
    expect_equal(v, t(v), tolerance = 1e-12)
    expect_equal(rownames(v), names(coef(count_fit)))
    expect_equal(colnames(v), names(coef(count_fit)))

    eigenvalues <- eigen(v, symmetric = TRUE, only.values = TRUE)$values
    expect_true(all(eigenvalues >= -1e-10))
  }
})

test_that("vcov returns the expected stored matrices and df-adjusted scaling", {
  expect_identical(
    vcov(binary_fit, cov_type = "robust"),
    binary_fit$robust_covariance
  )
  expect_identical(
    vcov(binary_fit, cov_type = "naive"),
    binary_fit$naive_covariance
  )
  expect_identical(
    vcov(binary_fit, cov_type = "bias-corrected"),
    binary_fit$bias_corrected_covariance
  )

  scaling_factor <- count_fit$clusters_no / (count_fit$clusters_no - length(coef(count_fit)))
  expect_equal(
    vcov(count_fit, cov_type = "df-adjusted"),
    scaling_factor * vcov(count_fit, cov_type = "robust"),
    tolerance = 1e-10
  )
})

test_that("vcov validates cov_type and df-adjusted requirements", {
  expect_error(
    vcov(binary_fit, cov_type = "not-a-type"),
    "should be one of"
  )

  stub <- unclass(count_fit)
  stub$clusters_no <- length(coef(count_fit))
  bad_fit <- new_geer(stub)

  expect_error(
    vcov(bad_fit, cov_type = "df-adjusted"),
    "clusters_no must be > number of coefficients"
  )
})

test_that("residuals return finite vectors for supported types", {
  for (type in c("working", "pearson", "deviance")) {
    r <- residuals(count_fit, type = type)
    expect_length(r, count_fit$obs_no)
    expect_true(all(is.finite(r)))
  }

  r_bin <- residuals(binary_fit, type = "working")
  expect_length(r_bin, binary_fit$obs_no)
  expect_true(all(is.finite(r_bin)))
})

test_that("pearson residuals match the Poisson reference formula", {
  y <- count_fit$y
  mu <- fitted(count_fit)
  wt <- count_fit$prior.weights

  expected <- (y - mu) / sqrt(mu) * sqrt(wt)
  observed <- residuals(count_fit, type = "pearson")

  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("residuals validates type", {
  expect_error(
    residuals(count_fit, type = "raw"),
    "should be one of"
  )
})

test_that("coef and fitted accessors return the stored quantities", {
  b <- coef(count_fit)

  expect_type(b, "double")
  expect_named(b)
  expect_equal(length(b), ncol(count_fit$x))
  expect_identical(b, count_fit$coefficients)
  expect_identical(coefficients(count_fit), b)

  expect_type(fitted(count_fit), "double")
  expect_identical(fitted.values(count_fit), fitted(count_fit))
})

test_that("model.matrix returns the stored design matrix", {
  x <- model.matrix(count_fit)

  expect_equal(x, count_fit$x)
  expect_equal(nrow(x), count_fit$obs_no)
  expect_equal(ncol(x), length(coef(count_fit)))
})

test_that("model.matrix.geer rejects non-geer objects", {
  expect_error(
    model.matrix.geer(list(terms = terms(~ 1), data = data.frame(x = 1))),
    "must be a 'geer' object"
  )
})
