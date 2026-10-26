test_that("vcov returns well-formed covariance matrices for every cov_type", {
  coef_names <- names(coef(fit_geewa_pois_exch))
  standard <- c("robust", "naive", "bias-corrected", "df-adjusted")
  for (cov_type in standard) {
    v <- vcov(fit_geewa_pois_exch, cov_type = cov_type)
    expect_identical(dimnames(v), list(coef_names, coef_names))
    expect_equal(v, t(v), tolerance = 1e-12)
    if (cov_type %in% standard) {
      eigenvalues <- eigen(v, symmetric = TRUE, only.values = TRUE)$values
      expect_true(all(eigenvalues >= -1e-10))
    }
  }
})

test_that("vcov returns the expected stored matrices and df-adjusted scaling", {
  expect_identical(
    vcov(fit_geewa_bin_exch, cov_type = "robust"),
    fit_geewa_bin_exch$robust_covariance
  )
  expect_identical(
    vcov(fit_geewa_bin_exch, cov_type = "naive"),
    fit_geewa_bin_exch$naive_covariance
  )
  expect_identical(
    vcov(fit_geewa_bin_exch, cov_type = "bias-corrected"),
    fit_geewa_bin_exch$bias_corrected_covariance
  )

  clusters <- fit_geewa_pois_exch$clusters_no
  scaling_factor <- clusters / (clusters - length(coef(fit_geewa_pois_exch)))
  expect_equal(
    vcov(fit_geewa_pois_exch, cov_type = "df-adjusted"),
    scaling_factor * vcov(fit_geewa_pois_exch, cov_type = "robust"),
    tolerance = 1e-10
  )
})

test_that("coef and fitted accessors return the stored quantities", {
  b <- coef(fit_geewa_pois_exch)
  expect_type(b, "double")
  expect_named(b)
  expect_length(b, ncol(fit_geewa_pois_exch$x))
  expect_identical(b, fit_geewa_pois_exch$coefficients)
  expect_identical(coefficients(fit_geewa_pois_exch), b)
  expect_type(fitted(fit_geewa_pois_exch), "double")
  expect_identical(fitted.values(fit_geewa_pois_exch), fitted(fit_geewa_pois_exch))
})

test_that("confint returns a well-formed confidence interval matrix", {
  ci <- confint(fit_geewa_bin_exch)
  expect_true(is.matrix(ci))
  expect_equal(dim(ci)[2], 2L)
  expect_identical(colnames(ci), c(" 2.5%", "97.5%"))
  expect_identical(rownames(ci), names(coef(fit_geewa_bin_exch)))
  expect_true(all(is.finite(ci)))
})

test_that("confint supports selecting parameters by name and index", {
  coef_names <- names(coef(fit_geewa_bin_exch))[1:2]
  ci_by_name <- confint(fit_geewa_bin_exch, parm = coef_names)
  ci_by_index <- confint(fit_geewa_bin_exch, parm = 1:2)
  expect_equal(dim(ci_by_name), c(2L, 2L))
  expect_equal(dim(ci_by_index), c(2L, 2L))
  expect_identical(rownames(ci_by_name), coef_names)
  expect_identical(rownames(ci_by_index), coef_names)
  expect_equal(unname(ci_by_name), unname(ci_by_index))
})

test_that("confint validates level and parm", {
  expect_error(confint(fit_geewa_bin_exch, level = 1), "level")
  expect_error(confint(fit_geewa_bin_exch, parm = "not_a_coef"), "parm")
})
