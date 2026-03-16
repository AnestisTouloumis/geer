data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ period + treatment,
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("vcov.geer returns expected covariance matrices and has correct structure", {
  robust_mat <- vcov(fit, cov_type = "robust")
  naive_mat <- vcov(fit, cov_type = "naive")
  bc_mat <- vcov(fit, cov_type = "bias-corrected")
  df_mat <- vcov(fit, cov_type = "df-adjusted")
  coef_names <- names(coef(fit))
  coef_dim <- length(coef_names)
  expect_true(is.matrix(robust_mat))
  expect_true(is.matrix(naive_mat))
  expect_true(is.matrix(bc_mat))
  expect_true(is.matrix(df_mat))
  expect_equal(dim(robust_mat), c(coef_dim, coef_dim))
  expect_equal(dim(naive_mat), c(coef_dim, coef_dim))
  expect_equal(dim(bc_mat), c(coef_dim, coef_dim))
  expect_equal(dim(df_mat), c(coef_dim, coef_dim))
  expect_identical(robust_mat, fit$robust_covariance)
  expect_identical(naive_mat, fit$naive_covariance)
  expect_identical(bc_mat, fit$bias_corrected_covariance)
  cl_no <- fit$clusters_no
  cov_dim <- ncol(fit$robust_covariance)
  expect_true(cl_no > cov_dim)
  expected_df <- (cl_no / (cl_no - cov_dim)) * fit$robust_covariance
  expect_equal(df_mat, expected_df, tolerance = 1e-12)
  mats <- list(
    robust = robust_mat,
    naive = naive_mat,
    `bias-corrected` = bc_mat,
    `df-adjusted` = df_mat
  )
  for (chosen_type in names(mats)) {
    cov_mat <- mats[[chosen_type]]
    expect_identical(colnames(cov_mat), coef_names)
    expect_identical(rownames(cov_mat), coef_names)
  }
})


test_that("vcov.geer errors for invalid cov_type", {
  expect_error(vcov(fit, cov_type = "not-a-type"))
})


test_that("vcov.geer df-adjusted errors when clusters_no <= number of coefficients", {
  cov_dim <- ncol(fit$robust_covariance)
  fit_bad <- fit
  fit_bad$clusters_no <- cov_dim
  expect_error(
    vcov(fit_bad, cov_type = "df-adjusted"),
    "clusters_no must be > number of coefficients"
  )
})
