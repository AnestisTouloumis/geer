testthat::local_edition(3)

data("cerebrovascular", package = "geer")

fit_bin_null <- geewa(
  formula = ecg ~ 1,
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_bin_trt <- geewa(
  formula = ecg ~ treatment,
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_bin_period <- geewa(
  formula = ecg ~ factor(period),
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_bin_full <- geewa(
  formula = ecg ~ treatment + factor(period),
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_bin_full_exch <- geewa(
  formula = ecg ~ treatment + factor(period),
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = cerebrovascular,
  corstr = "exchangeable",
  method = "gee"
)

expect_test_result <- function(x) {
  expect_type(x, "list")
  expect_true(all(c("test_stat", "test_df", "test_p") %in% names(x)))
  expect_true(is.numeric(x$test_stat))
  expect_true(is.numeric(x$test_df))
  expect_true(is.numeric(x$test_p))
  expect_length(x$test_stat, 1L)
  expect_length(x$test_df, 1L)
  expect_length(x$test_p, 1L)
  expect_true(is.finite(x$test_stat))
  expect_true(is.finite(x$test_df))
  expect_true(is.finite(x$test_p))
  expect_gte(x$test_stat, 0)
  expect_gt(x$test_df, 0)
  expect_gte(x$test_p, 0)
  expect_lte(x$test_p, 1)
}

test_that("check_nested_models validates class and nesting requirements", {
  expect_error(
    check_nested_models(1, fit_bin_full),
    "object0"
  )
  expect_error(
    check_nested_models(fit_bin_full, 1),
    "object1"
  )

  fit_bad_obs <- fit_bin_full
  fit_bad_obs$obs_no <- fit_bad_obs$obs_no + 1L
  expect_error(
    check_nested_models(fit_bin_trt, fit_bad_obs),
    "same observations"
  )

  expect_error(
    check_nested_models(fit_bin_trt, fit_bin_trt),
    "different numbers of coefficients"
  )

  expect_error(
    check_nested_models(fit_bin_trt, fit_bin_period),
    "nested"
  )
})

test_that("check_nested_models returns the smaller and larger model in the right order", {
  res1 <- check_nested_models(fit_bin_trt, fit_bin_full)
  res2 <- check_nested_models(fit_bin_full, fit_bin_trt)

  added_terms <- setdiff(names(coef(fit_bin_full)), names(coef(fit_bin_trt)))
  expected_index <- match(added_terms, names(coef(fit_bin_full)))

  expect_identical(names(coef(res1$object0)), names(coef(fit_bin_trt)))
  expect_identical(names(coef(res1$object1)), names(coef(fit_bin_full)))
  expect_identical(sort(res1$index), sort(expected_index))

  expect_identical(names(coef(res2$object0)), names(coef(fit_bin_trt)))
  expect_identical(names(coef(res2$object1)), names(coef(fit_bin_full)))
  expect_identical(sort(res2$index), sort(expected_index))
})

test_that("compute_chisq_mixture computes Rao-Scott and Satterthwaite approximations correctly", {
  x <- c(1, 2, 3)
  test_stat <- 10

  rs <- compute_chisq_mixture(x, test_stat, pmethod = "rao-scott")
  expect_equal(rs$test_df, 3)
  expect_equal(rs$test_stat, test_stat / mean(x))
  expect_equal(rs$test_p, 1 - stats::pchisq(rs$test_stat, df = rs$test_df))

  sat <- compute_chisq_mixture(x, test_stat, pmethod = "satterthwaite")
  x_bar <- mean(x)
  cv2 <- sum((x - x_bar)^2) / (length(x) * x_bar^2)
  expected_df <- length(x) / (1 + cv2)
  expected_stat <- test_stat / ((1 + cv2) * x_bar)
  expected_p <- 1 - stats::pchisq(expected_stat, df = expected_df)

  expect_equal(sat$test_df, expected_df)
  expect_equal(sat$test_stat, expected_stat)
  expect_equal(sat$test_p, expected_p)
})

test_that("compute_chisq_mixture rejects invalid inputs", {
  expect_error(compute_chisq_mixture(c(1, 2), NA_real_), "test_stat")
  expect_error(compute_chisq_mixture(numeric(), 1), "non-empty")
  expect_error(compute_chisq_mixture(c(0, 0), 1), "invalid eigenvalues")
})

test_that("compute_score_components returns expected matrix components for nested geewa fits", {
  coeffs_test <- coef(fit_bin_full)
  coeffs_test[setdiff(names(coeffs_test), names(coef(fit_bin_trt)))] <- 0
  coeffs_test[names(coef(fit_bin_trt))] <- coef(fit_bin_trt)

  sc <- compute_score_components(fit_bin_trt, fit_bin_full, coeffs_test)

  expect_type(sc, "list")
  expect_true(all(c("score_vector", "naive_covariance", "robust_covariance", "bc_covariance") %in% names(sc)))
  expect_true(is.numeric(sc$score_vector))
  expect_true(is.matrix(sc$naive_covariance))
  expect_true(is.matrix(sc$robust_covariance))
  expect_true(is.matrix(sc$bc_covariance))
  expect_identical(dim(sc$naive_covariance), dim(sc$robust_covariance))
  expect_identical(dim(sc$naive_covariance), dim(sc$bc_covariance))
})

test_that("wald_test returns a valid result for nested models and is order-invariant", {
  res1 <- wald_test(fit_bin_trt, fit_bin_full, cov_type = "robust")
  res2 <- wald_test(fit_bin_full, fit_bin_trt, cov_type = "robust")

  expect_test_result(res1)
  expect_test_result(res2)
  expect_identical(
    res1$test_df,
    length(setdiff(names(coef(fit_bin_full)), names(coef(fit_bin_trt))))
  )
  expect_equal(res1$test_stat, res2$test_stat, tolerance = 1e-8)
  expect_equal(res1$test_p, res2$test_p, tolerance = 1e-8)

  expect_error(
    wald_test(fit_bin_trt, fit_bin_trt),
    "different numbers of coefficients"
  )
})

test_that("working_wald_test returns a valid result for nested models", {
  res <- working_wald_test(
    fit_bin_trt,
    fit_bin_full,
    cov_type = "robust",
    pmethod = "rao-scott"
  )

  expect_test_result(res)
})

test_that("working_lrt_test returns a valid result and checks phi consistency", {
  res <- working_lrt_test(
    fit_bin_trt,
    fit_bin_full,
    cov_type = "robust",
    pmethod = "satterthwaite"
  )

  expect_test_result(res)

  fit_bad_phi <- fit_bin_full
  fit_bad_phi$phi <- fit_bad_phi$phi + 1
  expect_error(
    working_lrt_test(fit_bin_trt, fit_bad_phi, cov_type = "robust"),
    "dispersion parameters differ"
  )
})

test_that("score_test returns a valid result for supported covariance types", {
  res_robust <- score_test(fit_bin_trt, fit_bin_full, cov_type = "robust")
  res_naive <- score_test(fit_bin_trt, fit_bin_full, cov_type = "naive")
  res_df <- score_test(fit_bin_trt, fit_bin_full, cov_type = "df-adjusted")

  expect_test_result(res_robust)
  expect_test_result(res_naive)
  expect_test_result(res_df)
})

test_that("working_score_test returns a valid result for supported covariance types", {
  res_robust <- working_score_test(
    fit_bin_trt,
    fit_bin_full,
    cov_type = "robust",
    pmethod = "rao-scott"
  )
  res_bc <- working_score_test(
    fit_bin_trt,
    fit_bin_full,
    cov_type = "bias-corrected",
    pmethod = "satterthwaite"
  )

  expect_test_result(res_robust)
  expect_test_result(res_bc)
})

test_that("compute_anova_geer_list returns an anova table for multiple nested models", {
  out <- compute_anova_geer_list(
    list(fit_bin_null, fit_bin_trt, fit_bin_full),
    test = "wald",
    cov_type = "robust",
    pmethod = "rao-scott"
  )

  expect_s3_class(out, "anova")
  expect_true(is.data.frame(out))
  expect_identical(nrow(out), 3L)
  expect_true(all(c("Resid. Df", "Df", "Chi", "Pr(>Chi)") %in% names(out)))
})

test_that("compute_anova_geer_list filters non-independence models for working-lrt and still returns a result", {
  expect_warning(
    out <- compute_anova_geer_list(
      list(fit_bin_trt, fit_bin_full_exch),
      test = "working-lrt",
      cov_type = "robust",
      pmethod = "rao-scott"
    ),
    "association structure differs from independence"
  )

  expect_s3_class(out, "anova")
  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1L)
})
