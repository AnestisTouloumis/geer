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
    "models were not fit on the same number of observations"
  )
  expect_error(
    check_nested_models(fit_bin_trt, fit_bin_trt),
    "models must be nested and have different numbers of coefficients"
  )
  expect_error(
    check_nested_models(fit_bin_trt, fit_bin_period),
    "models must be nested and have different numbers of coefficients"
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
  expect_equal(
    rs$test_p,
    stats::pchisq(rs$test_stat, df = rs$test_df, lower.tail = FALSE)
  )
  sat <- compute_chisq_mixture(x, test_stat, pmethod = "satterthwaite")
  x_bar <- mean(x)
  cv2 <- sum((x - x_bar)^2) / (length(x) * x_bar^2)
  expected_df <- length(x) / (1 + cv2)
  expected_stat <- test_stat / ((1 + cv2) * x_bar)
  expected_p <- stats::pchisq(expected_stat, df = expected_df, lower.tail = FALSE)
  expect_equal(sat$test_df, expected_df)
  expect_equal(sat$test_stat, expected_stat)
  expect_equal(sat$test_p, expected_p)
})


test_that("compute_chisq_mixture keeps extreme p-values above zero", {
  ## 1 - pchisq() would return exactly 0 here; the upper tail must not.
  x <- c(1, 1)
  rs <- compute_chisq_mixture(x, 200, pmethod = "rao-scott")
  expect_gt(rs$test_p, 0)
  expect_equal(
    rs$test_p,
    stats::pchisq(200, df = 2, lower.tail = FALSE)
  )
  sat <- compute_chisq_mixture(x, 200, pmethod = "satterthwaite")
  expect_gt(sat$test_p, 0)
})


test_that("compute_chisq_mixture rejects invalid inputs", {
  expect_error(compute_chisq_mixture(c(1, 2), NA_real_), "test_stat")
  expect_error(compute_chisq_mixture(numeric(), 1), "non-empty")
  expect_error(compute_chisq_mixture(c(0, 0), 1), "invalid eigenvalues")
})


test_that("compute_score_components returns expected matrix components for nested geewa fits", {
  sc <- compute_score_components(fit_bin_trt, fit_bin_full)
  expect_type(sc, "list")
  expect_true(all(c("score_vector", "naive_covariance", "robust_covariance", "bc_covariance") %in% names(sc)))
  expect_true(is.numeric(sc$score_vector))
  expect_true(is.matrix(sc$naive_covariance))
  expect_true(is.matrix(sc$robust_covariance))
  expect_true(is.matrix(sc$bc_covariance))
  expect_identical(dim(sc$naive_covariance), dim(sc$robust_covariance))
  expect_identical(dim(sc$naive_covariance), dim(sc$bc_covariance))
})


test_that("compute_wald_test returns a valid result for nested models and is order-invariant", {
  res1 <- compute_wald_test(fit_bin_trt, fit_bin_full, cov_type = "robust")
  res2 <- compute_wald_test(fit_bin_full, fit_bin_trt, cov_type = "robust")
  expect_test_result(res1)
  expect_test_result(res2)
  expect_identical(
    res1$test_df,
    length(setdiff(names(coef(fit_bin_full)), names(coef(fit_bin_trt))))
  )
  expect_equal(res1$test_stat, res2$test_stat, tolerance = 1e-8)
  expect_equal(res1$test_p, res2$test_p, tolerance = 1e-8)
  expect_error(
    compute_wald_test(fit_bin_trt, fit_bin_trt),
    "different numbers of coefficients"
  )
})


test_that("compute_working_wald_test returns a valid result for nested models", {
  res <- compute_working_wald_test(
    fit_bin_trt,
    fit_bin_full,
    cov_type = "robust",
    pmethod = "rao-scott"
  )
  expect_test_result(res)
})


test_that("compute_working_lrt_test returns a valid result and scales each model by its own phi", {
  res <- compute_working_lrt_test(
    fit_bin_trt,
    fit_bin_full,
    cov_type = "robust",
    pmethod = "satterthwaite"
  )
  expect_test_result(res)
  fit_bad_phi <- fit_bin_full
  fit_bad_phi$phi <- fit_bad_phi$phi + 1
  fit_bad_phi0 <- fit_bin_trt
  fit_bad_phi0$phi <- fit_bad_phi0$phi + 1
  res_phi <- compute_working_lrt_test(
    fit_bad_phi0,
    fit_bad_phi,
    cov_type = "robust",
    pmethod = "satterthwaite"
  )
  expect_test_result(res_phi)
  ## The rescaled dispersions change the working log-likelihoods, so the
  ## statistic must differ from the one based on the estimated dispersions.
  expect_false(isTRUE(all.equal(res_phi$test_stat, res$test_stat)))
})


test_that("compute_score_test returns a valid result for supported covariance types", {
  res_robust <- compute_score_test(fit_bin_trt, fit_bin_full, cov_type = "robust")
  res_naive <- compute_score_test(fit_bin_trt, fit_bin_full, cov_type = "naive")
  res_df <- compute_score_test(fit_bin_trt, fit_bin_full, cov_type = "df-adjusted")
  expect_test_result(res_robust)
  expect_test_result(res_naive)
  expect_test_result(res_df)
})


test_that("compute_working_score_test returns a valid result for supported covariance types", {
  res_robust <- compute_working_score_test(
    fit_bin_trt,
    fit_bin_full,
    cov_type = "robust",
    pmethod = "rao-scott"
  )
  res_bc <- compute_working_score_test(
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


test_that("compute_anova_geer_list rejects non-independence models for working-lrt", {
  expect_error(
    compute_anova_geer_list(
      list(fit_bin_trt, fit_bin_full_exch),
      test = "working-lrt",
      cov_type = "robust",
      pmethod = "rao-scott"
    ),
    "the modified working LRT requires all models to use an independence working structure"
  )
})
test_that("get_geer_fit_function prefers the stored fit function over the call", {
  fit <- fit_bin_full
  expect_identical(get_geer_fit_function(fit), "geewa")
  fit$call[[1L]] <- quote(geer::geewa)
  expect_identical(get_geer_fit_function(fit), "geewa")
  expect_true(is_geewa_fit(fit))
  fit$call[[1L]] <- quote(some_wrapper)
  expect_true(is_geewa_fit(fit))
  expect_false(is_geewa_fit(fit_geewa_bin_exch))
})


test_that("the fit function is recovered from the call when it is not stored", {
  fit <- fit_bin_full
  fit$fit_function <- NULL
  expect_identical(get_geer_fit_function(fit), "geewa")
  fit$call[[1L]] <- quote(geer::geewa)
  expect_identical(get_geer_fit_function(fit), "geewa")
  fit$call[[1L]] <- quote(geer::geewa_binary)
  expect_identical(get_geer_fit_function(fit), "geewa_binary")
  fit$call[[1L]] <- quote(some_wrapper)
  expect_true(is.na(get_geer_fit_function(fit)))
  expect_error(is_geewa_fit(fit), "cannot determine whether the model was fitted")
})


test_that("score tests do not depend on how geewa was called", {
  reference <- compute_score_test(fit_bin_trt, fit_bin_full)
  fit0 <- fit_bin_trt
  fit1 <- fit_bin_full
  fit0$fit_function <- NULL
  fit1$fit_function <- NULL
  fit0$call[[1L]] <- quote(geer::geewa)
  fit1$call[[1L]] <- quote(geer::geewa)
  expect_equal(compute_score_test(fit0, fit1), reference)
})


test_that("is_phi_fixed uses the stored flag rather than parsing the call", {
  fixed_flag <- TRUE
  fit_flag <- geewa(
    formula = ecg ~ treatment,
    id = id,
    family = binomial(link = "logit"),
    phi_fixed = fixed_flag,
    phi_value = 1,
    data = cerebrovascular,
    corstr = "independence",
    method = "gee"
  )
  expect_true(is_phi_fixed(fit_flag))
  expect_true(is_phi_fixed(fit_bin_trt))
  expect_false(is_phi_fixed(fit_geewa_pois_indep))
  expect_true(is_phi_fixed(fit_geewa_bin_exch))
})

test_that("check_nested_models rejects fits with different settings", {
  fit_other_method <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "logit"),
    phi_fixed = TRUE,
    phi_value = 1,
    data = cerebrovascular,
    corstr = "independence",
    method = "brgee-robust"
  )
  expect_error(
    check_nested_models(fit_bin_trt, fit_other_method),
    "models differ in the estimation method"
  )
  expect_error(
    check_nested_models(fit_bin_trt, fit_bin_full_exch),
    "models differ in the working association structure"
  )
  fit_free_phi <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "logit"),
    data = cerebrovascular,
    corstr = "independence",
    method = "gee"
  )
  expect_error(
    check_nested_models(fit_bin_trt, fit_free_phi),
    "models differ in the treatment of the dispersion parameter"
  )
  expect_error(
    check_comparable_fit_settings(fit_bin_trt, fit_geewa_bin_exch),
    "models differ in the fitting function"
  )
  expect_silent(check_comparable_fit_settings(fit_bin_trt, fit_bin_full))
})


test_that("check_test_statistic reports a negative statistic as NA with a warning", {
  expect_equal(check_test_statistic(2.5, "Wald test"), 2.5)
  expect_equal(check_test_statistic(-1e-12, "Wald test"), 0)
  expect_warning(
    out <- check_test_statistic(-0.5, "Wald test"),
    "Wald test: the test statistic is negative"
  )
  expect_true(is.na(out))
  expect_error(
    check_test_statistic(Inf, "Wald test"),
    "non-finite test statistic"
  )
})


test_that("missing_test_result has an NA statistic and p-value", {
  res <- missing_test_result(2L)
  expect_true(is.na(res$test_stat))
  expect_true(is.na(res$test_p))
  expect_equal(res$test_df, 2L)
})


test_that("Wald and working Wald statistics match independent quadratic forms", {
  full <- fit_bin_full
  index <- match(
    setdiff(names(coef(full)), names(coef(fit_bin_trt))),
    names(coef(full))
  )
  b <- unname(coef(full)[index])
  k <- length(index)

  robust <- vcov(full, cov_type = "robust")[index, index, drop = FALSE]
  wald <- as.numeric(crossprod(b, solve(robust, b)))
  res <- compute_wald_test(fit_bin_trt, full, cov_type = "robust")
  expect_equal(res$test_stat, wald, tolerance = 1e-8)
  expect_equal(res$test_df, k)
  expect_equal(res$test_p, stats::pchisq(wald, k, lower.tail = FALSE),
               tolerance = 1e-8)

  naive <- vcov(full, cov_type = "naive")[index, index, drop = FALSE]
  working <- as.numeric(crossprod(b, solve(naive, b)))
  lambda <- Re(eigen(solve(naive, robust), only.values = TRUE)$values)
  lambda_bar <- mean(lambda)

  res_rs <- compute_working_wald_test(
    fit_bin_trt, full, cov_type = "robust", pmethod = "rao-scott"
  )
  expect_equal(res_rs$test_stat, working / lambda_bar, tolerance = 1e-8)
  expect_equal(
    res_rs$test_p,
    stats::pchisq(working / lambda_bar, k, lower.tail = FALSE),
    tolerance = 1e-8
  )

  cv2 <- sum((lambda - lambda_bar)^2) / (k * lambda_bar^2)
  res_sa <- compute_working_wald_test(
    fit_bin_trt, full, cov_type = "robust", pmethod = "satterthwaite"
  )
  expect_equal(res_sa$test_df, k / (1 + cv2), tolerance = 1e-8)
  expect_equal(res_sa$test_stat, working / ((1 + cv2) * lambda_bar),
               tolerance = 1e-8)
})
