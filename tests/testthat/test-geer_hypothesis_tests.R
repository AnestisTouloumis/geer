cerebrovascular <- test_data$cerebrovascular

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

ctl_tight <- geer_control(tolerance = 1e-12, maxiter = 500)

fixed_null_fit <- function(d, formula) {
  geewa(formula,
    data = d$df, id = id, repeated = visit, family = d$family, corstr = "fixed",
    alpha_vector = ref_corr(d$struct, d$alpha, d$tmax)[lower.tri(diag(d$tmax))],
    phi_fixed = TRUE, phi_value = d$phi, weights = w, control = ctl_tight
  )
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
    "models must be fitted to the same observations"
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
  expect_named(coef(res1$object0), names(coef(fit_bin_trt)))
  expect_named(coef(res1$object1), names(coef(fit_bin_full)))
  expect_identical(sort(res1$index), sort(expected_index))
  expect_named(coef(res2$object0), names(coef(fit_bin_trt)))
  expect_named(coef(res2$object1), names(coef(fit_bin_full)))
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
    pchisq(rs$test_stat, df = rs$test_df, lower.tail = FALSE)
  )
  sat <- compute_chisq_mixture(x, test_stat, pmethod = "satterthwaite")
  x_bar <- mean(x)
  cv2 <- sum((x - x_bar)^2) / (length(x) * x_bar^2)
  expected_df <- length(x) / (1 + cv2)
  expected_stat <- test_stat / ((1 + cv2) * x_bar)
  expected_p <- pchisq(expected_stat, df = expected_df, lower.tail = FALSE)
  expect_equal(sat$test_df, expected_df)
  expect_equal(sat$test_stat, expected_stat)
  expect_equal(sat$test_p, expected_p)
})

test_that("compute_chisq_mixture keeps extreme p-values above zero", {
  x <- c(1, 1)
  rs <- compute_chisq_mixture(x, 200, pmethod = "rao-scott")
  expect_gt(rs$test_p, 0)
  expect_equal(
    rs$test_p,
    pchisq(200, df = 2, lower.tail = FALSE)
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
  expect_true(all(c(
    "score_vector", "naive_covariance", "robust_covariance",
    "bc_covariance"
  ) %in% names(sc)))
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

test_that("compute_working_lrt_test returns a valid result and scales both models by the larger model's phi", {
  est_trt <- geewa(ecg ~ treatment,
    id = id, family = binomial(link = "logit"),
    data = cerebrovascular, corstr = "independence", method = "gee"
  )
  est_full <- geewa(ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "logit"), data = cerebrovascular,
    corstr = "independence", method = "gee"
  )
  res <- compute_working_lrt_test(est_trt, est_full,
    cov_type = "robust",
    pmethod = "satterthwaite"
  )
  expect_test_result(res)
  small_phi <- est_trt
  small_phi$phi <- small_phi$phi + 1
  res_small <- compute_working_lrt_test(small_phi, est_full,
    cov_type = "robust",
    pmethod = "satterthwaite"
  )
  expect_equal(res_small$test_stat, res$test_stat, tolerance = 1e-12)
  large_phi <- est_full
  large_phi$phi <- large_phi$phi + 1
  res_large <- compute_working_lrt_test(est_trt, large_phi,
    cov_type = "robust",
    pmethod = "satterthwaite"
  )
  expect_equal(res_large$test_stat, res$test_stat * est_full$phi / (est_full$phi + 1),
    tolerance = 1e-10
  )
  expect_test_result(compute_working_lrt_test(fit_bin_trt, fit_bin_full,
    cov_type = "robust", pmethod = "satterthwaite"
  ))
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
  expect_s3_class(out, "data.frame")
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
    "the modified working LRT requires an independence working structure"
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

test_that("check_nested_models rejects models fitted to different data", {
  n <- fit_bin_full$obs_no
  changes <- list(
    "models differ in the response" = function(x) {
      x$y[1] <- 1 - x$y[1]
      x
    },
    "models differ in the cluster identifiers" = function(x) {
      x$id <- rev(x$id)
      x
    },
    "models differ in the repeated indices" = function(x) {
      x$repeated[1] <- x$repeated[1] + 10
      x
    },
    "models differ in the prior weights" = function(x) {
      x$prior.weights[1] <- 2
      x
    },
    "models differ in the offset" = function(x) {
      x$offset <- rep(1, n)
      x
    },
    "models differ in the family or link" = function(x) {
      x$family <- binomial(link = "probit")
      x
    },
    "named coefficients" = function(x) {
      names(x$coefficients) <- NULL
      x
    },
    "models must be nested" = function(x) {
      names(x$coefficients)[2] <- "other"
      x
    }
  )
  for (msg in names(changes)) {
    expect_error(
      check_nested_models(fit_bin_trt, changes[[msg]](fit_bin_full)),
      msg,
      fixed = TRUE
    )
  }
})

test_that("anova omits models with a different response or association structure", {
  reference <- anova(fit_bin_trt, fit_bin_full, test = "wald")
  expect_warning(
    out <- anova(fit_bin_trt, fit_bin_full, fit_bin_full_exch, test = "wald"),
    "association structure differs from that of the first model were omitted: exchangeable"
  )
  expect_equal(out, reference)
  fit_other_response <- update(fit_bin_trt, I(1 - ecg) ~ .)
  expect_warning(
    out <- anova(fit_bin_trt, fit_other_response, fit_bin_full, test = "wald"),
    "response differs from that of the first model were omitted: I(1 - ecg)",
    fixed = TRUE
  )
  expect_equal(out, reference)
  expect_warning(
    out <- anova(fit_bin_full, fit_bin_full_exch, test = "wald"),
    "association structure differs from that of the first model"
  )
  expect_equal(out, anova(fit_bin_full, test = "wald"))
  fit_subset <- update(fit_bin_full, subset = as.numeric(id) > 10)
  expect_error(
    anova(fit_bin_trt, fit_subset, test = "wald"),
    "models must be fitted to the same observations"
  )
})

test_that("Wald and working Wald tests report singular and negative statistics", {
  singular <- fit_bin_full
  singular$robust_covariance[] <- 0
  expect_error(
    compute_wald_test(fit_bin_trt, singular, cov_type = "robust"),
    "Wald test failed: the covariance matrix is singular or invalid"
  )
  singular$naive_covariance[] <- 0
  expect_error(
    compute_working_wald_test(fit_bin_trt, singular, cov_type = "robust"),
    "naive covariance matrix is singular or invalid"
  )
  negative <- fit_bin_full
  negative$naive_covariance <- -negative$naive_covariance
  expect_warning(
    res <- compute_working_wald_test(fit_bin_trt, negative, cov_type = "robust"),
    "the test statistic is negative"
  )
  expect_true(is.na(res$test_stat))
  expect_true(is.na(res$test_p))
})

test_that("covariance helpers of the tests validate their inputs", {
  expect_error(
    compute_df_adjusted_covariance(diag(2), NA_real_, 2L, "Test"),
    "Test failed: invalid clusters_no"
  )
  expect_error(
    compute_df_adjusted_covariance(diag(2), 2, 2L, "Test"),
    "requires more clusters than coefficients"
  )
  expect_error(
    compute_mixture_eigenvalues(matrix(0, 2, 2), diag(2), "Test"),
    "Test failed: could not compute eigenvalues"
  )
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

test_that("build_missing_test_result has an NA statistic and p-value", {
  res <- build_missing_test_result(2L)
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
  expect_equal(res$test_p, pchisq(wald, k, lower.tail = FALSE),
    tolerance = 1e-8
  )

  naive <- vcov(full, cov_type = "naive")[index, index, drop = FALSE]
  working <- as.numeric(crossprod(b, solve(naive, b)))
  lambda <- Re(eigen(solve(naive, robust), only.values = TRUE)$values)
  lambda_bar <- mean(lambda)

  res_rs <- compute_working_wald_test(
    fit_bin_trt, full,
    cov_type = "robust", pmethod = "rao-scott"
  )
  expect_equal(res_rs$test_stat, working / lambda_bar, tolerance = 1e-8)
  expect_equal(
    res_rs$test_p,
    pchisq(working / lambda_bar, k, lower.tail = FALSE),
    tolerance = 1e-8
  )

  cv2 <- sum((lambda - lambda_bar)^2) / (k * lambda_bar^2)
  res_sa <- compute_working_wald_test(
    fit_bin_trt, full,
    cov_type = "robust", pmethod = "satterthwaite"
  )
  expect_equal(res_sa$test_df, k / (1 + cv2), tolerance = 1e-8)
  expect_equal(res_sa$test_stat, working / ((1 + cv2) * lambda_bar),
    tolerance = 1e-8
  )
})

test_that("p-values of large statistics are computed from the upper tail", {
  large <- fit_bin_full
  large$robust_covariance <- large$robust_covariance / 1e3
  res <- compute_wald_test(fit_bin_trt, large, cov_type = "robust")
  expect_gt(res$test_stat, 100)
  expect_gt(res$test_p, 0)
  expect_equal(res$test_p, pchisq(res$test_stat, res$test_df, lower.tail = FALSE),
    tolerance = 1e-10
  )
})

test_that("score and working-score tests match an independent evaluation at the null fit", {
  d <- ref_sim("cc", ref_times_unbalanced, poisson(), "log", "exchangeable", 0.4,
    phi = 1.3, weights = c(0.6, 1.2, 2, 1, 1.5),
    beta = c(0.3, 0.2, -0.1), seed = 31
  )
  f1 <- ref_fit(d, "gee", ctl_tight, beta_start = d$beta)
  f0 <- fixed_null_fit(d, y ~ x1)
  b0 <- c(unname(coef(f0)), 0)
  cv0 <- ref_covariances(d, b0)
  v <- solve(ref_Sigma0(d, b0), ref_U(d, b0))[3]
  K <- d$N
  for (ct in c("robust", "naive", "bias-corrected", "df-adjusted")) {
    C <- switch(ct,
      robust = cv0$robust,
      naive = cv0$naive,
      `bias-corrected` = cv0$bc,
      `df-adjusted` = K / (K - 3) * cv0$robust
    )[3, 3]
    a <- anova(f0, f1, test = "score", cov_type = ct)
    expect_equal(a$Chi[2], v^2 / C, tolerance = 1e-6)
  }
  st <- v^2 / cv0$naive[3, 3]
  lam <- cv0$robust[3, 3] / cv0$naive[3, 3]
  a <- anova(f0, f1, test = "working-score", cov_type = "robust", pmethod = "rao-scott")
  expect_equal(a$Chi[2], st / lam, tolerance = 1e-6)
  expect_equal(a[["Pr(>Chi)"]][2], pchisq(st / lam, 1, lower.tail = FALSE), tolerance = 1e-6)
})

test_that("score tests use the association and dispersion of the larger model", {
  d <- ref_sim("cc", ref_times_unbalanced, poisson(), "log", "exchangeable", 0.4,
    phi = 1.3, weights = c(0.6, 1.2, 2, 1, 1.5), beta = c(0.3, 0.2, -0.1), seed = 31
  )
  f1 <- geewa(y ~ x1 + x2,
    data = d$df, id = id, repeated = visit, family = poisson(),
    corstr = "exchangeable", weights = w, control = ctl_tight
  )
  f0 <- update(f1, . ~ . - x2)
  b0 <- c(unname(coef(f0)), 0)
  nuisance <- function(fit) {
    out <- d
    out$alpha <- as.numeric(fit$alpha)
    out$phi <- fit$phi
    out
  }
  d1 <- nuisance(f1)
  cv1 <- ref_covariances(d1, b0)
  v1 <- solve(ref_Sigma0(d1, b0), ref_U(d1, b0))[3]
  for (ct in c("robust", "naive")) {
    expected <- v1^2 / switch(ct,
      robust = cv1$robust,
      naive = cv1$naive
    )[3, 3]
    expect_equal(anova(f0, f1, test = "score", cov_type = ct)$Chi[2], expected,
      tolerance = 1e-6, label = ct
    )
  }
  d0 <- nuisance(f0)
  v0 <- solve(ref_Sigma0(d0, b0), ref_U(d0, b0))[3]
  with_null_nuisance <- v0^2 / ref_covariances(d0, b0)$robust[3, 3]
  expect_gt(
    abs(anova(f0, f1, test = "score", cov_type = "robust")$Chi[2] - with_null_nuisance),
    0.01
  )
})

test_that("working LRT with estimated dispersion does not depend on y-only terms of the working log-likelihood", {
  d <- ref_sim("cc", ref_times_unbalanced, poisson(), "log", "independence", 0,
    phi = 1.3,
    beta = c(0.3, 0.2, -0.1), seed = 52
  )
  g1 <- geewa(y ~ x1 + x2,
    data = d$df, id = id, repeated = visit, family = poisson(),
    corstr = "independence"
  )
  g0 <- geewa(y ~ x1,
    data = d$df, id = id, repeated = visit, family = poisson(),
    corstr = "independence"
  )
  out <- anova(g0, g1, test = "working-lrt", cov_type = "naive")
  expect_gte(out$Chi[2], 0)
  r1 <- glm(y ~ x1 + x2, family = poisson(), data = d$df)
  r0 <- glm(y ~ x1, family = poisson(), data = d$df)
  phi1 <- sum(residuals(r1, "pearson")^2) / df.residual(r1)
  expect_equal(g1$phi, phi1, tolerance = 1e-6)
  expect_equal(out$Chi[2], (deviance(r0) - deviance(r1)) / phi1, tolerance = 1e-6)
})

test_that("working LRT has power for a Gaussian response and matches the reference statistic", {
  withr::local_seed(9)
  K <- 60
  id <- rep(seq_len(K), each = 3)
  x1 <- rnorm(3 * K)
  x2 <- rnorm(3 * K)
  u <- rnorm(K)[id]
  df <- data.frame(
    id = id, t = rep(1:3, K), x1 = x1, x2 = x2,
    y = 1 + 0.5 * x1 + 3 * x2 + u + rnorm(3 * K)
  )
  g1 <- geewa(y ~ x1 + x2, data = df, id = id, repeated = t, corstr = "independence")
  g0 <- update(g1, . ~ . - x2)
  expect_lt(anova(g0, g1, test = "working-lrt")[2, 4], 1e-6)
  out <- anova(g0, g1, test = "working-lrt", cov_type = "naive")
  r1 <- lm(y ~ x1 + x2, data = df)
  r0 <- lm(y ~ x1, data = df)
  rss1 <- sum(resid(r1)^2)
  rss0 <- sum(resid(r0)^2)
  expect_equal(out$Chi[2], (rss0 - rss1) / (rss1 / df.residual(r1)), tolerance = 1e-8)
  df$y <- 1 + 0.5 * df$x1 + u + rnorm(3 * K)
  g1n <- update(g1, data = df)
  g0n <- update(g0, data = df)
  expect_gt(anova(g0n, g1n, test = "working-lrt")[2, 4], 0.001)
})

test_that("working LRT with a fixed dispersion agrees with the Wald and score tests asymptotically", {
  withr::local_seed(10)
  K <- 400
  id <- rep(seq_len(K), each = 3)
  x1 <- rnorm(3 * K)
  x2 <- rnorm(3 * K)
  df <- data.frame(
    id = id, t = rep(1:3, K), x1 = x1, x2 = x2,
    y = rbinom(3 * K, 1, plogis(-0.2 + 0.4 * x1 + 0.15 * x2))
  )
  g1 <- geewa(y ~ x1 + x2,
    family = binomial(), data = df, id = id, repeated = t,
    corstr = "independence", phi_fixed = TRUE
  )
  g0 <- update(g1, . ~ . - x2)
  lrt <- anova(g0, g1, test = "working-lrt", cov_type = "naive")$Chi[2]
  wald <- anova(g0, g1, test = "wald", cov_type = "naive")$Chi[2]
  score <- anova(g0, g1, test = "score", cov_type = "naive")$Chi[2]
  expect_lt(abs(lrt - wald), 0.25 * max(1, wald))
  expect_lt(abs(lrt - score), 0.25 * max(1, wald))
})
