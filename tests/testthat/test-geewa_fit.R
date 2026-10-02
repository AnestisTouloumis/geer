testthat::local_edition(3)

fit_geewa_gauss_indep <- geewa(
  formula = score ~ treatment + baseline + time,
  data = test_data$rinse,
  id = id,
  family = gaussian(link = "identity"),
  corstr = "independence",
  method = "gee"
)


test_that("geewa fits a coherent Gaussian independence model", {
  fit <- fit_geewa_gauss_indep
  expect_s3_class(fit, "geer")
  expect_true(fit$converged)
  expect_equal(fit$family$family, "gaussian")
  expect_equal(fit$family$link, "identity")
  expect_equal(names(coef(fit)), colnames(fit$x))
  expect_true(all(is.finite(fitted(fit))))
  expect_equal(length(fitted(fit)), fit$obs_no)
  expect_equal(
    residuals(fit, type = "working"),
    fit$y - fitted(fit),
    tolerance = 1e-10
  )
  expect_gt(fit$phi, 0)
})


test_that("geewa fits coherent Poisson independence and exchangeable models", {
  expect_true(fit_geewa_pois_indep$converged)
  expect_true(all(fitted(fit_geewa_pois_indep) > 0))
  expect_true(fit_geewa_pois_exch$converged)
  expect_true(all(fitted(fit_geewa_pois_exch) > 0))
  expect_length(fit_geewa_pois_exch$alpha, 1L)
  expect_gt(fit_geewa_pois_exch$alpha, -1)
  expect_lt(fit_geewa_pois_exch$alpha, 1)
  expect_equal(
    fit_geewa_pois_indep$clusters_no,
    length(unique(test_data$epilepsy$id))
  )
  expect_equal(
    fit_geewa_pois_indep$obs_no,
    nrow(test_data$epilepsy)
  )
})


test_that("geewa fits representative non-independence correlation structures", {
  fit_ar1 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "ar1"
  )
  fit_unstr <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "unstructured"
  )
  fit_mdep <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "m-dependent",
    Mv = 1
  )
  fit_fixed <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "fixed",
    alpha_vector = rep(0.2, choose(4, 2))
  )
  expect_s3_class(fit_ar1, "geer")
  expect_true(fit_ar1$converged)
  expect_length(fit_ar1$alpha, 1L)
  expect_s3_class(fit_unstr, "geer")
  expect_true(fit_unstr$converged)
  expect_length(fit_unstr$alpha, choose(4, 2))
  expect_s3_class(fit_mdep, "geer")
  expect_true(fit_mdep$converged)
  expect_s3_class(fit_fixed, "geer")
  expect_true(fit_fixed$converged)
})


test_that("geewa respects phi_fixed and phi_value", {
  fit <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "independence",
    phi_fixed = TRUE,
    phi_value = 2
  )
  expect_s3_class(fit, "geer")
  expect_equal(fit$phi, 2, tolerance = 1e-10)
})


test_that("geewa converges for representative alternative methods", {
  methods <- c("brgee-robust", "bcgee-robust", "pgee-jeffreys")
  for (method_name in methods) {
    fit <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "exchangeable",
      method = method_name
    )
    expect_s3_class(fit, "geer")
    expect_true(fit$converged)
    expect_true(all(is.finite(coef(fit))))
  }
})


test_that("bias-reduced estimates differ from plain GEE estimates", {
  fit_br <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "exchangeable",
    method = "brgee-robust"
  )
  expect_false(isTRUE(all.equal(
    coef(fit_geewa_pois_exch),
    coef(fit_br),
    tolerance = 1e-12
  )))
  expect_equal(
    coef(fit_geewa_pois_exch),
    coef(fit_br),
    tolerance = 0.5
  )
})


test_that("geewa is invariant to row order via internal sorting", {
  set.seed(1)
  fit_1 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  epilepsy_shuffled <- test_data$epilepsy[
    sample.int(nrow(test_data$epilepsy)),
    ,
    drop = FALSE
  ]
  fit_2 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy_shuffled,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  expect_equal(
    unname(fit_1$coefficients),
    unname(fit_2$coefficients),
    tolerance = 1e-8
  )
  expect_equal(names(fit_1$coefficients), names(fit_2$coefficients))
})


testthat::test_that("grouped binomial fit does not fail from an n - p^2 phi denominator", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1),
    id = 1:4
  )
  testthat::expect_no_error(
    fit <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = stats::binomial(link = "logit"),
      corstr = "independence",
      use_p = TRUE
    )
  )
  testthat::expect_s3_class(fit, "geer")
})


testthat::test_that("use_p = TRUE and use_p = FALSE both fit for a small grouped binomial model", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1),
    id = 1:4
  )
  testthat::expect_no_error(
    fit_true <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = stats::binomial(link = "logit"),
      corstr = "independence",
      use_p = TRUE
    )
  )
  testthat::expect_no_error(
    fit_false <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = stats::binomial(link = "logit"),
      corstr = "independence",
      use_p = FALSE
    )
  )
  testthat::expect_s3_class(fit_true, "geer")
  testthat::expect_s3_class(fit_false, "geer")
})
