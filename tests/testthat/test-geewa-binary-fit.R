testthat::local_edition(3)


test_that("geewa_binary fits a coherent exchangeable logit model", {
  fit <- fit_geewa_bin_exch
  expect_s3_class(fit, "geer")
  expect_true(fit$converged)
  expect_equal(fit$family$family, "binomial")
  expect_equal(fit$family$link, "logit")
  expect_true(all(fitted(fit) > 0 & fitted(fit) < 1))
  expect_equal(fit$phi, 1)
  expect_length(fit$alpha, 1L)
  expect_gt(fit$alpha, 0)
})


test_that("geewa_binary returns a well-formed fitted object on a simple independence fit", {
  fit <- geewa_binary(
    formula = ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = test_data$cerebrovascular,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )
  p <- length(fit$coefficients)
  expect_s3_class(fit, "geer")
  expect_named(fit$coefficients)
  expect_true(is.matrix(fit$x))
  expect_equal(length(fit$fitted.values), nrow(fit$x))
  expect_equal(dim(fit$robust_covariance), c(p, p))
})


test_that("geewa_binary fits representative alternative links and association structures", {
  fit_probit <- geewa_binary(
    ecg ~ period + treatment,
    id = id,
    data = test_data$cerebrovascular,
    link = "probit",
    orstr = "independence"
  )
  respiratory_c2 <- test_data$respiratory[test_data$respiratory$center == "C2", ]
  fit_unstr <- geewa_binary(
    status ~ treatment + baseline,
    id = id,
    repeated = visit,
    data = respiratory_c2,
    link = "logit",
    orstr = "unstructured"
  )
  expect_s3_class(fit_probit, "geer")
  expect_true(fit_probit$converged)
  expect_equal(fit_probit$family$link, "probit")
  expect_s3_class(fit_unstr, "geer")
  expect_true(fit_unstr$converged)
  expect_length(fit_unstr$alpha, choose(max(respiratory_c2$visit), 2))
})


test_that("geewa_binary converges for representative alternative methods", {
  methods <- c("brgee-robust", "bcgee-robust", "pgee-jeffreys")
  for (method_name in methods) {
    fit <- geewa_binary(
      ecg ~ period + treatment,
      id = id,
      data = test_data$cerebrovascular,
      link = "logit",
      orstr = "exchangeable",
      method = method_name
    )
    expect_s3_class(fit, "geer")
    expect_true(fit$converged)
    expect_true(all(is.finite(coef(fit))))
  }
})


test_that("geewa_binary is invariant to row order via internal sorting", {
  set.seed(1)
  fit_1 <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = test_data$cerebrovascular,
    link = "logit",
    orstr = "exchangeable",
    method = "gee"
  )
  cerebrovascular_shuffled <- test_data$cerebrovascular[
    sample.int(nrow(test_data$cerebrovascular)),
    ,
    drop = FALSE
  ]
  fit_2 <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cerebrovascular_shuffled,
    link = "logit",
    orstr = "exchangeable",
    method = "gee"
  )
  expect_equal(
    unname(fit_1$coefficients),
    unname(fit_2$coefficients),
    tolerance = 1e-8
  )
  expect_equal(names(fit_1$coefficients), names(fit_2$coefficients))
})
