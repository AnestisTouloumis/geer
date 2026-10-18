cerebrovascular <- test_data$cerebrovascular
respiratory2 <- test_data$respiratory2

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
    data = cerebrovascular,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )
  p <- length(fit$coefficients)
  expect_s3_class(fit, "geer")
  expect_named(fit$coefficients)
  expect_true(is.matrix(fit$x))
  expect_length(fit$fitted.values, nrow(fit$x))
  expect_equal(dim(fit$robust_covariance), c(p, p))
})

test_that("geewa_binary fits representative alternative links and association structures", {
  fit_probit <- geewa_binary(
    ecg ~ period + treatment,
    id = id,
    data = cerebrovascular,
    link = "probit",
    orstr = "independence"
  )
  fit_unstr <- geewa_binary(
    status ~ treatment + baseline,
    id = id,
    repeated = visit,
    data = respiratory2,
    link = "logit",
    orstr = "unstructured"
  )
  expect_s3_class(fit_probit, "geer")
  expect_true(fit_probit$converged)
  expect_equal(fit_probit$family$link, "probit")
  expect_s3_class(fit_unstr, "geer")
  expect_true(fit_unstr$converged)
  expect_length(fit_unstr$alpha, choose(max(respiratory2$visit), 2))
})

test_that("geewa_binary is invariant to row order via internal sorting", {
  withr::local_seed(1)
  fit_1 <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cerebrovascular,
    link = "logit",
    orstr = "exchangeable",
    method = "gee"
  )
  cerebrovascular_shuffled <- cerebrovascular[
    sample.int(nrow(cerebrovascular)), ,
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
  expect_named(fit_1$coefficients, names(fit_2$coefficients))
})

test_that("geewa_binary honors subset and na.action and rejects stray arguments", {
  cer <- cerebrovascular
  cer_sub <- cer[as.numeric(cer$id) <= 40, , drop = FALSE]
  fit_subset <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cer,
    subset = as.numeric(id) <= 40
  )
  fit_manual <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cer_sub
  )
  expect_equal(fit_subset$obs_no, nrow(cer_sub))
  expect_equal(coef(fit_subset), coef(fit_manual), tolerance = 1e-8)
  cer_na <- cer
  cer_na$treatment[1L] <- NA
  fit_na <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cer_na
  )
  expect_equal(fit_na$obs_no, nrow(cer) - 1L)
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      repeated = period,
      data = cer_na,
      na.action = na.fail
    )
  )
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      repeated = period,
      data = cer,
      control = geer_control(),
      control_glm = list(),
      not_an_argument = TRUE
    ),
    "does not use the argument"
  )
})
