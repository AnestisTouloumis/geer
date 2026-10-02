testthat::local_edition(3)


test_that("geewa rejects unsupported correlation structures and methods", {
  expect_error(
    geewa(
      seizures ~ treatment,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "compound"
    ),
    "'corstr' must be one of"
  )
  expect_error(
    geewa(
      seizures ~ treatment,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      method = "em-gee"
    ),
    "'method' must be one of"
  )
})


test_that("geewa validates m-dependent and fixed-correlation arguments", {
  expect_error(
    geewa(
      seizures ~ treatment,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "m-dependent",
      Mv = 0
    ),
    "'Mv' must be a positive integer"
  )
  expect_error(
    geewa(
      seizures ~ treatment,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "fixed"
    ),
    "'alpha_vector' must be provided"
  )
})


test_that("geewa validates beta_start and phi_value", {
  expect_error(
    geewa(
      seizures ~ treatment,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      beta_start = c(0, 0, 0)
    ),
    "'beta_start' must be a numeric vector of length"
  )
  expect_error(
    geewa(
      seizures ~ treatment,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      phi_fixed = TRUE,
      phi_value = -1
    ),
    "'phi_value' must be a single positive number"
  )
})


test_that("geewa rejects non-positive weights", {
  bad_weights <- rep(1, nrow(test_data$epilepsy))
  bad_weights[1] <- -1
  expect_error(
    geewa(
      seizures ~ treatment,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      weights = bad_weights
    ),
    "'weights' must be strictly positive"
  )
})


test_that("geewa rejects duplicated repeated values within cluster", {
  epilepsy_bad <- test_data$epilepsy
  epilepsy_bad$rep2 <- epilepsy_bad$visit
  idx <- which(epilepsy_bad$id == epilepsy_bad$id[1])[1:2]
  epilepsy_bad$rep2[idx] <- 1
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = epilepsy_bad,
      id = id,
      repeated = rep2,
      family = poisson(link = "log"),
      corstr = "independence"
    ),
    "unique",
    ignore.case = TRUE
  )
})


test_that("geewa accepts row-level prior weights for two-column binomial responses", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1),
    id = 1:4,
    w = c(1, 2, 1, 3)
  )
  expect_no_error(
    fit <- geewa(
      cbind(success, failure) ~ x,
      data = dat,
      id = id,
      family = binomial(link = "logit"),
      corstr = "independence",
      weights = w,
      phi_fixed = TRUE,
      phi_value = 1
    )
  )
  expect_length(fit$prior.weights, nrow(dat))
  expect_equal(
    unname(fit$prior.weights),
    unname(dat$w * (dat$success + dat$failure))
  )
})


test_that("extract_geer_response_weights treats quasibinomial grouped responses like binomial", {
  dat <- data.frame(
    success = c(1, 3, 2, 4),
    failure = c(4, 2, 3, 1),
    x = c(0, 1, 0, 1)
  )
  mf <- stats::model.frame(
    stats::as.formula(cbind(success, failure) ~ x),
    data = dat
  )
  out_binomial <- extract_geer_response_weights(
    mf,
    stats::binomial(link = "logit")
  )
  out_quasibinomial <- extract_geer_response_weights(
    mf,
    stats::quasibinomial(link = "logit")
  )
  expect_equal(out_quasibinomial$y, out_binomial$y)
  expect_equal(out_quasibinomial$weights, out_binomial$weights)
})


test_that("extract_geer_response_weights treats quasibinomial factor responses like binomial", {
  dat <- data.frame(
    y = factor(c("no", "yes", "no", "yes")),
    x = c(0, 1, 0, 1)
  )
  mf <- stats::model.frame(
    stats::as.formula(y ~ x),
    data = dat
  )
  out_binomial <- extract_geer_response_weights(
    mf,
    stats::binomial(link = "logit")
  )
  out_quasibinomial <- extract_geer_response_weights(
    mf,
    stats::quasibinomial(link = "logit")
  )
  expect_equal(out_quasibinomial$y, out_binomial$y)
  expect_equal(out_quasibinomial$weights, out_binomial$weights)
})
