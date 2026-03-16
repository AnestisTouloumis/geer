data("cerebrovascular", package = "geer")
test_that("geewa_binary validates link", {
  skip_if_not(exists("geewa_binary", mode = "function"))

  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      data = cerebrovascular,
      link = "not-a-link"
    ),
    "link"
  )
})


test_that("geewa_binary validates weights length and positivity", {
  skip_if_not(exists("geewa_binary", mode = "function"))
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      data = cerebrovascular,
      link = "logit",
      weights = rep.int(1, nrow(cerebrovascular) - 1)
    ),
    "weights"
  )
  bad_w <- rep.int(1, nrow(cerebrovascular))
  bad_w[1] <- 0
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      data = cerebrovascular,
      link = "logit",
      weights = bad_w
    ),
    "positive|strictly",
    ignore.case = TRUE
  )
})


test_that("geewa_binary rejects duplicated repeated values within cluster", {
  skip_if_not(exists("geewa_binary", mode = "function"))
  dat2 <- cerebrovascular
  dat2$rep2 <- dat2$period
  first_id <- dat2$id[1]
  idx <- which(dat2$id == first_id)[1:2]
  dat2$rep2[idx] <- 1
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      repeated = rep2,
      data = dat2,
      link = "logit",
      orstr = "independence"
    ),
    "unique",
    ignore.case = TRUE
  )
})


test_that("geewa_binary validates orstr and alpha_vector for fixed", {
  skip_if_not(exists("geewa_binary", mode = "function"))
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      data = cerebrovascular,
      link = "logit",
      orstr = "not-a-structure"
    ),
    "orstr"
  )
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      repeated = period,
      data = cerebrovascular,
      link = "logit",
      orstr = "fixed"
    ),
    "alpha_vector",
    ignore.case = TRUE
  )
})


test_that("geewa_binary fits a basic model and returns expected structure", {
  skip_if_not(exists("geewa_binary", mode = "function"))
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  fit <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cerebrovascular,
    link = "logit",
    orstr = "exchangeable",
    method = "gee"
  )
  expect_s3_class(fit, "geer")
  expect_true(is.numeric(fit$coefficients))
  expect_named(fit$coefficients)
  expect_true(is.logical(fit$converged) && length(fit$converged) == 1L)
  p <- length(fit$coefficients)
  expect_true(is.matrix(fit$naive_covariance))
  expect_equal(dim(fit$naive_covariance), c(p, p))
  expect_true(is.matrix(fit$robust_covariance))
  expect_equal(dim(fit$robust_covariance), c(p, p))
  expect_true(is.matrix(fit$bias_corrected_covariance))
  expect_equal(dim(fit$bias_corrected_covariance), c(p, p))
  expect_true(is.character(fit$association_structure))
  expect_true(is.numeric(fit$phi) && length(fit$phi) == 1L)
  expect_false(is.null(attr(fit$x, "assign")))
  expect_equal(length(attr(fit$x, "assign")), ncol(fit$x))
})


test_that("geewa_binary is invariant to row order via internal sorting", {
  skip_if_not(exists("geewa_binary", mode = "function"))
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  fit1 <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cerebrovascular,
    link = "logit",
    orstr = "exchangeable",
    method = "gee"
  )
  set.seed(1)
  dat2 <- cerebrovascular[sample.int(nrow(cerebrovascular)), , drop = FALSE]
  fit2 <- geewa_binary(
    ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = dat2,
    link = "logit",
    orstr = "exchangeable",
    method = "gee"
  )
  expect_equal(unname(fit1$coefficients), unname(fit2$coefficients), tolerance = 1e-8)
  expect_equal(names(fit1$coefficients), names(fit2$coefficients))
})
