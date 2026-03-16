test_that("geewa fits a basic model and returns a geer object", {
  skip_if_not(exists("geewa", mode = "function"))
  # geewa may use brglmFit for starting values in many non-identity link cases
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  data("epilepsy", package = "geer")
  fit <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
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
  expect_equal(fit$association_structure, "exchangeable")
  expect_true(is.numeric(fit$phi) && length(fit$phi) == 1L)
  expect_true(is.numeric(fit$alpha))
})


test_that("geewa is invariant to row order via internal sorting", {
  skip_if_not(exists("geewa", mode = "function"))
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  data("epilepsy", package = "geer")
  fit1 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  set.seed(1)
  epi2 <- epilepsy[sample.int(nrow(epilepsy)), , drop = FALSE]
  fit2 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = epi2,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  expect_equal(unname(fit1$coefficients), unname(fit2$coefficients), tolerance = 1e-8)
  expect_equal(names(fit1$coefficients), names(fit2$coefficients))
})


test_that("geewa rejects invalid phi_value when phi_fixed = TRUE", {
  skip_if_not(exists("geewa", mode = "function"))
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  data("epilepsy", package = "geer")
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = epilepsy,
      id = id,
      family = poisson(link = "log"),
      corstr = "independence",
      method = "gee",
      phi_fixed = TRUE,
      phi_value = 0
    ),
    "phi_value"
  )
})


test_that("geewa validates weights length and positivity", {
  skip_if_not(exists("geewa", mode = "function"))
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  data("epilepsy", package = "geer")
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = epilepsy,
      id = id,
      family = poisson(link = "log"),
      corstr = "independence",
      weights = rep(1, nrow(epilepsy) - 1)
    ),
    "weights"
  )
  bad_w <- rep(1, nrow(epilepsy))
  bad_w[1] <- 0
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = epilepsy,
      id = id,
      family = poisson(link = "log"),
      corstr = "independence",
      weights = bad_w
    ),
    "positive|strictly",
    ignore.case = TRUE
  )
})


test_that("geewa rejects duplicated repeated values within cluster", {
  skip_if_not(exists("geewa", mode = "function"))
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  data("epilepsy", package = "geer")
  epi2 <- epilepsy
  epi2$rep2 <- epi2$visit

  first_id <- epi2$id[1]
  idx <- which(epi2$id == first_id)[1:2]
  epi2$rep2[idx] <- 1
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = epi2,
      id = id,
      repeated = rep2,
      family = poisson(link = "log"),
      corstr = "independence"
    ),
    "unique",
    ignore.case = TRUE
  )
})


test_that("geewa validates alpha_vector length when corstr = fixed", {
  skip_if_not(exists("geewa", mode = "function"))
  if (!exists("brglmFit", mode = "function")) {
    skip("brglmFit not available (needed for starting values in this test)")
  }
  data("epilepsy", package = "geer")
  Tmax <- max(epilepsy$visit)
  need <- choose(Tmax, 2)
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = epilepsy,
      id = id,
      repeated = visit,
      family = poisson(link = "log"),
      corstr = "fixed",
      alpha_vector = rep(0, max(1, need - 1))
    ),
    "alpha_vector",
    ignore.case = TRUE
  )
})
