testthat::local_edition(3)


start_model_matrix <- function() {
  stats::model.matrix(seizures ~ treatment + lnbaseline + lnage, data = epilepsy)
}


start_values <- function(beta_start = NULL,
                         family = stats::poisson(),
                         link = "log",
                         method = "gee",
                         y = epilepsy$seizures,
                         offset = rep(0, nrow(epilepsy))) {
  geer:::compute_geer_start_values(
    model_matrix = start_model_matrix(),
    y = y,
    family = family,
    weights = rep(1, nrow(epilepsy)),
    offset = offset,
    method = method,
    link = link,
    beta_start = beta_start,
    control = geer_control(),
    control_glm = list()
  )
}


binary_start_values <- function(beta_start = NULL,
                                method = "gee",
                                offset = rep(0, nrow(cerebrovascular))) {
  geer:::compute_geer_binary_start_values(
    model_matrix = stats::model.matrix(ecg ~ period + treatment, data = cerebrovascular),
    y = cerebrovascular$ecg,
    family = stats::binomial(),
    weights = rep(1, nrow(cerebrovascular)),
    offset = offset,
    method = method,
    beta_start = beta_start,
    control_glm = list(),
    tolerance = 1e-8,
    jeffreys_power = 0.5
  )
}


## ---------------------------------------------------------------- beta_start ----

test_that("a supplied beta_start is returned as a numeric vector", {
  expect_identical(start_values(beta_start = c(0L, 1L, 2L, 3L)), c(0, 1, 2, 3))
  expect_identical(binary_start_values(beta_start = c(0L, 1L, 2L)), c(0, 1, 2))
})


test_that("beta_start must have one value per model-matrix column", {
  expect_error(
    start_values(beta_start = c(0, 0, 0)),
    "'beta_start' must be a numeric vector of length 4"
  )
  expect_error(
    binary_start_values(beta_start = c(0, 0)),
    "'beta_start' must be a numeric vector of length 3"
  )
})


## ---------------------------------------------------------- derived starts ----

test_that("starting values for the GEE method are the maximum-likelihood fit", {
  reference <- stats::glm(
    seizures ~ treatment + lnbaseline + lnage,
    family = stats::poisson(),
    data = epilepsy
  )

  out <- start_values()

  expect_equal(unname(out), unname(stats::coef(reference)), tolerance = 1e-5)
})


test_that("quasi families and the identity link start from stats::glm.fit", {
  reference <- stats::lm(
    seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy
  )

  out <- start_values(family = stats::gaussian(), link = "identity")
  expect_equal(unname(out), unname(stats::coef(reference)), tolerance = 1e-6)

  quasi_out <- start_values(family = stats::quasipoisson())
  expect_equal(unname(quasi_out), unname(start_values()), tolerance = 1e-5)
})


test_that("every estimation method yields finite starting values", {
  methods <- c(
    "gee", "bcgee-naive", "brgee-naive", "brgee-robust", "brgee-empirical",
    "pgee-jeffreys", "opgee-jeffreys", "hpgee-jeffreys"
  )
  for (method in methods) {
    expect_true(all(is.finite(start_values(method = method))), info = method)
    expect_true(all(is.finite(binary_start_values(method = method))), info = method)
  }
})


test_that("failure to compute starting values asks for beta_start", {
  bad_offset <- rep(Inf, nrow(epilepsy))
  expect_error(
    start_values(offset = bad_offset),
    "cannot compute starting values; please supply 'beta_start'"
  )
  expect_error(
    start_values(
      family = stats::gaussian(),
      link = "identity",
      offset = bad_offset
    ),
    "cannot compute starting values; please supply 'beta_start'"
  )
  expect_error(
    binary_start_values(offset = rep(Inf, nrow(cerebrovascular))),
    "cannot compute starting values; please supply 'beta_start'"
  )
})
