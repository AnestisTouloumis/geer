testthat::local_edition(3)


finalize_fit <- function(criterion = 1e-8,
                         fitted = c(0.2, 0.5, 0.8, 0.3, 0.6, 0.7),
                         method = "gee",
                         family = stats::binomial(),
                         structure = "exchangeable",
                         repeated = c(1, 2, 3, 1, 2, 3),
                         alpha = 0.3,
                         fit_function = "geewa",
                         tolerance = 1e-6,
                         iter = 1L) {
  geer:::finalize_geer_fit(
    fit = list(iter = iter, fitted.values = fitted),
    geesolver_fit = list(criterion = criterion, alpha = alpha),
    tolerance = tolerance,
    method = method,
    family = family,
    association_structure = structure,
    repeated = repeated,
    fit_function = fit_function
  )
}


test_that("a converged fit is flagged as such and raises no warning", {
  expect_no_warning(out <- finalize_fit())

  expect_true(out$converged)
  expect_identical(out$fit_function, "geewa")
  expect_equal(out$alpha, 0.3)
})


test_that("convergence is judged at the final iteration", {
  criterion <- c(1, 1e-9)

  expect_no_warning(out <- finalize_fit(criterion = criterion, iter = 2L))
  expect_true(out$converged)

  expect_warning(
    out <- finalize_fit(criterion = criterion, iter = 1L),
    "geewa: algorithm did not converge"
  )
  expect_false(out$converged)
})


test_that("one-step methods are reported as converged", {
  for (method in c("bcgee-naive", "bcgee-robust", "bcgee-empirical",
                   "opgee-jeffreys", "hpgee-jeffreys")) {
    expect_no_warning(out <- finalize_fit(criterion = 1, method = method))
    expect_true(out$converged)
  }
})


test_that("the warning names the calling fit function", {
  expect_warning(
    finalize_fit(criterion = 1, fit_function = "geewa_binary"),
    "geewa_binary: algorithm did not converge"
  )
})


test_that("independence is stored as 0 for geewa and 1 for geewa_binary", {
  expect_identical(
    finalize_fit(structure = "independence", alpha = numeric(0))$alpha,
    0
  )
  expect_identical(
    finalize_fit(
      structure = "independence",
      alpha = numeric(0),
      fit_function = "geewa_binary"
    )$alpha,
    1
  )
})


test_that("association parameters are named by occasion pair or lag", {
  unstructured <- finalize_fit(
    structure = "unstructured",
    alpha = c(0.1, 0.2, 0.3)
  )
  expect_identical(names(unstructured$alpha), c("alpha_1.2", "alpha_1.3", "alpha_2.3"))

  fixed <- finalize_fit(structure = "fixed", alpha = c(0.1, 0.2, 0.3))
  expect_identical(names(fixed$alpha), c("alpha_1.2", "alpha_1.3", "alpha_2.3"))

  toeplitz <- finalize_fit(structure = "toeplitz", alpha = c(0.4, 0.2))
  expect_identical(names(toeplitz$alpha), c("alpha_lag1", "alpha_lag2"))

  m_dependent <- finalize_fit(structure = "m-dependent", alpha = 0.4)
  expect_identical(names(m_dependent$alpha), "alpha_lag1")

  exchangeable <- finalize_fit(structure = "exchangeable", alpha = 0.3)
  expect_null(names(exchangeable$alpha))
})


test_that("binomial fits warn about fitted probabilities at the boundary", {
  expect_warning(
    finalize_fit(fitted = c(0.2, 0.5, 1, 0.3, 0.6, 0.7)),
    "geewa: fitted probabilities numerically 0 or 1 occurred"
  )
  expect_warning(
    finalize_fit(fitted = c(0.2, 0.5, 0.8, 0, 0.6, 0.7)),
    "fitted probabilities numerically 0 or 1 occurred"
  )
  expect_no_warning(finalize_fit(fitted = c(0.01, 0.5, 0.99, 0.3, 0.6, 0.7)))
})


test_that("poisson fits warn about fitted rates at the boundary", {
  expect_warning(
    finalize_fit(
      family = stats::poisson(),
      fitted = c(2, 0, 3, 1, 4, 2)
    ),
    "geewa: fitted rates numerically 0 occurred"
  )
  expect_no_warning(
    finalize_fit(family = stats::poisson(), fitted = c(2, 0.001, 3, 1, 4, 2))
  )
})


test_that("other families are not checked for boundary fitted values", {
  expect_no_warning(
    finalize_fit(
      family = stats::gaussian(),
      fitted = c(0, 0, 0, 1, 1, 1)
    )
  )
})


test_that("a non-converging fit warns through the public interface", {
  ## Starting from zero with a single iteration cannot reach the solution.
  expect_warning(
    geewa(
      formula = seizures ~ treatment + lnbaseline + lnage,
      data = epilepsy,
      id = id,
      family = poisson(link = "log"),
      corstr = "exchangeable",
      method = "gee",
      beta_start = rep(0, 4),
      control = geer_control(maxiter = 1)
    ),
    "geewa: algorithm did not converge"
  )
  expect_warning(
    geewa_binary(
      formula = ecg ~ period + treatment,
      data = cerebrovascular,
      id = id,
      link = "logit",
      orstr = "exchangeable",
      method = "gee",
      beta_start = rep(0, 3),
      control = geer_control(maxiter = 1)
    ),
    "geewa_binary: algorithm did not converge"
  )
})
