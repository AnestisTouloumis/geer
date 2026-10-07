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
                         iter = 1L,
                         failure = NULL) {
  geer:::finalize_geer_fit(
    fit = list(iter = iter, fitted.values = fitted),
    geesolver_fit = list(criterion = criterion, alpha = alpha,
                         failure = failure),
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


test_that("the bias-corrected covariance is NA and warned about when clusters do not exceed parameters", {
  few_clusters <- data.frame(
    id = rep(1:2, each = 5),
    y = c(1, 2, 0, 3, 1, 2, 1, 4, 0, 2),
    x1 = rep(c(-2, -1, 0, 1, 2), times = 2),
    x2 = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
  )
  expect_warning(
    fit <- geewa(
      formula = y ~ x1 + x2,
      family = poisson(link = "log"),
      id = id,
      data = few_clusters,
      corstr = "independence",
      method = "gee"
    ),
    "geewa: the bias-corrected covariance matrix is undefined"
  )
  expect_true(all(is.na(fit$bias_corrected_covariance)))
  expect_true(all(is.na(vcov(fit))))
  expect_true(all(is.finite(fit$robust_covariance)))
  expect_true(all(is.finite(vcov(fit, cov_type = "robust"))))
})


test_that("the bias-corrected covariance is available when clusters exceed parameters", {
  fit <- geewa(
    formula = ecg ~ period + treatment,
    family = binomial(link = "logit"),
    id = id,
    data = test_data$cerebrovascular,
    corstr = "exchangeable",
    method = "gee"
  )
  expect_true(all(is.finite(fit$bias_corrected_covariance)))
})


collect_warnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = messages)
}


test_that("a solver failure marks the fit as not converged and warns with the reason", {
  reason <- "update_beta_gee_cc: cluster 3 (id = 3): singular"
  res <- collect_warnings(
    finalize_fit(criterion = c(1, Inf), iter = 2L, failure = reason)
  )
  expect_false(res$value$converged)
  expect_length(res$messages, 2L)
  expect_match(
    res$messages[[1L]],
    paste0("geewa: the fitting algorithm stopped early because of a numerical ",
           "failure (", reason, "); returning the estimates from the last ",
           "accepted iteration"),
    fixed = TRUE
  )
  expect_match(res$messages[[2L]], "geewa: algorithm did not converge",
               fixed = TRUE)
})


test_that("a solver failure overrides the converged flag of one-step methods", {
  res <- collect_warnings(
    finalize_fit(criterion = Inf, method = "bcgee-robust",
                 failure = "singular", fit_function = "geewa_binary")
  )
  expect_false(res$value$converged)
  expect_match(res$messages[[1L]],
               "geewa_binary: the fitting algorithm stopped early")
  expect_match(res$messages[[1L]],
               "returning the preliminary estimates of the first stage",
               fixed = TRUE)
})


test_that("an empty failure message leaves the convergence logic unchanged", {
  expect_no_warning(out <- finalize_fit(failure = ""))
  expect_true(out$converged)
})


test_that("a damped step_multiplier reaches the same solution", {
  fit_default <- geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee",
    control = geer_control(tolerance = 1e-8)
  )
  fit_damped <- geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data = epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee",
    control = geer_control(tolerance = 1e-8, step_multiplier = 0.5)
  )
  expect_true(fit_damped$converged)
  expect_equal(fit_damped$coefficients, fit_default$coefficients,
               tolerance = 1e-6)
})
