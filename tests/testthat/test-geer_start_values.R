start_model_matrix <- function() {
  model.matrix(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy
  )
}

start_values <- function(beta_start = NULL,
                         family = poisson(),
                         link = "log",
                         method = "gee",
                         y = test_data$epilepsy$seizures,
                         offset = rep(0, nrow(test_data$epilepsy))) {
  compute_geer_start_values(
    model_matrix = start_model_matrix(),
    y = y,
    family = family,
    weights = rep(1, nrow(test_data$epilepsy)),
    offset = offset,
    method = method,
    link = link,
    beta_start = beta_start,
    control = geer_control(),
    control_glm = list()
  )
}

binary_start_values <- function(
    beta_start = NULL,
    method = "gee",
    offset = rep(0, nrow(test_data$cerebrovascular))) {
  cer <- test_data$cerebrovascular
  compute_binary_start_values(
    model_matrix = model.matrix(ecg ~ period + treatment, data = cer),
    y = cer$ecg,
    family = binomial(),
    weights = rep(1, nrow(cer)),
    offset = offset,
    method = method,
    beta_start = beta_start,
    control_glm = list(),
    tolerance = 1e-8,
    jeffreys_power = 0.5
  )
}

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

test_that("starting values for the GEE method are the maximum-likelihood fit", {
  reference <- glm(
    seizures ~ treatment + lnbaseline + lnage,
    family = poisson(),
    data = test_data$epilepsy
  )

  out <- start_values()

  expect_equal(unname(out), unname(coef(reference)), tolerance = 1e-5)
})

test_that("quasi families and the identity link start from stats::glm.fit", {
  reference <- lm(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy
  )

  out <- start_values(family = gaussian(), link = "identity")
  expect_equal(unname(out), unname(coef(reference)), tolerance = 1e-6)

  quasi_out <- start_values(family = quasipoisson())
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
  bad_offset <- rep(Inf, nrow(test_data$epilepsy))
  expect_error(
    start_values(offset = bad_offset),
    "starting values could not be computed; supply 'beta_start'"
  )
  expect_error(
    start_values(
      family = gaussian(),
      link = "identity",
      offset = bad_offset
    ),
    "starting values could not be computed; supply 'beta_start'"
  )
  expect_error(
    binary_start_values(offset = rep(Inf, nrow(test_data$cerebrovascular))),
    "starting values could not be computed; supply 'beta_start'"
  )
})

test_that("the starting-value fit does not emit per-observation non-integer warnings", {
  pois <- collect_warnings(
    geewa(seizures ~ treatment + lnbaseline,
      family = poisson(), data = test_data$epilepsy,
      id = id, corstr = "exchangeable"
    )
  )
  expect_false(any(grepl("non-integer", pois$warnings, fixed = TRUE)))
  expect_length(pois$warnings, 0L)
  bin <- collect_warnings(
    geewa(status ~ treatment + baseline,
      family = binomial(), data = test_data$respiratory,
      id = id, repeated = visit, corstr = "exchangeable"
    )
  )
  expect_false(any(grepl("non-integer", bin$warnings, fixed = TRUE)))
  expect_length(bin$warnings, 0L)
  bin_or <- collect_warnings(
    geewa_binary(status ~ treatment + baseline,
      data = test_data$respiratory,
      id = id, repeated = visit, orstr = "exchangeable"
    )
  )
  expect_length(bin_or$warnings, 0L)
})

test_that("evaluate_noninteger_quiet keeps every other warning", {
  seen <- character()
  withCallingHandlers(
    evaluate_noninteger_quiet({
      warning("non-integer x = 5.5")
      warning("non-integer #successes in a binomial glm!")
      warning("fitted probabilities numerically 0 or 1 occurred")
      warning("algorithm did not converge")
      1
    }),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(seen, c(
    "fitted probabilities numerically 0 or 1 occurred",
    "algorithm did not converge"
  ))
})

test_that("build_family_no_aic only replaces the AIC of the family used for starting values", {
  fam <- poisson()
  fam2 <- build_family_no_aic(fam)
  expect_identical(fam2$family, fam$family)
  expect_identical(fam2$link, fam$link)
  expect_identical(fam2$variance, fam$variance)
  expect_identical(fam2$linkinv, fam$linkinv)
  expect_true(is.na(fam2$aic(1, 1, 1, 1, 0)))
  x <- model.matrix(~ treatment + lnbaseline, test_data$epilepsy)
  y <- test_data$epilepsy$seizures
  w <- rep(1, length(y))
  b1 <- suppressWarnings(brglm2::brglmFit(
    x = x, y = y, family = fam, weights = w,
    offset = rep(0, length(y)),
    control = list(type = "AS_mean")
  ))$coefficients
  b2 <- brglm2::brglmFit(
    x = x, y = y, family = fam2, weights = w,
    offset = rep(0, length(y)),
    control = list(type = "AS_mean")
  )$coefficients
  expect_identical(b1, b2)
})

test_that("start value helpers validate beta_start and control_glm", {
  expect_error(
    check_beta_start(c("a", "b"), 2L),
    "'beta_start' must be a numeric vector of length 2"
  )
  expect_error(
    check_beta_start(c(1, Inf), 2L),
    "'beta_start' must contain only finite values"
  )
  expect_identical(check_beta_start(1:2, 2L), c(1, 2))
  expect_error(normalize_control_glm(5), "'control_glm' must be a list")
  expect_identical(get_start_type("gee"), "ML")
  expect_identical(get_start_type("bcgee-robust"), "ML")
  expect_identical(get_start_type("brgee-naive"), "AS_mean")
  expect_identical(get_start_type("pgee-jeffreys"), "MPL_Jeffreys")
})

test_that("supplied starts bypass control_glm validation", {
  fit <- fit_geewa_pois_exch
  expect_identical(
    compute_geer_start_values(
      model_matrix = fit$x, y = fit$y, family = fit$family,
      weights = fit$prior.weights, offset = fit$offset, method = fit$method,
      link = fit$family$link, beta_start = coef(fit), control = fit$control,
      control_glm = 17
    ),
    as.numeric(coef(fit))
  )
  fit <- fit_geewa_bin_exch
  expect_identical(
    compute_binary_start_values(
      model_matrix = fit$x, y = fit$y, family = fit$family,
      weights = fit$prior.weights, offset = fit$offset, method = fit$method,
      beta_start = coef(fit), control_glm = 17,
      tolerance = fit$control$tolerance,
      jeffreys_power = fit$control$jeffreys_power
    ),
    as.numeric(coef(fit))
  )
})
