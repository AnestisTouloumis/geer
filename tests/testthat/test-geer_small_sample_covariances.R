expect_small_sample_reference <- function(fit) {
  reference <- ref_small_sample_covariances(fit)
  for (cov_type in geer_small_sample_cov_types) {
    expect_equal(
      unname(vcov(fit, cov_type = cov_type)),
      unname(reference[[cov_type]]),
      tolerance = 1e-10,
      label = cov_type
    )
  }
}

test_that("small-sample covariances agree with a direct implementation for geewa fits", {
  expect_small_sample_reference(geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    repeated = visit,
    family = poisson(link = "log"),
    corstr = "ar1"
  ))
  expect_small_sample_reference(geewa(
    formula = seizures ~ treatment + lnbaseline,
    data = test_data$epilepsy,
    id = id,
    repeated = visit,
    family = poisson(link = "sqrt"),
    corstr = "exchangeable",
    method = "brgee-robust"
  ))
  d <- ref_sim("cc", rep(list(1:4), 30), family = Gamma(), link = "log",
               weights = c(0.5, 2, 1, 1.5, 0.8))
  expect_small_sample_reference(geewa(
    formula = y ~ x1 + x2,
    data = d$df,
    id = id,
    repeated = visit,
    family = Gamma(link = "log"),
    corstr = "unstructured",
    weights = w
  ))
})

test_that("small-sample covariances agree with a direct implementation for OR fits", {
  expect_small_sample_reference(fit_geewa_bin_exch)
  expect_small_sample_reference(geewa_binary(
    formula = ecg ~ period + treatment,
    data = test_data$cerebrovascular,
    id = id,
    link = "probit",
    orstr = "independence",
    method = "pgee-jeffreys"
  ))
})

test_that("pooled small-sample covariances stop with an unbalanced design", {
  d <- ref_sim("cc", rep(list(1:4), 30), family = gaussian(),
               link = "identity")
  fit <- geewa(
    formula = y ~ x1 + x2,
    data = d$df[-c(3, 10, 50), ],
    id = id,
    repeated = visit,
    family = gaussian(link = "identity"),
    corstr = "ar1"
  )
  for (cov_type in geer_pooled_cov_types) {
    expect_error(
      vcov(fit, cov_type = cov_type),
      sprintf(
        paste0(
          "vcov failed: the %s covariance requires a balanced design, in ",
          "which every cluster observes each repeated position once"
        ),
        cov_type
      ),
      fixed = TRUE
    )
  }
  reference <- ref_small_sample_covariances(fit)
  for (cov_type in setdiff(geer_small_sample_cov_types, geer_pooled_cov_types)) {
    expect_equal(
      unname(vcov(fit, cov_type = cov_type)),
      unname(reference[[cov_type]]),
      tolerance = 1e-10,
      label = cov_type
    )
  }
})

test_that("gosho-sato-takeuchi stops when clusters do not exceed coefficients", {
  d <- ref_sim("cc", rep(list(1:4), 3), family = gaussian(),
               link = "identity")
  fit <- suppressWarnings(geewa(
    formula = y ~ x1 + x2,
    data = d$df,
    id = id,
    repeated = visit,
    family = gaussian(link = "identity"),
    corstr = "exchangeable"
  ))
  expect_error(
    vcov(fit, cov_type = "gosho-sato-takeuchi"),
    paste0(
      "vcov failed: the gosho-sato-takeuchi covariance requires more ",
      "clusters than coefficients"
    ),
    fixed = TRUE
  )
  expect_true(all(is.finite(vcov(fit, cov_type = "pan"))))
})
