testthat::local_edition(3)

independence_formula <- status ~ baseline + treatment + gender + visit + age + center

fit_independence_pair <- function(method, link = "probit", control = NULL) {
  list(
    or = geewa_binary(
      formula = independence_formula,
      id = id,
      repeated = visit,
      data = test_data$respiratory,
      link = link,
      orstr = "independence",
      method = method,
      control = control
    ),
    cor = geewa(
      formula = independence_formula,
      id = id,
      repeated = visit,
      family = stats::binomial(link = link),
      data = test_data$respiratory,
      phi_fixed = TRUE,
      phi_value = 1,
      corstr = "independence",
      method = method,
      control = control
    )
  )
}

test_that("geewa_binary and geewa agree under independence for representative methods", {
  methods <- c(
    "gee",
    "brgee-robust",
    "bcgee-robust",
    "brgee-naive",
    "bcgee-naive",
    "brgee-empirical",
    "bcgee-empirical",
    "pgee-jeffreys"
  )

  for (method in methods) {
    pair <- fit_independence_pair(method = method, link = "probit")

    expect_s3_class(pair$or, "geer")
    expect_s3_class(pair$cor, "geer")
    expect_equal(
      coef(pair$or),
      coef(pair$cor),
      tolerance = 1e-5,
      label = method
    )
  }
})

test_that("geewa_binary and geewa agree on naive covariance under independence for gee", {
  pair <- fit_independence_pair(method = "gee", link = "probit")

  expect_equal(
    vcov(pair$or, cov_type = "naive"),
    vcov(pair$cor, cov_type = "naive"),
    tolerance = 1e-5
  )
})

test_that("geewa_binary and geewa pgee-jeffreys agree under independence for logit", {
  pair <- fit_independence_pair(
    method = "pgee-jeffreys",
    link = "logit",
    control = list(jeffreys_power = 0.5, tolerance = 1e-12)
  )

  expect_equal(coef(pair$or), coef(pair$cor), tolerance = 1e-5)
})
