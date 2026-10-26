list2env(get_cached_test_fixtures("geewa_binary_theory"), environment())

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
    expect_equal(
      coef(pair$or),
      coef(pair$cor),
      tolerance = 1e-8,
      label = method
    )
  }
})

test_that("geewa_binary and geewa agree on naive covariance under independence for gee", {
  pair <- fit_independence_pair(method = "gee", link = "probit")
  expect_equal(
    vcov(pair$or, cov_type = "naive"),
    vcov(pair$cor, cov_type = "naive"),
    tolerance = 1e-8
  )
})

test_that("geewa_binary and geewa pgee-jeffreys agree under independence for logit", {
  pair <- fit_independence_pair(
    method = "pgee-jeffreys",
    link = "logit",
    control = list(jeffreys_power = 0.5, tolerance = 1e-12)
  )
  expect_equal(coef(pair$or), coef(pair$cor), tolerance = 1e-8)
})
