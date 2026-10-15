cerebrovascular <- test_data$cerebrovascular
respiratory <- test_data$respiratory

fit_gee_probit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  link = "probit",
  data = cerebrovascular,
  id = id,
  orstr = "independence",
  method = "gee"
)

fit_glm_probit <- glm(
  formula = ecg ~ treatment + factor(period),
  family = binomial(link = "probit"),
  data = cerebrovascular
)

fit_gee_logit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  link = "logit",
  data = cerebrovascular,
  orstr = "independence",
  id = id,
  method = "gee"
)

fit_glm_logit <- glm(
  formula = ecg ~ treatment + factor(period),
  family = binomial(link = "logit"),
  data = cerebrovascular
)

fit_brgee_naive_logit <- update(
  fit_gee_logit,
  method = "brgee-naive"
)

fit_pgee_jeffreys_logit <- update(
  fit_brgee_naive_logit,
  method = "pgee-jeffreys"
)

fit_brgee_naive_probit <- update(
  fit_gee_probit,
  method = "brgee-naive"
)

fit_pgee_jeffreys_probit <- update(
  fit_brgee_naive_probit,
  method = "pgee-jeffreys"
)

independence_formula <- status ~ baseline + treatment + gender + visit + age + center

fit_independence_pair <- function(method, link = "probit", control = NULL) {
  list(
    or = geewa_binary(
      formula = independence_formula,
      id = id,
      repeated = visit,
      data = respiratory,
      link = link,
      orstr = "independence",
      method = method,
      control = control
    ),
    cor = geewa(
      formula = independence_formula,
      id = id,
      repeated = visit,
      family = binomial(link = link),
      data = respiratory,
      phi_fixed = TRUE,
      phi_value = 1,
      corstr = "independence",
      method = method,
      control = control
    )
  )
}

test_that("geewa_binary under independence matches the corresponding GLM", {
  skip_on_cran()
  expect_equal(coef(fit_gee_probit), coef(fit_glm_probit), tolerance = 1e-5)
  expect_equal(
    vcov(fit_gee_probit, cov_type = "naive"),
    vcov(fit_glm_probit),
    tolerance = 1e-5
  )
  expect_equal(coef(fit_gee_logit), coef(fit_glm_logit), tolerance = 1e-5)
  expect_equal(
    vcov(fit_gee_logit, cov_type = "naive"),
    vcov(fit_glm_logit),
    tolerance = 1e-5
  )
})

test_that("geewa_binary brgee-naive matches brglm2 AS_mean for representative links", {
  skip_on_cran()
  fit_glm_as_mean_probit <- update(
    fit_glm_probit,
    method = brglm2::brglmFit,
    type = "AS_mean"
  )
  fit_glm_as_mean_logit <- update(
    fit_glm_logit,
    method = brglm2::brglmFit,
    type = "AS_mean"
  )
  expect_equal(
    coef(fit_brgee_naive_probit),
    coef(fit_glm_as_mean_probit),
    tolerance = 1e-5
  )
  expect_equal(
    vcov(fit_brgee_naive_probit, cov_type = "naive"),
    vcov(fit_glm_as_mean_probit),
    tolerance = 1e-5
  )
  expect_equal(
    coef(fit_brgee_naive_logit),
    coef(fit_glm_as_mean_logit),
    tolerance = 1e-5
  )
  expect_equal(
    vcov(fit_brgee_naive_logit, cov_type = "naive"),
    vcov(fit_glm_as_mean_logit),
    tolerance = 1e-5
  )
})

test_that("pgee-jeffreys matches brgee-naive for logit but not for probit in geewa_binary", {
  skip_on_cran()
  expect_equal(
    coef(fit_brgee_naive_logit),
    coef(fit_pgee_jeffreys_logit),
    tolerance = 1e-10
  )
  expect_false(isTRUE(all.equal(
    coef(fit_brgee_naive_probit),
    coef(fit_pgee_jeffreys_probit)
  )))
})

test_that("pgee-jeffreys matches brglm2 Jeffreys fits for representative links in geewa_binary", {
  skip_on_cran()
  fit_glm_jeffreys_logit <- update(
    fit_glm_logit,
    method = brglm2::brglmFit,
    type = "MPL_Jeffreys"
  )
  fit_glm_jeffreys_probit <- update(
    fit_glm_probit,
    method = brglm2::brglmFit,
    type = "MPL_Jeffreys"
  )
  expect_equal(
    coef(fit_pgee_jeffreys_logit),
    coef(fit_glm_jeffreys_logit),
    tolerance = 1e-5
  )
  expect_equal(
    vcov(fit_pgee_jeffreys_logit, cov_type = "naive"),
    vcov(fit_glm_jeffreys_logit),
    tolerance = 1e-5
  )
  expect_equal(
    coef(fit_pgee_jeffreys_probit),
    coef(fit_glm_jeffreys_probit),
    tolerance = 1e-5
  )
  expect_equal(
    vcov(fit_pgee_jeffreys_probit, cov_type = "naive"),
    vcov(fit_glm_jeffreys_probit),
    tolerance = 1e-5
  )
})

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
