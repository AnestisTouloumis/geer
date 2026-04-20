testthat::local_edition(3)

fit_gee_probit <- geewa(
  formula = ecg ~ treatment + factor(period),
  id = id,
  family = binomial(link = "probit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = test_data$cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_glm_probit <- glm(
  formula = ecg ~ treatment + factor(period),
  family = binomial(link = "probit"),
  data = test_data$cerebrovascular
)

fit_gee_logit <- geewa(
  formula = ecg ~ treatment + factor(period),
  id = id,
  family = binomial(link = "logit"),
  phi_fixed = TRUE,
  phi_value = 1,
  data = test_data$cerebrovascular,
  corstr = "independence",
  method = "gee"
)

fit_glm_logit <- glm(
  formula = ecg ~ treatment + factor(period),
  family = binomial(link = "logit"),
  data = test_data$cerebrovascular
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
  fit_gee_logit,
  family = binomial(link = "probit"),
  method = "brgee-naive"
)

fit_pgee_jeffreys_probit <- update(
  fit_brgee_naive_probit,
  method = "pgee-jeffreys"
)


test_that("geewa under independence matches the corresponding GLM", {
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


test_that("brgee-naive matches brglm2 AS_mean for representative links", {
  skip_if_not_installed("brglm2")
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


test_that("pgee-jeffreys matches brgee-naive for logit but not for probit", {
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


test_that("pgee-jeffreys matches brglm2 Jeffreys fits for representative links", {
  skip_if_not_installed("brglm2")
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
