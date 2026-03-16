local({
  data("cerebrovascular", package = "geer")
  fit_gee <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "logit"),
    phi_fixed = TRUE,
    phi_value = 1,
    data = cerebrovascular,
    corstr = "independence",
    method = "gee"
  )
  fit_glm <- glm(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "logit"),
    data = cerebrovascular
  )
  fit_brgee_naive <- update(fit_gee, method = "brgee-naive")
  fit_pgee_jeffreys <- update(fit_brgee_naive, method = "pgee-jeffreys")
  fit_brgee_naive_probit <- update(
    fit_gee,
    family = binomial(link = "probit"),
    method = "brgee-naive"
  )
  fit_pgee_jeffreys_probit <- update(fit_brgee_naive_probit, method = "pgee-jeffreys")


  test_that("gee matches glm (logit): coefficients and naive vcov", {
    expect_equal(coef(fit_gee), coef(fit_glm), tolerance = 1e-5)
    expect_equal(vcov(fit_gee, cov_type = "naive"), vcov(fit_glm), tolerance = 1e-5)
  })


  test_that("Jeffreys equals naive for logit (pgee-jeffreys vs brgee-naive)", {
    expect_equal(coef(fit_brgee_naive), coef(fit_pgee_jeffreys), tolerance = 1e-10)
  })


  test_that("Jeffreys differs from naive for probit", {
    expect_false(
      isTRUE(all.equal(coef(fit_brgee_naive_probit), coef(fit_pgee_jeffreys_probit)))
    )
  })


  test_that("Jeffreys matches brglm2 GLM Jeffreys (logit): coef and vcov", {
    skip_if_not_installed("brglm2")
    fit_glm_jeffreys <- update(
      fit_glm,
      method = brglmFit,
      type = "MPL_Jeffreys"
    )
    expect_equal(coef(fit_glm_jeffreys), coef(fit_pgee_jeffreys), tolerance = 1e-5)
    expect_equal(
      vcov(fit_glm_jeffreys),
      vcov(fit_pgee_jeffreys, cov_type = "naive"),
      tolerance = 1e-5
    )
  })


  test_that("Jeffreys matches brglm2 GLM Jeffreys (probit): coef and vcov", {
    skip_if_not_installed("brglm2")
    fit_glm_jeffreys_probit <- update(
      fit_glm,
      family = binomial(link = "probit"),
      method = brglmFit,
      type = "MPL_Jeffreys"
    )
    expect_equal(coef(fit_glm_jeffreys_probit), coef(fit_pgee_jeffreys_probit), tolerance = 1e-5)
    expect_equal(
      vcov(fit_glm_jeffreys_probit),
      vcov(fit_pgee_jeffreys_probit, cov_type = "naive"),
      tolerance = 1e-5
    )
  })
})
