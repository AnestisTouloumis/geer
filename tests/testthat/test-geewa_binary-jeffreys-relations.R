local({
  data("cerebrovascular", package = "geer")
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
  fit_brgee_naive_logit <- update(fit_gee_logit, method = "brgee-naive")
  fit_pgee_jeffreys_logit <- update(fit_brgee_naive_logit, method = "pgee-jeffreys")
  fit_brgee_naive_probit <- update(fit_gee_logit, link = "probit", method = "brgee-naive")
  fit_pgee_jeffreys_probit <- update(fit_brgee_naive_probit, method = "pgee-jeffreys")


  test_that("GEE independence matches GLM (logit): coefficients and naive vcov", {
    expect_equal(coef(fit_gee_logit), coef(fit_glm_logit), tolerance = 1e-5)
    expect_equal(vcov(fit_gee_logit, cov_type = "naive"), vcov(fit_glm_logit), tolerance = 1e-5)
  })


  test_that("Jeffreys penalty: pgee-jeffreys equals brgee-naive for logit", {
    expect_equal(coef(fit_brgee_naive_logit), coef(fit_pgee_jeffreys_logit), tolerance = 1e-10)
  })


  test_that("Jeffreys penalty: pgee-jeffreys differs from brgee-naive for probit", {
    expect_false(
      isTRUE(all.equal(coef(fit_brgee_naive_probit), coef(fit_pgee_jeffreys_probit)))
    )
  })


  test_that("Jeffreys penalty matches brglm2 GLM Jeffreys (logit): coef and vcov", {
    skip_if_not_installed("brglm2")

    fit_glm_jeffreys_logit <- update(
      fit_glm_logit,
      method = brglmFit,
      type = "MPL_Jeffreys"
    )
    expect_equal(coef(fit_glm_jeffreys_logit), coef(fit_pgee_jeffreys_logit), tolerance = 1e-5)
    expect_equal(
      vcov(fit_glm_jeffreys_logit),
      vcov(fit_pgee_jeffreys_logit, cov_type = "naive"),
      tolerance = 1e-5
    )
  })


  test_that("Jeffreys penalty matches brglm2 GLM Jeffreys (probit): coef and vcov", {
    skip_if_not_installed("brglm2")
    fit_glm_jeffreys_probit <- update(
      fit_glm_logit,
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
