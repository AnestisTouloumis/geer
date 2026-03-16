local({
  data("cerebrovascular", package = "geer")
  fit_gee_probit <- geewa(
    formula = ecg ~ treatment + factor(period),
    id = id,
    family = binomial(link = "probit"),
    phi_fixed = TRUE,
    phi_value = 1,
    data = cerebrovascular,
    corstr = "independence",
    method = "gee"
  )
  fit_glm_probit <- glm(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "probit"),
    data = cerebrovascular
  )
  fit_brgee_naive_probit <- update(fit_gee_probit, method = "brgee-naive")
  fit_brgee_naive_logit <- update(
    fit_gee_probit,
    family = binomial(link = "logit"),
    method = "brgee-naive"
  )


  test_that("gee matches glm (probit): coef and naive vcov", {
    expect_equal(coef(fit_gee_probit), coef(fit_glm_probit), tolerance = 1e-5)
    expect_equal(vcov(fit_gee_probit, cov_type = "naive"), vcov(fit_glm_probit), tolerance = 1e-5)
  })


  test_that("brgee-naive matches brglm2 Firth (AS_mean) for probit: coef and vcov", {
    skip_if_not_installed("brglm2")
    fit_glm_firth_probit <- update(
      fit_glm_probit,
      method = brglmFit,
      type = "AS_mean"
    )
  expect_equal(coef(fit_brgee_naive_probit), coef(fit_glm_firth_probit), tolerance = 1e-5)
    expect_equal(
      vcov(fit_brgee_naive_probit, cov_type = "naive"),
      vcov(fit_glm_firth_probit),
      tolerance = 1e-5
    )
  })


  test_that("brgee-naive matches brglm2 Firth (AS_mean) for logit: coef and vcov", {
    skip_if_not_installed("brglm2")
    fit_glm_firth_logit <- update(
      fit_glm_probit,
      family = binomial(link = "logit"),
      method = brglmFit,
      type = "AS_mean"
    )
    expect_equal(coef(fit_brgee_naive_logit), coef(fit_glm_firth_logit), tolerance = 1e-5)
    expect_equal(
      vcov(fit_brgee_naive_logit, cov_type = "naive"),
      vcov(fit_glm_firth_logit),
      tolerance = 1e-5
    )
  })
})
