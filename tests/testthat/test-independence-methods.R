local({
  skip_if_not_installed("brglm2")
  data("respiratory", package = "geer")
  link <- "probit"
  jeffreys_power <- 0.5
  fmla <- status ~ visit + age + baseline + center + treatment + gender
  fit_model_or <- geewa_binary(
    formula = fmla,
    link = link,
    id = id,
    repeated = visit,
    data = respiratory,
    orstr = "independence",
    method = "gee"
  )
  fit_model_cor <- geewa(
    formula = fmla,
    id = id,
    repeated = visit,
    family = binomial(link = link),
    data = respiratory,
    phi_fixed = TRUE,
    phi_value = 1,
    corstr = "independence",
    method = "gee"
  )
  fit_model_glm <- glm(
    formula = fmla,
    family = binomial(link = link),
    data = respiratory,
    method = brglmFit,
    control = list(type = "ML", epsilon = 1e-12)
  )
  fit_or_brgee_robust <- update(fit_model_or, method = "brgee-robust")
  fit_cor_brgee_robust <- update(fit_model_cor, method = "brgee-robust")
  fit_or_brgee_naive <- update(fit_model_or, method = "brgee-naive")
  fit_cor_brgee_naive <- update(fit_model_cor, method = "brgee-naive")
  fit_or_brgee_emp <- update(fit_model_or, method = "brgee-empirical")
  fit_cor_brgee_emp <- update(fit_model_cor, method = "brgee-empirical")
  fit_or_bcgee_robust <- update(fit_model_or, method = "bcgee-robust")
  fit_cor_bcgee_robust <- update(fit_model_cor, method = "bcgee-robust")
  fit_or_bcgee_naive <- update(fit_model_or, method = "bcgee-naive")
  fit_cor_bcgee_naive <- update(fit_model_cor, method = "bcgee-naive")
  fit_or_bcgee_emp <- update(fit_model_or, method = "bcgee-empirical")
  fit_cor_bcgee_emp <- update(fit_model_cor, method = "bcgee-empirical")
  fit_or_pgee <- update(fit_model_or, method = "pgee-jeffreys")
  fit_cor_pgee <- update(fit_model_cor, method = "pgee-jeffreys")
  fit_glm_asmean <- update(fit_model_glm, control = list(type = "AS_mean"))
  fit_or_pgee_pow <- update(
    fit_model_or,
    method = "pgee-jeffreys",
    control = list(jeffreys_power = jeffreys_power, tolerance = 1e-12)
  )
  fit_glm_jeff <- update(
    fit_model_glm,
    control = list(type = "MPL_Jeffreys", a = jeffreys_power, epsilon = 1e-12)
  )
  fit_or_logit_naive <- update(fit_model_or, link = "logit", method = "brgee-naive")
  fit_or_logit_pgee <- update(
    fit_or_logit_naive,
    method = "pgee-jeffreys",
    control = list(jeffreys_power = 0.5, tolerance = 1e-12)
  )


  test_that("respiratory: method equivalences (independence), probit link", {
    expect_equal(coef(fit_model_or), coef(fit_model_cor), tolerance = 1e-5)
    expect_equal(coef(fit_model_or), coef(fit_model_glm), tolerance = 1e-5)
    expect_equal(coef(fit_or_brgee_robust), coef(fit_cor_brgee_robust), tolerance = 1e-5)
    expect_equal(coef(fit_or_brgee_naive), coef(fit_cor_brgee_naive), tolerance = 1e-5)
    expect_equal(coef(fit_or_brgee_naive), coef(fit_glm_asmean), tolerance = 1e-5)
    expect_equal(coef(fit_or_brgee_emp), coef(fit_cor_brgee_emp), tolerance = 1e-5)
    expect_equal(coef(fit_or_bcgee_robust), coef(fit_cor_bcgee_robust), tolerance = 1e-5)
    expect_equal(coef(fit_or_bcgee_naive), coef(fit_cor_bcgee_naive), tolerance = 1e-5)
    expect_equal(coef(fit_or_bcgee_emp), coef(fit_cor_bcgee_emp), tolerance = 1e-5)
    expect_equal(coef(fit_or_pgee), coef(fit_cor_pgee), tolerance = 1e-5)
    expect_equal(coef(fit_or_pgee_pow), coef(fit_glm_jeff), tolerance = 1e-5)
    expect_equal(coef(fit_or_logit_naive), coef(fit_or_logit_pgee), tolerance = 1e-5)
  })
})
