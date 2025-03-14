data("cerebrovascular")
fit_gee <- geewa(formula = I(ecg == "normal") ~ treatment + factor(period),
                 id = id,
                 family = binomial(link = "probit"),
                 phi_fixed = TRUE,
                 phi_value = 1,
                 data = cerebrovascular,
                 corstr = "independence",
                 method = "gee")
fit_glm <- glm(formula = I(ecg == "normal") ~ treatment + factor(period),
               family = binomial(link = "probit"),
               data = cerebrovascular)


test_that("gee-glm", {
  expect_equal(coef(fit_gee),
               coef(fit_glm),
               tolerance = 1e-5)
  expect_equal(vcov(fit_gee,
                    cov_type = "naive"),
               vcov(fit_glm),
               tolerance = 1e-5)
})


test_that("brgee_naive-brglm_firth - probit", {
  fit_brgee_naive <- update(fit_gee,
                            method = "brgee_naive")
  fit_glm_firth <- update(fit_glm,
                          method = brglmFit,
                          type = "AS_mean")

  expect_equal(coef(fit_brgee_naive, cov_type = "naive"),
               coef(fit_glm_firth),
               tolerance = 1e-5)
  expect_equal(vcov(fit_brgee_naive, cov_type = "naive"),
               vcov(fit_glm_firth),
               tolerance = 1e-5)
})


test_that("brgee_naive-brglm_firth - logit", {
  fit_brgee_naive <- update(fit_gee,
                            family = binomial(link = "logit"),
                            method = "brgee_naive")
  fit_glm_firth <- update(fit_glm,
                          family = binomial(link = "logit"),
                          method = brglmFit,
                          type = "AS_mean")

  expect_equal(coef(fit_brgee_naive, cov_type = "naive"),
               coef(fit_glm_firth),
               tolerance = 1e-5)
  expect_equal(vcov(fit_brgee_naive, cov_type = "naive"),
               vcov(fit_glm_firth),
               tolerance = 1e-5)
})
