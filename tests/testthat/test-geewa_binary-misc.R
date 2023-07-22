data("cerebrovascular")
fit_gee <- geewa_binary(formula = I(ecg == "normal") ~ treatment + factor(period),
                        link =  "probit",
                        data = cerebrovascular,
                        or_structure = "independence",
                        method = "gee")
fit_glm <- glm(formula = I(ecg == "normal") ~ treatment + factor(period),
               family = binomial(link = "probit"),
               data = cerebrovascular)


test_that("gee-glm", {
  expect_equal(coef(fit_gee),
               coef(fit_glm),
               tolerance = 1e-5)
  expect_equal(vcov(fit_gee,
                    type = "naive"),
               vcov(fit_glm),
               tolerance = 1e-5)
})


test_that("brgee_naive-brglm_firth - probit", {
  fit_brgee_naive <- update(fit_gee,
                            method = "brgee_naive")
  fit_glm_firth <- update(fit_glm,
                          method = brglmFit,
                          type = "AS_mean")

  expect_equal(coef(fit_brgee_naive),
               coef(fit_glm_firth),
               tolerance = 1e-5)
  expect_equal(vcov(fit_brgee_naive, type = "naive"),
               vcov(fit_glm_firth),
               tolerance = 1e-5)
})


test_that("brgee_naive-brglm_firth - logit", {
  fit_brgee_naive <- update(fit_gee,
                            link = "logit",
                            method = "brgee_naive")
  fit_glm_firth <- update(fit_glm,
                          family = binomial(link = "logit"),
                          method = brglmFit,
                          type = "AS_mean")

  expect_equal(coef(fit_brgee_naive),
               coef(fit_glm_firth),
               tolerance = 1e-5)
  expect_equal(vcov(fit_brgee_naive, type = "naive"),
               vcov(fit_glm_firth),
               tolerance = 1e-5)
})
