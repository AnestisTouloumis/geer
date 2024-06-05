data("cerebrovascular")
fit_gee <- geewa(formula = I(ecg == "normal") ~ treatment + factor(period),
                 id = id,
                 family = binomial(link = "logit"),
                 phi_fixed = TRUE,
                 phi_value = 1,
                 data = cerebrovascular,
                 correlation_structure = "independence",
                 method = "gee")

fit_glm <- glm(formula = I(ecg == "normal") ~ treatment + factor(period),
               family = binomial(link = "logit"),
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



test_that("jeffreys = naive - logit", {
  fit_brgee_naive <- update(fit_gee,
                            method = "brgee_naive")
  fit_pgee_jeffreys <- update(fit_brgee_naive,
                              method = "pgee_jeffreys")
  expect_equal(coef(fit_brgee_naive),
               coef(fit_pgee_jeffreys))
})


test_that("jeffreys != naive - probit", {
  fit_brgee_naive_probit <- update(fit_gee,
                                   family = binomial(link = "probit"),
                                   method = "brgee_naive")
  fit_pgee_jeffreys_probit <- update(fit_brgee_naive_probit,
                                     method = "pgee_jeffreys")
  expect_false(all((coef(fit_brgee_naive_probit) - coef(fit_pgee_jeffreys_probit)) == 0))
})


test_that("jeffreys = glm - logit", {
  fit_glm_jeffreys <-  update(fit_glm,
                              method = brglmFit,
                              type = "MPL_Jeffreys")
  fit_pgee_jeffreys <- update(fit_gee,
                              method = "pgee_jeffreys")
  expect_equal(coef(fit_glm_jeffreys),
               coef(fit_pgee_jeffreys),
               tolerance = 1e-5)
  expect_equal(vcov(fit_glm_jeffreys),
               vcov(fit_pgee_jeffreys, cov_type = "naive"),
               tolerance = 1e-5)
})

test_that("jeffreys = glm - probit", {
  fit_glm_jeffreys_probit <-  update(fit_glm,
                                     method = brglmFit,
                                     family = binomial(link = "probit"),
                                     type = "MPL_Jeffreys")
  fit_pgee_jeffreys_probit <- update(fit_gee,
                                     family = binomial(link = "probit"),
                                     method = "pgee_jeffreys")
  expect_equal(coef(fit_glm_jeffreys_probit),
               coef(fit_pgee_jeffreys_probit),
               tolerance = 1e-5)
  expect_equal(vcov(fit_glm_jeffreys_probit),
               vcov(fit_pgee_jeffreys_probit, cov_type = "naive"),
               tolerance = 1e-5)
})
