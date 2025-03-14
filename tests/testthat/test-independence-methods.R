link <- sample(c("logit", "probit", "cloglog", "cauchit"), 1)
jeffreys_power <- sample(seq(0.1, 0.6, 0.1), 1)
print(c(link, jeffreys_power))
data("respiratory")
library("brglm2")

fmla <- y ~  visit + age + baseline + center + treatment + gender

fit_model_or <-
  geewa_binary(formula = fmla,
               link =  link,
               id = id,
               repeated = visit,
               data = respiratory,
               orstr = "independence",
               method = "gee",
               control = list(jeffreys_power = jeffreys_power,
                              tolerance = 1e-12,
                              step_maxiter = 12))

fit_model_cor <- geewa(formula = fmla,
                       id = id,
                       repeated = visit,
                       family = binomial(link = link),
                       data = respiratory,
                       phi_fixed = TRUE,
                       phi_value = 1,
                       corstr = "independence",
                       method = "gee",
                       control = list(jeffreys_power = jeffreys_power,
                                      tolerance = 1e-12,
                                      step_maxiter = 12))

fit_model_glm <- glm(formula = fmla,
                     family = binomial(link = link),
                     data = respiratory,
                     method = "brglmFit",
                     control = list(type = "ML",
                                    epsilon = 1e-12))


test_that("gee", {
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})

test_that("gee - glm", {
  expect_equal(coef(fit_model_or),
               coef(fit_model_glm))
})

test_that("brgee_robust", {
  fit_model_or <- update(fit_model_or,
                         method = "brgee_robust")
  fit_model_cor <- update(fit_model_cor,
                          method = "brgee_robust")
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})


test_that("brgee_naive", {
  fit_model_or <- update(fit_model_or,
                         method = "brgee_naive")
  fit_model_cor <- update(fit_model_cor,
                          method = "brgee_naive")
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})

test_that("brgee_naive - firth", {
  fit_model_or <- update(fit_model_or,
                         method = "brgee_naive")
  fit_model_glm <- update(fit_model_glm,
                          control = list(type = "AS_mean",
                                         epsilon = 1e-12))
  expect_equal(coef(fit_model_or),
               coef(fit_model_glm))
})

test_that("brgee_empirical", {
  fit_model_or <- update(fit_model_or,
                         method = "brgee_empirical")
  fit_model_cor <- update(fit_model_cor,
                          method = "brgee_empirical")
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})

test_that("bcgee_robust", {
  fit_model_or <- update(fit_model_or,
                         method = "bcgee_robust")
  fit_model_cor <- update(fit_model_cor,
                          method = "bcgee_robust")
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})

test_that("bcgee_naive", {
  fit_model_or <- update(fit_model_or,
                         method = "bcgee_naive")
  fit_model_cor <- update(fit_model_cor,
                          method = "bcgee_naive")
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})


test_that("bcgee_empirical", {
  fit_model_or <- update(fit_model_or,
                         method = "bcgee_empirical")
  fit_model_cor <- update(fit_model_cor,
                          method = "bcgee_empirical")
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})

test_that("pgee_jeffreys", {
  fit_model_or <- update(fit_model_or,
                         method = "pgee_jeffreys")
  fit_model_cor <- update(fit_model_cor,
                          method = "pgee_jeffreys")
  expect_equal(coef(fit_model_or),
               coef(fit_model_cor))
})

test_that("pgee_jeffreys - glm", {
  fit_model_or <- update(fit_model_or,
                         method = "pgee_jeffreys",
                         control = list(jeffreys_power = jeffreys_power,
                                        tolerance = 1e-12))
  fit_model_glm <- update(fit_model_glm,
                          control = list(type = "MPL_Jeffreys",
                                         a = jeffreys_power,
                                         epsilon = 1e-12))
  expect_equal(coef(fit_model_or),
               coef(fit_model_glm))
})


test_that("pgee_jeffreys - brgee_naive - logit", {
  fit_model_or <- update(fit_model_or,
                         link = "logit",
                         method = "brgee_naive")
  fit_model_or2 <- update(fit_model_or,
                          method = "pgee_jeffreys",
                          control = list(jeffreys_power = 0.5,
                                         tolerance = 1e-12))
  expect_equal(coef(fit_model_or),
               coef(fit_model_or2))
})
