data("cerebrovascular", package = "geer")
fitor1 <- geewa_binary(formula = ecg ~ treatment + offset(period),
                       link =  "logit",
                       id = id,
                       data = cerebrovascular,
                       orstr = "independence",
                       method = "gee")
fitor2 <- geewa_binary(formula = ecg ~ treatment,
                       link =  "logit",
                       id = id,
                       offset = period,
                       data = cerebrovascular,
                       orstr = "independence",
                       method = "gee")
fitcc1 <- geewa(formula = ecg ~ treatment + offset(period),
                family = binomial(link =  "logit"),
                id = id,
                data = cerebrovascular,
                corstr = "independence",
                method = "gee",
                phi_fixed = TRUE,
                phi_value = 1)
fitcc2 <- geewa(formula = ecg ~ treatment,
                family = binomial(link =  "logit"),
                id = id,
                offset = period,
                data = cerebrovascular,
                corstr = "independence",
                method = "gee",
                phi_fixed = TRUE,
                phi_value = 1)


test_that("bcgee-offset", {
  expect_equal(coef(fitor1),
               coef(fitor2))
  expect_equal(coef(fitcc1),
               coef(fitcc2))
  expect_equal(coef(fitor1),
               coef(fitcc1))
})
