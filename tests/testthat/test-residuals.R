link <- sample(c("logit", "probit", "cloglog", "cauchit"),
               1)
association <- sample(c("independence", "exchangeable", "unstructured"),
                      1)
method_gee <-  sample(c("gee",
                        "brgee-naive", "brgee-robust", "brgee-empirical",
                        "bcgee-naive", "bcgee-robust", "bcgee-empirical",
                        "pgee-jeffreys"),
                      1)
print(c(link, association, method_gee))

data("cerebrovascular")

fit_geewa <- geewa(formula = I(ecg == "normal") ~ treatment + factor(period),
                   id = id,
                   family = binomial(link = link),
                   data = cerebrovascular,
                   corstr = association,
                   method = method_gee)
fit_geewa_binary <-
  geewa_binary(formula = I(ecg == "normal") ~ treatment + factor(period),
               id = id,
               link = link,
               data = cerebrovascular,
               orstr = association,
               method = method_gee)

response_vector <- I(cerebrovascular$ecg == "normal")
mu_vector_geewa <- fitted(fit_geewa)
mu_vector_geewa_binary <- fitted(fit_geewa_binary)
weights_geewa <- fit_geewa$prior.weights
weights_geewa_binary <- fit_geewa_binary$prior.weights

test_that("working residuals", {
  expect_equal(response_vector - mu_vector_geewa,
               residuals(fit_geewa))
  expect_equal(response_vector - mu_vector_geewa_binary,
               residuals(fit_geewa_binary))
})

test_that("pearson residuals", {
  expect_equal(c(get_pearson_residuals("binomial",
                                       response_vector,
                                       mu_vector_geewa,
                                       weights_geewa)),
               residuals(fit_geewa, type = "pearson"))
  expect_equal(c(get_pearson_residuals("binomial",
                                       response_vector,
                                       mu_vector_geewa_binary,
                                       weights_geewa_binary)),
               residuals(fit_geewa_binary, type = "pearson"))
})

test_that("deviance residuals", {
  dr_geewa1 <- sign(response_vector - mu_vector_geewa) * weights_geewa
  dr_geewa2 <- response_vector * log(response_vector / mu_vector_geewa)
  dr_geewa3 <- (1 - response_vector) * log((1 - response_vector) / (1 - mu_vector_geewa))
  deviance_residuals_geewa <-
    dr_geewa1 *
    sqrt(2 * ifelse(is.na(dr_geewa2), dr_geewa3, dr_geewa2))
  deviance_residuals_geewa[is.nan(deviance_residuals_geewa)] <- 0
  expect_equal(deviance_residuals_geewa,
               residuals(fit_geewa, type = "deviance"))
  dr_geewa_binary1 <- sign(response_vector - mu_vector_geewa_binary) * weights_geewa_binary
  dr_geewa_binary2 <- response_vector * log(response_vector / mu_vector_geewa_binary)
  dr_geewa_binary3 <- (1 - response_vector) * log((1 - response_vector) / (1 - mu_vector_geewa_binary))
  deviance_residuals_geewa_binary <-
    dr_geewa_binary1 *
    sqrt(2 * ifelse(is.na(dr_geewa_binary2), dr_geewa_binary3, dr_geewa_binary2))
  deviance_residuals_geewa_binary[is.nan(deviance_residuals_geewa_binary)] <- 0
  expect_equal(deviance_residuals_geewa_binary,
               residuals(fit_geewa_binary, type = "deviance"))

})
