data("cerebrovascular")
fit_geewa <- geewa(formula = I(ecg == "normal") ~ treatment + factor(period),
                   id = id,
                   family = binomial(link = "logit"),
                   phi_fixed = TRUE,
                   phi_value = 1,
                   data = cerebrovascular,
                   correlation_structure = "independence",
                   method = "gee")
fit_geewa_binary <-
  geewa_binary(formula = I(ecg == "normal") ~ treatment + factor(period),
               id = id,
               link = "logit",
               data = cerebrovascular,
               or_structure = "independence",
               method = "gee")


test_that("vcov-robust", {
  expect_equal(vcov(fit_geewa,
                    cov_type = "robust"),
               fit_geewa$robust_covariance)
  expect_equal(vcov(fit_geewa),
               fit_geewa$robust_covariance)
  expect_equal(vcov(fit_geewa_binary,
                    cov_type = "robust"),
               fit_geewa_binary$robust_covariance)
  expect_equal(vcov(fit_geewa_binary),
               fit_geewa_binary$robust_covariance)
})

test_that("vcov-naive", {
  expect_equal(vcov(fit_geewa,
                    cov_type = "naive"),
               fit_geewa$naive_covariance)
  expect_equal(vcov(fit_geewa_binary,
                    cov_type = "naive"),
               fit_geewa_binary$naive_covariance)
})

test_that("vcov-bias-corrected", {
  expect_equal(vcov(fit_geewa,
                    cov_type = "bias-corrected"),
               fit_geewa$bias_corrected_covariance)
  expect_equal(vcov(fit_geewa_binary,
                    cov_type = "bias-corrected"),
               fit_geewa_binary$bias_corrected_covariance)
})
