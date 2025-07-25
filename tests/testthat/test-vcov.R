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
fit_geewa <- geewa(formula = ecg ~ treatment + factor(period),
                   id = id,
                   family = binomial(link = link),
                   data = cerebrovascular,
                   corstr = association,
                   method = method_gee)
fit_geewa_binary <-
  geewa_binary(formula = ecg ~ treatment + factor(period),
               id = id,
               link = link,
               data = cerebrovascular,
               orstr = association,
               method = method_gee)

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

test_that("vcov-df-adjusted", {
  expect_equal(vcov(fit_geewa,
                    cov_type = "df-adjusted"),
               67/64 * fit_geewa$robust_covariance)
  expect_equal(vcov(fit_geewa_binary,
                    cov_type = "df-adjusted"),
               67/64 * fit_geewa_binary$robust_covariance)
})
