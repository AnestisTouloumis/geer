testthat::local_edition(3)

test_that("geewa bcgee methods return geer objects and one-step updates match converged fits", {
  fit_gee <- geewa(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "logit"),
    id = id,
    phi_fixed = TRUE,
    phi_value = 1,
    data = test_data$cerebrovascular,
    corstr = "independence",
    method = "gee"
  )
  beta_start_gee <- coef(fit_gee)
  for (variant in c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    fit_conv <- update(fit_gee, method = variant)
    fit_one  <- update(fit_gee, beta_start = beta_start_gee,
                       method = variant, maxiter = 1)
    expect_s3_class(fit_conv, "geer")
    expect_equal(coef(fit_conv), coef(fit_one),
                 tolerance = 1e-5, label = variant)
  }
})



test_that("geewa_binary bcgee methods return geer objects and one-step updates match converged fits", {
  fit_gee <- geewa_binary(
    formula = ecg ~ treatment + factor(period),
    link = "logit",
    id = id,
    data = test_data$cerebrovascular,
    orstr = "independence",
    method = "gee"
  )
  beta_start_gee <- coef(fit_gee)
  for (variant in c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    fit_conv <- update(fit_gee, method = variant)
    fit_one  <- update(fit_gee, beta_start = beta_start_gee,
                       method = variant, maxiter = 1)
    expect_s3_class(fit_conv, "geer")
    expect_equal(coef(fit_conv), coef(fit_one),
                 tolerance = 1e-6, label = variant)
  }
})
