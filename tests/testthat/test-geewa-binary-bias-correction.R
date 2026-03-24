testthat::local_edition(3)

fit_gee <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  link = "logit",
  id = id,
  data = test_data$cerebrovascular,
  orstr = "independence",
  method = "gee"
)

beta_start_gee <- coef(fit_gee)

fit_bcgee_naive <- update(
  fit_gee,
  method = "bcgee-naive"
)

fit_bcgee_robust <- update(
  fit_gee,
  method = "bcgee-robust"
)

fit_bcgee_empirical <- update(
  fit_gee,
  method = "bcgee-empirical"
)

test_that("bcgee methods return geer objects and one-step updates match converged fits for geewa_binary", {
  expect_s3_class(fit_bcgee_naive, "geer")
  expect_s3_class(fit_bcgee_robust, "geer")
  expect_s3_class(fit_bcgee_empirical, "geer")

  fit_one_step_naive <- update(
    fit_gee,
    beta_start = beta_start_gee,
    method = "bcgee-naive",
    maxiter = 1
  )

  fit_one_step_robust <- update(
    fit_gee,
    beta_start = beta_start_gee,
    method = "bcgee-robust",
    maxiter = 1
  )

  fit_one_step_empirical <- update(
    fit_gee,
    beta_start = beta_start_gee,
    method = "bcgee-empirical",
    maxiter = 1
  )

  expect_equal(coef(fit_bcgee_naive), coef(fit_one_step_naive), tolerance = 1e-6)
  expect_equal(coef(fit_bcgee_robust), coef(fit_one_step_robust), tolerance = 1e-6)
  expect_equal(coef(fit_bcgee_empirical), coef(fit_one_step_empirical), tolerance = 1e-6)
})
