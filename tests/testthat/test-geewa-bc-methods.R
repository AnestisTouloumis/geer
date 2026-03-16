local({
  data("cerebrovascular", package = "geer")
  fit_gee <- geewa(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "logit"),
    id = id,
    phi_fixed = TRUE,
    phi_value = 1,
    data = cerebrovascular,
    corstr = "independence",
    method = "gee"
  )
  beta0 <- coef(fit_gee)
  fit_full_naive <- update(fit_gee, method = "bcgee-naive")
  fit_full_robust <- update(fit_gee, method = "bcgee-robust")
  fit_full_emp <- update(fit_gee, method = "bcgee-empirical")


  test_that("bcgee-naive: one-step from GEE start matches converged fit", {
    fit_one_step <- update(
      fit_gee,
      beta_start = beta0,
      method = "bcgee-naive",
      maxiter = 1
    )
    expect_equal(coef(fit_full_naive), coef(fit_one_step), tolerance = 1e-5)
  })


  test_that("bcgee-robust: one-step from GEE start matches converged fit", {
    fit_one_step <- update(
      fit_gee,
      beta_start = beta0,
      method = "bcgee-robust",
      maxiter = 1
    )
    expect_equal(coef(fit_full_robust), coef(fit_one_step), tolerance = 1e-5)
  })


  test_that("bcgee-empirical: one-step from GEE start matches converged fit", {
    fit_one_step <- update(
      fit_gee,
      beta_start = beta0,
      method = "bcgee-empirical",
      maxiter = 1
    )
    expect_equal(coef(fit_full_emp), coef(fit_one_step), tolerance = 1e-5)
  })
})
