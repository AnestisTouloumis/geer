data("cerebrovascular")
fit_gee <- geewa(formula = ecg ~ treatment + factor(period),
                 family = binomial(link = "logit"),
                 id = id,
                 phi_fixed = TRUE,
                 phi_value = 1,
                 data = cerebrovascular,
                 corstr = "independence",
                 method = "gee")

test_that("bcgee_naive", {
  fit_bcgee_naive <- update(fit_gee,
                            method = "bcgee-naive")
  fit_bcgee_naive_one_step <- update(fit_gee,
                                     beta_start = coef(fit_gee),
                                     method = "bcgee-naive",
                                     maxiter = 1)
  expect_equal(coef(fit_bcgee_naive),
               coef(fit_bcgee_naive_one_step))
})


test_that("bcgee_robust", {
  fit_bcgee_robust <- update(fit_gee,
                             method = "bcgee-robust")
  fit_bcgee_robust_one_step <- update(fit_gee,
                                      beta_start = coef(fit_gee),
                                      method = "bcgee-robust",
                                      maxiter = 1)
  expect_equal(coef(fit_bcgee_robust),
               coef(fit_bcgee_robust_one_step))
})


test_that("bcgee_empirical", {
  fit_bcgee_empirical <- update(fit_gee,
                                method = "bcgee-empirical")
  fit_bcgee_empirical_one_step <- update(fit_gee,
                                         beta_start = coef(fit_gee),
                                         method = "bcgee-empirical",
                                         maxiter = 1)
  expect_equal(coef(fit_bcgee_empirical),
               coef(fit_bcgee_empirical_one_step))
})
