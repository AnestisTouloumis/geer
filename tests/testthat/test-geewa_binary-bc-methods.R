data("cerebrovascular")
fit_gee <- geewa_binary(formula = I(ecg == "normal") ~ treatment + factor(period),
                        link =  "logit",
                        id = id,
                        data = cerebrovascular,
                        orstr = "independence",
                        method = "gee")

test_that("bcgee-naive", {
  fit_bcgee_naive <- update(fit_gee,
                            method = "bcgee-naive")
  fit_bcgee_naive_one_step <- update(fit_gee,
                                     beta_start = coef(fit_gee),
                                     method = "bcgee-naive",
                                     maxiter = 1)
  expect_equal(coef(fit_bcgee_naive),
               coef(fit_bcgee_naive_one_step))
})


test_that("bcgee-robust", {
  fit_bcgee_robust <- update(fit_gee,
                             method = "bcgee-robust")
  fit_bcgee_robust_one_step <- update(fit_gee,
                                      beta_start = coef(fit_gee),
                                      method = "bcgee-robust",
                                      maxiter = 1)
  expect_equal(coef(fit_bcgee_robust),
               coef(fit_bcgee_robust_one_step))
})


test_that("bcgee-empirical", {
  fit_bcgee_empirical <- update(fit_gee,
                                method = "bcgee-empirical")
  fit_bcgee_empirical_one_step <- update(fit_gee,
                                         beta_start = coef(fit_gee),
                                         method = "bcgee-empirical",
                                         maxiter = 1)
  expect_equal(coef(fit_bcgee_empirical),
               coef(fit_bcgee_empirical_one_step))
})
