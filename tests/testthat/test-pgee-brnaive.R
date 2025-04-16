data("respiratory")
formula <- y ~ baseline + treatment + gender + visit + age + center
fitted_pgee_corr <- geewa(formula = formula,
                          id = id,
                          repeated = visit,
                          family = binomial(link = "logit"),
                          data = respiratory,
                          phi_fixed = TRUE,
                          phi_value = 1,
                          corstr = "independence",
                          method = "pgee-jeffreys")

fitted_brnaive_corr <- update(fitted_pgee_corr, method = "brgee-naive")


test_that("jeffreys = naive - independence - correlation", {
  expect_equal(coef(fitted_pgee_corr),
               coef(fitted_brnaive_corr))
})


fitted_pgee_or <- geewa_binary(formula = formula,
                                 id = id,
                                 repeated = visit,
                                 link = "logit",
                                 data = respiratory,
                                 orstr = "independence",
                                 method = "pgee-jeffreys")


fitted_brnaive_or <- update(fitted_pgee_or, method = "brgee-naive")


test_that("jeffreys = naive - independence1", {
  expect_equal(coef(fitted_pgee_or),
               coef(fitted_pgee_corr))
})

test_that("jeffreys = naive - independence2", {
  expect_equal(coef(fitted_pgee_or),
               coef(fitted_brnaive_or))
})

