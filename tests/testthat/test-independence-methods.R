link <- sample(c("logit", "probit", "cloglog", "cauchit"), 1)
print(link)

data("respiratory")
fit_gee_binary <- geewa_binary(formula = y ~ baseline + treatment*gender + visit*age + center,
                               link =  link,
                               id = id,
                               repeated = visit,
                               data = respiratory,
                               or_structure = "independence",
                               method = "gee",
                               tolerance = 1e-10)

fit_gee <- geewa(formula = y ~ baseline + treatment*gender + visit*age + center,
                 id = id,
                 repeated = visit,
                 family = binomial(link = link),
                 data = respiratory,
                 phi_fixed = TRUE,
                 phi_value = 1,
                 correlation_structure = "independence",
                 method = "gee",
                 tolerance = 1e-10)



test_that("gee", {
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})


test_that("brgee_robust", {
  fit_gee_binary <- update(fit_gee_binary,
                           method = "brgee_robust")
  fit_gee <- update(fit_gee,
                    method = "brgee_robust")
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})


test_that("brgee_naive", {
  fit_gee_binary <- update(fit_gee_binary,
                           method = "brgee_naive")
  fit_gee <- update(fit_gee,
                    method = "brgee_naive")
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})


test_that("brgee_empirical", {
  fit_gee_binary <- update(fit_gee_binary,
                           method = "brgee_empirical")
  fit_gee <- update(fit_gee,
                    method = "brgee_empirical")
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})


test_that("bcgee_robust", {
  fit_gee_binary <- update(fit_gee_binary,
                           method = "bcgee_robust")
  fit_gee <- update(fit_gee,
                    method = "bcgee_robust")
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})


test_that("bcgee_naive", {
  fit_gee_binary <- update(fit_gee_binary,
                           method = "bcgee_naive")
  fit_gee <- update(fit_gee,
                    method = "bcgee_naive")
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})


test_that("bcgee_empirical", {
  fit_gee_binary <- update(fit_gee_binary,
                           method = "bcgee_empirical")
  fit_gee <- update(fit_gee,
                    method = "bcgee_empirical")
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})

test_that("pgee_jeffreys", {
  fit_gee_binary <- update(fit_gee_binary,
                           method = "pgee_jeffreys")
  fit_gee <- update(fit_gee,
                    method = "pgee_jeffreys")
  expect_equal(coef(fit_gee_binary),
               coef(fit_gee))
})
