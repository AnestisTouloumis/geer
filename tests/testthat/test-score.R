link <- sample(c("logit", "probit", "cloglog", "cauchit"), 1)
association <- sample(c("independence", "exchangeable", "unstructured"), 1)
method_gee <-  sample(c("gee", "brgee_naive", "brgee_robust", "brgee_empirical",
                        "bcgee_naive", "bcgee_robust", "bcgee_empirical",
                        "pgee_jeffreys"), 1)
cov_type <- sample(c("robust", "naive", "bias-corrected", "df-adjusted"), 1)

print(c(link, association, method_gee, cov_type))

data("cerebrovascular")
fit_model <-
  geewa_binary(formula = I(ecg == "normal") ~ treatment + factor(period),
               id = id,
               link = link,
               data = cerebrovascular,
               or_structure = association,
               method = method_gee)
reduced_model <-
  update(fit_model, formula = I(ecg == "normal") ~ factor(period))


test_that("invariance", {
  expect_equal(
    score_test(fit_model, reduced_model, cov_type = cov_type),
    score_test(reduced_model, fit_model, cov_type = cov_type))
})


fit_model <-
  geewa(formula = I(ecg == "normal") ~ treatment + factor(period),
        id = id,
        family = binomial(link = link),
        data = cerebrovascular,
        correlation_structure = association,
        method = method_gee)
reduced_model <-
  update(fit_model, formula = I(ecg == "normal") ~ factor(period))


test_that("invariance", {
  expect_equal(
    score_test(fit_model, reduced_model, cov_type = cov_type),
    score_test(reduced_model, fit_model, cov_type = cov_type))
})
