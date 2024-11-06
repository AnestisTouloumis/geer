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


test_that("simple test", {
  placebo_pvalue <-
    fit_model |>
    summary(cov_type = cov_type) |>
    _$coefficients["treatmentplacebo", "Pr(>|z|)"]
  wald_pvalue <- wald_test(fit_model, reduced_model, cov_type = cov_type)$"P(>Chi)"
  expect_equal(placebo_pvalue,
               wald_pvalue)
})


test_that("invariance", {
  expect_equal(
    wald_test(fit_model, reduced_model, cov_type = cov_type),
    wald_test(reduced_model, fit_model, cov_type = cov_type))
})
