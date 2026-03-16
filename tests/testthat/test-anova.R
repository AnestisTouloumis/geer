data("cerebrovascular", package = "geer")
fit_full <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)
fit0 <- geewa_binary(
  formula = ecg ~ 1,
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)
fit1 <- geewa_binary(
  formula = ecg ~ treatment,
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("anova.geer (single model) returns an anova data frame with expected columns", {
  tab <- anova(fit_full, test = "wald", cov_type = "robust")
  expect_s3_class(tab, "anova")
  expect_true(is.data.frame(tab))
  expect_true(all(c("Df", "Resid. Df", "Chi", "Pr(>Chi)") %in% colnames(tab)))
  expect_true(nrow(tab) >= 1L)
  expect_true(rownames(tab)[1] %in% c("NULL", "(Intercept)", "1"))
  expect_true(is.numeric(tab[["Chi"]]))
  expect_true(all(is.na(tab[["Pr(>Chi)"]]) | (tab[["Pr(>Chi)"]] >= 0 & tab[["Pr(>Chi)"]] <= 1)))
})


test_that("anova.geer supports test = 'score' for a single model", {
  tab <- anova(fit_full, test = "score", cov_type = "robust")
  expect_s3_class(tab, "anova")
  expect_true(is.data.frame(tab))
  expect_true(all(c("Df", "Resid. Df", "Chi", "Pr(>Chi)") %in% colnames(tab)))
})


test_that("anova.geer errors for working-lrt when association_structure is not independence", {
  expect_error(
    anova(fit_full, test = "working-lrt"),
    "independence"
  )
})


test_that("anova.geer with multiple models returns sequential tests (wald)", {
  tab <- anova(fit0, fit1, test = "wald", cov_type = "robust")
  expect_s3_class(tab, "anova")
  expect_true(is.data.frame(tab))
  expect_true(all(c("Resid. Df", "Df", "Chi", "Pr(>Chi)") %in% colnames(tab)))
  expect_equal(nrow(tab), 2L)
  expect_true(is.finite(tab[2, "Chi"]))
  expect_true(is.finite(tab[2, "Pr(>Chi)"]))
  expect_true(tab[2, "Pr(>Chi)"] >= 0 && tab[2, "Pr(>Chi)"] <= 1)
})


test_that("anova.geer with multiple models returns sequential tests (score)", {
  tab <- anova(fit0, fit1, test = "score", cov_type = "robust")
  expect_s3_class(tab, "anova")
  expect_true(is.data.frame(tab))
  expect_true(all(c("Resid. Df", "Df", "Chi", "Pr(>Chi)") %in% colnames(tab)))
  expect_equal(nrow(tab), 2L)
  expect_true(is.finite(tab[2, "Chi"]))
  expect_true(is.finite(tab[2, "Pr(>Chi)"]))
  expect_true(tab[2, "Pr(>Chi)"] >= 0 && tab[2, "Pr(>Chi)"] <= 1)
})


test_that("anova.geer accepts cov_type options for wald/score", {
  cov_types <- c("robust", "naive", "bias-corrected", "df-adjusted")
  for (ct in cov_types) {
    tab_w <- anova(fit_full, test = "wald", cov_type = ct)
    tab_s <- anova(fit_full, test = "score", cov_type = ct)
    expect_s3_class(tab_w, "anova")
    expect_s3_class(tab_s, "anova")
    expect_true(is.data.frame(tab_w))
    expect_true(is.data.frame(tab_s))
  }
})


test_that("anova.geer ignores pmethod for wald/score (does not error when supplied)", {
  expect_s3_class(anova(fit_full, test = "wald", pmethod = "rao-scott"), "anova")
  expect_s3_class(anova(fit_full, test = "score", pmethod = "satterthwaite"), "anova")
})


test_that("anova.geer warns and drops named extra arguments", {
  expect_warning(
    anova(fit_full, bogus = 123),
    "invalid and dropped"
  )
})


test_that("anova.geer silently drops unnamed non-geer objects in ...", {
  tab <- anova(fit0, fit1, 123, test = "wald", cov_type = "robust")
  expect_s3_class(tab, "anova")
  expect_equal(nrow(tab), 2L)
})
