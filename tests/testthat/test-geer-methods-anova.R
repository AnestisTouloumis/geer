testthat::local_edition(3)

fit_null <- geewa(
  formula = seizures ~ 1,
  data = test_data$epilepsy,
  id = id,
  family = poisson(link = "log"),
  corstr = "independence",
  method = "gee"
)

fit_trt <- geewa(
  formula = seizures ~ treatment,
  data = test_data$epilepsy,
  id = id,
  family = poisson(link = "log"),
  corstr = "independence",
  method = "gee"
)

fit_full <- fit_geewa_pois_indep   # seizures ~ treatment + lnbaseline + lnage


test_that("anova.geer single-model returns an anova data frame with correct structure", {
  out <- anova(fit_full, test = "wald", cov_type = "robust")
  expect_s3_class(out, "anova")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("Df", "Resid. Df", "Chi", "Pr(>Chi)") %in% names(out)))
  expect_equal(nrow(out), 1L + length(attr(fit_full$terms, "term.labels")))
})


test_that("anova.geer single-model first row has NA test statistics", {
  out <- anova(fit_full, test = "wald", cov_type = "robust")
  expect_true(is.na(out$Chi[1L]))
  expect_true(is.na(out$`Pr(>Chi)`[1L]))
})


test_that("anova.geer single-model subsequent rows have finite test statistics", {
  out <- anova(fit_full, test = "wald", cov_type = "robust")
  expect_true(all(is.finite(out$Chi[-1L])))
  expect_true(all(is.finite(out$`Pr(>Chi)`[-1L])))
  expect_true(all(out$Chi[-1L] >= 0))
  expect_true(all(out$`Pr(>Chi)`[-1L] >= 0 & out$`Pr(>Chi)`[-1L] <= 1))
})


test_that("anova.geer single-model row names match term labels", {
  out <- anova(fit_full, test = "wald", cov_type = "robust")
  terms <- attr(fit_full$terms, "term.labels")
  expect_identical(rownames(out)[-1L], terms)
  expect_identical(rownames(out)[1L], "NULL")
})


test_that("anova.geer heading contains family and link information", {
  out <- anova(fit_full, test = "wald", cov_type = "robust")
  heading <- attr(out, "heading")
  expect_type(heading, "character")
  expect_true(grepl("poisson", heading, ignore.case = TRUE))
  expect_true(grepl("log", heading, ignore.case = TRUE))
})


test_that("anova.geer single-model works for wald and score test types", {
  for (tst in c("wald", "score")) {
    out <- anova(fit_full, test = tst, cov_type = "robust")
    expect_s3_class(out, "anova")
    expect_true(all(is.finite(out$Chi[-1L])))
  }
})


test_that("anova.geer multi-model returns an anova data frame with correct structure", {
  out <- anova(fit_null, fit_trt, fit_full, test = "wald", cov_type = "robust")
  expect_s3_class(out, "anova")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("Df", "Resid. Df", "Chi", "Pr(>Chi)") %in% names(out)))
  expect_equal(nrow(out), 3L)
})


test_that("anova.geer multi-model first row has NA test statistics", {
  out <- anova(fit_null, fit_trt, fit_full, test = "wald", cov_type = "robust")
  expect_true(is.na(out$Chi[1L]))
  expect_true(is.na(out$`Pr(>Chi)`[1L]))
})


test_that("anova.geer multi-model subsequent rows have finite test statistics", {
  out <- anova(fit_null, fit_trt, fit_full, test = "wald", cov_type = "robust")
  expect_true(all(is.finite(out$Chi[-1L])))
  expect_true(all(is.finite(out$`Pr(>Chi)`[-1L])))
  expect_true(all(out$Chi[-1L] >= 0))
  expect_true(all(out$`Pr(>Chi)`[-1L] >= 0 & out$`Pr(>Chi)`[-1L] <= 1))
})


test_that("anova.geer multi-model Resid. Df decreases monotonically", {
  out <- anova(fit_null, fit_trt, fit_full, test = "wald", cov_type = "robust")
  expect_true(all(diff(out$`Resid. Df`) < 0))
})


test_that("anova.geer rejects working-lrt for non-independence models", {
  expect_error(
    anova(fit_trt, fit_geewa_pois_exch, test = "working-lrt", cov_type = "robust"),
    "independence"
  )
})


test_that("anova.geer rejects non-geer first argument", {
  expect_error(anova.geer(list()), "'object' must be of 'geer' class")
})


test_that("anova.geer silently drops non-geer objects from dots", {
  out <- anova(fit_null, fit_trt, "not_a_model",
               test = "wald", cov_type = "robust")
  expect_s3_class(out, "anova")
  expect_equal(nrow(out), 2L)
})


test_that("anova.geer warns about named arguments in dots", {
  expect_warning(
    anova(fit_full, foo = "bar", test = "wald", cov_type = "robust"),
    "named arguments"
  )
})
