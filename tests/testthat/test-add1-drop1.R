data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("add1.geer returns an anova-like table with expected structure", {
  res <- add1(
    fit,
    scope = . ~ . + treatment:factor(period),
    test = "score",
    cov_type = "robust"
  )
  expect_s3_class(res, "anova")
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 2L)
  expect_true(all(c("Df", "CIC", "Chi", "Pr(>Chi)") %in% names(res)))
  expect_true("<none>" %in% rownames(res))
  expect_true(any(grepl("treatment:factor\\(period\\)", rownames(res))))
  expect_true(is.numeric(res$Df))
  expect_true(is.numeric(res$CIC))
  expect_true(is.numeric(res$Chi))
  expect_true(is.numeric(res$`Pr(>Chi)`))
  expect_true(all(is.finite(res$CIC)))
  expect_true(all(is.finite(res$Chi[-1])))
  expect_true(all(res$`Pr(>Chi)` >= 0 & res$`Pr(>Chi)` <= 1, na.rm = TRUE))
})


test_that("add1.geer errors when scope is missing or NULL", {
  expect_error(add1(fit), "no terms in scope")
  expect_error(add1(fit, scope = NULL), "no terms in scope")
})


test_that("add1.geer working-lrt requires independence association_structure", {
  expect_error(
    add1(
      fit,
      scope = . ~ . + treatment:factor(period),
      test = "working-lrt"
    ),
    "independence"
  )
})


test_that("drop1.geer returns an anova-like table with expected structure", {
  res <- drop1(
    fit,
    scope = "treatment",
    test = "score",
    cov_type = "robust"
  )
  expect_s3_class(res, "anova")
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 2L)
  expect_true(all(c("Df", "CIC", "Chi", "Pr(>Chi)") %in% names(res)))
  expect_true("<none>" %in% rownames(res))
  expect_true("treatment" %in% rownames(res))
  expect_true(is.numeric(res$Df))
  expect_true(is.numeric(res$CIC))
  expect_true(is.numeric(res$Chi))
  expect_true(is.numeric(res$`Pr(>Chi)`))
  expect_true(all(is.finite(res$CIC)))
  expect_true(all(is.finite(res$Chi[-1])))
  expect_true(all(res$`Pr(>Chi)` >= 0 & res$`Pr(>Chi)` <= 1, na.rm = TRUE))
})


test_that("drop1.geer uses default scope when scope is missing", {
  res <- drop1(fit, test = "score", cov_type = "robust")
  expect_s3_class(res, "anova")
  expect_true(is.data.frame(res))
  expect_true(nrow(res) >= 2L)
  expect_true("<none>" %in% rownames(res))
})


test_that("drop1.geer working-lrt requires independence association_structure", {
  expect_error(
    drop1(
      fit,
      scope = "treatment",
      test = "working-lrt"
    ),
    "independence"
  )
})
