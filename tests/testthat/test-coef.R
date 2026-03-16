data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)

test_that("coef.geer returns the stored coefficient vector", {
  co <- coef(fit)
  expect_type(co, "double")
  expect_true(!is.null(names(co)))
  expect_true(length(co) > 0L)
  expect_identical(co, fit$coefficients)
})

test_that("coefficients() is an alias for coef() on geer objects", {
  expect_identical(coefficients(fit), coef(fit))
})
