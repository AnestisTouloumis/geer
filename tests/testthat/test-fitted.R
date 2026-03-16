data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("fitted.geer returns the stored fitted values vector", {
  fv <- fitted(fit)
  expect_type(fv, "double")
  expect_true(length(fv) > 0L)
  expect_identical(fv, fit$fitted.values)
})


test_that("fitted.values() is an alias for fitted() on geer objects", {
  expect_identical(fitted.values(fit), fitted(fit))
})
