data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("residuals.geer returns working residuals identical to stored residuals", {
  r_work <- residuals(fit, type = "working")
  expect_type(r_work, "double")
  expect_equal(length(r_work), fit$obs_no)
  expect_identical(r_work, fit$residuals)
})


test_that("residuals.geer returns deviance residuals with correct length and finite values", {
  r_dev <- residuals(fit, type = "deviance")
  expect_type(r_dev, "double")
  expect_equal(length(r_dev), fit$obs_no)
  expect_true(all(is.finite(r_dev)))
  if (isTRUE(fit$df.residual > 0)) {
    expect_equal(sign(r_dev), sign(fit$y - fit$fitted.values))
  }
})


test_that("residuals.geer returns Pearson residuals consistent with family variance (binomial)", {
  r_p <- residuals(fit, type = "pearson")
  expect_type(r_p, "double")
  expect_equal(length(r_p), fit$obs_no)
  expect_true(all(is.finite(r_p)))
  y <- fit$y
  mu <- fit$fitted.values
  w <- fit$prior.weights
  v <- fit$family$variance(mu)
  denom <- sqrt(pmax(v, 0))
  denom[denom == 0] <- NA_real_
  expected <- (y - mu) / denom
  if (all(w == 1)) {
    expect_equal(r_p, expected, tolerance = 1e-10)
  } else {
    ok <- is.finite(expected) & is.finite(r_p)
    expect_gt(stats::cor(r_p[ok], expected[ok]), 0.9)
  }
})


test_that("residuals.geer errors for non-geer objects", {
  expect_error(
    residuals.geer(list()),
    "must be a 'geer' object"
  )
})
