# tests/testthat/test-predict-geer.R
#
# Uses the full cerebrovascular dataset.
# Keeps factor(period) inside the formula.
# Speed improvement: fit once per file and reuse across tests.

data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("predict.geer returns in-sample predictions consistent with stored values", {
  pr_link <- predict(fit, type = "link")
  pr_resp <- predict(fit, type = "response")
  expect_type(pr_link, "double")
  expect_type(pr_resp, "double")
  expect_equal(length(pr_link), fit$obs_no)
  expect_equal(length(pr_resp), fit$obs_no)
  expect_identical(pr_link, fit$linear.predictors)
  expect_identical(pr_resp, fit$fitted.values)
})


test_that("predict.geer works with newdata and returns correct lengths", {
  nd <- cerebrovascular[1:10, , drop = FALSE]
  pr1 <- predict(fit, newdata = nd, type = "link")
  pr2 <- predict(fit, newdata = nd, type = "response")
  expect_equal(length(pr1), nrow(nd))
  expect_equal(length(pr2), nrow(nd))
  expect_true(all(is.finite(pr1)))
  expect_true(all(is.finite(pr2)))
})


test_that("predict.geer with se.fit=TRUE returns list with nonnegative SEs", {
  pred_link <- predict(fit, type = "link", se.fit = TRUE, cov_type = "robust")
  pred_resp <- predict(fit, type = "response", se.fit = TRUE, cov_type = "robust")
  expect_type(pred_link, "list")
  expect_true(all(c("fit", "se.fit") %in% names(pred_link)))
  expect_equal(length(pred_link$fit), fit$obs_no)
  expect_equal(length(pred_link$se.fit), fit$obs_no)
  expect_true(all(is.finite(pred_link$se.fit)))
  expect_true(all(pred_link$se.fit >= 0))
  expect_type(pred_resp, "list")
  expect_true(all(c("fit", "se.fit") %in% names(pred_resp)))
  expect_equal(length(pred_resp$fit), fit$obs_no)
  expect_equal(length(pred_resp$se.fit), fit$obs_no)
  expect_true(all(is.finite(pred_resp$se.fit)))
  expect_true(all(pred_resp$se.fit >= 0))
})


test_that("predict.geer errors for non-geer objects", {
  expect_error(
    predict.geer(list()),
    "must be a 'geer' object"
  )
})
