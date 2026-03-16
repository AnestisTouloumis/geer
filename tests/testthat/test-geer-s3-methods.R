data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  repeated = period,
  data = cerebrovascular,
  link = "logit",
  orstr = "independence",
  method = "gee"
)


test_that("basic S3 methods behave consistently", {
  expect_s3_class(fit, "geer")
  expect_identical(coef(fit), fit$coefficients)
  expect_identical(fitted(fit), fit$fitted.values)
  expect_identical(predict(fit, type = "link"), fit$linear.predictors)
  expect_identical(predict(fit, type = "response"), fit$fitted.values)
  newdata <- cerebrovascular[1:5, , drop = FALSE]
  pred_link <- predict(fit, newdata = newdata, type = "link")
  pred_resp <- predict(fit, newdata = newdata, type = "response")
  expect_length(pred_link, 5L)
  expect_length(pred_resp, 5L)
  expect_true(all(is.finite(pred_link)))
  expect_true(all(pred_resp >= 0 & pred_resp <= 1))
  pred_se <- predict(fit, newdata = newdata, type = "response", se.fit = TRUE)
  expect_type(pred_se, "list")
  expect_named(pred_se, c("fit", "se.fit"))
  expect_length(pred_se$fit, 5L)
  expect_length(pred_se$se.fit, 5L)
    s <- summary(fit, cov_type = "robust")
  expect_s3_class(s, "summary.geer")
  expect_true(is.matrix(s$coefficients))
  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in% colnames(s$coefficients)))
})
