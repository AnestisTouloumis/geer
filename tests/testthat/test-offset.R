testthat::local_edition(3)

cerebrovascular <- test_data$cerebrovascular

fit_or_offset_formula <- geewa_binary(
  formula = ecg ~ treatment + offset(period),
  link = "logit",
  id = id,
  data = cerebrovascular,
  orstr = "independence",
  method = "gee"
)

fit_or_offset_argument <- geewa_binary(
  formula = ecg ~ treatment,
  link = "logit",
  id = id,
  offset = period,
  data = cerebrovascular,
  orstr = "independence",
  method = "gee"
)

fit_cc_offset_formula <- geewa(
  formula = ecg ~ treatment + offset(period),
  family = binomial(link = "logit"),
  id = id,
  data = cerebrovascular,
  corstr = "independence",
  method = "gee",
  phi_fixed = TRUE,
  phi_value = 1
)

test_that("geewa_binary treats formula and argument offsets equivalently", {
  expect_s3_class(fit_or_offset_formula, "geer")
  expect_s3_class(fit_or_offset_argument, "geer")
  expect_equal(coef(fit_or_offset_formula), coef(fit_or_offset_argument))
})


test_that("geewa and geewa_binary agree for the same offset specification under independence", {
  expect_s3_class(fit_cc_offset_formula, "geer")
  expect_equal(coef(fit_or_offset_formula), coef(fit_cc_offset_formula))
})


test_that("predict.geer uses offset supplied via the offset argument with newdata", {
  nd <- cerebrovascular[1:8, , drop = FALSE]
  pred_formula <- predict(
    fit_or_offset_formula,
    newdata = nd,
    type = "link"
  )
  pred_argument <- predict(
    fit_or_offset_argument,
    newdata = nd,
    type = "link"
  )
  expect_equal(pred_formula, pred_argument)
})


test_that("predict.geer gives the same response-scale predictions for formula and argument offsets", {
  nd <- cerebrovascular[1:8, , drop = FALSE]
  pred_formula <- predict(
    fit_or_offset_formula,
    newdata = nd,
    type = "response"
  )
  pred_argument <- predict(
    fit_or_offset_argument,
    newdata = nd,
    type = "response"
  )
  expect_equal(pred_formula, pred_argument)
})
