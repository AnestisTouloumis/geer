testthat::local_edition(3)

fit_or_offset_formula <- geewa_binary(
  formula = ecg ~ treatment + offset(period),
  link = "logit",
  id = id,
  data = test_data$cerebrovascular,
  orstr = "independence",
  method = "gee"
)

fit_or_offset_argument <- geewa_binary(
  formula = ecg ~ treatment,
  link = "logit",
  id = id,
  offset = period,
  data = test_data$cerebrovascular,
  orstr = "independence",
  method = "gee"
)

fit_cc_offset_formula <- geewa(
  formula = ecg ~ treatment + offset(period),
  family = binomial(link = "logit"),
  id = id,
  data = test_data$cerebrovascular,
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
