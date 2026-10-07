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


test_that("formula offsets are reported by geer_formula_offset_terms", {
  expect_identical(
    geer:::geer_formula_offset_terms(fit_or_offset_formula),
    "offset(period)"
  )
  expect_identical(
    geer:::geer_formula_offset_terms(fit_or_offset_argument),
    character(0)
  )
})


test_that("sequential anova keeps a formula offset in the null model", {
  fit_null <- geewa(
    formula = ecg ~ 1 + offset(period),
    family = binomial(link = "logit"),
    id = id,
    data = cerebrovascular,
    corstr = "independence",
    method = "gee",
    phi_fixed = TRUE,
    phi_value = 1
  )
  ## The score statistic depends on the null fit, so it exposes a null model
  ## that lost the offset.
  tab <- anova(fit_cc_offset_formula, test = "score", cov_type = "robust")
  expected <- geer:::compute_score_test(fit_null, fit_cc_offset_formula, "robust")
  expect_equal(tab["treatment", "Chi"], expected$test_stat)
})


test_that("a scalar 'offset' is recycled and equals the equivalent vector offset", {
  fit_scalar <- geewa(
    ecg ~ treatment,
    family = binomial(link = "logit"),
    id = id,
    data = cerebrovascular,
    offset = 0.5
  )
  fit_vector <- geewa(
    ecg ~ treatment,
    family = binomial(link = "logit"),
    id = id,
    data = cerebrovascular,
    offset = rep(0.5, nrow(cerebrovascular))
  )
  expect_equal(coef(fit_scalar), coef(fit_vector), tolerance = 1e-8)
  expect_true(all(fit_scalar$call$offset == 0.5))
  fit_subset <- geewa(
    ecg ~ treatment,
    family = binomial(link = "logit"),
    id = id,
    data = cerebrovascular,
    offset = 0.5,
    subset = as.numeric(id) <= 40
  )
  expect_s3_class(fit_subset, "geer")
})


test_that("refits for anova, add1, drop1 and step_p keep a formula offset", {
  refit_drop <- refit_geer(fit_cc_offset_formula, . ~ . - treatment)
  expect_identical(geer_formula_offset_terms(refit_drop), "offset(period)")
  refit_add <- refit_geer(fit_cc_offset_formula, . ~ . + factor(period))
  expect_identical(geer_formula_offset_terms(refit_add), "offset(period)")
  out_drop <- drop1(fit_cc_offset_formula, test = "wald", cov_type = "robust")
  expect_s3_class(out_drop, "anova")
  expect_true(all(is.finite(out_drop[["CIC"]])))
})