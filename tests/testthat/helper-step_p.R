fit_resp_full_indep <- geewa_binary(
  formula = status ~ baseline + treatment + gender + visit + age,
  id = id,
  repeated = visit,
  data = test_data$respiratory2,
  link = "probit",
  orstr = "independence",
  method = "pgee-jeffreys"
)

term_labels <- function(object) {
  attr(terms(formula(object)), "term.labels")
}

expect_step_p_result <- function(object) {
  expect_s3_class(object, "geer")
  expect_false(is.null(object$anova))
  expect_s3_class(object$anova, "anova")
  expect_s3_class(object$anova, "data.frame")
  expect_true(all(c("Step", "Df", "Chi", "Pr(>Chi)", "CIC") %in% names(object$anova)))
}

expect_no_step_taken <- function(object) {
  expect_step_p_result(object)
  expect_equal(nrow(object$anova), 1L)
  expect_identical(object$anova$Step[[1L]], NA_integer_)
}
