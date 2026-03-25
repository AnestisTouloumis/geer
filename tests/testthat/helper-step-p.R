term_labels <- function(object) {
  attr(stats::terms(stats::formula(object)), "term.labels")
}

expect_step_p_result <- function(object) {
  testthat::expect_s3_class(object, "geer")
  testthat::expect_true(!is.null(object$anova))
  testthat::expect_s3_class(object$anova, "anova")
  testthat::expect_true(is.data.frame(object$anova))
  testthat::expect_true(all(c("Step", "Df", "Chi", "Pr(>Chi)", "CIC") %in% names(object$anova)))
}


expect_no_step_taken <- function(object) {
  expect_step_p_result(object)
  testthat::expect_equal(nrow(object$anova), 1L)
  testthat::expect_identical(object$anova$Step[[1L]], "")
}
