term_labels <- function(object) {
  attr(stats::terms(stats::formula(object)), "term.labels")
}

expect_step_p_path <- function(x) {
  testthat::expect_s3_class(x, "geer")
  testthat::expect_true("anova" %in% names(x))
  testthat::expect_s3_class(x$anova, "anova")
  testthat::expect_s3_class(x$anova, "data.frame")

  needed_cols <- c("Step", "Df", "Chi", "Pr(>Chi)", "CIC")
  testthat::expect_true(all(needed_cols %in% names(x$anova)))

  testthat::expect_type(x$anova$Step, "character")
  testthat::expect_true(is.numeric(x$anova$Df))
  testthat::expect_true(is.numeric(x$anova$Chi))
  testthat::expect_true(is.numeric(x$anova$`Pr(>Chi)`))
  testthat::expect_true(is.numeric(x$anova$CIC))
  testthat::expect_true(nrow(x$anova) >= 1L)
}

expect_no_step_taken <- function(x) {
  expect_step_p_path(x)
  testthat::expect_equal(nrow(x$anova), 1L)
  testthat::expect_identical(x$anova$Step[1], "")
}
