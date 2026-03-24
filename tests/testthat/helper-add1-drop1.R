expect_anova_table <- function(out, candidate_pattern = NULL, n_rows = NULL) {
  testthat::expect_s3_class(out, "anova")
  testthat::expect_s3_class(out, "data.frame")

  needed_cols <- c("Df", "CIC", "Chi", "Pr(>Chi)")
  testthat::expect_true(all(needed_cols %in% names(out)))

  testthat::expect_true("<none>" %in% rownames(out))

  testthat::expect_true(is.numeric(out$Df))
  testthat::expect_true(is.numeric(out$CIC))
  testthat::expect_true(is.numeric(out$Chi))
  testthat::expect_true(is.numeric(out$`Pr(>Chi)`))

  testthat::expect_true(all(is.finite(out$CIC)))
  if (nrow(out) > 1L) {
    testthat::expect_true(all(is.finite(out$Chi[-1])))
  }
  testthat::expect_true(
    all(out$`Pr(>Chi)` >= 0 & out$`Pr(>Chi)` <= 1, na.rm = TRUE)
  )

  if (!is.null(candidate_pattern)) {
    testthat::expect_true(any(grepl(candidate_pattern, rownames(out))))
  }

  if (!is.null(n_rows)) {
    testthat::expect_equal(nrow(out), n_rows)
  }
}
