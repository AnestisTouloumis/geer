expect_anova_table <- function(out, candidate_pattern = NULL, n_rows = NULL) {
  expect_s3_class(out, "anova")
  expect_s3_class(out, "data.frame")
  needed_cols <- c("Df", "CIC", "Chi", "Pr(>Chi)")
  expect_true(all(needed_cols %in% names(out)))
  expect_true("<none>" %in% rownames(out))
  expect_true(is.numeric(out$Df))
  expect_true(is.numeric(out$CIC))
  expect_true(is.numeric(out$Chi))
  expect_true(is.numeric(out$`Pr(>Chi)`))
  expect_true(all(is.finite(out$CIC)))
  if (nrow(out) > 1L) {
    expect_true(all(is.finite(out$Chi[-1])))
  }
  expect_true(
    all(out$`Pr(>Chi)` >= 0 & out$`Pr(>Chi)` <= 1, na.rm = TRUE)
  )
  if (!is.null(candidate_pattern)) {
    expect_true(any(grepl(candidate_pattern, rownames(out))))
  }
  if (!is.null(n_rows)) {
    expect_equal(nrow(out), n_rows)
  }
}
