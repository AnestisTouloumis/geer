expected_geecriteria_cols <- c(
  "QIC", "CIC", "RJC", "QICu", "GESSC", "GPC", "Parameters"
)

expect_geecriteria_table <- function(out, n_rows = NULL, row_names = NULL) {
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_identical(names(out), expected_geecriteria_cols)

  if (!is.null(n_rows)) {
    testthat::expect_equal(nrow(out), n_rows)
  }

  if (!is.null(row_names)) {
    testthat::expect_identical(rownames(out), row_names)
  }

  testthat::expect_true(all(vapply(out[expected_geecriteria_cols], is.numeric, logical(1))))
  testthat::expect_true(all(is.finite(as.matrix(out[expected_geecriteria_cols]))))
}
