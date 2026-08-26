expected_geecriteria_cols <- c(
  "QIC", "QICHH", "QICC", "CIC", "RJC", "QICu", "EQIC", "GESSC", "GPC",
  "AGPC", "SGPC", "GHYC", "PAC", "Parameters"
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
  potentially_undefined <- c("QICC", "GHYC", "PAC")
  always_finite <- setdiff(expected_geecriteria_cols, potentially_undefined)
  testthat::expect_true(all(is.finite(as.matrix(out[always_finite]))))
  for (criterion in potentially_undefined) {
    testthat::expect_true(all(is.finite(out[[criterion]]) | is.na(out[[criterion]])))
  }
}
