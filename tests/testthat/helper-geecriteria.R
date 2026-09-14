## Resolved when a test runs rather than when this file is sourced, so that
## evaluation order between helper files and the loaded package cannot matter.
## Derived from the package constant rather than restated, so that adding a
## criterion cannot leave this helper behind.
geecriteria_cols <- function() {
  c(geer:::geer_criteria_columns, "Parameters")
}


## Criteria that are legitimately NA for some fits: QICC when the residual
## cluster degrees of freedom are exhausted, GHYC and PAC outside balanced
## designs, the eigenvalue criteria when the independence information matrix is
## not positive definite, and RJC, QICHH and EQIC when the quantities they need
## cannot be evaluated.
undefined_geecriteria_cols <- c(
  "QICC", "GHYC", "PAC", "PT", "WR", "RMR", "RJC", "QICHH", "EQIC"
)


expect_geecriteria_table <- function(out, n_rows = NULL, row_names = NULL) {
  expected_cols <- geecriteria_cols()
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_identical(names(out), expected_cols)
  if (!is.null(n_rows)) {
    testthat::expect_equal(nrow(out), n_rows)
  }
  if (!is.null(row_names)) {
    testthat::expect_identical(rownames(out), row_names)
  }
  testthat::expect_true(all(vapply(out[expected_cols], is.numeric, logical(1))))
  potentially_undefined <- intersect(undefined_geecriteria_cols, expected_cols)
  always_finite <- setdiff(expected_cols, potentially_undefined)
  testthat::expect_true(all(is.finite(as.matrix(out[always_finite]))))
  for (criterion in potentially_undefined) {
    testthat::expect_true(all(is.finite(out[[criterion]]) | is.na(out[[criterion]])))
  }
}
