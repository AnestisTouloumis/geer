geecriteria_cols <- function() {
  c(geer_criteria_columns, "Parameters")
}

undefined_geecriteria_cols <- c(
  "QICC", "GHYC", "PAC", "PT", "WR", "RMR", "RJC", "QICHH", "EQIC"
)

expect_geecriteria_table <- function(out, n_rows = NULL, row_names = NULL) {
  expected_cols <- geecriteria_cols()
  expect_s3_class(out, "data.frame")
  expect_named(out, expected_cols)
  if (!is.null(n_rows)) {
    expect_equal(nrow(out), n_rows)
  }
  if (!is.null(row_names)) {
    expect_identical(rownames(out), row_names)
  }
  expect_true(all(vapply(out[expected_cols], is.numeric, logical(1))))
  potentially_undefined <- intersect(undefined_geecriteria_cols, expected_cols)
  always_finite <- setdiff(expected_cols, potentially_undefined)
  expect_true(all(is.finite(as.matrix(out[always_finite]))))
  values <- as.matrix(out[potentially_undefined])
  expect_true(all(is.finite(values) | is.na(values)))
}
