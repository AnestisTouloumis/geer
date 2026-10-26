test_that("the sigma generic is imported for method registration", {
  imports <- parent.env(asNamespace("geer"))
  expect_identical(
    get0("sigma", envir = imports, inherits = FALSE), stats::sigma
  )
})
