testthat::local_edition(3)

test_that("new_geer and geer construct geer objects from valid lists", {
  stub <- unclass(fit_geewa_pois_indep)

  expect_s3_class(new_geer(stub), "geer")
  expect_s3_class(geer(stub), "geer")
})

test_that("geer rejects invalid inputs", {
  expect_error(geer(42), "'x' must be a list")
  expect_error(geer("a string"), "'x' must be a list")

  incomplete <- list(coefficients = 1, residuals = 1)

  expect_error(
    geer(incomplete),
    "missing required components"
  )
  expect_error(
    geer(incomplete),
    "fitted.values|rank|family"
  )
})
