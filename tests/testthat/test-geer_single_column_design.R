testthat::local_edition(3)


test_that("a one-column design keeps its column name and assign attribute", {
  fit_no_intercept <- geewa(
    seizures ~ 0 + lnbaseline,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "independence"
  )
  expect_identical(colnames(fit_no_intercept$x), "lnbaseline")
  expect_identical(attr(fit_no_intercept$x, "assign"), 1L)

  fit_intercept_only <- geewa(
    seizures ~ 1,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    corstr = "independence"
  )
  expect_identical(colnames(fit_intercept_only$x), "(Intercept)")
  expect_identical(attr(fit_intercept_only$x, "assign"), 0L)
})
