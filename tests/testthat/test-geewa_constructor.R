testthat::local_edition(3)

test_that("a real geewa() fit passes new_geer() and validate_geer()", {
  fit <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  expect_s3_class(fit, "geer")
  expect_s3_class(validate_geer(fit), "geer")
  expect_equal(coef(validate_geer(fit)), coef(fit))
})
