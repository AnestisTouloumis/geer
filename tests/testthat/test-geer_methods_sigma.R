test_that("sigma returns the square root of working dispersion", {
  fit <- fit_geewa_pois_exch
  expect_equal(sigma(fit), sqrt(fit$phi), tolerance = 0)
  expect_identical(sigma(fit_geewa_bin_exch), 1)
  expect_error(sigma(fit, unused = TRUE))
  expect_error(sigma.geer(1))
  fit$phi <- 0
  expect_error(sigma(fit), "'object$phi' must be positive", fixed = TRUE)
})
