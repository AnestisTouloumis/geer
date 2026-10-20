test_that("model.matrix returns the stored design matrix without the data", {
  x <- model.matrix(fit_geewa_pois_exch)
  expect_equal(x, fit_geewa_pois_exch$x)
  expect_identical(
    dim(x),
    c(fit_geewa_pois_exch$obs_no, length(coef(fit_geewa_pois_exch)))
  )
  fit <- fit_geewa_pois_exch
  fit$data <- NULL
  expect_equal(model.matrix(fit), fit_geewa_pois_exch$x)
  fit$data <- new.env()
  expect_equal(model.matrix(fit), fit_geewa_pois_exch$x)
})
