test_that("model.matrix returns the stored design matrix", {
  x <- model.matrix(fit_geewa_pois_exch)
  expect_equal(x, fit_geewa_pois_exch$x)
  expect_equal(nrow(x), fit_geewa_pois_exch$obs_no)
  expect_equal(ncol(x), length(coef(fit_geewa_pois_exch)))
})

test_that("model.matrix does not depend on the stored data component", {
  fit <- fit_geewa_pois_exch
  fit$data <- NULL
  expect_equal(model.matrix(fit), fit_geewa_pois_exch$x)
  fit$data <- new.env()
  expect_equal(model.matrix(fit), fit_geewa_pois_exch$x)
})
