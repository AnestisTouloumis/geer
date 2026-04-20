testthat::local_edition(3)
count_fit <- fit_geewa_pois_exch


test_that("summary.geer returns a coherent summary object", {
  out <- summary(count_fit)
  expect_s3_class(out, "summary.geer")
  expect_true(is.matrix(out$coefficients) || is.data.frame(out$coefficients))
  expect_equal(
    colnames(out$coefficients),
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  )
  expect_equal(nrow(out$coefficients), length(coef(count_fit)))
  expect_equal(out$family$family, "poisson")
  expect_equal(out$family$link, "log")
})


test_that("summary.geer respects the requested covariance type", {
  for (cov_type in c("robust", "naive", "bias-corrected", "df-adjusted")) {
    out <- summary(count_fit, cov_type = cov_type)

    expect_s3_class(out, "summary.geer")
    expect_identical(out$cov_type, cov_type)
    expect_equal(nrow(out$coefficients), length(coef(count_fit)))
  }
})


test_that("summary.geer handles zero standard errors without crashing", {
  fit <- fit_geewa_pois_exch
  fit$naive_covariance[, ] <- 0
  fit$robust_covariance[, ] <- 0
  fit$bias_corrected_covariance[, ] <- 0
  out <- summary(fit)
  expect_s3_class(out, "summary.geer")
  expect_true(all(
    is.na(out$coefficients[, "z value"]) |
      is.finite(out$coefficients[, "z value"])
  ))
})


test_that("print methods return invisibly", {
  expect_invisible(print(count_fit))
  expect_invisible(print(summary(count_fit)))
})
