testthat::local_edition(3)

skip_if_no_marginaleffects <- function() {
  testthat::skip_if_not_installed(
    "marginaleffects",
    minimum_version = "0.32.0"
  )
}

count_fit <- fit_geewa_pois_exch
binary_fit <- fit_geewa_bin_exch


test_that("marginaleffects low-level methods expose geer components", {
  skip_if_no_marginaleffects()
  expect_equal(insight::get_data(count_fit), count_fit$data)

  expect_equal(
    marginaleffects::get_coef(count_fit),
    coef(count_fit)
  )

  expect_equal(
    marginaleffects::get_vcov(count_fit),
    vcov(count_fit, cov_type = "bias-corrected")
  )
  expect_equal(
    marginaleffects::get_vcov(count_fit, vcov = "robust"),
    vcov(count_fit, cov_type = "robust")
  )

  custom_vcov <- vcov(count_fit, cov_type = "naive")
  expect_equal(
    marginaleffects::get_vcov(count_fit, vcov = custom_vcov),
    custom_vcov
  )

  shifted <- coef(count_fit)
  shifted[[1L]] <- shifted[[1L]] + 0.1
  modified <- marginaleffects::set_coef(count_fit, shifted)
  expect_equal(coef(modified), shifted)
  expect_equal(coef(count_fit), marginaleffects::get_coef(count_fit))

  newdata <- test_data$epilepsy[1:4, , drop = FALSE]
  original_prediction <- marginaleffects::get_predict(
    count_fit,
    newdata = newdata,
    type = "response"
  )$estimate
  modified_prediction <- marginaleffects::get_predict(
    modified,
    newdata = newdata,
    type = "response"
  )$estimate
  expect_false(isTRUE(all.equal(original_prediction, modified_prediction)))
})


test_that("marginaleffects get_predict agrees with predict.geer", {
  skip_if_no_marginaleffects()
  newdata <- test_data$epilepsy[1:6, , drop = FALSE]

  response <- marginaleffects::get_predict(
    count_fit,
    newdata = newdata,
    type = "response"
  )
  link <- marginaleffects::get_predict(
    count_fit,
    newdata = newdata,
    type = "link"
  )

  expect_named(response, c("rowid", "estimate"))
  expect_equal(response$rowid, seq_len(nrow(newdata)))
  expect_equal(
    response$estimate,
    unname(predict(count_fit, newdata = newdata, type = "response"))
  )
  expect_equal(
    link$estimate,
    unname(predict(count_fit, newdata = newdata, type = "link"))
  )
})


test_that("marginaleffects predictions work for geer models", {
  skip_if_no_marginaleffects()
  newdata <- test_data$cerebrovascular[1:8, , drop = FALSE]
  out <- marginaleffects::predictions(
    binary_fit,
    newdata = newdata,
    type = "response"
  )

  expect_s3_class(out, "predictions")
  expect_length(out$estimate, nrow(newdata))
  expect_equal(
    out$estimate,
    unname(predict(binary_fit, newdata = newdata, type = "response")),
    tolerance = 1e-7
  )
  expect_true(all(is.finite(out$std.error)))
})


test_that("marginaleffects comparisons and slopes work for geer models", {
  skip_if_no_marginaleffects()
  comparison <- marginaleffects::avg_comparisons(
    binary_fit,
    variables = "treatment",
    type = "response"
  )
  slope <- marginaleffects::avg_slopes(
    count_fit,
    variables = "lnage",
    type = "response"
  )

  expect_s3_class(comparison, "comparisons")
  expect_s3_class(slope, "slopes")
  expect_true(all(is.finite(comparison$estimate)))
  expect_true(all(is.finite(comparison$std.error)))
  expect_true(all(is.finite(slope$estimate)))
  expect_true(all(is.finite(slope$std.error)))
})


test_that("marginaleffects robust covariance is available for geer models", {
  skip_if_no_marginaleffects()
  default <- marginaleffects::avg_slopes(
    count_fit,
    variables = "lnage",
    type = "response"
  )
  robust <- marginaleffects::avg_slopes(
    count_fit,
    variables = "lnage",
    type = "response",
    vcov = "robust"
  )

  expect_equal(default$estimate, robust$estimate, tolerance = 1e-10)
  expect_true(all(is.finite(default$std.error)))
  expect_true(all(is.finite(robust$std.error)))
})
