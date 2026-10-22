skip_if_no_marginaleffects <- function() {
  skip_if_not_installed(
    "marginaleffects",
    minimum_version = "0.32.0"
  )
}

test_that("marginaleffects low-level methods expose geer components", {
  skip_if_no_marginaleffects()
  expect_equal(insight::get_data(fit_geewa_pois_exch), fit_geewa_pois_exch$data)

  expect_equal(
    marginaleffects::get_coef(fit_geewa_pois_exch),
    coef(fit_geewa_pois_exch)
  )

  expect_equal(
    marginaleffects::get_vcov(fit_geewa_pois_exch),
    vcov(fit_geewa_pois_exch, cov_type = "bias-corrected")
  )
  expect_equal(
    marginaleffects::get_vcov(fit_geewa_pois_exch, vcov = "robust"),
    vcov(fit_geewa_pois_exch, cov_type = "robust")
  )

  custom_vcov <- vcov(fit_geewa_pois_exch, cov_type = "naive")
  expect_equal(
    marginaleffects::get_vcov(fit_geewa_pois_exch, vcov = custom_vcov),
    custom_vcov
  )

  expect_equal(
    marginaleffects::get_vcov(
      fit_geewa_pois_exch,
      vcov = function(x) vcov(x, cov_type = "robust")
    ),
    vcov(fit_geewa_pois_exch, cov_type = "robust")
  )
  expect_error(
    marginaleffects::get_vcov(fit_geewa_pois_exch, vcov = function(x) 1)
  )

  shifted <- coef(fit_geewa_pois_exch)
  shifted[[1L]] <- shifted[[1L]] + 0.1
  modified <- marginaleffects::set_coef(fit_geewa_pois_exch, shifted)
  expect_equal(coef(modified), shifted)
  expect_equal(coef(fit_geewa_pois_exch), marginaleffects::get_coef(fit_geewa_pois_exch))

  newdata <- test_data$epilepsy[1:4, , drop = FALSE]
  original_prediction <- marginaleffects::get_predict(
    fit_geewa_pois_exch,
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

test_that("marginaleffects get_predict agrees with predict", {
  skip_if_no_marginaleffects()
  newdata <- test_data$epilepsy[1:6, , drop = FALSE]

  response <- marginaleffects::get_predict(
    fit_geewa_pois_exch,
    newdata = newdata,
    type = "response"
  )
  link <- marginaleffects::get_predict(
    fit_geewa_pois_exch,
    newdata = newdata,
    type = "link"
  )

  expect_named(response, c("rowid", "estimate"))
  expect_equal(response$rowid, seq_len(nrow(newdata)))
  expect_equal(
    response$estimate,
    unname(predict(fit_geewa_pois_exch, newdata = newdata, type = "response"))
  )
  expect_equal(
    link$estimate,
    unname(predict(fit_geewa_pois_exch, newdata = newdata, type = "link"))
  )
})

test_that("marginaleffects predictions work for geer models", {
  skip_if_no_marginaleffects()
  newdata <- test_data$cerebrovascular[1:8, , drop = FALSE]
  out <- marginaleffects::predictions(
    fit_geewa_bin_exch,
    newdata = newdata,
    type = "response"
  )

  expect_s3_class(out, "predictions")
  expect_length(out$estimate, nrow(newdata))
  expect_equal(
    out$estimate,
    unname(predict(fit_geewa_bin_exch, newdata = newdata, type = "response")),
    tolerance = 1e-7
  )
  expect_true(all(is.finite(out$std.error)))
})

test_that("marginaleffects comparisons and slopes work for geer models", {
  skip_if_no_marginaleffects()
  comparison <- marginaleffects::avg_comparisons(
    fit_geewa_bin_exch,
    variables = "treatment",
    type = "response"
  )
  slope <- marginaleffects::avg_slopes(
    fit_geewa_pois_exch,
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
    fit_geewa_pois_exch,
    variables = "lnage",
    type = "response"
  )
  robust <- marginaleffects::avg_slopes(
    fit_geewa_pois_exch,
    variables = "lnage",
    type = "response",
    vcov = "robust"
  )

  expect_equal(default$estimate, robust$estimate, tolerance = 1e-10)
  expect_true(all(is.finite(default$std.error)))
  expect_true(all(is.finite(robust$std.error)))
})
