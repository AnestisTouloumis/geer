list2env(
  get_cached_test_fixtures("geer_methods_marginaleffects"),
  environment()
)

test_that("the marginaleffects example displays both result tables", {
  skip_if_no_marginaleffects()
  output <- compute_example_output("marginaleffects-support")
  expect_equal(sum(grepl("Estimate", output, fixed = TRUE)), 2L)
})

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

test_that("marginaleffects covariance names are matched exactly", {
  expect_equal(
    get_vcov.geer(fit_geewa_pois_exch, vcov = "robust"),
    vcov(fit_geewa_pois_exch, cov_type = "robust")
  )
  expect_error(
    get_vcov.geer(fit_geewa_pois_exch, vcov = "Robust"),
    "'vcov' must be one of: ",
    fixed = TRUE
  )
})

test_that("marginaleffects replacement coefficients must be finite", {
  fit <- fit_geewa_pois_exch
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    coefs <- coef(fit)
    coefs[1L] <- value
    expect_error(
      set_coef.geer(fit, coefs),
      "'coefs' must contain only finite values", fixed = TRUE
    )
  }
})

test_that("marginaleffects covariance inputs must be finite", {
  fit <- fit_geewa_pois_exch
  for (value in c(NA_real_, NaN, Inf, -Inf)) {
    covariance <- diag(length(coef(fit)))
    covariance[1L, 1L] <- value
    expect_error(
      get_vcov.geer(fit, vcov = covariance),
      "'vcov' must contain only finite values", fixed = TRUE
    )
    expect_error(get_vcov.geer(fit, vcov = function(object) covariance))
  }
})

test_that("marginaleffects adapters retain their integration dots contract", {
  fit <- fit_geewa_pois_exch
  expect_identical(get_coef.geer(fit, framework_argument = TRUE), coef(fit))
  expect_identical(get_data.geer(fit, framework_argument = TRUE), fit$data)
  expect_identical(
    coef(set_coef.geer(fit, coef(fit), framework_argument = TRUE)), coef(fit)
  )
  expect_identical(
    get_vcov.geer(fit, vcov = "robust", framework_argument = TRUE),
    vcov(fit, cov_type = "robust")
  )
  newdata <- head(test_data$epilepsy)
  expect_identical(
    get_predict.geer(fit, newdata, framework_argument = TRUE),
    get_predict.geer(fit, newdata)
  )
})

test_that("marginaleffects adapters support recovery and coercion paths", {
  fit <- fit_geewa_pois_exch
  expect_identical(get_predict.geer(fit), get_predict.geer(fit, fit$data))
  newdata <- head(fit$data)
  expect_identical(
    get_predict.geer(fit, as.list(newdata)), get_predict.geer(fit, newdata)
  )
  covariance <- unname(vcov(fit, cov_type = "robust"))
  expected <- list(names(coef(fit)), names(coef(fit)))
  expect_identical(dimnames(get_vcov.geer(fit, vcov = covariance)), expected)
  fit$data <- NULL
  expect_error(get_predict.geer(fit), "supply 'newdata' explicitly")
})
