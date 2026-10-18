cerebrovascular <- test_data$cerebrovascular
respiratory2 <- test_data$respiratory2

fit_bin_exch_T2 <- geewa(
  formula = ecg ~ treatment + factor(period),
  family = binomial(link = "logit"),
  data = cerebrovascular,
  id = id,
  corstr = "exchangeable",
  method = "gee"
)

fit_bin_exch_T4 <- geewa(
  formula = status ~ baseline + treatment,
  family = binomial(link = "logit"),
  data = respiratory2,
  id = id,
  repeated = visit,
  corstr = "exchangeable",
  method = "gee"
)

fit_bin_unstr_T4 <- geewa(
  formula = status ~ baseline + treatment,
  family = binomial(link = "logit"),
  data = respiratory2,
  id = id,
  repeated = visit,
  corstr = "unstructured",
  method = "gee"
)

fit_bin_ar1_T4 <- geewa(
  formula = status ~ baseline + treatment,
  family = binomial(link = "logit"),
  data = respiratory2,
  id = id,
  repeated = visit,
  corstr = "ar1",
  method = "gee"
)

test_that("frechet_bounds_cor returns a data frame with the documented columns", {
  out <- frechet_bounds_cor(fit_bin_exch_T4)
  expect_s3_class(out, "data.frame")
  expect_named(
    out,
    c(
      "alpha_name", "alpha_value", "lower_max", "upper_min",
      "n_clusters", "n_violated"
    )
  )
  expect_type(out$alpha_name, "character")
  expect_type(out$alpha_value, "double")
  expect_type(out$lower_max, "double")
  expect_type(out$upper_min, "double")
  expect_type(out$n_clusters, "integer")
  expect_type(out$n_violated, "integer")
  expect_true(all(grepl("^alpha_\\d+\\.\\d+$", out$alpha_name)))
})

test_that("frechet_bounds_cor returns one row per unique time pair", {
  out_T2 <- frechet_bounds_cor(fit_bin_exch_T2)
  expect_equal(nrow(out_T2), 1L)

  out_T4 <- frechet_bounds_cor(fit_bin_exch_T4)
  expect_equal(nrow(out_T4), choose(4L, 2L))
})

test_that("frechet_bounds_cor bounds are numerically valid", {
  out <- frechet_bounds_cor(fit_bin_exch_T4)
  expect_true(all(is.finite(out$alpha_value)))
  expect_true(all(is.finite(out$lower_max)))
  expect_true(all(is.finite(out$upper_min)))
  expect_true(all(out$lower_max < out$upper_min))
  expect_true(all(out$lower_max >= -1))
  expect_true(all(out$upper_min <= 1))
  expect_true(all(out$n_violated >= 0L))
  expect_true(all(out$n_violated <= out$n_clusters))
  expect_true(all(out$n_clusters > 0L))
})

test_that("frechet_bounds_cor alpha_value matches the fitted working correlation", {
  out <- frechet_bounds_cor(fit_bin_exch_T2)
  expect_equal(out$alpha_value, fit_bin_exch_T2$alpha, tolerance = 1e-10)
})

test_that("frechet_bounds_cor works for ar1 association structure", {
  out <- frechet_bounds_cor(fit_bin_ar1_T4)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), choose(4L, 2L))
  expect_true(all(is.finite(out$lower_max)))
  expect_true(all(is.finite(out$upper_min)))
})

test_that("frechet_bounds_cor works for an unstructured association structure", {
  out <- frechet_bounds_cor(fit_bin_unstr_T4)
  expect_equal(nrow(out), choose(4L, 2L))
  expect_true(all(out$lower_max < out$upper_min))
  expect_identical(out$alpha_name, names(fit_bin_unstr_T4$alpha))
  expect_equal(
    out$alpha_value,
    unname(fit_bin_unstr_T4$alpha),
    tolerance = 1e-10
  )
})

test_that("frechet_bounds_cor errors on a geewa_binary fit", {
  expect_identical(fit_geewa_bin_exch$family$family, "binomial")
  expect_false(
    identical(fit_geewa_bin_exch$association_structure, "independence")
  )
  expect_identical(fit_geewa_bin_exch$fit_function, "geewa_binary")
  expect_error(
    frechet_bounds_cor(fit_geewa_bin_exch),
    "must be a geewa() fit",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor recovers the fit function from the call when the tag is absent", {
  fit_no_tag <- fit_bin_exch_T2
  fit_no_tag$fit_function <- NULL
  expect_identical(
    frechet_bounds_cor(fit_no_tag),
    frechet_bounds_cor(fit_bin_exch_T2)
  )
})

test_that("frechet_bounds_cor errors when the fit function cannot be determined", {
  fit_unknown <- fit_bin_exch_T2
  fit_unknown$fit_function <- NULL
  fit_unknown$call[[1L]] <- quote(some_wrapper)
  expect_error(
    frechet_bounds_cor(fit_unknown),
    "must be a geewa() fit",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor errors on non-binomial family", {
  expect_error(
    frechet_bounds_cor(fit_geewa_pois_exch),
    "must be a binomial fit",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor errors on independence association structure", {
  fit_indep <- geewa(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "logit"),
    data = cerebrovascular,
    id = id,
    corstr = "independence",
    method = "gee"
  )
  expect_error(
    frechet_bounds_cor(fit_indep),
    "must not have an 'independence' association structure",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor errors on degenerate fitted probabilities", {
  fit_degenerate <- fit_bin_exch_T4
  fit_degenerate$fitted.values[1L] <- 1
  expect_error(
    frechet_bounds_cor(fit_degenerate),
    "the Frechet bounds are undefined",
    fixed = TRUE
  )

  fit_zero <- fit_bin_exch_T4
  fit_zero$fitted.values[3L] <- 0
  expect_error(
    frechet_bounds_cor(fit_zero),
    "the Frechet bounds are undefined",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor errors on a non-binary response", {
  fit_grouped <- fit_bin_exch_T4
  fit_grouped$y[2L] <- 0.5
  expect_error(
    frechet_bounds_cor(fit_grouped),
    "must have a binary (0/1) response",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor accepts a quasibinomial fit", {
  fit_quasi <- geewa(
    formula = ecg ~ treatment + factor(period),
    family = quasibinomial(link = "logit"),
    data = cerebrovascular,
    id = id,
    corstr = "exchangeable",
    method = "gee"
  )
  out <- frechet_bounds_cor(fit_quasi)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
})
