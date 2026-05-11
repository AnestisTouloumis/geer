testthat::local_edition(3)

# ── Fixtures ──────────────────────────────────────────────────────────────────

# cerebrovascular: T = 2, so choose(2, 2) = 1 time pair
fit_bin_exch_T2 <- geewa(
  formula = ecg ~ treatment + factor(period),
  family  = binomial(link = "logit"),
  data    = test_data$cerebrovascular,
  id      = id,
  corstr  = "exchangeable",
  method  = "gee"
)

# respiratory (center C2): T = 4, so choose(4, 2) = 6 time pairs
fit_bin_exch_T4 <- geewa(
  formula = status ~ baseline + treatment,
  family  = binomial(link = "logit"),
  data    = test_data$respiratory2,
  id      = id,
  repeated = visit,
  corstr  = "exchangeable",
  method  = "gee"
)

fit_bin_ar1_T4 <- geewa(
  formula = status ~ baseline + treatment,
  family  = binomial(link = "logit"),
  data    = test_data$respiratory2,
  id      = id,
  repeated = visit,
  corstr  = "ar1",
  method  = "gee"
)

# Poisson fit — used to test the non-binomial error path
fit_pois_exch <- fit_geewa_pois_exch


# ── Output structure ──────────────────────────────────────────────────────────

test_that("frechet_bounds_cor returns a data frame with the correct columns", {
  out <- frechet_bounds_cor(fit_bin_exch_T2)
  expect_s3_class(out, "data.frame")
  expect_identical(
    names(out),
    c("alpha_name", "alpha_value", "lower_max", "upper_min", "n_violated")
  )
})

test_that("frechet_bounds_cor returns one row per unique time pair", {
  out_T2 <- frechet_bounds_cor(fit_bin_exch_T2)
  expect_equal(nrow(out_T2), 1L)

  out_T4 <- frechet_bounds_cor(fit_bin_exch_T4)
  expect_equal(nrow(out_T4), choose(4L, 2L))
})

test_that("frechet_bounds_cor column types are correct", {
  out <- frechet_bounds_cor(fit_bin_exch_T4)
  expect_type(out$alpha_name,  "character")
  expect_type(out$alpha_value, "double")
  expect_type(out$lower_max,   "double")
  expect_type(out$upper_min,   "double")
  expect_type(out$n_violated,  "integer")
})

test_that("frechet_bounds_cor alpha_name labels follow alpha_j.k convention", {
  out <- frechet_bounds_cor(fit_bin_exch_T2)
  expect_match(out$alpha_name, "^alpha_\\d+\\.\\d+$")

  out4 <- frechet_bounds_cor(fit_bin_exch_T4)
  expect_true(all(grepl("^alpha_\\d+\\.\\d+$", out4$alpha_name)))
})

test_that("frechet_bounds_cor bounds are numerically valid", {
  out <- frechet_bounds_cor(fit_bin_exch_T4)
  # All values must be finite
  expect_true(all(is.finite(out$alpha_value)))
  expect_true(all(is.finite(out$lower_max)))
  expect_true(all(is.finite(out$upper_min)))
  # lower_max < upper_min (Frechet bounds must form a proper interval)
  expect_true(all(out$lower_max < out$upper_min))
  # Bounds lie within [-1, 1]
  expect_true(all(out$lower_max >= -1))
  expect_true(all(out$upper_min <=  1))
  # n_violated is non-negative
  expect_true(all(out$n_violated >= 0L))
})

test_that("frechet_bounds_cor alpha_value matches the fitted working correlation", {
  out <- frechet_bounds_cor(fit_bin_exch_T2)
  # For exchangeable with T = 2 there is one alpha; it should equal the stored alpha
  expect_equal(out$alpha_value, fit_bin_exch_T2$alpha, tolerance = 1e-10)
})

test_that("frechet_bounds_cor works for ar1 association structure", {
  out <- frechet_bounds_cor(fit_bin_ar1_T4)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), choose(4L, 2L))
  expect_true(all(is.finite(out$lower_max)))
  expect_true(all(is.finite(out$upper_min)))
})


# ── Error paths ───────────────────────────────────────────────────────────────

test_that("frechet_bounds_cor errors on non-geer input", {
  expect_error(
    frechet_bounds_cor(list()),
    "'object' must be of class \"geer\"",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor errors on non-binomial family", {
  expect_error(
    frechet_bounds_cor(fit_pois_exch),
    "family must be \"binomial\"",
    fixed = TRUE
  )
})

test_that("frechet_bounds_cor errors on independence association structure", {
  fit_indep <- geewa(
    formula  = ecg ~ treatment + factor(period),
    family   = binomial(link = "logit"),
    data     = test_data$cerebrovascular,
    id       = id,
    corstr   = "independence",
    method   = "gee"
  )
  expect_error(
    frechet_bounds_cor(fit_indep),
    "association structure must not be \"independence\"",
    fixed = TRUE
  )
})
