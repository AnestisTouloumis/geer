testthat::local_edition(3)


test_that("format_percent returns a character vector with percent suffixes", {
  out <- format_percent(c(0.025, 0.975), digits = 3)
  expect_type(out, "character")
  expect_length(out, 2L)
  expect_true(all(grepl("%$", out)))
})


test_that("format_percent warns when probabilities are outside [0, 1] and check = TRUE", {
  expect_warning(
    format_percent(c(-0.1, 0.5, 1.1), check = TRUE),
    "outside \\[0, 1\\]"
  )
})


test_that("format_percent validates its inputs", {
  expect_error(format_percent("0.5"), "'probs' must be numeric")
  expect_error(format_percent(0.5, digits = 0), "'digits' must be a positive integer")
  expect_error(format_percent(0.5, digits = 1.5), "'digits' must be a positive integer")
  expect_error(format_percent(0.5, check = NA), "'check' must be a single logical value")
  expect_error(format_percent(0.5, check = c(TRUE, FALSE)), "'check' must be a single logical value")
})


test_that("format_test_label returns the expected display labels", {
  expect_identical(format_test_label("wald"), "Wald")
  expect_identical(format_test_label("score"), "Score")
  expect_identical(format_test_label("working-wald"), "Modified Working Wald")
  expect_identical(format_test_label("working-score"), "Modified Working Score")
  expect_identical(format_test_label("working-lrt"), "Modified Working LRT")
})


test_that("scalar check helpers accept valid inputs and reject invalid ones", {
  expect_true(is_positive_scalar(1))
  expect_true(is_positive_scalar(1.5))
  expect_false(is_positive_scalar(0))
  expect_false(is_positive_scalar(-1))
  expect_false(is_positive_scalar(NA_real_))
  expect_false(is_positive_scalar(Inf))
  expect_false(is_positive_scalar(c(1, 2)))
  expect_true(is_positive_integer_scalar(1))
  expect_true(is_positive_integer_scalar(2 + .Machine$double.eps^0.5 / 2))
  expect_false(is_positive_integer_scalar(0))
  expect_false(is_positive_integer_scalar(1.2))
  expect_false(is_positive_integer_scalar(NA_real_))
  expect_no_error(check_single_numeric(1, "x"))
  expect_error(check_single_numeric("1", "x"), "'x' must be a single finite numeric value")
  expect_error(check_single_numeric(c(1, 2), "x"), "'x' must be a single finite numeric value")
  expect_error(check_single_numeric(NA_real_, "x"), "'x' must be a single finite numeric value")
  expect_error(check_single_numeric(Inf, "x"), "'x' must be a single finite numeric value")
  expect_no_error(check_probability_open(0.5, "p"))
  expect_error(check_probability_open(0, "p"), "'p' must be strictly between 0 and 1")
  expect_error(check_probability_open(1.1, "p"), "'p' must be strictly between 0 and 1")
  expect_no_error(check_choice("robust", c("robust", "naive"), "cov_type"))
  expect_error(
    check_choice("bad", c("robust", "naive"), "cov_type"),
    "'cov_type' must be one of: robust, naive"
  )
  expect_error(
    check_choice(NA_character_, c("robust", "naive"), "cov_type"),
    "'cov_type' must be a single character value"
  )
})

test_that("summary and tidy warn about a negative variance and give NA", {
  fit <- fit_geewa_pois_exch
  fit$bias_corrected_covariance[2L, 2L] <- -1
  expect_warning(
    out <- summary(fit),
    "variance is negative or non-finite for coefficient"
  )
  expect_true(is.na(out$coefficients[2L, "Std. Error"]))
  expect_false(anyNA(out$coefficients[-2L, "Std. Error"]))
  expect_warning(td <- tidy(fit), "variance is negative or non-finite")
  expect_true(is.na(td$std.error[2L]))
  expect_error(confint(fit), "negative variance")
})

test_that("predict warns when a prediction variance is negative", {
  fit <- fit_geewa_pois_exch
  p <- ncol(fit$x)
  fit$bias_corrected_covariance[] <- -diag(p)
  expect_warning(
    out <- predict(fit, se.fit = TRUE),
    "variance is negative or non-finite"
  )
  expect_true(anyNA(out$se.fit))
})

test_that("standard_errors_or_na leaves valid variances unchanged", {
  expect_silent(se <- standard_errors_or_na(c(a = 4, b = 0, c = 9), "x"))
  expect_equal(se, c(a = 2, b = 0, c = 3))
})

test_that("jackknife covariance is cached on the fit and invalidated by coefficient changes", {
  fit <- fit_geewa_pois_indep
  expect_true(is.environment(fit$cache))
  expect_null(fit$cache$jackknife)
  v1 <- vcov(fit, cov_type = "jackknife")
  expect_false(is.null(fit$cache$jackknife))
  ## a second call must come from the cache, not from refitting
  fit$cache$jackknife$covariance[1L, 1L] <- 12345
  expect_equal(unname(vcov(fit, cov_type = "jackknife")[1L, 1L]), 12345)
  ## a changed coefficient vector bypasses the stale entry
  fit2 <- fit
  fit2$coefficients <- fit2$coefficients + 1e-3
  v2 <- vcov(fit2, cov_type = "jackknife")
  expect_false(isTRUE(all.equal(unname(v2[1L, 1L]), 12345)))
  expect_equal(dim(v1), dim(v2))
})

test_that("frechet_bounds_cor uses the shared fit-function detection", {
  fit <- fit_geewa_bin_exch
  expect_error(frechet_bounds_cor(fit), "must be fitted by 'geewa'")
})
