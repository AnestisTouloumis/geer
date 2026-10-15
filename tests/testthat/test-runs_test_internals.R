epilepsy <- test_data$epilepsy

runs_order_data <- function() {
  data <- epilepsy
  data$date <- as.Date("2020-01-01") + 7 * data$visit
  data$label <- c("first", "second", "third", "fourth")[data$visit]
  data$late <- data$visit > 2
  data$stage <- factor(data$label, levels = c("first", "second", "third", "fourth"))
  data$complex_key <- complex(real = 1, imaginary = data$visit)
  data
}

fit_runs_order <- function(data) {
  geewa(
    formula = seizures ~ treatment + lnbaseline + lnage,
    data = data,
    id = id,
    family = poisson(link = "log"),
    corstr = "independence",
    method = "gee"
  )
}

test_that("resolve_runs_order requires a single non-empty character value", {
  for (bad in list(character(0), NA_character_, "", c("natural", "fitted"))) {
    expect_error(
      resolve_runs_order(fit_geewa_pois_exch, bad),
      "'order_by' must be a single non-empty character string or a numeric vector"
    )
  }
})

test_that("resolve_runs_order rejects values that are neither text nor numeric", {
  expect_error(
    resolve_runs_order(fit_geewa_pois_exch, TRUE),
    "'order_by' must be a single non-empty character string or a numeric vector"
  )
  expect_error(
    resolve_runs_order(fit_geewa_pois_exch, list(1, 2)),
    "'order_by' must be a single non-empty character string or a numeric vector"
  )
})

test_that("resolve_runs_order flags the natural order", {
  natural <- resolve_runs_order(fit_geewa_pois_exch, "natural")
  expect_true(natural$natural)
  expect_identical(natural$index, order(fit_geewa_pois_exch$id, fit_geewa_pois_exch$repeated))

  tied <- resolve_runs_order(fit_geewa_pois_exch, rep(1, fit_geewa_pois_exch$obs_no))
  expect_true(tied$natural)
  expect_identical(tied$label, "supplied ordering vector")
})

test_that("extract_runs_order_variable needs the call and formula of the fit", {
  fit <- fit_geewa_pois_exch
  fit$call <- NULL

  expect_error(
    extract_runs_order_variable(fit, "visit"),
    "unknown 'order_by' value 'visit'"
  )
})

test_that("extract_runs_order_variable falls back to the supplied environment", {
  fit <- fit_geewa_pois_exch
  formula_without_env <- fit$formula
  attr(formula_without_env, ".Environment") <- NULL
  fit$formula <- formula_without_env

  key <- extract_runs_order_variable(fit, "visit", env = environment())

  expect_equal(key, as.numeric(epilepsy$visit))
})

test_that("extract_runs_order_variable returns numeric keys for common types", {
  data <- runs_order_data()
  fit <- fit_runs_order(data)

  expect_equal(
    extract_runs_order_variable(fit, "visit"),
    as.numeric(data$visit)
  )
  expect_equal(
    extract_runs_order_variable(fit, "date"),
    as.numeric(data$date)
  )
  expect_equal(
    extract_runs_order_variable(fit, "stage"),
    as.numeric(data$stage)
  )
  expect_equal(
    extract_runs_order_variable(fit, "late"),
    as.numeric(factor(data$late))
  )
  expect_equal(
    extract_runs_order_variable(fit, "label"),
    as.numeric(factor(data$label))
  )
})

test_that("runs_test accepts date, factor and character ordering variables", {
  data <- runs_order_data()
  fit <- fit_runs_order(data)

  by_date <- runs_test(fit, order_by = "date")
  by_visit <- runs_test(fit, order_by = "visit")

  expect_identical(by_date$order_by, "variable 'date'")
  expect_identical(by_date$runs, by_visit$runs)
  expect_identical(
    runs_test(fit, order_by = "stage")$runs,
    by_visit$runs
  )
  expect_identical(
    runs_test(fit, order_by = "stage")$order_by,
    "variable 'stage'"
  )
})

test_that("extract_runs_order_variable rejects keys that are not orderable", {
  fit <- fit_runs_order(runs_order_data())

  expect_error(
    extract_runs_order_variable(fit, "complex_key"),
    "'order_by' variable 'complex_key' must be numeric, a factor, or date-like"
  )
  expect_error(
    extract_runs_order_variable(fit, "no_such_column"),
    "unknown 'order_by' value 'no_such_column'"
  )
})

test_that("an ordering variable that drops observations cannot be aligned", {
  data <- runs_order_data()
  data$partial <- seq_len(nrow(data))
  data$partial[5] <- NA_real_
  fit <- fit_runs_order(data)

  expect_error(
    runs_test(fit, order_by = "partial"),
    "including it changes the number of complete cases"
  )
})

test_that("an ordering variable must align with the stored observations", {
  fit <- fit_geewa_pois_exch
  fit$id <- rev(fit$id)

  expect_error(
    extract_runs_order_variable(fit, "visit"),
    "cannot be aligned with the observations stored in the fitted model"
  )
})

test_that("compute_runs_statistics validates the residuals", {
  expect_error(
    compute_runs_statistics(c("a", "b"), "two.sided"),
    "residuals must be a non-empty numeric vector"
  )
  expect_error(
    compute_runs_statistics(numeric(0), "two.sided"),
    "residuals must be a non-empty numeric vector"
  )
  expect_error(
    compute_runs_statistics(c(1, NA, -1), "two.sided"),
    "residuals must all be finite and non-missing"
  )
  expect_error(
    compute_runs_statistics(c(1, Inf, -1), "two.sided"),
    "residuals must all be finite and non-missing"
  )
})

test_that("compute_runs_statistics validates the alternative", {
  residual_values <- c(1, -1, 1, 1, -1, -1, 1)

  expect_error(
    compute_runs_statistics(residual_values, 1),
    "'alternative' must be a character string"
  )
  expect_error(
    compute_runs_statistics(residual_values, NA_character_),
    "'alternative' must be a character string"
  )
  expect_error(
    compute_runs_statistics(residual_values, c("less", "greater")),
    "'alternative' must be a character string"
  )
  expect_error(
    compute_runs_statistics(residual_values, "sideways"),
    "'alternative' must be one of: two.sided, less, greater"
  )
})

test_that("compute_runs_statistics rejects a degenerate sign sequence", {
  expect_error(
    compute_runs_statistics(c(1, -1), "two.sided"),
    "the runs-test variance is not positive"
  )
  expect_error(
    compute_runs_statistics(c(1, 2, 3), "two.sided"),
    "at least one positive and one negative residual"
  )
})
