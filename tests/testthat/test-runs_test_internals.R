testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch


## epilepsy is sorted by id and visit, so a column of the data used to fit the
## model is already in the natural cluster/repeated order.
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


## ----------------------------------------------------- resolve_runs_order ----

test_that("resolve_runs_order requires a single non-empty character value", {
  for (bad in list(character(0), NA_character_, "", c("natural", "fitted"))) {
    expect_error(
      geer:::resolve_runs_order(count_fit, bad),
      "'order_by' must be a single non-empty character value or a numeric vector"
    )
  }
})


test_that("resolve_runs_order rejects values that are neither text nor numeric", {
  expect_error(
    geer:::resolve_runs_order(count_fit, TRUE),
    "'order_by' must be a character value or a numeric vector"
  )
  expect_error(
    geer:::resolve_runs_order(count_fit, list(1, 2)),
    "'order_by' must be a character value or a numeric vector"
  )
})


test_that("resolve_runs_order flags the natural order", {
  natural <- geer:::resolve_runs_order(count_fit, "natural")
  expect_true(natural$natural)
  expect_identical(natural$index, order(count_fit$id, count_fit$repeated))

  ## A numeric key that is constant leaves the natural order unchanged.
  tied <- geer:::resolve_runs_order(count_fit, rep(1, count_fit$obs_no))
  expect_true(tied$natural)
  expect_identical(tied$label, "supplied ordering vector")
})


## ------------------------------------------------ extract_runs_order_variable ----

test_that("extract_runs_order_variable needs the call and formula of the fit", {
  fit <- count_fit
  fit$call <- NULL

  expect_error(
    geer:::extract_runs_order_variable(fit, "visit"),
    "unknown 'order_by' value 'visit'"
  )
})


test_that("extract_runs_order_variable falls back to the supplied environment", {
  fit <- count_fit
  formula_without_env <- fit$formula
  attr(formula_without_env, ".Environment") <- NULL
  fit$formula <- formula_without_env

  key <- geer:::extract_runs_order_variable(fit, "visit", env = environment())

  expect_equal(key, as.numeric(epilepsy$visit))
})


test_that("extract_runs_order_variable returns numeric keys for common types", {
  data <- runs_order_data()
  fit <- fit_runs_order(data)

  expect_equal(
    geer:::extract_runs_order_variable(fit, "visit"),
    as.numeric(data$visit)
  )
  expect_equal(
    geer:::extract_runs_order_variable(fit, "date"),
    as.numeric(data$date)
  )
  expect_equal(
    geer:::extract_runs_order_variable(fit, "stage"),
    as.numeric(data$stage)
  )
  expect_equal(
    geer:::extract_runs_order_variable(fit, "late"),
    as.numeric(factor(data$late))
  )
  expect_equal(
    geer:::extract_runs_order_variable(fit, "label"),
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
    geer:::extract_runs_order_variable(fit, "complex_key"),
    "'order_by' variable 'complex_key' must be numeric, a factor, or date-like"
  )
  expect_error(
    geer:::extract_runs_order_variable(fit, "no_such_column"),
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
  fit <- count_fit
  fit$id <- rev(fit$id)

  expect_error(
    geer:::extract_runs_order_variable(fit, "visit"),
    "cannot be aligned with the observations stored in the fitted model"
  )
})


## ------------------------------------------------- compute_runs_statistics ----

test_that("compute_runs_statistics validates the residuals", {
  expect_error(
    geer:::compute_runs_statistics(c("a", "b"), "two.sided"),
    "residuals must be a non-empty numeric vector"
  )
  expect_error(
    geer:::compute_runs_statistics(numeric(0), "two.sided"),
    "residuals must be a non-empty numeric vector"
  )
  expect_error(
    geer:::compute_runs_statistics(c(1, NA, -1), "two.sided"),
    "residuals must all be finite and non-missing"
  )
  expect_error(
    geer:::compute_runs_statistics(c(1, Inf, -1), "two.sided"),
    "residuals must all be finite and non-missing"
  )
})


test_that("compute_runs_statistics validates the alternative", {
  residual_values <- c(1, -1, 1, 1, -1, -1, 1)

  expect_error(
    geer:::compute_runs_statistics(residual_values, 1),
    "'alternative' must be a single character value"
  )
  expect_error(
    geer:::compute_runs_statistics(residual_values, NA_character_),
    "'alternative' must be a single character value"
  )
  expect_error(
    geer:::compute_runs_statistics(residual_values, c("less", "greater")),
    "'alternative' must be a single character value"
  )
  expect_error(
    geer:::compute_runs_statistics(residual_values, "sideways"),
    "'alternative' must be one of 'two.sided', 'less' or 'greater'"
  )
})


test_that("one-sided p-values are complementary", {
  residual_values <- c(1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1)

  less <- geer:::compute_runs_statistics(residual_values, "less")
  greater <- geer:::compute_runs_statistics(residual_values, "greater")
  two_sided <- geer:::compute_runs_statistics(residual_values, "two.sided")

  expect_equal(less$p_value + greater$p_value, 1, tolerance = 1e-12)
  expect_equal(
    two_sided$p_value,
    2 * min(less$p_value, greater$p_value),
    tolerance = 1e-12
  )
})


test_that("compute_runs_statistics rejects a degenerate sign sequence", {
  ## One positive and one negative residual give a zero variance.
  expect_error(
    geer:::compute_runs_statistics(c(1, -1), "two.sided"),
    "the runs-test variance is not positive"
  )
  expect_error(
    geer:::compute_runs_statistics(c(1, 2, 3), "two.sided"),
    "at least one positive and one negative residual"
  )
})
