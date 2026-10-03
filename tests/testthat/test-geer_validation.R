testthat::local_edition(3)


test_that("check_step_thresholds accepts valid thresholds and rejects invalid ones", {
  out <- check_step_thresholds(0.10, 0.20)
  expect_type(out, "list")
  expect_identical(out$p_enter, 0.10)
  expect_identical(out$p_remove, 0.20)
  expect_error(check_step_thresholds(0, 0.20),    "'p_enter' must be strictly between 0 and 1")
  expect_error(check_step_thresholds("0.10", 0.20), "'p_enter' must be a single finite numeric value")
  expect_error(check_step_thresholds(0.10, 1.10),  "'p_remove' must be strictly between 0 and 1")
  expect_error(check_step_thresholds(0.10, "0.20"), "'p_remove' must be a single finite numeric value")
})


test_that("check_step_count accepts valid counts and rejects invalid ones", {
  expect_identical(check_step_count(10), 10L)
  expect_identical(check_step_count(0),  0L)
  expect_error(check_step_count(-1),    "'steps' must be a single nonnegative integer")
  expect_error(check_step_count(1.5),   "'steps' must be a single nonnegative integer")
  expect_error(check_step_count(c(1, 2)), "'steps' must be a single finite numeric value")
})


test_that("normalize_geer_test_options drops pmethod for non-working tests", {
  out <- normalize_geer_test_options("wald", "robust", "rao-scott")
  expect_type(out, "list")
  expect_identical(out$test, "wald")
  expect_identical(out$cov_type, "robust")
  expect_null(out$pmethod)
})


test_that("normalize_geer_test_options keeps pmethod for working tests", {
  out <- normalize_geer_test_options("working-score", "robust", "satterthwaite")
  expect_type(out, "list")
  expect_identical(out$test, "working-score")
  expect_identical(out$cov_type, "robust")
  expect_identical(out$pmethod, "satterthwaite")
})


test_that("jackknife is a valid package-wide covariance choice", {
  out <- normalize_geer_test_options("wald", "jackknife", "rao-scott")
  expect_identical(out$cov_type, "jackknife")
  expect_true("jackknife" %in% geer_cov_type_choices)
})


test_that("literal choice defaults in public signatures match the constants", {
  default_of <- function(fn, arg) eval(formals(fn)[[arg]])
  cov_fns <- list(
    summary = geer:::summary.geer,
    vcov = geer:::vcov.geer,
    confint = geer:::confint.geer,
    predict = geer:::predict.geer,
    tidy = geer:::tidy.geer,
    add1 = geer:::add1.geer,
    drop1 = geer:::drop1.geer,
    anova = geer:::anova.geer,
    step_p = step_p,
    mcar_logistic_test = mcar_logistic_test
  )
  for (name in names(cov_fns)) {
    expect_identical(
      default_of(cov_fns[[name]], "cov_type"),
      geer:::geer_cov_type_choices,
      label = paste(name, "cov_type")
    )
  }
  for (fn in list(geer:::add1.geer, geer:::drop1.geer, geer:::anova.geer,
                  step_p, mcar_logistic_test)) {
    expect_identical(default_of(fn, "test"), geer:::geer_test_choices)
    expect_identical(default_of(fn, "pmethod"), geer:::geer_pmethod_choices)
  }
  expect_identical(
    default_of(geecriteria, "cov_type"),
    geer:::geer_criteria_cov_type_choices
  )
  expect_identical(
    default_of(step_p, "direction"),
    geer:::geer_direction_choices
  )
  expect_identical(
    default_of(mcar_little_test, "reference"),
    geer:::geer_mcar_reference_choices
  )
  expect_identical(
    default_of(mcar_logistic_test, "orstr"),
    geer:::geer_mcar_orstr_choices
  )
  expect_identical(
    default_of(mcar_homoscedasticity_test, "method"),
    geer:::geer_mcar_homoscedasticity_method_choices
  )
  expect_identical(
    default_of(mcar_homoscedasticity_test, "imputation"),
    geer:::geer_mcar_imputation_choices
  )
})
