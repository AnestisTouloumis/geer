test_that("check_geer_object accepts geer fits and rejects anything else", {
  expect_identical(check_geer_object(fit_geewa_pois_exch), fit_geewa_pois_exch)
  expect_error(check_geer_object(list()), "'object' must be a 'geer' object")
  expect_error(
    check_geer_object(list(), name = "object1"),
    "'object1' must be a 'geer' object"
  )
})

test_that("check_step_thresholds accepts valid thresholds and rejects invalid ones", {
  out <- check_step_thresholds(0.10, 0.20)
  expect_type(out, "list")
  expect_identical(out$p_enter, 0.10)
  expect_identical(out$p_remove, 0.20)
  expect_error(check_step_thresholds(0, 0.20), "'p_enter' must be strictly between 0 and 1")
  expect_error(
    check_step_thresholds("0.10", 0.20),
    "'p_enter' must be a finite number"
  )
  expect_error(check_step_thresholds(0.10, 1.10), "'p_remove' must be strictly between 0 and 1")
  expect_error(
    check_step_thresholds(0.10, "0.20"),
    "'p_remove' must be a finite number"
  )
})

test_that("check_step_count accepts valid counts and rejects invalid ones", {
  expect_identical(check_step_count(10), 10L)
  expect_identical(check_step_count(0), 0L)
  expect_error(check_step_count(-1), "'steps' must be a non-negative integer")
  expect_error(check_step_count(1.5), "'steps' must be a non-negative integer")
  expect_error(check_step_count(c(1, 2)), "'steps' must be a finite number")
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

test_that("match_cov_type accepts the exact names only", {
  for (cov_type in geer_cov_type_choices) {
    expect_identical(match_cov_type(cov_type), cov_type)
  }
  expect_identical(match_cov_type(NULL), "bias-corrected")
  expect_identical(match_cov_type(NULL, default = "robust"), "robust")
  expect_error(match_cov_type("rob"), "'cov_type' must be one of")
  expect_error(match_cov_type("kauermann"), "'cov_type' must be one of")
})

test_that("literal choice defaults in public signatures match the constants", {
  default_of <- function(fn, arg) eval(formals(fn)[[arg]])
  cov_fns <- list(
    summary = summary.geer,
    vcov = vcov.geer,
    confint = confint.geer,
    predict = predict.geer,
    tidy = tidy.geer,
    add1 = add1.geer,
    drop1 = drop1.geer,
    anova = anova.geer,
    step_p = step_p,
    mcar_logistic_test = mcar_logistic_test
  )
  for (name in names(cov_fns)) {
    expect_identical(
      default_of(cov_fns[[name]], "cov_type"),
      geer_cov_type_choices[1L],
      label = paste(name, "cov_type")
    )
  }
  expect_identical(default_of(geecriteria, "cov_type"), "robust")
  expect_true("robust" %in% geer_cov_type_choices)
  for (fn in list(
    add1.geer, drop1.geer, anova.geer,
    step_p, mcar_logistic_test
  )) {
    expect_identical(default_of(fn, "test"), geer_test_choices)
    expect_identical(default_of(fn, "pmethod"), geer_pmethod_choices)
  }
  expect_identical(
    default_of(step_p, "direction"),
    geer_direction_choices
  )
})
