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
  expect_error(check_step_count(-1), "'steps' must be a nonnegative integer")
  expect_error(check_step_count(1.5), "'steps' must be a nonnegative integer")
  expect_error(check_step_count(c(1, 2)), "'steps' must be a finite number")
})

test_that("normalize_geer_test_options drops pmethod for non-working tests", {
  out <- normalize_geer_test_options("wald", "robust", "ignored-value")
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

test_that("match helpers for the new choice sets require exact names", {
  expect_identical(match_alternative("less"), "less")
  expect_identical(match_mcar_reference("asymptotic"), "asymptotic")
  expect_identical(match_mcar_method("hawkins"), "hawkins")
  expect_identical(match_mcar_imputation("normal"), "normal")
  expect_identical(match_mcar_orstr("unstructured"), "unstructured")
  expect_identical(match_influence_unit("observations"), "observations")
  expect_error(match_alternative("two"), "'alternative' must be one of")
  expect_error(match_mcar_reference("asym"), "'reference' must be one of")
  expect_error(match_mcar_method("haw"), "'method' must be one of")
  expect_error(match_mcar_imputation("norm"), "'imputation' must be one of")
  expect_error(match_mcar_orstr("fixed"), "'orstr' must be one of")
  expect_error(match_influence_unit("clus"), "'unit' must be one of")
})

test_that("check_working_lrt_allowed requires independence", {
  expect_true(check_working_lrt_allowed("working-lrt", "independence"))
  expect_true(check_working_lrt_allowed("wald", "exchangeable"))
  expect_error(
    check_working_lrt_allowed("working-lrt", "exchangeable"),
    "requires an independence working structure"
  )
})

test_that("check_integer_at_least accepts integers and rejects anything else", {
  expect_identical(check_integer_at_least(3, "n"), 3L)
  for (x in list(0, "10", 1e12)) {
    expect_error(
      check_integer_at_least(x, "n"),
      "'n' must be a positive integer"
    )
  }
  expect_error(
    check_integer_at_least(1, "n", lower = 2L),
    "'n' must be an integer greater than or equal to 2"
  )
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
  for (x in list("1", c(1, 2), NA_real_, Inf)) {
    expect_error(check_single_numeric(x, "x"), "'x' must be a finite number")
  }
  expect_no_error(check_probability_open(0.5, "p"))
  for (p in c(0, 1.1)) {
    expect_error(
      check_probability_open(p, "p"),
      "'p' must be strictly between 0 and 1"
    )
  }
  expect_no_error(check_choice("robust", c("robust", "naive"), "cov_type"))
  expect_error(
    check_choice("bad", c("robust", "naive"), "cov_type"),
    "'cov_type' must be one of: robust, naive"
  )
  expect_error(
    check_choice(NA_character_, c("robust", "naive"), "cov_type"),
    "'cov_type' must be a character string"
  )
})

test_that("integer count validators reject overflow without coercion", {
  expect_identical(
    check_nonnegative_integerish(.Machine$integer.max, "n"),
    .Machine$integer.max
  )
  expect_identical(check_nonnegative_integerish(0, "n"), 0L)
  expect_error(
    check_nonnegative_integerish(2^31, "n"),
    "'n' must be a nonnegative integer", fixed = TRUE
  )
  expect_error(check_step_count(2^31))
})

test_that("fitting choice matchers are exact and retain strict NULL handling", {
  matchers <- list(
    method = match_method, family = match_family, link = match_link,
    corstr = match_corstr, orstr = match_orstr, action = match_step_action
  )
  choices <- list(
    method = geer_method_choices, family = geer_family_choices,
    link = geer_link_choices, corstr = geer_corstr_choices,
    orstr = geer_orstr_choices, action = geer_step_actions
  )
  for (name in names(matchers)) {
    matcher <- matchers[[name]]
    for (value in choices[[name]]) {
      expect_identical(matcher(value), value)
    }
    expect_error(matcher(NULL))
    expect_error(matcher(factor(choices[[name]][1L])))
    expect_error(matcher(paste0(choices[[name]][1L], "-invalid")))
  }
})
