testthat::local_edition(3)

fit_bin <- fit_geewa_bin_exch

test_that("add1.geer requires a non-empty scope", {
  expect_error(add1(fit_bin), "no terms in scope")
  expect_error(add1(fit_bin, scope = NULL), "no terms in scope")
})

test_that("add1.geer returns an anova-like table with the expected structure", {
  out <- add1(
    fit_bin,
    scope = . ~ . + treatment:period,
    test = "score",
    cov_type = "robust"
  )

  expect_anova_table(
    out,
    candidate_pattern = "period:treatment"
  )
  expect_true(nrow(out) >= 2L)
})

test_that("drop1.geer returns an anova-like table for explicit and default scope", {
  out_explicit <- drop1(
    fit_bin,
    scope = "treatment",
    test = "score",
    cov_type = "robust"
  )

  expect_anova_table(
    out_explicit,
    candidate_pattern = "^treatment$",
    n_rows = 2L
  )

  out_default <- drop1(
    fit_bin,
    test = "score",
    cov_type = "robust"
  )

  expect_anova_table(out_default)
  expect_true(nrow(out_default) >= 2L)
})

test_that("add1.geer and drop1.geer require independence for working-lrt", {
  expect_error(
    add1(
      fit_bin,
      scope = . ~ . + treatment:period,
      test = "working-lrt"
    ),
    "independence"
  )

  expect_error(
    drop1(
      fit_bin,
      scope = "treatment",
      test = "working-lrt"
    ),
    "independence"
  )
})
