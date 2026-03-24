testthat::local_edition(3)

fit_binary_indep <- geewa_binary(
  formula = ecg ~ period * treatment,
  id = id,
  data = test_data$cerebrovascular,
  link = "logit",
  orstr = "independence",
  method = "gee"
)

fit_binary_exch <- update(
  fit_binary_indep,
  orstr = "exchangeable"
)

cerebrovascular_small <- test_data$cerebrovascular[
  seq_len(nrow(test_data$cerebrovascular) - 1L),
  ,
  drop = FALSE
]

fit_binary_small <- geewa_binary(
  formula = ecg ~ period * treatment,
  id = id,
  data = cerebrovascular_small,
  link = "logit",
  orstr = "independence",
  method = "gee"
)

test_that("geecriteria returns the expected structure for single and multiple models", {
  out_single <- geecriteria(fit_geewa_pois_exch)
  expect_geecriteria_table(out_single, n_rows = 1L)
  expect_equal(out_single$Parameters, length(coef(fit_geewa_pois_exch)))

  out_multi <- geecriteria(
    fit_binary_indep,
    fit_binary_exch,
    cov_type = "robust",
    digits = 3
  )

  expect_geecriteria_table(
    out_multi,
    n_rows = 2L,
    row_names = c("fit_binary_indep", "fit_binary_exch")
  )
})

test_that("geecriteria warns when models differ in number of observations", {
  expect_warning(
    geecriteria(fit_binary_indep, fit_binary_small),
    regexp = "same number of observations",
    ignore.case = TRUE
  )
})

test_that("geecriteria rejects invalid inputs", {
  lm_fit <- lm(seizures ~ treatment, data = test_data$epilepsy)

  expect_error(
    geecriteria(1),
    regexp = "geer",
    ignore.case = TRUE
  )

  expect_error(
    geecriteria(fit_geewa_pois_exch, lm_fit),
    regexp = "geer|Only 'geer' objects are supported",
    ignore.case = TRUE
  )

  expect_error(
    geecriteria(fit_binary_indep, cov_type = "not-a-type"),
    regexp = "cov_type|arg",
    ignore.case = TRUE
  )

  expect_error(
    geecriteria(fit_binary_indep, digits = 1.2),
    regexp = "digits",
    ignore.case = TRUE
  )
})

test_that("geecriteria works for representative model and cov_type variants", {
  expect_geecriteria_table(geecriteria(fit_geewa_bin_exch), n_rows = 1L)

  for (cov_type in c("robust", "naive", "bias-corrected", "df-adjusted")) {
    out <- geecriteria(fit_geewa_pois_exch, cov_type = cov_type)
    expect_geecriteria_table(out, n_rows = 1L)
  }
})
