test_that("geecriteria: setup works", {
  expect_true(TRUE)
})


local({
  data("cerebrovascular", package = "geer")
  fit_ind <- geewa_binary(
    formula = ecg ~ period * treatment,
    id = id,
    data = cerebrovascular,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )
  fit_exch <- update(fit_ind, orstr = "exchangeable")
  cereb_small <- cerebrovascular[seq_len(nrow(cerebrovascular) - 1L), , drop = FALSE]
  fit_small <- geewa_binary(
    formula = ecg ~ period * treatment,
    id = id,
    data = cereb_small,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )
  expected_core_cols <- c("QIC", "CIC", "RJC", "QICu", "GESSC", "GPC", "Parameters")
  has_param_col <- function(x) {
    any(c("p", "npar", "n_params", "Parameters") %in% names(x))
  }


  test_that("geecriteria: single model returns expected structure", {
    out <- geecriteria(fit_ind)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 1L)
    expect_true(all(expected_core_cols %in% names(out)))
    expect_true(has_param_col(out))
    expect_true(all(vapply(out[expected_core_cols], is.numeric, logical(1))))
  })


  test_that("geecriteria: cov_type and digits validation", {
    expect_error(
      geecriteria(fit_ind, cov_type = "not-a-type"),
      regexp = "cov_type|arg",
      ignore.case = TRUE
    )
    expect_error(
      geecriteria(fit_ind, digits = -1),
      regexp = "digits",
      ignore.case = TRUE
    )
    expect_error(
      geecriteria(fit_ind, digits = 1.2),
      regexp = "digits",
      ignore.case = TRUE
    )
    expect_error(
      geecriteria(fit_ind, digits = NA_real_),
      regexp = "digits",
      ignore.case = TRUE
    )
  })


  test_that("geecriteria: errors on non-geer inputs", {
    expect_error(
      geecriteria(1),
      regexp = "geer",
      ignore.case = TRUE
    )
    expect_error(
      geecriteria(fit_ind, 2),
      regexp = "geer",
      ignore.case = TRUE
    )
  })


  test_that("geecriteria: multiple models bind rows and set rownames", {
    out <- geecriteria(fit_ind, fit_exch, cov_type = "robust", digits = 3)
    expect_s3_class(out, "data.frame")
    expect_equal(nrow(out), 2L)
    expect_true(all(expected_core_cols %in% names(out)))
    expect_true(has_param_col(out))
    expect_equal(rownames(out), c("fit_ind", "fit_exch"))
  })


  test_that("geecriteria: warns when models differ in number of observations", {
    expect_warning(
      geecriteria(fit_ind, fit_small),
      regexp = "same number of observations",
      ignore.case = TRUE
    )
  })
})
