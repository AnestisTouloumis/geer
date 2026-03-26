testthat::local_edition(3)


test_that("format_perc returns a character vector with percent suffixes", {
  out <- format_perc(c(0.025, 0.975), digits = 3)
  expect_type(out, "character")
  expect_length(out, 2L)
  expect_true(all(grepl("%$", out)))
  expect_false(all(grepl(" ^\\s|\\s$", out)))
})


test_that("format_perc warns when probabilities are outside [0, 1] and check = TRUE", {
  expect_warning(
    format_perc(c(-0.1, 0.5, 1.1), check = TRUE),
    "outside \\[0, 1\\]"
  )
})


test_that("format_perc validates its inputs", {
  expect_error(format_perc("0.5"), "'probs' must be numeric")
  expect_error(format_perc(0.5, digits = 0), "'digits' must be a positive integer")
  expect_error(format_perc(0.5, digits = 1.5), "'digits' must be a positive integer")
  expect_error(format_perc(0.5, check = NA), "'check' must be a single logical value")
  expect_error(format_perc(0.5, check = c(TRUE, FALSE)), "'check' must be a single logical value")
})


test_that("is_pos_scalar identifies positive finite numeric scalars", {
  expect_true(is_pos_scalar(1))
  expect_true(is_pos_scalar(1.5))
  expect_false(is_pos_scalar(0))
  expect_false(is_pos_scalar(-1))
  expect_false(is_pos_scalar(NA_real_))
  expect_false(is_pos_scalar(Inf))
  expect_false(is_pos_scalar(c(1, 2)))
  expect_false(is_pos_scalar("1"))
})


test_that("is_pos_int_scalar identifies positive integer-like scalars", {
  expect_true(is_pos_int_scalar(1))
  expect_true(is_pos_int_scalar(2 + .Machine$double.eps^0.5 / 2))
  expect_false(is_pos_int_scalar(0))
  expect_false(is_pos_int_scalar(-1))
  expect_false(is_pos_int_scalar(1.2))
  expect_false(is_pos_int_scalar(NA_real_))
  expect_false(is_pos_int_scalar(c(1, 2)))
})


test_that("check_single_numeric accepts a single finite numeric value", {
  expect_no_error(check_single_numeric(1, "x"))
  expect_no_error(check_single_numeric(pi, "x"))
})


test_that("check_single_numeric rejects invalid inputs", {
  expect_error(check_single_numeric("1", "x"), "'x' must be a single finite numeric value")
  expect_error(check_single_numeric(c(1, 2), "x"), "'x' must be a single finite numeric value")
  expect_error(check_single_numeric(NA_real_, "x"), "'x' must be a single finite numeric value")
  expect_error(check_single_numeric(Inf, "x"), "'x' must be a single finite numeric value")
})


test_that("check_positive_numeric accepts positive values and rejects non-positive ones", {
  expect_no_error(check_positive_numeric(0.1, "x"))
  expect_error(check_positive_numeric(0, "x"), "'x' must be positive")
  expect_error(check_positive_numeric(-1, "x"), "'x' must be positive")
})


test_that("check_probability validates probabilities on open and closed intervals", {
  expect_no_error(check_probability(0.5, "p"))
  expect_no_error(check_probability(0, "p", open = FALSE))
  expect_no_error(check_probability(1, "p", open = FALSE))
  expect_error(check_probability(0, "p"), "'p' must be strictly between 0 and 1")
  expect_error(check_probability(1, "p"), "'p' must be strictly between 0 and 1")
  expect_error(check_probability(-0.1, "p", open = FALSE), "'p' must be between 0 and 1")
  expect_error(check_probability(1.1, "p", open = FALSE), "'p' must be between 0 and 1")
})


test_that("check_choice accepts valid choices and rejects invalid ones", {
  expect_no_error(check_choice("robust", c("robust", "naive"), "cov_type"))
  expect_error(
    check_choice("bad", c("robust", "naive"), "cov_type"),
    "'cov_type' should be one of: robust, naive"
  )
  expect_error(
    check_choice(c("robust", "naive"), c("robust", "naive"), "cov_type"),
    "'cov_type' must be a single character value"
  )
  expect_error(
    check_choice(NA_character_, c("robust", "naive"), "cov_type"),
    "'cov_type' must be a single character value"
  )
})


test_that("test_label returns the expected display labels", {
  expect_identical(test_label("wald"), "Wald")
  expect_identical(test_label("score"), "Score")
  expect_identical(test_label("working-wald"), "Modified Working Wald")
  expect_identical(test_label("working-score"), "Modified Working Score")
  expect_identical(test_label("working-lrt"), "Modified Working LRT")
})
