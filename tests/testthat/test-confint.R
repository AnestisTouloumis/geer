testthat::local_edition(3)

fit <- fit_geewa_bin_exch

test_that("confint.geer returns a well-formed confidence interval matrix", {
  ci <- confint(fit)

  expect_true(is.matrix(ci))
  expect_equal(dim(ci)[2], 2L)
  expect_identical(colnames(ci), c(" 2.5%", "97.5%"))
  expect_identical(rownames(ci), names(coef(fit)))
  expect_true(all(is.finite(ci)))
})

test_that("confint.geer supports selecting parameters by name and index", {
  coef_names <- names(coef(fit))[1:2]

  ci_by_name <- confint(fit, parm = coef_names)
  ci_by_index <- confint(fit, parm = 1:2)

  expect_equal(dim(ci_by_name), c(2L, 2L))
  expect_equal(dim(ci_by_index), c(2L, 2L))
  expect_identical(rownames(ci_by_name), coef_names)
  expect_identical(rownames(ci_by_index), coef_names)
  expect_equal(unname(ci_by_name), unname(ci_by_index))
})

test_that("confint.geer validates level and parm", {
  expect_error(confint(fit, level = 1), "level")
  expect_error(confint(fit, parm = "not_a_coef"), "parm")
})
