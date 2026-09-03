testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch
binary_fit <- fit_geewa_bin_exch


test_that("update can change the fitting method", {
  fit2 <- update(count_fit, method = "brgee-robust")
  expect_s3_class(fit2, "geer")
  expect_identical(fit2$method, "brgee-robust")
  expect_identical(fit2$association_structure, count_fit$association_structure)
  expect_identical(fit2$family$family, count_fit$family$family)
  expect_identical(fit2$family$link, count_fit$family$link)
})


test_that("update can change association structure for geewa fits", {
  fit2 <- update(count_fit, corstr = "ar1")
  expect_s3_class(fit2, "geer")
  expect_identical(fit2$association_structure, "ar1")
})


test_that("update can change association structure for geewa_binary fits", {
  fit2 <- update(binary_fit, orstr = "independence")
  expect_s3_class(fit2, "geer")
  expect_identical(fit2$association_structure, "independence")
})


test_that("update can add a term to the model formula", {
  fit2 <- update(count_fit, . ~ . + treatment:lnbaseline)
  expect_s3_class(fit2, "geer")
  expect_true("treatmentprogabide:lnbaseline" %in% names(coef(fit2)))
  expect_gt(length(coef(fit2)), length(coef(count_fit)))
})


test_that("update can remove a term from the model formula", {
  fit2 <- update(count_fit, . ~ . - lnage)
  expect_s3_class(fit2, "geer")
  expect_false("lnage" %in% names(coef(fit2)))
  expect_lt(length(coef(fit2)), length(coef(count_fit)))
})


test_that("formula is stored as a formula object when passed as a symbol", {
  # Regression guard. fit$formula used to be copied from fit$call$formula, so a
  # formula reaching geewa() through a variable was stored as the symbol bound
  # to it rather than the formula itself, breaking formula(), update.formula()
  # and therefore add1()/drop1() scope handling.
  user_formula <- seizures ~ treatment + lnbaseline
  fit <- geewa(
    formula = user_formula,
    data = epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  expect_s3_class(fit$formula, "formula")
  expect_identical(fit$call$formula, quote(user_formula))
  expect_equal(fit$formula, user_formula, ignore_attr = TRUE)
})


test_that("add1 accepts a scope when the formula was passed as a symbol", {
  user_formula <- seizures ~ treatment
  fit <- geewa(
    formula = user_formula,
    data = epilepsy,
    id = id,
    family = poisson(link = "log"),
    corstr = "exchangeable",
    method = "gee"
  )
  out <- add1(fit, scope = ~ . + lnbaseline)
  expect_s3_class(out, "data.frame")
  expect_true("lnbaseline" %in% rownames(out))
})


test_that("update preserves the stored formula as a formula object", {
  fit2 <- update(count_fit, . ~ . - lnage)
  expect_s3_class(fit2$formula, "formula")
  expect_false("lnage" %in% all.vars(fit2$formula))
})
