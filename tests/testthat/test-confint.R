data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("confint.geer returns a 2-column matrix with correct row names", {
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2L)
  expect_true(all(rownames(ci) %in% names(coef(fit))))
  expect_true(all(is.finite(ci)))
})


test_that("confint.geer supports parm by name and by index", {
  nm <- names(coef(fit))[1]
  ci1 <- confint(fit, parm = nm)
  ci2 <- confint(fit, parm = 1)
  expect_equal(dim(ci1), c(1L, 2L))
  expect_equal(dim(ci2), c(1L, 2L))
  expect_identical(rownames(ci1), nm)
  expect_identical(rownames(ci2), nm)
})


test_that("confint.geer errors on invalid level and invalid parm", {
  expect_error(confint(fit, level = 1), "level")
  expect_error(confint(fit, level = 0), "level")
  expect_error(confint(fit, parm = "not_a_coef"), "parm")
})
