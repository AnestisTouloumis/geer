data("cerebrovascular", package = "geer")
fit <- geewa_binary(
  formula = ecg ~ treatment + factor(period),
  id = id,
  data = cerebrovascular,
  link = "logit",
  orstr = "exchangeable"
)


test_that("model.matrix.geer returns a matrix matching the stored design matrix", {
  X <- model.matrix(fit)
  expect_true(is.matrix(X))
  expect_true(is.numeric(X))
  expect_equal(nrow(X), nrow(fit$x))
  expect_equal(ncol(X), ncol(fit$x))
  expect_identical(colnames(X), colnames(fit$x))
  expect_true(!is.null(attr(X, "assign")))
})


test_that("model.matrix.geer errors for non-geer objects", {
  expect_error(
    model.matrix.geer(list(terms = terms(~ 1), data = data.frame(x = 1))),
    "must be a 'geer' object"
  )
})
