testthat::local_edition(3)

binary_start_data <- data.frame(
  y = c(0, 1, 0, 1, 1, 0, 1, 0, 0, 1,
        1, 0, 1, 0, 1, 0, 0, 1, 0, 1,
        1, 0, 0, 1, 0, 1, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 0, 1, 0, 0, 1),
  x1 = c(-0.63, 0.18, -0.84, 1.60, 0.33, -0.82, 0.49, 0.74, 0.58, -0.31,
         1.51, 0.39, -0.62, -2.21, 1.12, -0.04, -0.02, 0.94, 0.82, 0.59,
         -0.12, -0.03, -0.69, -0.45, 1.10, -1.12, 0.73, -0.09, -0.35, -1.47,
         0.13, 1.72, -0.69, 0.53, -0.18, -1.40, -0.18, 1.30, -0.28, -0.78),
  x2 = c(1, 0, 1, 0, 1, 1, 0, 1, 0, 0,
         1, 1, 0, 1, 0, 1, 0, 1, 1, 0,
         1, 0, 1, 0, 1, 0, 1, 1, 0, 0,
         1, 0, 1, 1, 0, 1, 0, 1, 0, 1),
  id = rep(1:10, each = 4),
  time = rep(1:4, times = 10)
)

test_that("geewa computes automatic starting values and returns a valid fit", {
  fit <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    repeated = visit
  )

  expect_s3_class(fit, "geer")
  expect_equal(
    length(coef(fit)),
    ncol(model.matrix(
      seizures ~ treatment + lnbaseline + lnage,
      data = test_data$epilepsy
    ))
  )
  expect_true(is.numeric(fit$coefficients))
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("geewa accepts a valid user-supplied beta_start", {
  fit0 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    repeated = visit
  )

  fit1 <- geewa(
    seizures ~ treatment + lnbaseline + lnage,
    data = test_data$epilepsy,
    id = id,
    repeated = visit,
    beta_start = coef(fit0)
  )

  expect_s3_class(fit1, "geer")
  expect_equal(length(coef(fit1)), length(coef(fit0)))
  expect_true(all(is.finite(coef(fit1))))
})

test_that("geewa rejects beta_start of the wrong length", {
  expect_error(
    geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = test_data$epilepsy,
      id = id,
      repeated = visit,
      beta_start = c(0, 0)
    ),
    "beta_start"
  )
})

test_that("geewa_binary computes automatic starting values and returns a valid fit", {
  fit <- geewa_binary(
    y ~ x1 + x2,
    data = binary_start_data,
    id = id,
    repeated = time
  )

  expect_s3_class(fit, "geer")
  expect_equal(
    length(coef(fit)),
    ncol(model.matrix(y ~ x1 + x2, data = binary_start_data))
  )
  expect_true(is.numeric(fit$coefficients))
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("geewa_binary accepts a valid user-supplied beta_start", {
  fit0 <- geewa_binary(
    y ~ x1 + x2,
    data = binary_start_data,
    id = id,
    repeated = time
  )

  fit1 <- geewa_binary(
    y ~ x1 + x2,
    data = binary_start_data,
    id = id,
    repeated = time,
    beta_start = coef(fit0)
  )

  expect_s3_class(fit1, "geer")
  expect_equal(length(coef(fit1)), length(coef(fit0)))
  expect_true(all(is.finite(coef(fit1))))
})

test_that("geewa_binary rejects beta_start of the wrong length", {
  expect_error(
    geewa_binary(
      y ~ x1 + x2,
      data = binary_start_data,
      id = id,
      repeated = time,
      beta_start = c(0, 0)
    ),
    "beta_start"
  )
})
