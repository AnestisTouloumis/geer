test_that("geewa_binary validates key inputs", {
  df <- data.frame(
    y = c(0, 1, 0, 1, 0, 1),
    x = c(0.2, -0.1, 0.4, 0.3, 0.1, 0.1),
    id = c(1, 1, 2, 2, 4, 4),
    rep = c(1, 2, 1, 2, 1, 2)
  )
  expect_error(
    geewa_binary(y ~ x, data = df, id = id, repeated = rep, link = "not-a-link"),
    "link"
  )
  expect_error(
    geewa_binary(y ~ x, data = df, id = id, repeated = rep, orstr = "not-a-orstr"),
    "orstr"
  )
  expect_error(
    geewa_binary(y ~ x, data = df, id = id, repeated = rep, orstr = "fixed"),
    "alpha_vector"
  )
  expect_error(
    geewa_binary(
      y ~ x, data = df, id = id, repeated = rep, orstr = "fixed",
      alpha_vector = c(1, 1)
    ),
    "length"
  )
  expect_error(
    geewa_binary(
      y ~ x, data = df, id = id, repeated = rep,
      beta_start = c(0, 1, 1)
    ),
    "beta_start"
  )
})


local({
  data("cerebrovascular", package = "geer")
  fit <- geewa_binary(
    formula = ecg ~ treatment + factor(period),
    id = id,
    repeated = period,
    data = cerebrovascular,
    link = "logit",
    orstr = "independence",
    method = "gee"
  )


  test_that("geewa_binary returns a well-formed geer object on a simple fit", {
    expect_s3_class(fit, "geer")
    expect_named(fit$coefficients)
    expect_true(is.matrix(fit$x))
    expect_equal(length(fit$fitted.values), nrow(fit$x))
    p <- length(fit$coefficients)
    expect_equal(dim(fit$robust_covariance), c(p, p))
  })
})
