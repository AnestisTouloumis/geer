test_that("geewa_binary rejects unsupported link and odds-ratio structure", {
  expect_error(
    geewa_binary(
      ecg ~ treatment,
      id = id,
      data = test_data$cerebrovascular,
      link = "log-log"
    ),
    "'link' must be one of"
  )
  expect_error(
    geewa_binary(
      ecg ~ treatment,
      id = id,
      data = test_data$cerebrovascular,
      link = "logit",
      orstr = "ar1"
    ),
    "'orstr' must be one of"
  )
})

test_that("geewa_binary validates weights", {
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      data = test_data$cerebrovascular,
      link = "logit",
      weights = rep.int(1, nrow(test_data$cerebrovascular) - 1)
    ),
    "weights"
  )
  bad_weights <- rep.int(1, nrow(test_data$cerebrovascular))
  bad_weights[1] <- 0
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      data = test_data$cerebrovascular,
      link = "logit",
      weights = bad_weights
    ),
    "positive|strictly",
    ignore.case = TRUE
  )
})

test_that("geewa_binary requires alpha_vector for fixed odds-ratio structure", {
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      repeated = period,
      data = test_data$cerebrovascular,
      link = "logit",
      orstr = "fixed"
    ),
    "alpha_vector",
    ignore.case = TRUE
  )
})

test_that("geewa_binary validates fixed odds ratios", {
  fit_fixed <- function(alpha_vector) {
    geewa_binary(
      formula = ecg ~ treatment,
      data = test_data$cerebrovascular,
      id = id,
      repeated = period,
      orstr = "fixed",
      alpha_vector = alpha_vector
    )
  }
  expect_error(
    fit_fixed(c(1, 2)),
    "'alpha_vector' must be a numeric vector of length 1"
  )
  expect_error(
    fit_fixed(-1),
    "'alpha_vector' must be finite and strictly positive"
  )
})
