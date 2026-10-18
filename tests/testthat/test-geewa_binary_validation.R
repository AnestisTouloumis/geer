cerebrovascular <- test_data$cerebrovascular

test_that("geewa_binary rejects unsupported link and odds-ratio structure", {
  expect_error(
    geewa_binary(
      ecg ~ treatment,
      id = id,
      data = cerebrovascular,
      link = "log-log"
    ),
    "'log-log' link not recognised"
  )
  expect_error(
    geewa_binary(
      ecg ~ treatment,
      id = id,
      data = cerebrovascular,
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
      data = cerebrovascular,
      link = "logit",
      weights = rep.int(1, nrow(cerebrovascular) - 1)
    ),
    "weights"
  )
  bad_weights <- rep.int(1, nrow(cerebrovascular))
  bad_weights[1] <- 0
  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      data = cerebrovascular,
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
      data = cerebrovascular,
      link = "logit",
      orstr = "fixed"
    ),
    "alpha_vector",
    ignore.case = TRUE
  )
})

test_that("geewa_binary rejects a non-finite beta_start", {
  cer <- cerebrovascular
  expect_error(
    geewa_binary(ecg ~ treatment,
      data = cer, id = id, repeated = period,
      beta_start = c(NA, 0)
    ),
    "'beta_start' must contain only finite values"
  )
})
