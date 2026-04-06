testthat::local_edition(3)

test_that("geewa_binary rejects unsupported link and odds-ratio structure", {
  expect_error(
    geewa_binary(
      ecg ~ treatment,
      id = id,
      data = test_data$cerebrovascular,
      link = "log-log"
    ),
    "'log-log' link not recognised"
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

test_that("geewa_binary rejects duplicated repeated values within cluster", {
  cerebrovascular_bad <- test_data$cerebrovascular
  cerebrovascular_bad$rep2 <- cerebrovascular_bad$period

  idx <- which(cerebrovascular_bad$id == cerebrovascular_bad$id[1])[1:2]
  cerebrovascular_bad$rep2[idx] <- 1

  expect_error(
    geewa_binary(
      ecg ~ treatment + factor(period),
      id = id,
      repeated = rep2,
      data = cerebrovascular_bad,
      link = "logit",
      orstr = "independence"
    ),
    "unique",
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
