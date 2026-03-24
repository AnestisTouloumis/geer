testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch
binary_fit <- fit_geewa_bin_exch

test_that("predict.geer returns the expected in-sample predictions", {
  expect_equal(
    predict(count_fit, type = "response"),
    fitted(count_fit),
    tolerance = 1e-10
  )

  expect_equal(
    predict(count_fit, type = "link"),
    count_fit$linear.predictors,
    tolerance = 1e-10
  )

  pred_bin <- predict(binary_fit, type = "response")
  expect_true(all(is.finite(pred_bin)))
  expect_true(all(pred_bin >= 0 & pred_bin <= 1))
})

test_that("predict.geer returns fit and standard errors when requested", {
  pred_in <- predict(count_fit, type = "response", se.fit = TRUE)

  expect_type(pred_in, "list")
  expect_named(pred_in, c("fit", "se.fit"))
  expect_equal(pred_in$fit, predict(count_fit, type = "response"))
  expect_length(pred_in$fit, count_fit$obs_no)
  expect_length(pred_in$se.fit, count_fit$obs_no)
  expect_true(all(is.finite(pred_in$se.fit)))
  expect_true(all(pred_in$se.fit >= 0))

  newdata_bin <- test_data$cerebrovascular[1:8, , drop = FALSE]
  pred_bin <- predict(
    binary_fit,
    newdata = newdata_bin,
    type = "response",
    se.fit = TRUE
  )

  expect_type(pred_bin, "list")
  expect_named(pred_bin, c("fit", "se.fit"))
  expect_equal(
    pred_bin$fit,
    predict(binary_fit, newdata = newdata_bin, type = "response")
  )
  expect_length(pred_bin$fit, nrow(newdata_bin))
  expect_length(pred_bin$se.fit, nrow(newdata_bin))
  expect_true(all(is.finite(pred_bin$fit)))
  expect_true(all(is.finite(pred_bin$se.fit)))
  expect_true(all(pred_bin$fit >= 0 & pred_bin$fit <= 1))
  expect_true(all(pred_bin$se.fit >= 0))
})

test_that("predict.geer returns coherent predictions on new count data", {
  newdata <- test_data$epilepsy[1:5, , drop = FALSE]

  pred_response <- predict(count_fit, newdata = newdata, type = "response")
  pred_link <- predict(count_fit, newdata = newdata, type = "link")
  pred_link_se <- predict(count_fit, newdata = newdata, type = "link", se.fit = TRUE)

  expect_length(pred_response, nrow(newdata))
  expect_true(all(is.finite(pred_response)))
  expect_true(all(pred_response > 0))

  expect_equal(log(pred_response), pred_link, tolerance = 1e-10)

  expect_type(pred_link_se, "list")
  expect_named(pred_link_se, c("fit", "se.fit"))
  expect_equal(pred_link_se$fit, pred_link)
  expect_length(pred_link_se$fit, nrow(newdata))
  expect_length(pred_link_se$se.fit, nrow(newdata))
  expect_true(all(is.finite(pred_link_se$fit)))
  expect_true(all(is.finite(pred_link_se$se.fit)))
  expect_true(all(pred_link_se$se.fit >= 0))
})

test_that("predict.geer handles reordered and extra columns in newdata", {
  newdata_1 <- test_data$epilepsy[
    1:5,
    c("seizures", "treatment", "lnbaseline", "lnage", "id", "visit"),
    drop = FALSE
  ]
  newdata_2 <- test_data$epilepsy[
    1:5,
    c("visit", "id", "lnage", "lnbaseline", "treatment", "seizures"),
    drop = FALSE
  ]

  pred_1 <- predict(count_fit, newdata = newdata_1, type = "link")
  pred_2 <- predict(count_fit, newdata = newdata_2, type = "link")

  expect_equal(pred_1, pred_2)

  newdata_extra <- newdata_1
  newdata_extra$junk <- seq_len(nrow(newdata_extra))

  expect_equal(
    predict(count_fit, newdata = newdata_extra, type = "link"),
    pred_1
  )
})

test_that("predict.geer validates newdata and prediction type", {
  newdata_missing <- subset(
    test_data$epilepsy[1:5, , drop = FALSE],
    select = -lnage
  )

  expect_error(
    predict(count_fit, newdata = newdata_missing),
    "object 'lnage' not found",
    fixed = FALSE
  )

  expect_error(
    predict(count_fit, type = "invalid-type"),
    "'arg' should be one of"
  )
})
