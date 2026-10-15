cerebrovascular <- test_data$cerebrovascular
epilepsy <- test_data$epilepsy

test_that("predict returns the expected in-sample predictions", {
  expect_equal(
    predict(fit_geewa_pois_exch, type = "response"),
    fitted(fit_geewa_pois_exch),
    tolerance = 1e-10
  )
  expect_equal(
    predict(fit_geewa_pois_exch, type = "link"),
    fit_geewa_pois_exch$linear.predictors,
    tolerance = 1e-10
  )
  pred_bin <- predict(fit_geewa_bin_exch, type = "response")
  expect_true(all(is.finite(pred_bin)))
  expect_true(all(pred_bin >= 0 & pred_bin <= 1))
})

test_that("predict returns fit and standard errors when requested", {
  pred_in <- predict(fit_geewa_pois_exch, type = "response", se.fit = TRUE)
  expect_type(pred_in, "list")
  expect_named(pred_in, c("fit", "se.fit"))
  expect_equal(pred_in$fit, predict(fit_geewa_pois_exch, type = "response"))
  expect_length(pred_in$fit, fit_geewa_pois_exch$obs_no)
  expect_length(pred_in$se.fit, fit_geewa_pois_exch$obs_no)
  expect_true(all(is.finite(pred_in$se.fit)))
  expect_true(all(pred_in$se.fit >= 0))
  newdata_bin <- cerebrovascular[1:8, , drop = FALSE]
  pred_bin <- predict(
    fit_geewa_bin_exch,
    newdata = newdata_bin,
    type = "response",
    se.fit = TRUE
  )
  expect_type(pred_bin, "list")
  expect_named(pred_bin, c("fit", "se.fit"))
  expect_equal(
    pred_bin$fit,
    predict(fit_geewa_bin_exch, newdata = newdata_bin, type = "response")
  )
  expect_length(pred_bin$fit, nrow(newdata_bin))
  expect_length(pred_bin$se.fit, nrow(newdata_bin))
  expect_true(all(is.finite(pred_bin$fit)))
  expect_true(all(is.finite(pred_bin$se.fit)))
  expect_true(all(pred_bin$fit >= 0 & pred_bin$fit <= 1))
  expect_true(all(pred_bin$se.fit >= 0))
})

test_that("predict returns coherent predictions on new count data", {
  newdata <- epilepsy[1:5, , drop = FALSE]
  pred_response <- predict(fit_geewa_pois_exch, newdata = newdata, type = "response")
  pred_link <- predict(fit_geewa_pois_exch, newdata = newdata, type = "link")
  pred_link_se <- predict(fit_geewa_pois_exch, newdata = newdata, type = "link", se.fit = TRUE)
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

test_that("predict handles reordered and extra columns in newdata", {
  newdata_1 <- epilepsy[
    1:5,
    c("seizures", "treatment", "lnbaseline", "lnage", "id", "visit"),
    drop = FALSE
  ]
  newdata_2 <- epilepsy[
    1:5,
    c("visit", "id", "lnage", "lnbaseline", "treatment", "seizures"),
    drop = FALSE
  ]
  pred_1 <- predict(fit_geewa_pois_exch, newdata = newdata_1, type = "link")
  pred_2 <- predict(fit_geewa_pois_exch, newdata = newdata_2, type = "link")
  expect_equal(pred_1, pred_2)
  newdata_extra <- newdata_1
  newdata_extra$junk <- seq_len(nrow(newdata_extra))
  expect_equal(
    predict(fit_geewa_pois_exch, newdata = newdata_extra, type = "link"),
    pred_1
  )
})

test_that("predict validates newdata and prediction type", {
  newdata_missing <- subset(
    epilepsy[1:5, , drop = FALSE],
    select = -lnage
  )
  expect_error(
    predict(fit_geewa_pois_exch, newdata = newdata_missing),
    "object 'lnage' not found",
    fixed = FALSE
  )
  expect_error(
    predict(fit_geewa_pois_exch, type = "invalid-type"),
    "'arg' should be one of"
  )
})

test_that("predict validates se.fit and the stored coefficients", {
  for (bad in list(NA, c(TRUE, FALSE), "yes")) {
    expect_error(
      predict(fit_geewa_pois_exch, se.fit = bad),
      "'se.fit' must be TRUE or FALSE"
    )
  }
  unnamed <- fit_geewa_pois_exch
  names(unnamed$coefficients) <- NULL
  expect_error(predict(unnamed), "'coefficients' must be named")
})

test_that("predict accepts newdata given as a list", {
  newdata <- epilepsy[1:6, ]
  expect_equal(
    predict(fit_geewa_pois_exch, newdata = as.list(newdata)),
    predict(fit_geewa_pois_exch, newdata = newdata)
  )
})

test_that("prediction offsets are taken from newdata and must be finite", {
  log_weeks <- rep(log(2), nrow(epilepsy))
  fit_formula <- geewa(
    seizures ~ treatment + offset(log_weeks),
    data = epilepsy, id = id, family = poisson()
  )
  expect_error(
    predict(fit_formula, newdata = epilepsy),
    "supplied via the model formula are not present in 'newdata'"
  )
  newdata <- epilepsy[1:4, ]
  newdata$log_weeks <- c(log(2), Inf, log(2), log(2))
  expect_error(
    predict(fit_formula, newdata = newdata),
    "formula-based prediction offset must be a finite numeric vector"
  )
  fit_argument <- geewa(
    seizures ~ treatment,
    data = epilepsy, id = id, family = poisson(), offset = log_weeks
  )
  expect_error(
    predict(fit_argument, newdata = epilepsy),
    "supplied via the 'offset' argument are not present in 'newdata'"
  )
  newdata$log_weeks <- c(NA, log(2), log(2), log(2))
  expect_error(
    predict(fit_argument, newdata = newdata),
    "prediction offset supplied via the 'offset' argument must be a finite"
  )
  fit_scalar <- geewa(
    seizures ~ treatment,
    data = epilepsy, id = id, family = poisson(), offset = 0.5
  )
  newdata <- epilepsy[1:4, ]
  expected <- drop(model.matrix(~treatment, newdata) %*% coef(fit_scalar)) + 0.5
  expect_equal(unname(predict(fit_scalar, newdata = newdata)), unname(expected))
})

test_that("residuals return finite vectors for supported types", {
  for (type in c("working", "pearson", "deviance")) {
    r <- residuals(fit_geewa_pois_exch, type = type)
    expect_length(r, fit_geewa_pois_exch$obs_no)
    expect_true(all(is.finite(r)))
  }
  r_mahalanobis <- residuals(fit_geewa_pois_exch, type = "mahalanobis")
  expect_length(r_mahalanobis, fit_geewa_pois_exch$clusters_no)
  expect_true(all(is.finite(r_mahalanobis)))
  expect_true(all(r_mahalanobis >= 0))
  expect_named(r_mahalanobis, as.character(unique(fit_geewa_pois_exch$id)))

  r_bin <- residuals(fit_geewa_bin_exch, type = "mahalanobis")
  expect_length(r_bin, fit_geewa_bin_exch$clusters_no)
  expect_true(all(is.finite(r_bin)))
  expect_true(all(r_bin >= 0))
})

test_that("pearson residuals match the Poisson reference formula", {
  y <- fit_geewa_pois_exch$y
  mu <- fitted(fit_geewa_pois_exch)
  wt <- fit_geewa_pois_exch$prior.weights
  expected <- ((y - mu) / sqrt(mu) * sqrt(wt)) / sqrt(fit_geewa_pois_exch$phi)
  observed <- residuals(fit_geewa_pois_exch, type = "pearson")
  expect_equal(observed, expected, tolerance = 1e-6)
})

test_that("deviance residuals are scaled by the fitted dispersion", {
  y <- fit_geewa_pois_exch$y
  mu <- fitted(fit_geewa_pois_exch)
  wt <- fit_geewa_pois_exch$prior.weights
  sign_term <- sign(y - mu)
  expected <- sign_term * sqrt(
    pmax(fit_geewa_pois_exch$family$dev.resids(y, mu, wt), 0) / fit_geewa_pois_exch$phi
  )
  observed <- residuals(fit_geewa_pois_exch, type = "deviance")
  expect_equal(observed, expected, tolerance = 1e-10)
})

test_that("Mahalanobis residuals match their cluster-level definition", {
  id_labels <- as.character(fit_geewa_pois_exch$id)
  cluster_ids <- unique(id_labels)
  alpha <- as.numeric(fit_geewa_pois_exch$alpha[[1L]])
  expected <- vapply(
    cluster_ids,
    FUN.VALUE = numeric(1),
    FUN = function(cluster_id) {
      indices <- which(id_labels == cluster_id)
      mu <- fit_geewa_pois_exch$fitted.values[indices]
      weights <- fit_geewa_pois_exch$prior.weights[indices]
      n_i <- length(indices)
      correlation <- matrix(alpha, n_i, n_i)
      diag(correlation) <- 1
      marginal_sd <- sqrt(fit_geewa_pois_exch$phi * mu / weights)
      working_covariance <- correlation * tcrossprod(marginal_sd)
      raw_residual <- fit_geewa_pois_exch$y[indices] - mu
      as.numeric(
        crossprod(raw_residual, solve(working_covariance, raw_residual)) / n_i
      )
    }
  )
  names(expected) <- cluster_ids
  observed <- residuals(fit_geewa_pois_exch, type = "mahalanobis")
  expect_equal(observed, expected, tolerance = 1e-10)
})

test_that("deviance residuals do not depend on the residual degrees of freedom", {
  reference <- residuals(fit_geewa_pois_exch, type = "deviance")
  fit <- fit_geewa_pois_exch
  fit$df.residual <- 0
  expect_identical(residuals(fit, type = "deviance"), reference)
  expect_true(any(reference != 0))
})

test_that("deviance residual signs follow the raw residual", {
  fit <- fit_geewa_pois_exch
  dev <- residuals(fit, type = "deviance")
  raw <- residuals(fit, type = "working")
  expect_identical(sign(dev), sign(raw))
  expect_true(any(dev > 0))
  expect_true(any(dev < 0))
})

test_that("scaled residuals require a usable dispersion estimate", {
  fit <- fit_geewa_pois_exch
  fit$phi <- 0
  for (type in c("pearson", "deviance", "mahalanobis")) {
    expect_error(
      residuals(fit, type = type),
      "finite positive dispersion estimate"
    )
  }
  expect_identical(
    residuals(fit, type = "working"),
    residuals(fit_geewa_pois_exch, type = "working")
  )
})

test_that("Mahalanobis residuals reject a non-positive-definite covariance", {
  fit <- fit_geewa_pois_exch
  first_cluster <- which(fit$id == fit$id[[1L]])
  fit$fitted.values[first_cluster[[1L]]] <- 0
  expect_error(
    residuals(fit, type = "mahalanobis"),
    "not positive definite"
  )
  expect_error(
    residuals(fit, type = "mahalanobis"),
    as.character(fit$id[[1L]])
  )
})

test_that("residuals validates type", {
  expect_error(
    residuals(fit_geewa_pois_exch, type = "raw"),
    "should be one of"
  )
})
