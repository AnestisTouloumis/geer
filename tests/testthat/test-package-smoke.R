testthat::local_edition(3)

count_fit <- fit_geewa_pois_exch
binary_fit <- fit_geewa_bin_exch


test_that("package exports the core public API", {
  exported <- getNamespaceExports(asNamespace("geer"))

  expect_true(all(c(
    "geecriteria",
    "geer_control",
    "geewa",
    "geewa_binary",
    "glance",
    "step_p",
    "tidy"
  ) %in% exported))
})


test_that("package datasets load and contain essential columns", {
  expected_columns <- list(
    cerebrovascular = c("id", "period", "ecg", "treatment"),
    epilepsy = c("id", "visit", "seizures", "treatment", "lnbaseline", "lnage"),
    leprosy = c("id", "period", "bacilli", "treatment"),
    respiratory = c("id", "visit", "status", "treatment"),
    rinse = c("id", "time", "score", "treatment", "baseline")
  )
  for (nm in names(expected_columns)) {
    dat <- test_data[[nm]]
    expect_s3_class(dat, "data.frame")
    expect_true(all(expected_columns[[nm]] %in% names(dat)))
    expect_gt(nrow(dat), 0L)
  }
})


test_that("count-response workflow returns a coherent fitted object", {
  fit <- count_fit
  p <- length(coef(fit))
  expect_s3_class(fit, "geer")
  expect_true(isTRUE(fit$converged))
  expect_identical(fit$association_structure, "exchangeable")
  expect_identical(fit$family$family, "poisson")
  expect_identical(fit$family$link, "log")
  expect_gt(p, 0L)
  expect_equal(length(fitted(fit)), fit$obs_no)
  expect_equal(length(residuals(fit)), fit$obs_no)
  expect_equal(length(predict(fit, type = "link")), fit$obs_no)
  expect_equal(length(predict(fit, type = "response")), fit$obs_no)
  expect_equal(dim(vcov(fit, cov_type = "robust")), c(p, p))
  crit <- geecriteria(fit)
  expect_true(all(c("QIC", "QICu", "CIC") %in% names(crit)))
  expect_true(all(is.finite(unlist(crit[c("QIC", "QICu", "CIC")]))))
})


test_that("binary-response workflow returns a coherent fitted object", {
  fit <- binary_fit
  predicted <- predict(fit, type = "response")
  expect_s3_class(fit, "geer")
  expect_true(isTRUE(fit$converged))
  expect_identical(fit$association_structure, "exchangeable")
  expect_identical(fit$family$family, "binomial")
  expect_identical(fit$family$link, "logit")
  expect_equal(fit$phi, 1, tolerance = 1e-10)
  expect_equal(length(predicted), fit$obs_no)
  expect_true(all(is.finite(predicted)))
  expect_true(all(predicted >= 0 & predicted <= 1))
  expect_equal(nrow(tidy(fit)), length(coef(fit)))
  expect_equal(glance(fit)$nobs, fit$obs_no)
  expect_equal(glance(fit)$nclusters, fit$clusters_no)
})
