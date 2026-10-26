test_that("build_percent_labels returns percent-suffixed labels", {
  out <- build_percent_labels(c(0.025, 0.975), digits = 3)
  expect_type(out, "character")
  expect_length(out, 2L)
  expect_true(all(grepl("%$", out)))
})

test_that("checked probability labels reject values outside [0, 1]", {
  expect_error(
    build_percent_labels(c(-0.1, 0.5, 1.1), check = TRUE),
    "'probs' must lie in [0, 1] when 'check = TRUE'", fixed = TRUE
  )
})

test_that("build_percent_labels validates its inputs", {
  expect_error(build_percent_labels("0.5"), "'probs' must be numeric")
  expect_error(build_percent_labels(0.5, digits = 0), "'digits' must be a positive integer")
  expect_error(build_percent_labels(0.5, digits = 1.5), "'digits' must be a positive integer")
  expect_error(
    build_percent_labels(0.5, check = NA),
    "'check' must be TRUE or FALSE"
  )
  expect_error(
    build_percent_labels(0.5, check = c(TRUE, FALSE)),
    "'check' must be TRUE or FALSE"
  )
})

test_that("build_test_label returns the expected display labels", {
  expect_identical(build_test_label("wald"), "Wald")
  expect_identical(build_test_label("score"), "Score")
  expect_identical(build_test_label("working-wald"), "Modified Working Wald")
  expect_identical(build_test_label("working-score"), "Modified Working Score")
  expect_identical(build_test_label("working-lrt"), "Modified Working LRT")
})

test_that("summary, tidy, and confint stop on a negative variance", {
  fit <- fit_geewa_pois_exch
  fit$bias_corrected_covariance[2L, 2L] <- -1
  message_text <- paste0(
    " failed: negative variance for treatmentprogabide; the covariance ",
    "estimator is not positive semi-definite, so no standard error is ",
    "available"
  )
  expect_error(summary(fit), paste0("summary.geer", message_text), fixed = TRUE)
  expect_error(tidy(fit), paste0("tidy.geer", message_text), fixed = TRUE)
  expect_error(confint(fit), paste0("confint.geer", message_text), fixed = TRUE)
})

test_that("predict stops when a prediction variance is negative", {
  fit <- fit_geewa_pois_exch
  p <- ncol(fit$x)
  fit$bias_corrected_covariance[] <- -diag(p)
  expect_error(
    predict(fit, se.fit = TRUE),
    "predict.geer failed: negative variance for prediction 1, ",
    fixed = TRUE
  )
  fit$bias_corrected_covariance[] <- NA_real_
  expect_error(
    predict(fit, newdata = test_data$epilepsy[1:2, ], se.fit = TRUE),
    "predict.geer failed: non-finite variance for prediction 1, prediction 2",
    fixed = TRUE
  )
})

test_that("compute_se_or_stop leaves valid variances unchanged", {
  expect_silent(se <- compute_se_or_stop(c(a = 4, b = 0, c = 9), "vcov.geer"))
  expect_equal(se, c(a = 2, b = 0, c = 3))
})

test_that("jackknife covariance is cached on the fit and invalidated by coefficient changes", {
  fit <- fit_geewa_pois_indep
  expect_type(fit$cache, "environment")
  expect_null(fit$cache$jackknife)
  v1 <- vcov(fit, cov_type = "jackknife")
  expect_false(is.null(fit$cache$jackknife))
  fit$cache$jackknife$value[1L, 1L] <- 12345
  expect_equal(unname(vcov(fit, cov_type = "jackknife")[1L, 1L]), 12345)
  fit2 <- fit
  fit2$coefficients <- fit2$coefficients + 1e-3
  v2 <- vcov(fit2, cov_type = "jackknife")
  expect_false(isTRUE(all.equal(unname(v2[1L, 1L]), 12345)))
  expect_equal(dim(v1), dim(v2))
})

test_that("small-sample covariances are cached on the fit and invalidated by coefficient changes", {
  fit <- fit_geewa_pois_indep
  expect_null(fit$cache$small_sample)
  v1 <- vcov(fit, cov_type = "mancl-derouen")
  expect_false(is.null(fit$cache$small_sample))
  fit$cache$small_sample$value[["mancl-derouen"]][1L, 1L] <- 12345
  expect_equal(unname(vcov(fit, cov_type = "mancl-derouen")[1L, 1L]), 12345)
  fit2 <- fit
  fit2$coefficients <- fit2$coefficients + 1e-3
  v2 <- vcov(fit2, cov_type = "mancl-derouen")
  expect_false(isTRUE(all.equal(unname(v2[1L, 1L]), 12345)))
  expect_equal(dim(v1), dim(v2))
})

test_that("is_balanced_design requires every cluster to observe each position once", {
  expect_true(is_balanced_design(fit_geewa_pois_exch))
  fit <- fit_geewa_pois_exch
  fit$repeated[2L] <- 1L
  expect_false(is_balanced_design(fit))
  fit <- fit_geewa_pois_exch
  keep <- -2L
  fit$id <- fit$id[keep]
  fit$repeated <- fit$repeated[keep]
  expect_false(is_balanced_design(fit))
})

test_that("refit_geer stops when the refit did not converge", {
  fit <- geewa(
    seizures ~ treatment,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log")
  )
  fit_ok <- refit_geer(fit, . ~ . + lnbaseline, "add1.geer")
  expect_s3_class(fit_ok, "geer")
  fit_limited <- suppressWarnings(
    geewa(
      seizures ~ treatment + visit,
      data = test_data$epilepsy,
      id = id,
      family = poisson("log"),
      corstr = "exchangeable",
      beta_start = rep(0, 3),
      control = geer_control(maxiter = 1L, tolerance = 1e-12)
    )
  )
  expect_false(fit_limited$converged)
  expect_error(
    suppressWarnings(refit_geer(fit_limited, . ~ . + lnbaseline, "add1.geer")),
    paste0(
      "add1.geer failed: the fit of 'seizures ~ treatment + visit + ",
      "lnbaseline' did not converge"
    ),
    fixed = TRUE
  )
})

test_that("refit_geer ignores a fixed-length beta_start of the original fit", {
  fit <- geewa(
    seizures ~ treatment,
    data = test_data$epilepsy,
    id = id,
    family = poisson("log"),
    beta_start = c(0, 0)
  )
  refit <- refit_geer(fit, . ~ . + lnbaseline, "add1.geer")
  expect_length(coef(refit), 3L)
})

test_that("data variable checks respect data and formula environments", {
  data <- test_data$epilepsy
  expect_null(check_data_variables(
    quote(seizures + treatment), data, "data", "in the model formula"
  ))
  expect_error(
    check_data_variables(quote(absent_variable), data, "data", "in a formula"),
    "'data' must contain the variable(s) 'absent_variable'", fixed = TRUE
  )
  env <- new.env(parent = emptyenv())
  env$constant <- 2
  expect_null(check_data_variables(
    quote(seizures + constant), data, "data", "in a formula", enclos = env
  ))
  expect_null(check_data_variables(
    quote(constant), env, "data", "in a formula"
  ))
})
