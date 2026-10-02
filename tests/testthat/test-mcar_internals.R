testthat::local_edition(3)


mcar_internals_matrix <- function() {
  cbind(
    y1 = c(1.0, 2.1, 2.9, 4.2, 5.1, 5.8, 7.2, 8.0, 9.1, 9.9),
    y2 = c(2.2, 1.7, 3.9, 3.1, 6.3, 5.2, 7.9, 6.8, 9.6, 8.7),
    y3 = c(0.5, 1.9, 1.1, 3.4, 2.8, 4.1, 3.6, 5.9, 5.2, 6.4)
  )
}


## ----------------------------------------------- compute_mcar_normal_initial_parameters ----

test_that("initial parameters use available-case means and a valid covariance", {
  x <- mcar_internals_matrix()
  x[c(2, 5), 2] <- NA

  out <- geer:::compute_mcar_normal_initial_parameters(x)

  expect_equal(out$mu, colMeans(x, na.rm = TRUE))
  expect_true(isSymmetric(out$sigma))
  expect_gt(min(eigen(out$sigma, symmetric = TRUE, only.values = TRUE)$values), 0)
  expect_equal(
    diag(out$sigma),
    apply(x, 2L, function(v) mean((v - mean(v, na.rm = TRUE))^2, na.rm = TRUE)),
    ignore_attr = TRUE
  )
})


test_that("initial parameters add a ridge to a singular covariance", {
  x <- cbind(a = as.numeric(1:6), b = 2 * as.numeric(1:6))

  out <- geer:::compute_mcar_normal_initial_parameters(x)

  expect_true(all(is.finite(out$sigma)))
  expect_no_error(chol(out$sigma))
})


test_that("initial parameters fail when a variable has no observed values", {
  x <- cbind(a = c(1, 2, 3, 4), b = NA_real_)

  expect_error(
    geer:::compute_mcar_normal_initial_parameters(x),
    "initial covariance matrix for normal-theory EM estimation is singular"
  )
})


## ------------------------------------------------------------ fit_mcar_normal_em ----

test_that("EM reproduces the complete-data maximum-likelihood estimates", {
  x <- mcar_internals_matrix()

  out <- geer:::fit_mcar_normal_em(x, maxit = 100L, tol = 1e-10)

  expect_true(out$converged)
  expect_equal(out$mu, colMeans(x), tolerance = 1e-8)
  expect_equal(
    out$sigma,
    stats::cov(x) * (nrow(x) - 1) / nrow(x),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})


test_that("EM ignores rows with no observed values", {
  x <- mcar_internals_matrix()
  x[c(3, 7), 2] <- NA
  x[9, 3] <- NA
  with_empty_row <- rbind(x, c(NA_real_, NA_real_, NA_real_))

  reference <- geer:::fit_mcar_normal_em(x, maxit = 5000L, tol = 1e-12)
  out <- geer:::fit_mcar_normal_em(with_empty_row, maxit = 5000L, tol = 1e-12)

  expect_true(out$converged)
  expect_equal(out$mu, reference$mu, tolerance = 1e-6)
  expect_equal(out$sigma, reference$sigma, tolerance = 1e-6)
})


test_that("EM warns when the iteration limit is reached", {
  x <- mcar_internals_matrix()
  x[c(3, 7), 2] <- NA
  x[9, 3] <- NA

  expect_warning(
    out <- geer:::fit_mcar_normal_em(x, maxit = 1L, tol = 1e-14),
    "did not converge within 1 iterations"
  )
  expect_false(out$converged)
  expect_identical(out$iterations, 1L)
})


test_that("EM rejects a singular maximum-likelihood covariance", {
  x <- mcar_internals_matrix()
  x[, 3] <- x[, 1] + x[, 2]

  expect_error(
    geer:::fit_mcar_normal_em(x, maxit = 100L, tol = 1e-10),
    "singular or nearly singular"
  )
})


## ----------------------------------------------- extract_geer_mcar_matrix ----

mcar_fit <- function(data = epilepsy, repeated = TRUE, ...) {
  if (repeated) {
    geewa(
      formula = seizures ~ treatment + lnbaseline + lnage,
      data = data,
      id = id,
      repeated = visit,
      family = poisson(link = "log"),
      corstr = "independence",
      method = "gee",
      ...
    )
  } else {
    geewa(
      formula = seizures ~ treatment + lnbaseline + lnage,
      data = data,
      id = id,
      family = poisson(link = "log"),
      corstr = "independence",
      method = "gee",
      ...
    )
  }
}


test_that("extract_geer_mcar_matrix returns one row per cluster", {
  fit <- mcar_fit()

  out <- geer:::extract_geer_mcar_matrix(fit)

  expect_identical(dim(out), c(length(unique(epilepsy$id)), 4L))
  expect_identical(colnames(out), as.character(1:4))
  expect_identical(
    unname(out[cbind(match(epilepsy$id, rownames(out)), epilepsy$visit)]),
    as.numeric(epilepsy$seizures)
  )
})


test_that("extract_geer_mcar_matrix numbers occasions when repeated is absent", {
  fit <- mcar_fit(repeated = FALSE)

  out <- geer:::extract_geer_mcar_matrix(fit)

  expect_identical(colnames(out), as.character(1:4))
  expect_false(anyNA(out))
})


test_that("extract_geer_mcar_matrix marks absent occasions as missing", {
  fit <- mcar_fit()
  incomplete <- epilepsy[-c(2, 3, 10), ]

  out <- geer:::extract_geer_mcar_matrix(fit, data = incomplete)

  expect_identical(sum(is.na(out)), 3L)
  expect_true(is.na(out[1L, 2L]))
  expect_true(is.na(out[1L, 3L]))
})


test_that("extract_geer_mcar_matrix works when the formula has no environment", {
  fit <- mcar_fit()
  reference <- geer:::extract_geer_mcar_matrix(fit)
  attr(fit$terms, ".Environment") <- NULL

  expect_identical(geer:::extract_geer_mcar_matrix(fit), reference)
})


test_that("extract_geer_mcar_matrix needs the original data", {
  fit <- mcar_fit()
  fit$data <- NULL

  expect_error(
    geer:::extract_geer_mcar_matrix(fit),
    "the original data are not available in the fitted object"
  )
  expect_error(
    geer:::extract_geer_mcar_matrix(fit, data = data.frame(z = 1)),
    "could not reconstruct the original data for the MCAR diagnostic"
  )
})


test_that("extract_geer_mcar_matrix needs a univariate numeric response", {
  data <- epilepsy
  data$succ <- pmin(data$seizures, 20)
  data$fail <- 20 - data$succ
  data$high <- factor(data$seizures > 3)

  matrix_fit <- geewa(
    formula = cbind(succ, fail) ~ treatment,
    data = data,
    id = id,
    family = binomial(link = "logit"),
    corstr = "independence",
    method = "gee"
  )
  expect_error(
    geer:::extract_geer_mcar_matrix(matrix_fit),
    "require a univariate response at each repeated measurement"
  )

  factor_fit <- geewa(
    formula = high ~ treatment,
    data = data,
    id = id,
    family = binomial(link = "logit"),
    corstr = "independence",
    method = "gee"
  )
  expect_error(
    geer:::extract_geer_mcar_matrix(factor_fit),
    "the response used for the MCAR diagnostic must be numeric"
  )
})


test_that("extract_geer_mcar_matrix validates id and repeated", {
  fit <- mcar_fit()

  no_id <- fit
  no_id$call$id <- NULL
  expect_error(
    geer:::extract_geer_mcar_matrix(no_id),
    "'id' could not be recovered from the fitted model"
  )

  na_id <- epilepsy
  na_id$id[3] <- NA
  expect_error(
    geer:::extract_geer_mcar_matrix(fit, data = na_id),
    "'id' cannot contain missing values for the MCAR diagnostic"
  )

  na_repeated <- epilepsy
  na_repeated$visit[3] <- NA
  expect_error(
    geer:::extract_geer_mcar_matrix(fit, data = na_repeated),
    "'repeated' cannot contain missing values for the MCAR diagnostic"
  )

  duplicated_repeated <- epilepsy
  duplicated_repeated$visit <- 1L
  expect_error(
    geer:::extract_geer_mcar_matrix(fit, data = duplicated_repeated),
    "'repeated' must identify unique measurements within each 'id'"
  )
})
