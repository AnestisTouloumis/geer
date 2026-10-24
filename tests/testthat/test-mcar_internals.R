mcar_internals_matrix <- function() {
  cbind(
    y1 = c(1.0, 2.1, 2.9, 4.2, 5.1, 5.8, 7.2, 8.0, 9.1, 9.9),
    y2 = c(2.2, 1.7, 3.9, 3.1, 6.3, 5.2, 7.9, 6.8, 9.6, 8.7),
    y3 = c(0.5, 1.9, 1.1, 3.4, 2.8, 4.1, 3.6, 5.9, 5.2, 6.4)
  )
}

mcar_fit <- function(data = test_data$epilepsy, repeated = TRUE, ...) {
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

test_that("initial parameters use available-case means and a valid covariance", {
  x <- mcar_internals_matrix()
  x[c(2, 5), 2] <- NA

  out <- compute_mcar_em_start(x, "mcar_little_test")

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

  out <- compute_mcar_em_start(x, "mcar_little_test")

  expect_true(all(is.finite(out$sigma)))
  expect_no_error(chol(out$sigma))
})

test_that("initial parameters fail when a variable has no observed values", {
  x <- cbind(a = c(1, 2, 3, 4), b = NA_real_)

  expect_error(
    compute_mcar_em_start(x, "mcar_little_test"),
    paste0(
      "mcar_little_test failed: the initial covariance matrix for ",
      "normal-theory EM estimation is singular"
    ),
    fixed = TRUE
  )
})

test_that("EM reproduces the complete-data maximum-likelihood estimates", {
  x <- mcar_internals_matrix()

  out <- fit_mcar_normal_em(
    x,
    maxiter = 100L,
    tolerance = 1e-10,
    context = "mcar_little_test"
  )

  expect_true(out$converged)
  expect_equal(out$mu, colMeans(x), tolerance = 1e-8)
  expect_equal(
    out$sigma,
    cov(x) * (nrow(x) - 1) / nrow(x),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
})

test_that("EM ignores rows with no observed values", {
  x <- mcar_internals_matrix()
  x[c(3, 7), 2] <- NA
  x[9, 3] <- NA
  with_empty_row <- rbind(x, c(NA_real_, NA_real_, NA_real_))

  reference <- fit_mcar_normal_em(
    x,
    maxiter = 5000L,
    tolerance = 1e-12,
    context = "mcar_little_test"
  )
  out <- fit_mcar_normal_em(
    with_empty_row,
    maxiter = 5000L,
    tolerance = 1e-12,
    context = "mcar_little_test"
  )

  expect_true(out$converged)
  expect_equal(out$mu, reference$mu, tolerance = 1e-6)
  expect_equal(out$sigma, reference$sigma, tolerance = 1e-6)
})

test_that("EM stops when the iteration limit is reached", {
  x <- mcar_internals_matrix()
  x[c(3, 7), 2] <- NA
  x[9, 3] <- NA

  expect_error(
    fit_mcar_normal_em(
      x,
      maxiter = 1L,
      tolerance = 1e-14,
      context = "mcar_homoscedasticity_test"
    ),
    paste0(
      "mcar_homoscedasticity_test failed: the EM algorithm did not converge ",
      "within 1 iterations; increase 'maxiter'"
    ),
    fixed = TRUE
  )
})

test_that("EM rejects a singular maximum-likelihood covariance", {
  x <- mcar_internals_matrix()
  x[, 3] <- x[, 1] + x[, 2]

  expect_error(
    fit_mcar_normal_em(
      x,
      maxiter = 100L,
      tolerance = 1e-10,
      context = "mcar_little_test"
    ),
    "singular or nearly singular"
  )
})

test_that("extract_mcar_matrix returns one row per cluster", {
  fit <- mcar_fit()

  out <- extract_mcar_matrix(fit)

  expect_identical(dim(out), c(length(unique(test_data$epilepsy$id)), 4L))
  expect_identical(colnames(out), as.character(1:4))
  expect_identical(
    unname(out[cbind(
      match(test_data$epilepsy$id, rownames(out)),
      test_data$epilepsy$visit
    )]),
    as.numeric(test_data$epilepsy$seizures)
  )
})

test_that("extract_mcar_matrix numbers occasions when repeated is absent", {
  fit <- mcar_fit(repeated = FALSE)

  out <- extract_mcar_matrix(fit)

  expect_identical(colnames(out), as.character(1:4))
  expect_false(anyNA(out))
})

test_that("extract_mcar_matrix marks absent occasions as missing", {
  fit <- mcar_fit()
  incomplete <- test_data$epilepsy[-c(2, 3, 10), ]

  out <- extract_mcar_matrix(fit, data = incomplete)

  expect_identical(sum(is.na(out)), 3L)
  expect_true(is.na(out[1L, 2L]))
  expect_true(is.na(out[1L, 3L]))
})

test_that("extract_mcar_matrix works when the formula has no environment", {
  fit <- mcar_fit()
  reference <- extract_mcar_matrix(fit)
  attr(fit$terms, ".Environment") <- NULL

  expect_identical(extract_mcar_matrix(fit), reference)
})

test_that("extract_mcar_matrix needs the original data", {
  fit <- mcar_fit()
  fit$data <- NULL

  expect_error(
    extract_mcar_matrix(fit),
    "the original data are not available in the fitted object"
  )
  expect_error(
    extract_mcar_matrix(fit, data = data.frame(z = 1))
  )
})

test_that("extract_mcar_matrix needs a univariate numeric response", {
  data <- test_data$epilepsy
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
    extract_mcar_matrix(matrix_fit),
    "require a univariate response"
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
    extract_mcar_matrix(factor_fit),
    "the MCAR tests require a numeric response"
  )
})

test_that("EM agrees with the monotone normal MLE in different units", {
  x <- test_data$mcar_monotone
  expected <- compute_reference_monotone(x)
  for (scales in list(c(1, 1), c(1e-5, 1e-5), c(1e-4, 1e4), c(-2, 0.5))) {
    centers_list <- if (all(scales == 1)) {
      list(c(0, 0), c(1e8, -1e8))
    } else {
      list(c(0, 0))
    }
    for (centers in centers_list) {
      transformed <- sweep(sweep(x, 2L, scales, `*`), 2L, centers, `+`)
      out <- fit_mcar_normal_em(transformed, 5000L, 1e-12,
                                "mcar_little_test")
      recovered <- (out$mu - centers) / scales
      expected_mu <- expected$mu
      if (any(centers != 0)) {
        reference <- compute_reference_monotone(
          sweep(sweep(transformed, 2L, centers, `-`), 2L, scales, `/`)
        )
        expected_mu <- reference$mu
      } else {
        reference <- expected
      }
      expect_equal(out$sigma / tcrossprod(scales), reference$sigma,
                   tolerance = 1e-6, ignore_attr = TRUE)
      expect_equal(recovered, expected_mu, tolerance = 1e-6,
                   ignore_attr = TRUE)
    }
  }
})

test_that("MCAR reconstruction errors identify the calling diagnostic", {
  fit <- fit_geewa_pois_exch
  fit$call$subset <- quote(stop("subset evaluation failed", call. = FALSE))
  error <- tryCatch(mcar_little_test(fit), error = identity)
  expect_s3_class(error, "error")
  expect_null(conditionCall(error))
  expect_match(conditionMessage(error), "^mcar_little_test failed:")
})
