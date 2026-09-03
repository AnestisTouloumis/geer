testthat::local_edition(3)

expect_jackknife_test_result <- function(x) {
  expect_type(x, "list")
  expect_true(all(c("test_stat", "test_df", "test_p") %in% names(x)))
  expect_true(is.numeric(x$test_stat))
  expect_true(is.numeric(x$test_df))
  expect_true(is.numeric(x$test_p))
  expect_length(x$test_stat, 1L)
  expect_length(x$test_df, 1L)
  expect_length(x$test_p, 1L)
  expect_true(is.finite(x$test_stat))
  expect_true(is.finite(x$test_df))
  expect_true(is.finite(x$test_p))
  expect_gte(x$test_stat, 0)
  expect_gt(x$test_df, 0)
  expect_gte(x$test_p, 0)
  expect_lte(x$test_p, 1)
}

jackknife_gaussian_data <- data.frame(
  id = rep(seq_len(8), each = 3),
  time = rep(seq_len(3), times = 8),
  group = rep(rep(c(0, 1), each = 4), each = 3),
  y = c(
    1.2, 1.7, 2.4,
    2.0, 2.2, 2.9,
    1.5, 2.1, 2.5,
    2.3, 2.7, 3.2,
    2.1, 2.8, 3.4,
    2.7, 3.0, 3.8,
    2.4, 2.6, 3.5,
    3.0, 3.5, 3.9
  )
)

jackknife_binary_data <- data.frame(
  id = rep(seq_len(10), each = 3),
  time = rep(seq_len(3), times = 10),
  group = rep(rep(c(0, 1), each = 5), each = 3),
  y = c(
    0, 0, 1,
    0, 1, 1,
    1, 0, 1,
    0, 1, 0,
    1, 1, 0,
    0, 1, 1,
    1, 1, 1,
    0, 0, 1,
    1, 0, 0,
    1, 1, 0
  )
)


test_that("jackknife covariance uses full leave-one-cluster refits with fixed correlation", {
  fit <- geewa(
    y ~ time + group,
    family = gaussian(),
    data = jackknife_gaussian_data,
    id = id,
    repeated = time,
    corstr = "exchangeable",
    method = "gee"
  )

  ids <- unique(jackknife_gaussian_data$id)
  fixed_alpha <- rep(as.numeric(fit$alpha), choose(max(jackknife_gaussian_data$time), 2L))
  delete_estimates <- t(vapply(
    ids,
    function(id_out) {
      fit_i <- geewa(
        y ~ time + group,
        family = gaussian(),
        data = jackknife_gaussian_data[jackknife_gaussian_data$id != id_out, ],
        id = id,
        repeated = time,
        corstr = "fixed",
        alpha_vector = fixed_alpha,
        method = "gee",
        beta_start = coef(fit)
      )
      coef(fit_i)
    },
    numeric(length(coef(fit)))
  ))
  dimnames(delete_estimates) <- list(as.character(ids), names(coef(fit)))
  centered <- sweep(delete_estimates, 2L, colMeans(delete_estimates), `-`)
  expected <- ((nrow(delete_estimates) - 1L) / nrow(delete_estimates)) *
    crossprod(centered)
  dimnames(expected) <- list(names(coef(fit)), names(coef(fit)))

  observed_estimates <- compute_jackknife_delete_estimates(fit)
  observed <- vcov(fit, cov_type = "jackknife")

  expect_equal(observed_estimates, delete_estimates, tolerance = 1e-7)
  expect_equal(observed, expected, tolerance = 1e-7)
  expect_equal(observed, t(observed), tolerance = 1e-12)
})


test_that("jackknife covariance fixes the full-data odds-ratio vector for geewa_binary", {
  fit <- geewa_binary(
    y ~ time + group,
    data = jackknife_binary_data,
    id = id,
    repeated = time,
    orstr = "exchangeable",
    method = "gee"
  )

  ids <- unique(jackknife_binary_data$id)
  delete_estimates <- t(vapply(
    ids,
    function(id_out) {
      fit_i <- geewa_binary(
        y ~ time + group,
        data = jackknife_binary_data[jackknife_binary_data$id != id_out, ],
        id = id,
        repeated = time,
        orstr = "fixed",
        alpha_vector = as.numeric(fit$alpha),
        method = "gee",
        beta_start = coef(fit)
      )
      coef(fit_i)
    },
    numeric(length(coef(fit)))
  ))
  dimnames(delete_estimates) <- list(as.character(ids), names(coef(fit)))
  centered <- sweep(delete_estimates, 2L, colMeans(delete_estimates), `-`)
  expected <- ((nrow(delete_estimates) - 1L) / nrow(delete_estimates)) *
    crossprod(centered)
  dimnames(expected) <- list(names(coef(fit)), names(coef(fit)))

  expect_equal(
    compute_jackknife_delete_estimates(fit),
    delete_estimates,
    tolerance = 1e-7
  )
  expect_equal(vcov(fit, cov_type = "jackknife"), expected, tolerance = 1e-7)
})


test_that("jackknife covariance works for an independence odds-ratio fit", {
  # Regression guard. fit_bingee_or() indexes alpha_vector by pair position, so
  # it needs length choose(max(repeated), 2), but an independence geewa_binary()
  # fit stores the scalar 1. Passing that scalar straight through read past the
  # end of the vector, which Armadillo does not bounds-check.
  fit <- geewa_binary(
    y ~ time + group,
    data = jackknife_binary_data,
    id = id,
    repeated = time,
    orstr = "independence",
    method = "gee"
  )
  expect_length(fit$alpha, 1L)

  ids <- unique(jackknife_binary_data$id)
  delete_estimates <- t(vapply(
    ids,
    function(id_out) {
      fit_i <- geewa_binary(
        y ~ time + group,
        data = jackknife_binary_data[jackknife_binary_data$id != id_out, ],
        id = id,
        repeated = time,
        orstr = "independence",
        method = "gee",
        beta_start = coef(fit)
      )
      coef(fit_i)
    },
    numeric(length(coef(fit)))
  ))
  dimnames(delete_estimates) <- list(as.character(ids), names(coef(fit)))
  centered <- sweep(delete_estimates, 2L, colMeans(delete_estimates), `-`)
  expected <- ((nrow(delete_estimates) - 1L) / nrow(delete_estimates)) *
    crossprod(centered)
  dimnames(expected) <- list(names(coef(fit)), names(coef(fit)))

  expect_equal(
    compute_jackknife_delete_estimates(fit),
    delete_estimates,
    tolerance = 1e-7
  )
  expect_equal(vcov(fit, cov_type = "jackknife"), expected, tolerance = 1e-7)
})


test_that("jackknife_or_alpha expands the independence odds-ratio vector", {
  fit <- geewa_binary(
    y ~ time + group,
    data = jackknife_binary_data,
    id = id,
    repeated = time,
    orstr = "independence",
    method = "gee"
  )
  alpha <- jackknife_or_alpha(fit, fit$repeated)
  expect_length(alpha, choose(max(fit$repeated), 2L))
  expect_true(all(alpha == 1))
})


test_that("jackknife_pair_subset rejects an alpha of the wrong length", {
  expect_error(
    jackknife_pair_subset(1, full_max = 4L, subset_max = 3L),
    "failed to map the fitted association parameters"
  )
})


test_that("jackknife covariance supports adjusted and penalized estimation methods", {
  methods <- c(
    "brgee-robust",
    "bcgee-robust",
    "pgee-jeffreys",
    "opgee-jeffreys",
    "hpgee-jeffreys"
  )

  for (method in methods) {
    fit <- geewa(
      y ~ time + group,
      family = gaussian(),
      data = jackknife_gaussian_data,
      id = id,
      repeated = time,
      corstr = "independence",
      method = method
    )
    v <- vcov(fit, cov_type = "jackknife")
    expect_true(is.matrix(v), info = method)
    expect_equal(dim(v), c(length(coef(fit)), length(coef(fit))), info = method)
    expect_true(all(is.finite(v)), info = method)
    expect_equal(v, t(v), tolerance = 1e-12, info = method)
  }
})


test_that("jackknife covariance applies the (K - 1) / K finite-sample factor", {
  fit <- geewa(
    y ~ time + group,
    family = gaussian(),
    data = jackknife_gaussian_data,
    id = id,
    repeated = time,
    corstr = "exchangeable"
  )
  delete_estimates <- compute_jackknife_delete_estimates(fit)
  k <- nrow(delete_estimates)
  expect_identical(k, fit$clusters_no)
  centered <- sweep(delete_estimates, 2L, colMeans(delete_estimates), `-`)
  unscaled <- crossprod(centered)
  dimnames(unscaled) <- list(names(coef(fit)), names(coef(fit)))
  expect_equal(
    vcov(fit, cov_type = "jackknife"),
    ((k - 1L) / k) * unscaled,
    tolerance = 1e-10
  )
})


test_that("jackknife association mapping preserves pair identities when maximum occasion drops", {
  alpha <- c(12, 13, 14, 23, 24, 34)
  expect_equal(
    jackknife_pair_subset(alpha, full_max = 4L, subset_max = 3L),
    c(12, 13, 23)
  )
})


test_that("summary accepts jackknife covariance", {
  fit <- geewa(
    y ~ time + group,
    family = gaussian(),
    data = jackknife_gaussian_data,
    id = id,
    repeated = time,
    corstr = "independence"
  )
  out <- summary(fit, cov_type = "jackknife")
  expect_s3_class(out, "summary.geer")
  expect_identical(out$cov_type, "jackknife")
})


test_that("all public cov_type interfaces advertise jackknife", {
  functions <- list(
    vcov.geer = vcov.geer,
    confint.geer = confint.geer,
    summary.geer = summary.geer,
    predict.geer = predict.geer,
    tidy.geer = tidy.geer,
    add1.geer = add1.geer,
    drop1.geer = drop1.geer,
    anova.geer = anova.geer,
    step_p = step_p,
    geecriteria = geecriteria,
    mcar_logistic_test = mcar_logistic_test
  )
  for (name in names(functions)) {
    choices <- eval(formals(functions[[name]])$cov_type)
    expect_true("jackknife" %in% choices, info = name)
  }
})


test_that("all hypothesis-test helpers accept jackknife covariance", {
  fit0 <- geewa(
    y ~ time,
    family = gaussian(),
    data = jackknife_gaussian_data,
    id = id,
    repeated = time,
    corstr = "independence",
    method = "gee"
  )
  fit1 <- geewa(
    y ~ time + group,
    family = gaussian(),
    data = jackknife_gaussian_data,
    id = id,
    repeated = time,
    corstr = "independence",
    method = "gee"
  )

  expect_jackknife_test_result(wald_test(fit0, fit1, cov_type = "jackknife"))
  expect_jackknife_test_result(score_test(fit0, fit1, cov_type = "jackknife"))
  expect_jackknife_test_result(
    working_wald_test(
      fit0, fit1, cov_type = "jackknife", pmethod = "rao-scott"
    )
  )
  expect_jackknife_test_result(
    working_score_test(
      fit0, fit1, cov_type = "jackknife", pmethod = "rao-scott"
    )
  )
  expect_jackknife_test_result(
    working_lrt_test(
      fit0, fit1, cov_type = "jackknife", pmethod = "rao-scott"
    )
  )
})


test_that("criteria and model-comparison interfaces accept jackknife", {
  fit <- geewa(
    y ~ time + group,
    family = gaussian(),
    data = jackknife_gaussian_data,
    id = id,
    repeated = time,
    corstr = "independence",
    method = "gee"
  )

  expect_s3_class(
    anova(fit, test = "wald", cov_type = "jackknife"),
    "anova"
  )
  expect_true(is.data.frame(geecriteria(fit, cov_type = "jackknife")))
  expect_s3_class(
    drop1(fit, scope = "group", test = "wald", cov_type = "jackknife"),
    "anova"
  )
})
