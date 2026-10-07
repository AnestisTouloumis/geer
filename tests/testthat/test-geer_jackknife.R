testthat::local_edition(3)

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
  fixed_alpha <- rep(as.numeric(fit$alpha), choose(max(jackknife_gaussian_data$time), 2L))
  delete_estimates <- jackknife_delete_estimates_by_refit(
    fit,
    jackknife_gaussian_data,
    function(d) {
      geewa(
        y ~ time + group,
        family = gaussian(),
        data = d,
        id = id,
        repeated = time,
        corstr = "fixed",
        alpha_vector = fixed_alpha,
        method = "gee",
        beta_start = coef(fit)
      )
    }
  )
  expect_jackknife_vcov(fit, delete_estimates)
})


test_that("jackknife covariance refits geewa_binary with the correct odds-ratio vector", {
  # The independence case is a regression guard. fit_geesolver_or() indexes
  # alpha_vector by pair position, so it needs length choose(max(repeated), 2),
  # but an independence geewa_binary() fit stores the scalar 1. Passing that
  # scalar straight through read past the end of the vector, which Armadillo
  # does not bounds-check.
  for (orstr in c("exchangeable", "independence")) {
    fit <- geewa_binary(
      y ~ time + group,
      data = jackknife_binary_data,
      id = id,
      repeated = time,
      orstr = orstr,
      method = "gee"
    )
    if (orstr == "independence") {
      expect_length(fit$alpha, 1L)
    }
    delete_estimates <- jackknife_delete_estimates_by_refit(
      fit,
      jackknife_binary_data,
      function(d) {
        if (orstr == "independence") {
          geewa_binary(
            y ~ time + group,
            data = d,
            id = id,
            repeated = time,
            orstr = "independence",
            method = "gee",
            beta_start = coef(fit)
          )
        } else {
          geewa_binary(
            y ~ time + group,
            data = d,
            id = id,
            repeated = time,
            orstr = "fixed",
            alpha_vector = as.numeric(fit$alpha),
            method = "gee",
            beta_start = coef(fit)
          )
        }
      }
    )
    expect_jackknife_vcov(fit, delete_estimates)
  }
})


test_that("compute_jackknife_alpha_or expands the independence odds-ratio vector", {
  fit <- geewa_binary(
    y ~ time + group,
    data = jackknife_binary_data,
    id = id,
    repeated = time,
    orstr = "independence",
    method = "gee"
  )
  alpha <- compute_jackknife_alpha_or(fit, fit$repeated)
  expect_length(alpha, choose(max(fit$repeated), 2L))
  expect_true(all(alpha == 1))
})


test_that("select_jackknife_pair_subset rejects an alpha of the wrong length", {
  expect_error(
    select_jackknife_pair_subset(1, full_max = 4L, subset_max = 3L),
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
    select_jackknife_pair_subset(alpha, full_max = 4L, subset_max = 3L),
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

  expect_test_result(compute_wald_test(fit0, fit1, cov_type = "jackknife"))
  expect_test_result(compute_score_test(fit0, fit1, cov_type = "jackknife"))
  expect_test_result(
    compute_working_wald_test(
      fit0, fit1, cov_type = "jackknife", pmethod = "rao-scott"
    )
  )
  expect_test_result(
    compute_working_score_test(
      fit0, fit1, cov_type = "jackknife", pmethod = "rao-scott"
    )
  )
  expect_test_result(
    compute_working_lrt_test(
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


test_that("jackknife convergence checks report the solver failure reason", {
  bad <- list(criterion = c(1, Inf), beta_mat = matrix(0, 2, 3),
              failure = "singular matrix")
  expect_error(
    geer:::check_jackknife_convergence(bad, 1e-6, "7"),
    "leave-one-cluster fit did not converge (singular matrix)",
    fixed = TRUE
  )
  expect_error(
    geer:::check_jackknife_solver_failure(bad, "7"),
    "jackknife covariance failed for cluster '7': singular matrix",
    fixed = TRUE
  )
  ok <- list(criterion = c(1, 0), beta_mat = matrix(0, 2, 3))
  expect_silent(geer:::check_jackknife_solver_failure(ok, "7"))
})
