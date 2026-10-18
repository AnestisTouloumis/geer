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

test_that("jackknife covariance equals the aggregation of independently refitted deletion estimates", {
  ctl_tight <- geer_control(tolerance = 1e-12, maxiter = 500)
  d <- ref_sim("cc", ref_times_unbalanced, binomial(), "logit", "ar1", 0.5,
    phi = 1.2,
    beta = c(0.3, 0.2, -0.1), seed = 71
  )
  fit <- ref_fit(d, "gee", ctl_tight, beta_start = d$beta)
  ids <- unique(d$id)
  B <- sapply(ids, function(i) {
    e <- d
    keep <- d$id != i
    for (nm in c("y", "id", "tt", "w")) e[[nm]] <- d[[nm]][keep]
    e$X <- d$X[keep, , drop = FALSE]
    b <- unname(coef(fit))
    for (it in 1:100) {
      st <- solve(ref_Sigma0(e, b), ref_U(e, b))
      b <- b + st
      if (sqrt(sum(st^2)) < 1e-13) break
    }
    b
  })
  N <- length(ids)
  Bm <- rowMeans(B)
  S <- (N - 1) / N * Reduce(`+`, lapply(seq_len(N), function(i) tcrossprod(B[, i] - Bm)))
  expect_equal(unname(vcov(fit, cov_type = "jackknife")), S, tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("jackknife covariance refits geewa_binary with the correct odds-ratio vector", {
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

test_that("jackknife covariance holds the odds ratios fixed for one-step methods", {
  for (method in c("bcgee-robust", "bcgee-naive", "opgee-jeffreys", "hpgee-jeffreys")) {
    fit <- geewa_binary(
      y ~ time + group,
      data = jackknife_binary_data,
      id = id,
      repeated = time,
      orstr = "exchangeable",
      method = method
    )
    delete_estimates <- jackknife_delete_estimates_by_refit(
      fit,
      jackknife_binary_data,
      function(d) {
        geewa_binary(
          y ~ time + group,
          data = d,
          id = id,
          repeated = time,
          orstr = "fixed",
          alpha_vector = as.numeric(fit$alpha),
          method = method,
          beta_start = coef(fit)
        )
      }
    )
    expect_jackknife_vcov(fit, delete_estimates)
  }
})

test_that("jackknife covariance supports toeplitz and m-dependent correlations", {
  dat <- withr::with_seed(4, {
    out <- data.frame(id = rep(1:20, each = 3), time = rep(1:3, 20), x = rnorm(60))
    out$y <- 1 + 0.5 * out$x + rep(rnorm(20, sd = 0.7), each = 3) + rnorm(60)
    out
  })
  for (corstr in c("toeplitz", "m-dependent")) {
    fit <- geewa(y ~ x, data = dat, id = id, repeated = time, corstr = corstr, Mv = 1)
    alpha <- as.numeric(fit$alpha)
    correlation <- toeplitz(c(1, alpha, rep(0, 2 - length(alpha))))
    delete_estimates <- jackknife_delete_estimates_by_refit(
      fit,
      dat,
      function(d) {
        geewa(y ~ x,
          data = d, id = id, repeated = time, corstr = "fixed",
          alpha_vector = correlation[lower.tri(correlation)], beta_start = coef(fit)
        )
      }
    )
    expect_jackknife_vcov(fit, delete_estimates)
  }
})

test_that("jackknife covariance needs two clusters and a known fitting function", {
  one_cluster <- fit_geewa_pois_indep
  one_cluster$id[] <- one_cluster$id[1L]
  expect_error(
    compute_jackknife_delete_estimates(one_cluster),
    "jackknife covariance requires at least two clusters"
  )
  unknown <- fit_geewa_pois_indep
  unknown$fit_function <- NULL
  unknown$call[[1L]] <- quote(some_wrapper)
  expect_error(
    compute_jackknife_delete_estimates(unknown),
    "could not determine whether the model was fitted by geewa\\(\\) or geewa_binary\\(\\)"
  )
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
    "brgee-naive",
    "brgee-empirical",
    "bcgee-robust",
    "bcgee-naive",
    "bcgee-empirical",
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

test_that("jackknife association mapping preserves pair identities when maximum occasion drops", {
  alpha <- c(12, 13, 14, 23, 24, 34)
  expect_equal(
    select_jackknife_pair_subset(alpha, full_max = 4L, subset_max = 3L),
    c(12, 13, 23)
  )
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
      fit0, fit1,
      cov_type = "jackknife", pmethod = "rao-scott"
    )
  )
  expect_test_result(
    compute_working_score_test(
      fit0, fit1,
      cov_type = "jackknife", pmethod = "rao-scott"
    )
  )
  expect_test_result(
    compute_working_lrt_test(
      fit0, fit1,
      cov_type = "jackknife", pmethod = "rao-scott"
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
  expect_s3_class(geecriteria(fit, cov_type = "jackknife"), "data.frame")
  expect_s3_class(
    drop1(fit, scope = "group", test = "wald", cov_type = "jackknife"),
    "anova"
  )
})

test_that("confint, predict and tidy use the jackknife covariance", {
  fit <- geewa(
    y ~ time + group,
    family = gaussian(),
    data = jackknife_gaussian_data,
    id = id,
    repeated = time,
    corstr = "exchangeable"
  )
  v <- vcov(fit, cov_type = "jackknife")
  se <- sqrt(diag(v))
  ci <- confint(fit, cov_type = "jackknife")
  expect_equal(unname(ci[, 2]), unname(coef(fit) + qnorm(0.975) * se), tolerance = 1e-10)
  expect_equal(tidy(fit, cov_type = "jackknife")$std.error, unname(se), tolerance = 1e-10)
  pred <- predict(fit, se.fit = TRUE, cov_type = "jackknife")
  x <- model.matrix(fit)
  expect_equal(pred$se.fit, unname(sqrt(rowSums((x %*% v) * x))), tolerance = 1e-10)
})

test_that("jackknife convergence checks report the solver failure reason", {
  bad <- list(
    criterion = c(1, Inf), beta_mat = matrix(0, 2, 3),
    failure = "singular matrix"
  )
  expect_error(
    check_jackknife_convergence(bad, 1e-6, "7"),
    "leave-one-cluster fit did not converge (singular matrix)",
    fixed = TRUE
  )
  expect_error(
    check_jackknife_solver_failure(bad, "7"),
    "jackknife covariance failed for cluster '7': singular matrix",
    fixed = TRUE
  )
  ok <- list(criterion = c(1, 0), beta_mat = matrix(0, 2, 3))
  expect_silent(check_jackknife_solver_failure(ok, "7"))
})
