list2env(get_cached_test_fixtures("estimation_convergence"), environment())

test_that("the solution does not depend on the starting value", {
  ctl <- geer_control(tolerance = 1e-11, maxiter = 500)
  for (m in c("gee", "brgee-naive", "brgee-robust", "pgee-jeffreys")) {
    ref <- geewa(y ~ x,
      data = conv_df, id = id, repeated = visit, family = binomial(), corstr = "exchangeable",
      method = m, beta_start = c(0, 0), control = ctl
    )
    for (b0 in list(c(1, 1), c(-1, 0), c(0.5, 0.5))) {
      fit <- geewa(y ~ x,
        data = conv_df, id = id, repeated = visit, family = binomial(), corstr = "exchangeable",
        method = m, beta_start = b0, control = ctl
      )
      expect_true(fit$converged, label = m)
      expect_equal(coef(fit), coef(ref), tolerance = 1e-7, label = m)
    }
  }
})

test_that("the stopping tolerance bounds the error of the returned estimate", {
  exact <- coef(geewa(y ~ x,
    data = conv_df, id = id, repeated = visit, family = binomial(),
    corstr = "exchangeable", beta_start = c(1, 1), control = geer_control(tolerance = 1e-13)
  ))
  for (tol in c(1e-3, 1e-5, 1e-8)) {
    fit <- geewa(y ~ x,
      data = conv_df, id = id, repeated = visit, family = binomial(),
      corstr = "exchangeable", beta_start = c(1, 1), control = geer_control(tolerance = tol)
    )
    expect_true(fit$converged)
    expect_lte(max(abs(coef(fit) - exact)), tol)
  }
})

test_that("a divergent Newton path from a poor start is flagged, and step_multiplier rescues it", {
  good <- geewa(y ~ x,
    data = conv_df, id = id, repeated = visit, family = binomial(), beta_start = c(0, 0),
    control = geer_control(tolerance = 1e-11)
  )
  expect_warning(
    bad <- geewa(y ~ x,
      data = conv_df, id = id, repeated = visit, family = binomial(),
      beta_start = c(2, -2), control = geer_control(tolerance = 1e-11)
    ),
    "converge"
  )
  expect_false(bad$converged)
  rescued <- geewa(y ~ x,
    data = conv_df, id = id, repeated = visit, family = binomial(), beta_start = c(2, -2),
    control = geer_control(tolerance = 1e-11, step_multiplier = 0.5)
  )
  expect_true(rescued$converged)
  expect_equal(coef(rescued), coef(good), tolerance = 1e-7)
})

test_that("plain GEE does not claim convergence under complete separation, Jeffreys-type methods stay finite", {
  sep <- data.frame(id = rep(1:12, each = 2), visit = rep(1:2, 12))
  withr::local_seed(2)
  sep$x <- rnorm(24)
  sep$y <- as.numeric(sep$x > 0)
  expect_warning(
    g <- geewa(y ~ x,
      data = sep, id = id, repeated = visit, family = binomial(), phi_fixed = TRUE,
      beta_start = c(0, 0), control = geer_control(maxiter = 200)
    ),
    "converge"
  )
  expect_false(g$converged)
  for (m in c("brgee-naive", "pgee-jeffreys")) {
    fit <- geewa(y ~ x,
      data = sep, id = id, repeated = visit, family = binomial(), phi_fixed = TRUE,
      method = m, beta_start = c(0, 0)
    )
    expect_true(fit$converged, label = m)
    expect_true(all(is.finite(coef(fit)) & abs(coef(fit)) < 50), label = m)
  }
  fit_or <- suppressWarnings(geewa_binary(y ~ x,
    data = sep, id = id, repeated = visit, beta_start = c(0, 0),
    control = geer_control(maxiter = 200)
  ))
  expect_false(fit_or$converged)
  fo <- geewa_binary(y ~ x,
    data = sep, id = id, repeated = visit, method = "pgee-jeffreys",
    beta_start = c(0, 0)
  )
  expect_true(fo$converged)
  expect_true(all(abs(coef(fo)) < 50))
})

test_that("an estimated dispersion removes the protection against separation", {
  sep <- withr::with_seed(2, {
    out <- data.frame(id = rep(1:12, each = 2), visit = rep(1:2, 12), x = rnorm(24))
    out$y <- as.numeric(out$x > 0)
    out
  })
  for (m in c("brgee-naive", "pgee-jeffreys")) {
    fixed <- geewa(y ~ x,
      data = sep, id = id, repeated = visit, family = binomial(), phi_fixed = TRUE,
      method = m, beta_start = c(0, 0)
    )
    estimated <- suppressWarnings(geewa(y ~ x,
      data = sep, id = id, repeated = visit, family = binomial(), method = m,
      beta_start = c(0, 0), control = geer_control(maxiter = 200)
    ))
    expect_true(fixed$converged, label = m)
    expect_lt(abs(coef(fixed)[["x"]]), 50, label = m)
    expect_lt(estimated$phi, 1e-8, label = m)
    expect_gt(abs(coef(estimated)[["x"]]), 100, label = m)
  }
  fit_or <- geewa_binary(y ~ x,
    data = sep, id = id, repeated = visit, method = "opgee-jeffreys",
    beta_start = c(0, 0)
  )
  expect_true(fit_or$converged)
  expect_lt(abs(coef(fit_or)[["x"]]), 50)
})

test_that("near-separation: bias-reduced fits stay finite from zero starting values when phi is fixed at 1", {
  withr::local_seed(1)
  x <- rnorm(90)
  y <- as.numeric(x > 0)
  i <- which.min(abs(x))
  y[i] <- 1 - y[i]
  df <- data.frame(id = rep(1:30, each = 3), visit = rep(1:3, 30), x = x, y = y)
  f <- suppressWarnings(
    geewa(y ~ x,
      data = df, id = id, repeated = visit, family = binomial(), corstr = "independence",
      method = "brgee-naive", phi_fixed = TRUE, phi_value = 1, beta_start = c(0, 0)
    )
  )
  expect_true(f$converged)
  expect_true(all(is.finite(coef(f))))
})

test_that("extreme linear predictors give finite estimating equations", {
  eta <- c(-800, -40, 0, 40, 800)
  for (link in c("logit", "probit", "cloglog")) {
    mu <- make.link(link)$linkinv(eta)
    expect_no_error(
      u <- get_estimating_equations_cc(
        c(0, 0, 1, 1, 1), cbind(1, seq_along(eta)), rep(1:5, each = 1),
        rep(1, 5), rep(1, 5), link, "binomial", mu, eta,
        "independence", numeric(0), 1
      )
    )
    expect_true(all(is.finite(u)), label = link)
  }
})

test_that("a single cluster, K = p and K = p + 1 produce estimates and well-defined covariances", {
  df <- conv_df
  one <- suppressWarnings(geewa(y ~ x,
    data = df[df$id == 1, ], id = id, repeated = visit, family = binomial(),
    beta_start = c(0, 0), phi_fixed = TRUE
  ))
  expect_true(all(is.na(one$bias_corrected_covariance)))
  expect_error(
    vcov(one, "bias-corrected"),
    "the bias-corrected covariance requires more clusters than coefficients"
  )
  two <- suppressWarnings(geewa(y ~ x,
    data = df[df$id <= 2, ], id = id, repeated = visit, family = binomial(),
    beta_start = c(0, 0), phi_fixed = TRUE
  ))
  expect_true(all(is.na(two$bias_corrected_covariance)))
  three <- suppressWarnings(geewa(y ~ x,
    data = df[df$id <= 3, ], id = id, repeated = visit, family = binomial(),
    beta_start = c(0, 0), phi_fixed = TRUE
  ))
  expect_true(all(is.finite(vcov(three, "bias-corrected"))))
})

test_that("non-finite covariates are rejected before reaching C++", {
  df <- conv_df
  d1 <- df
  d1$x[2] <- Inf
  expect_error(geewa(y ~ x, data = d1, id = id, family = binomial()))
  expect_error(geewa_binary(y ~ x, data = d1, id = id))
  d3 <- df
  d3$y[5] <- 0.5
  expect_no_error(geewa_binary(y ~ x, data = d3, id = id, beta_start = c(0, 0)))
})

test_that("cloglog fitted values keep full relative accuracy near zero", {
  withr::local_seed(3)
  d <- data.frame(id = rep(1:25, each = 4), x = seq(-4, 4, length.out = 100))
  d$y <- rbinom(100, 1, binomial("cloglog")$linkinv(-6 + 3.5 * d$x))
  fit <- suppressWarnings(geewa(
    y ~ x,
    family = binomial("cloglog"), data = d, id = id, phi_fixed = TRUE
  ))
  reference <- binomial("cloglog")$linkinv(fit$linear.predictors)
  expect_lt(min(reference), 1e-7)
  expect_equal(fit$fitted.values / reference, rep(1, 100), tolerance = 1e-12)
})

test_that("the solver keeps the start when every step-halving is invalid", {
  control <- geer_control(step_maxiter = 1L, step_multiplier = 50)
  res <- collect_warnings(geewa(
    formula = seizures ~ lnbaseline,
    family = poisson(link = "identity"),
    data = test_data$epilepsy,
    id = id,
    repeated = visit,
    beta_start = c(0.05, 0.001),
    control = control
  ))
  expect_false(res$value$converged)
  expect_match(
    res$warnings[[1L]],
    "every step-halving attempt produced invalid fitted values",
    fixed = TRUE
  )
  expect_equal(unname(coef(res$value)), c(0.05, 0.001))
  res <- collect_warnings(geewa_binary(
    formula = ecg ~ treatment,
    link = "identity",
    data = test_data$cerebrovascular,
    id = id,
    repeated = period,
    beta_start = c(0.5, 0.01),
    control = control
  ))
  expect_false(res$value$converged)
  expect_match(
    res$warnings[[1L]],
    "every step-halving attempt produced",
    fixed = TRUE
  )
  expect_equal(unname(coef(res$value)), c(0.5, 0.01))
})

test_that("a failure after a valid step reverts to the last iterate", {
  withr::local_seed(11)
  d <- data.frame(id = rep(1:6, each = 4), t = rep(1:4, 6))
  d$y <- rnorm(24) + rep(c(2, -2, 2, -2), 6) * rnorm(24, sd = 0.2)
  res <- collect_warnings(geewa(
    formula = y ~ 1,
    data = d,
    id = id,
    repeated = t,
    corstr = "unstructured"
  ))
  expect_false(res$value$converged)
  expect_match(
    res$warnings[[1L]],
    "unstructured working correlation matrix is not positive definite",
    fixed = TRUE
  )
})
