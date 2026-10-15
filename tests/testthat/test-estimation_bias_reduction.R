cerebrovascular <- test_data$cerebrovascular

adjusted_residual <- function(d, method, ...) {
  fit <- ref_fit(d, method, ...)
  b <- unname(coef(fit))
  lam <- ref_lambda(d, method, b)
  list(converged = fit$converged, res = ref_U(d, b) + lam, lam = lam)
}

expect_adjusted_equation <- function(d, method, label, ...) {
  r <- adjusted_residual(d, method, ...)
  expect_true(r$converged, label = paste(label, "converged"))
  expect_lt(max(abs(r$res)), 1e-5 * max(1, max(abs(r$lam))))
}

times8 <- list(1:3, c(1, 2, 4), 2:4, 1:4, c(1, 4), 3:4, 1:4, 1:3)

cc_methods <- c("brgee-robust", "brgee-naive", "brgee-empirical", "pgee-jeffreys")

joint_binary <- function(mu, lam) {
  m <- length(mu)
  Y <- as.matrix(expand.grid(rep(list(0:1), m)))
  pr <- combn(m, 2)
  pair <- rowSums(Y[, pr[1, ], drop = FALSE] * Y[, pr[2, ], drop = FALSE])
  th <- qlogis(mu)
  for (it in 1:500) {
    lp <- drop(Y %*% th) + lam * pair
    P <- exp(lp - max(lp))
    P <- P / sum(P)
    marg <- colSums(Y * P)
    th <- th + 0.8 * (qlogis(mu) - qlogis(pmin(pmax(marg, 1e-12), 1 - 1e-12)))
  }
  list(Y = Y, P = drop(P))
}

ctl_tight <- geer_control(tolerance = 1e-12, maxiter = 500)

test_that("correlation-based adjustments satisfy the preprint definitions (non-canonical links, weights)", {
  cases <- list(
    list(binomial(), "logit", "exchangeable", 0.3, 1, NULL),
    list(binomial(), "probit", "ar1", 0.4, 1, c(0.5, 1, 2.5, 1.5)),
    list(binomial(), "cloglog", "exchangeable", 0.2, 1, NULL),
    list(poisson(), "log", "exchangeable", 0.3, 1.7, NULL),
    list(poisson(), "sqrt", "ar1", 0.5, 1, c(0.5, 1, 2.5, 1.5)),
    list(Gamma(), "log", "ar1", 0.4, 0.4, NULL),
    list(Gamma(), "inverse", "exchangeable", 0.3, 0.4, NULL)
  )
  for (cf in cases) {
    beta <- if (cf[[2]] == "inverse") {
      c(1, 0.1, -0.1)
    } else if (cf[[2]] == "sqrt") {
      c(1.5, 0.2, 0.1)
    } else {
      c(0.2, 0.3, -0.4)
    }
    d <- ref_sim("cc", times8, cf[[1]], cf[[2]], cf[[3]], cf[[4]], cf[[5]],
      weights = cf[[6]], beta = beta, seed = 21
    )
    for (m in cc_methods) {
      expect_adjusted_equation(d, m, paste(cf[[1]]$family, cf[[2]], cf[[3]], m))
    }
  }
})

test_that("odds-ratio adjustments satisfy the preprint definitions", {
  d1 <- ref_sim("or", times8, link = "logit", alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2), seed = 21)
  for (m in cc_methods) expect_adjusted_equation(d1, m, paste("or logit", m))
  d2 <- ref_sim("or", times8,
    link = "probit", alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2),
    weights = c(0.5, 1, 2.5, 1.5), seed = 21
  )
  for (m in setdiff(cc_methods, "brgee-empirical")) {
    expect_adjusted_equation(
      d2, m,
      paste("or probit w", m)
    )
  }
})

test_that("adjustments satisfy the definitions for the cauchit and 1/mu^2 links", {
  damped <- geer_control(tolerance = 1e-10, maxiter = 500, step_multiplier = 0.5)
  cases <- list(
    list(binomial(), "cauchit", "exchangeable", 0.3, 1, NULL, c(0.2, 0.3, -0.4)),
    list(binomial(), "cauchit", "ar1", 0.4, 1, c(0.5, 1, 2.5, 1.5), c(0.2, 0.3, -0.4)),
    list(inverse.gaussian(), "1/mu^2", "ar1", 0.4, 0.3, c(0.5, 1, 2.5, 1.5), c(1, 0.1, -0.1))
  )
  for (cf in cases) {
    d <- ref_sim("cc", times8, cf[[1]], cf[[2]], cf[[3]], cf[[4]], cf[[5]],
      weights = cf[[6]], beta = cf[[7]], seed = 21
    )
    start <- unname(coef(ref_fit(d, "gee", damped)))
    for (m in cc_methods) {
      expect_adjusted_equation(d, m, paste(cf[[1]]$family, cf[[2]], cf[[3]], m),
        control = damped, beta_start = start
      )
    }
  }
  d <- ref_sim("or", times8, link = "cauchit", alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2), seed = 21)
  start <- unname(coef(ref_fit(d, "gee", damped)))
  for (m in cc_methods) {
    expect_adjusted_equation(d, m, paste("or cauchit", m),
      control = damped, beta_start = start
    )
  }
})

test_that("brgee-empirical with the odds-ratio structure is correct for non-uniform weights", {
  d <- ref_sim("or", times8,
    link = "probit", alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2),
    weights = c(0.5, 1, 2.5, 1.5), seed = 21
  )
  expect_adjusted_equation(d, "brgee-empirical", "or probit w brgee-empirical")
})

test_that("odds-ratio adjustments are correct with clusters of size one", {
  times <- list(1:3, 2, c(1, 2, 4), 4, 2:4, 1:4, 1, c(1, 4), 3:4, 1:4, 3)
  for (link in c("logit", "probit")) {
    d <- ref_sim("or", times,
      link = link, alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2),
      weights = c(0.5, 1, 2.5, 1.5), seed = 21
    )
    for (m in cc_methods) expect_adjusted_equation(d, m, paste("or singletons", link, m))
  }
})

test_that("Gaussian identity-link adjustments vanish (analytic: G is constant in beta)", {
  d <- ref_sim("cc", times8, gaussian(), "identity", "exchangeable", 0.4, 2, seed = 4)
  gee <- ref_fit(d, "gee")
  for (m in c("brgee-robust", "brgee-naive", "pgee-jeffreys")) {
    expect_equal(coef(ref_fit(d, m)), coef(gee), tolerance = 1e-8, label = m)
  }
})

test_that("variance-stabilizing links give Sigma_0 free of beta, so the Jeffreys penalty vanishes", {
  d <- ref_sim("cc", times8, Gamma(), "log", "ar1", 0.4, 0.4, seed = 21)
  expect_equal(coef(ref_fit(d, "pgee-jeffreys")), coef(ref_fit(d, "gee")), tolerance = 1e-8)
  d <- ref_sim("cc", times8, poisson(), "sqrt", "ar1", 0.5, 1, beta = c(1.5, 0.2, 0.1), seed = 21)
  expect_equal(coef(ref_fit(d, "pgee-jeffreys")), coef(ref_fit(d, "gee")), tolerance = 1e-8)
})

test_that("jeffreys_power scales the penalty linearly (a * d log|Sigma_0| / d beta)", {
  d <- ref_sim("cc", times8, poisson(), "log", "exchangeable", 0.3, 1.7, seed = 21)
  for (a in c(0.25, 1)) {
    fit <- ref_fit(d, "pgee-jeffreys",
      control = geer_control(tolerance = 1e-10, maxiter = 200, jeffreys_power = a)
    )
    b <- unname(coef(fit))
    expect_lt(max(abs(ref_U(d, b) + ref_lambda(d, "pgee-jeffreys", b, a = a))), 1e-6)
  }
})

test_that("geewa bcgee methods return geer objects and one-step updates match converged fits", {
  fit_gee <- geewa(
    formula = ecg ~ treatment + factor(period),
    family = binomial(link = "logit"),
    id = id,
    phi_fixed = TRUE,
    phi_value = 1,
    data = cerebrovascular,
    corstr = "independence",
    method = "gee"
  )
  beta_start_gee <- coef(fit_gee)
  for (variant in c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    fit_conv <- update(fit_gee, method = variant)
    fit_one <- update(fit_gee,
      beta_start = beta_start_gee,
      method = variant, maxiter = 1
    )
    expect_s3_class(fit_conv, "geer")
    expect_equal(coef(fit_conv), coef(fit_one),
      tolerance = 1e-5, label = variant
    )
  }
})

test_that("geewa_binary bcgee methods return geer objects and one-step updates match converged fits", {
  fit_gee <- geewa_binary(
    formula = ecg ~ treatment + factor(period),
    link = "logit",
    id = id,
    data = cerebrovascular,
    orstr = "independence",
    method = "gee"
  )
  beta_start_gee <- coef(fit_gee)
  for (variant in c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    fit_conv <- update(fit_gee, method = variant)
    fit_one <- update(fit_gee,
      beta_start = beta_start_gee,
      method = variant, maxiter = 1
    )
    expect_s3_class(fit_conv, "geer")
    expect_equal(coef(fit_conv), coef(fit_one),
      tolerance = 1e-6, label = variant
    )
  }
})

test_that("bcgee equals the one-step closed form built from the reference Lambda", {
  times8 <- list(1:3, c(1, 2, 4), 2:4, 1:4, c(1, 4), 3:4, 1:4, 1:3)
  ctl <- geer_control(tolerance = 1e-12, maxiter = 200)
  for (kind in c("cc", "or")) {
    d <- if (kind == "cc") {
      ref_sim("cc", times8, poisson(), "log", "ar1", 0.4, 1.5, beta = c(1.5, 0.2, 0.1), seed = 8)
    } else {
      ref_sim("or", times8, link = "logit", alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2), seed = 8)
    }
    b0 <- unname(coef(ref_fit(d, "gee", ctl)))
    for (v in c("robust", "naive", "empirical")) {
      info <- if (v == "empirical") -ref_jac(function(b) ref_U(d, b), b0) else ref_Sigma0(d, b0)
      expected <- b0 + solve(info, ref_lambda(d, paste0("brgee-", v), b0))
      expect_equal(unname(coef(ref_fit(d, paste0("bcgee-", v), ctl))), drop(expected),
        tolerance = 1e-6,
        label = paste(kind, v)
      )
    }
  }
})

test_that("E[robust Lambda] over the exact outcome distribution equals the M-estimation bias adjustment", {
  d0 <- ref_sim("cc", rep(list(1:3), 3), binomial(), "logit", "exchangeable", 0.25, 1, seed = 11)
  b <- d0$beta
  p <- length(b)
  cl <- ref_clusters(d0)
  laws <- lapply(cl, function(idx) joint_binary(ref_cluster(d0, b, idx)$mu, 1.0))
  expect_lt(max(abs(unlist(Map(
    function(J, idx) colSums(J$Y * J$P) - ref_cluster(d0, b, idx)$mu,
    laws, cl
  )))), 1e-12)
  S0 <- 0
  S1 <- 0
  for (c in seq_along(cl)) {
    g <- ref_cluster(d0, b, cl[[c]])
    S0 <- S0 + g$G %*% g$D
    for (k in seq_along(laws[[c]]$P)) {
      u <- drop(g$G %*% (laws[[c]]$Y[k, ] - g$mu))
      S1 <- S1 + laws[[c]]$P[k] * tcrossprod(u)
    }
  }
  Si <- solve(S0)
  SG <- Si %*% S1 %*% Si
  t1 <- numeric(p)
  t2 <- numeric(p)
  for (c in seq_along(cl)) {
    idx <- cl[[c]]
    for (k in seq_along(laws[[c]]$P)) {
      y <- laws[[c]]$Y[k, ]
      Uf <- function(bb) {
        gg <- ref_cluster(d0, bb, idx)
        drop(gg$G %*% (y - gg$mu))
      }
      Jk <- ref_jac(Uf, b)
      Hk <- ref_tensor(Uf, b)
      t1 <- t1 + laws[[c]]$P[k] * drop(Jk %*% Si %*% Uf(b))
      for (r in 1:p) t2[r] <- t2[r] + laws[[c]]$P[k] * sum(SG * Hk[r, , ])
    }
  }
  lambda_ideal <- -(t1 + 0.5 * t2)
  step <- function(y, method) {
    out <- fit_geesolver_cc(
      y, d0$X, as.numeric(d0$id), as.numeric(d0$tt), d0$w, "logit", "binomial", b,
      rep(0, length(y)), 1L, 1e-8, 1L, 1, 0.5, method, 1L, d0$alpha, 1L,
      "exchangeable", 1L, 1, 1L, 1L
    )
    out$beta_mat[, 2] - out$beta_mat[, 1]
  }
  combos <- as.matrix(expand.grid(lapply(laws, function(J) seq_along(J$P))))
  e_lambda <- numeric(p)
  for (r in seq_len(nrow(combos))) {
    y <- unlist(lapply(seq_along(laws), function(c) laws[[c]]$Y[combos[r, c], ]))
    w <- prod(vapply(seq_along(laws), function(c) laws[[c]]$P[combos[r, c]], 0))
    e_lambda <- e_lambda + w * drop(S0 %*% (step(y, "brgee-robust") - step(y, "gee")))
  }
  expect_lt(max(abs(e_lambda - lambda_ideal)), 1e-6)
  expect_gt(max(abs(lambda_ideal)), 0.01)
})

test_that("independence logistic brgee-naive and pgee-jeffreys equal Firth's estimator for a covariate model (phi fixed at 1)", {
  d <- ref_sim("cc", ref_times_unbalanced, binomial(), "logit", "independence", 0,
    phi = 1,
    beta = c(0.3, 0.2, -0.1), seed = 21
  )
  X <- d$X
  y <- d$y
  b <- rep(0, 3)
  for (it in 1:200) {
    mu <- plogis(drop(X %*% b))
    W <- mu * (1 - mu)
    I <- crossprod(X, W * X)
    h <- W * rowSums((X %*% solve(I)) * X)
    st <- solve(I, crossprod(X, y - mu + h * (0.5 - mu)))
    b <- b + drop(st)
    if (max(abs(st)) < 1e-14) break
  }
  for (m in c("brgee-naive", "pgee-jeffreys")) {
    f <- geewa(y ~ x1 + x2,
      data = d$df, id = id, repeated = visit, family = binomial(), corstr = "independence",
      method = m, phi_fixed = TRUE, phi_value = 1, control = ctl_tight
    )
    expect_equal(unname(coef(f)), b, tolerance = 1e-8)
  }
})

test_that("pgee-jeffreys under a non-canonical link maximizes l + 1/2 log|X'WX| (independent optimizer)", {
  for (lk in c("probit", "cloglog")) {
    d <- ref_sim("cc", ref_times_unbalanced, binomial(), lk, "independence", 0,
      phi = 1,
      beta = c(0.3, 0.2, -0.1), seed = 33
    )
    X <- d$X
    y <- d$y
    mk <- make.link(lk)
    pen <- function(b) {
      eta <- drop(X %*% b)
      mu <- pmin(pmax(mk$linkinv(eta), 1e-12), 1 - 1e-12)
      dm <- mk$mu.eta(eta)
      sum(y * log(mu) + (1 - y) * log(1 - mu)) + 0.5 * determinant(crossprod(
        X,
        dm^2 / (mu * (1 - mu)) * X
      ))$modulus[1]
    }
    o <- optim(rep(0, 3), function(b) -pen(b),
      method = "BFGS",
      control = list(reltol = 1e-15, maxit = 2000)
    )
    o <- optim(o$par, function(b) -pen(b),
      method = "Nelder-Mead",
      control = list(reltol = 1e-16, maxit = 5000)
    )
    f <- geewa(y ~ x1 + x2,
      data = d$df, id = id, repeated = visit, family = binomial(link = lk),
      corstr = "independence", method = "pgee-jeffreys", phi_fixed = TRUE, phi_value = 1,
      control = ctl_tight
    )
    expect_equal(unname(coef(f)), o$par, tolerance = 1e-5)
  }
})
