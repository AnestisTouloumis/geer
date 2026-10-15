test_that("Gaussian identity GEE with fixed working correlation is generalized least squares", {
  d <- ref_sim("cc", list(1:3, c(1, 3), 1:3, c(2, 3), 1:3, 1:2, 1:3, 1:3), gaussian(), "identity",
    "exchangeable", 0.35, 1.9,
    weights = c(1, 2, 0.5), seed = 12
  )
  A <- 0
  B <- 0
  for (idx in ref_clusters(d)) {
    V <- ref_V(d, rep(0, length(idx)), idx)
    A <- A + t(d$X[idx, ]) %*% solve(V, d$X[idx, ])
    B <- B + t(d$X[idx, ]) %*% solve(V, d$y[idx])
  }
  beta_gls <- drop(solve(A, B))
  fit <- ref_fit(d, "gee")
  expect_equal(unname(coef(fit)), beta_gls, tolerance = 1e-9)
  expect_lte(fit$iter, 3L)
  M <- 0
  for (idx in ref_clusters(d)) {
    V <- ref_V(d, rep(0, length(idx)), idx)
    u <- t(d$X[idx, ]) %*% solve(V, d$y[idx] - d$X[idx, ] %*% beta_gls)
    M <- M + tcrossprod(u)
  }
  expect_equal(unname(vcov(fit, "naive")), solve(A), tolerance = 1e-9)
  expect_equal(unname(vcov(fit, "robust")), solve(A) %*% M %*% solve(A), tolerance = 1e-9)
})

test_that("intercept-only logistic regression has closed-form GEE, Firth and Jeffreys estimates", {
  for (kind in c("cc", "or")) {
    n_i <- c(3, 2, 4, 1, 3, 2, 3)
    withr::local_seed(5)
    df <- data.frame(id = rep(seq_along(n_i), n_i), visit = unlist(lapply(n_i, seq_len)))
    df$y <- rbinom(nrow(df), 1, 0.3)
    n <- nrow(df)
    s <- sum(df$y)
    fit <- function(method, control = geer_control(tolerance = 1e-12)) {
      if (kind == "cc") {
        geewa(y ~ 1,
          data = df, id = id, repeated = visit, family = binomial(),
          phi_fixed = TRUE, phi_value = 1, method = method, control = control
        )
      } else {
        geewa_binary(y ~ 1,
          data = df, id = id, repeated = visit, method = method,
          control = control
        )
      }
    }
    expect_equal(unname(coef(fit("gee"))), qlogis(s / n), tolerance = 1e-9)
    expect_equal(unname(coef(fit("pgee-jeffreys"))), log((s + 0.5) / (n - s + 0.5)),
      tolerance = 1e-9
    )
    expect_equal(unname(coef(fit("brgee-naive"))), log((s + 0.5) / (n - s + 0.5)), tolerance = 1e-9)
    a <- 0.2
    expect_equal(
      unname(coef(fit(
        "pgee-jeffreys",
        geer_control(tolerance = 1e-12, jeffreys_power = a)
      ))),
      log((s + a) / (n - s + a)),
      tolerance = 1e-9
    )
  }
})

test_that("pgee-jeffreys has a closed form for an all-zero binary response", {
  df <- data.frame(id = rep(1:6, each = 2), y = 0)
  fit <- geewa(y ~ 1,
    data = df, id = id, family = binomial(), phi_fixed = TRUE, phi_value = 1,
    method = "pgee-jeffreys", beta_start = 0, control = geer_control(tolerance = 1e-12)
  )
  expect_true(fit$converged)
  expect_equal(unname(coef(fit)), log(0.5 / (12 + 0.5)), tolerance = 1e-9)
})

test_that("Poisson log intercept-only estimates have closed forms (canonical link)", {
  withr::local_seed(7)
  n_i <- c(4, 4, 3, 4, 2, 4)
  df <- data.frame(id = rep(seq_along(n_i), n_i), visit = unlist(lapply(n_i, seq_len)))
  df$y <- rpois(nrow(df), 3)
  n <- nrow(df)
  s <- sum(df$y)
  ctl <- geer_control(tolerance = 1e-12)
  f <- function(method, phi = 1.7) {
    geewa(y ~ 1,
      data = df, id = id, repeated = visit, family = poisson(),
      phi_fixed = TRUE, phi_value = phi, method = method, control = ctl
    )
  }
  expect_equal(unname(coef(f("gee"))), log(s / n), tolerance = 1e-9)
  expect_equal(unname(coef(f("pgee-jeffreys"))), log((s + 0.5 * 1.7) / n), tolerance = 1e-9)
  expect_equal(unname(coef(f("brgee-naive"))), log((s + 0.5 * 1.7) / n), tolerance = 1e-9)
})

test_that("with one observation per cluster every estimable association structure reduces to independence", {
  withr::local_seed(3)
  df <- data.frame(id = 1:40, visit = 1L, x = rnorm(40))
  df$y <- rpois(40, exp(0.4 + 0.3 * df$x))
  ctl <- geer_control(tolerance = 1e-12)
  for (m in c("gee", "brgee-robust", "brgee-naive", "brgee-empirical", "pgee-jeffreys")) {
    ind <- geewa(y ~ x,
      data = df, id = id, repeated = visit, family = poisson(), phi_fixed = TRUE,
      phi_value = 1, method = m, control = ctl
    )
    fixed <- geewa(y ~ x,
      data = df, id = id, repeated = visit, family = poisson(), phi_fixed = TRUE,
      phi_value = 1, method = m, corstr = "fixed", alpha_vector = numeric(0), control = ctl
    )
    unstr <- geewa(y ~ x,
      data = df, id = id, repeated = visit, family = poisson(), phi_fixed = TRUE,
      phi_value = 1, method = m, corstr = "unstructured", control = ctl
    )
    expect_equal(coef(fixed), coef(ind), tolerance = 1e-9, label = paste(m, "fixed"))
    expect_equal(coef(unstr), coef(ind), tolerance = 1e-9, label = paste(m, "unstructured"))
    expect_equal(vcov(fixed, "robust"), vcov(ind, "robust"), tolerance = 1e-9)
  }
})

test_that("balanced intercept-only GEE returns the sample mean for every working correlation", {
  withr::local_seed(9)
  df <- data.frame(id = rep(1:15, each = 4), visit = rep(1:4, 15))
  df$y <- rnorm(60, 2)
  for (alpha in c(-0.2, 0, 0.5, 0.9)) {
    R <- ref_corr("exchangeable", alpha, 4)
    fit <- geewa(y ~ 1,
      data = df, id = id, repeated = visit, corstr = "fixed",
      alpha_vector = R[lower.tri(R)], phi_fixed = TRUE, phi_value = 1
    )
    expect_equal(unname(coef(fit)), mean(df$y), tolerance = 1e-10)
  }
})
