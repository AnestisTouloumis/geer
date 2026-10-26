test_native_fixtures <- list(cpp_contract = function() {
  cc_args <- function(...) {
    withr::local_seed(1)
    N <- 10
    a <- list(
      y_vector = rbinom(30, 1, 0.5), model_matrix = cbind(
        1,
        rnorm(30)
      ), id_vector = rep(1:N, each = 3), repeated_vector = rep(
        1:3,
        N
      ), weights_vector = rep(1, 30), link = "logit", family = "binomial",
      beta_vector = c(0, 0), offset = rep(0, 30), maxiter = 50L,
      tolerance = 1e-10,
      step_maxiter = 10L, step_multiplier = 1, jeffreys_power = 0.5,
      method = "gee", use_params = 1L, alpha_vector = numeric(0),
      alpha_fixed = 1L,
      correlation_structure = "independence", mdependence = 1L, phi = 1,
      phi_fixed = 1L, hold_nuisance = 0L
    )
    mod <- list(...)
    a[names(mod)] <- mod
    a
  }
  solve_cc <- function(...) do.call(fit_geesolver_cc, cc_args(...))
  or_solver_x <- c(
    -1.2, 0.3, 0.8, -0.5, 1.1, -0.2, 0.9, -1.4, 0.1, 0.4,
    -0.7, 1.3, -0.9, 0.6, -0.1, 1.5, -0.3, 0.2, -1.1, 0.7, -0.6, 0, 1,
    -0.8
  )
  or_solver_y <- c(
    0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1,
    0, 1, 0, 0, 1, 0
  )
  call_or_solver <- function(y) {
    fit_geesolver_or(
      y, cbind(1, or_solver_x), rep(1:8, each = 3), rep(1:3,
        times = 8
      ), rep(1, 24), "logit", c(0, 0), rep(0, 24), 100L, 1e-06,
      10L, 1, 0.5, "gee", rep(1.5, 3)
    )
  }
  as.list(environment(), all.names = TRUE)
}, cpp_estimating_equations = function() {
  cc_ee <- function(d, b) {
    eta <- drop(d$X %*% b)
    mu <- d$lk$linkinv(eta)
    drop(get_estimating_equations_cc(
      d$y, d$X, as.numeric(d$id),
      as.numeric(d$tt),
      d$w, d$link, d$family$family, mu, eta, d$struct, ref_corr_alpha(d),
      d$phi
    ))
  }
  or_ee <- function(d, b) {
    eta <- drop(d$X %*% b)
    mu <- d$lk$linkinv(eta)
    drop(get_estimating_equations_or(
      d$y, d$X, as.numeric(d$id),
      as.numeric(d$tt),
      d$w, d$link, mu, eta, d$alpha
    ))
  }
  as.list(environment(), all.names = TRUE)
}, cpp_odds_ratio_model = function() {
  or_V <- function(mu, psi, w = c(1, 1)) {
    ee <- function(r) {
      get_estimating_equations_or(
        mu + r, diag(2), c(1, 1), c(1, 2),
        w, "identity", mu, mu, psi
      )
    }
    solve(cbind(ee(c(1, 0)), ee(c(0, 1))))
  }
  as.list(environment(), all.names = TRUE)
})
