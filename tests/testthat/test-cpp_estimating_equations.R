cc_ee <- function(d, b) {
  eta <- drop(d$X %*% b)
  mu <- d$lk$linkinv(eta)
  drop(get_estimating_equations_cc(
    d$y, d$X, as.numeric(d$id), as.numeric(d$tt), d$w, d$link, d$family$family,
    mu, eta, d$struct, ref_corr_alpha(d), d$phi
  ))
}

or_ee <- function(d, b) {
  eta <- drop(d$X %*% b)
  mu <- d$lk$linkinv(eta)
  drop(get_estimating_equations_or(
    d$y, d$X, as.numeric(d$id), as.numeric(d$tt), d$w, d$link, mu, eta, d$alpha
  ))
}

test_that("get_estimating_equations_cc equals the textbook sum of D' V^-1 (y - mu)", {
  cfg <- list(
    list(binomial(), "logit", "exchangeable", 0.3, 1),
    list(binomial(), "probit", "ar1", 0.4, 1),
    list(binomial(), "cloglog", "exchangeable", 0.2, 1),
    list(binomial(), "cauchit", "ar1", 0.2, 1),
    list(poisson(), "log", "exchangeable", 0.3, 1.7),
    list(poisson(), "sqrt", "ar1", 0.5, 1),
    list(poisson(), "identity", "exchangeable", 0.3, 1.2),
    list(gaussian(), "identity", "exchangeable", 0.5, 2),
    list(gaussian(), "log", "ar1", 0.3, 1.4),
    list(Gamma(), "inverse", "exchangeable", 0.3, 0.4),
    list(Gamma(), "log", "ar1", 0.4, 0.4),
    list(inverse.gaussian(), "1/mu^2", "exchangeable", 0.3, 0.3),
    list(inverse.gaussian(), "log", "exchangeable", 0.3, 0.3)
  )
  times <- list(1:3, c(1, 2, 4), 2:4, 1:4, c(1, 4), 3:4)
  for (cf in cfg) {
    beta <- if (cf[[2]] %in% c("inverse", "1/mu^2")) {
      c(1, 0.1, -0.1)
    } else if (cf[[1]]$family == "poisson" && cf[[2]] %in% c("sqrt", "identity")) {
      c(1.5, 0.2, 0.1)
    } else {
      c(0.2, 0.3, -0.4)
    }
    for (wts in list(NULL, c(0.5, 1, 2.5, 1.5))) {
      d <- ref_sim("cc", times, cf[[1]], cf[[2]], cf[[3]], cf[[4]], cf[[5]],
        weights = wts, beta = beta, seed = 11
      )
      expect_equal(cc_ee(d, beta), ref_U(d, beta),
        tolerance = 1e-10,
        label = paste(cf[[1]]$family, cf[[2]], cf[[3]], if (is.null(wts)) "w=1" else "w!=1")
      )
    }
  }
})

test_that("get_estimating_equations_cc handles every correlation structure", {
  times <- list(1:4, c(1, 2, 4), 2:4, 1:4)
  a_un <- c(0.1, 0.2, 0.15, 0.25, 0.3, 0.05)
  d <- ref_sim("cc", times, poisson(), "log", "unstructured", a_un, 1.3,
    beta = c(1, 0.2, 0.1),
    seed = 3
  )
  expect_equal(cc_ee(d, d$beta), ref_U(d), tolerance = 1e-10)
  d$struct <- "independence"
  d$alpha <- numeric(0)
  expect_equal(cc_ee(d, d$beta), ref_U(d), tolerance = 1e-10)
})

test_that("get_estimating_equations_or equals the textbook sum under the Plackett working model", {
  tmax <- 4
  alpha <- c(1.5, 2.2, 0.8, 3, 1.1, 2)
  times <- list(1:4, c(1, 2, 4), 2:4, 1:4, c(1, 3), 3:4)
  for (link in c("logit", "probit", "cloglog")) {
    for (wts in list(NULL, c(0.5, 1, 2.5, 1.5))) {
      d <- ref_sim("or", times, link = link, alpha = alpha, weights = wts, seed = 5)
      expect_equal(or_ee(d, d$beta), ref_U(d),
        tolerance = 1e-10,
        label = paste(link, if (is.null(wts)) "w=1" else "w!=1")
      )
    }
  }
})

test_that("the C++ solution solves the independent reference estimating equation", {
  for (kind in c("cc", "or")) {
    d <- if (kind == "cc") {
      ref_sim("cc", list(1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3), binomial(), "logit",
        "exchangeable", 0.3, 1,
        seed = 2
      )
    } else {
      ref_sim("or", rep(list(1:3), 8), link = "logit", alpha = c(2, 1.5, 3), seed = 2)
    }
    fit <- ref_fit(d, "gee")
    expect_true(fit$converged)
    expect_lt(max(abs(ref_U(d, unname(coef(fit))))), 1e-7)
  }
})

test_that("naive, robust and bias-corrected covariances match the reference sandwich", {
  times <- list(1:3, c(1, 2, 4), 2:4, 1:4, c(1, 4), 3:4, 1:4, 1:3)
  for (kind in c("cc", "or")) {
    d <- if (kind == "cc") {
      ref_sim("cc", times, poisson(), "log", "ar1", 0.4, 1.5, beta = c(1.5, 0.2, 0.1), seed = 8)
    } else {
      ref_sim("or", times, link = "logit", alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2), seed = 8)
    }
    fit <- ref_fit(d, "gee")
    b <- unname(coef(fit))
    ref <- ref_covariances(d, b)
    expect_equal(unname(vcov(fit, cov_type = "naive")), ref$naive, tolerance = 1e-7)
    expect_equal(unname(vcov(fit, cov_type = "robust")), ref$robust, tolerance = 1e-7)
    expect_equal(unname(vcov(fit, cov_type = "bias-corrected")), ref$bc, tolerance = 1e-7)
    expect_equal(unname(vcov(fit, cov_type = "robust")), t(unname(vcov(fit, cov_type = "robust"))),
      tolerance = 1e-12
    )
  }
})

test_that("get_information_matrix_independence equals X' diag(w d^2 / (phi V)) X", {
  for (cf in list(
    list(poisson(), "log"), list(binomial(), "probit"), list(Gamma(), "inverse"),
    list(inverse.gaussian(), "1/mu^2"), list(gaussian(), "identity")
  )) {
    beta <- if (cf[[2]] %in% c("inverse", "1/mu^2")) c(1, 0.1, -0.1) else c(0.2, 0.3, -0.4)
    d <- ref_sim("cc", list(1:3, 1:2, 1:3, 1:3, 1:4), cf[[1]], cf[[2]], "independence", numeric(0),
      1.4,
      weights = c(1, 2, 0.5), beta = beta, seed = 71
    )
    eta <- drop(d$X %*% beta)
    mu <- d$lk$linkinv(eta)
    dmu <- d$lk$mu.eta(eta)
    ref <- t(d$X) %*% (d$X * (d$w * dmu^2 / (d$phi * d$family$variance(mu))))
    got <- get_information_matrix_independence(
      d$X, as.numeric(d$id), d$link, d$family$family,
      mu, eta, d$phi, d$w
    )
    expect_equal(got, ref, tolerance = 1e-12, label = paste(cf[[1]]$family, cf[[2]]))
  }
})

test_that("working-covariance criteria kernels equal sum s'V^-1 s and -(sum + sum logdet V)/2", {
  times <- list(1:3, c(1, 2, 4), 2:4, 1:4, c(1, 4))
  d <- ref_sim("cc", times, poisson(), "log", "ar1", 0.4, 1.4,
    weights = c(1, 2, 0.5),
    beta = c(1, 0.2, 0.1), seed = 72
  )
  mu <- d$lk$linkinv(drop(d$X %*% d$beta))
  sc <- 0
  ld <- 0
  for (idx in ref_clusters(d)) {
    V <- ref_V(d, mu[idx], idx)
    s <- d$y[idx] - mu[idx]
    sc <- sc + drop(t(s) %*% solve(V, s))
    ld <- ld + as.numeric(determinant(V)$modulus)
  }
  got <- get_working_covariance_criteria_cc(
    d$y, as.numeric(d$id), as.numeric(d$tt), "poisson", mu,
    "ar1", 0.4, 1.4, d$w
  )
  expect_equal(got$sc, sc, tolerance = 1e-10)
  expect_equal(got$gp, -0.5 * (sc + ld), tolerance = 1e-10)
  d <- ref_sim("or", times,
    link = "logit", alpha = c(1.5, 2.2, 0.8, 3, 1.1, 2),
    weights = c(1, 2, 0.5), seed = 72
  )
  mu <- d$lk$linkinv(drop(d$X %*% d$beta))
  sc <- 0
  ld <- 0
  for (idx in ref_clusters(d)) {
    V <- ref_V(d, mu[idx], idx)
    s <- d$y[idx] - mu[idx]
    sc <- sc + drop(t(s) %*% solve(V, s))
    ld <- ld + as.numeric(determinant(V)$modulus)
  }
  got <- get_working_covariance_criteria_or(
    d$y, as.numeric(d$id), as.numeric(d$tt), mu, d$alpha,
    d$w
  )
  expect_equal(got$sc, sc, tolerance = 1e-10)
  expect_equal(got$gp, -0.5 * (sc + ld), tolerance = 1e-10)
})

test_that("criteria kernels fail with cluster context when V_i is not positive definite", {
  y <- c(1, 0, 1)
  mu <- rep(0.5, 3)
  expect_error(
    get_working_covariance_criteria_or(y, c(7, 7, 7), 1:3, mu, c(1e4, 1e4, 1e-4), rep(1, 3)),
    "not positive definite"
  )
  expect_error(
    get_working_covariance_criteria_or(y, c(7, 7, 7), 1:3, mu, c(1e4, 1e4, 1e-4), rep(1, 3)),
    "id = 7"
  )
})

test_that("a converged odds-ratio fit has positive-definite working covariances in every cluster", {
  times <- list(1:3, 1:3, c(1, 2, 3), 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3)
  for (a in list(c(1.5, 2.2, 3), c(0.4, 5, 0.9), c(8, 8, 8))) {
    d <- ref_sim("or", times, link = "logit", alpha = a, seed = 73)
    fit <- suppressWarnings(ref_fit(d, "gee"))
    if (fit$converged) {
      for (idx in split(seq_along(fit$id), fit$id)) {
        V <- get_working_covariance(fit, idx)
        expect_gt(min(eigen(V, symmetric = TRUE, only.values = TRUE)$values), 0)
      }
    }
  }
})
