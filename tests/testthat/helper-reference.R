ref_jac <- function(f, x, h = 1e-5) {
  p <- length(x)
  f0 <- f(x)
  out <- matrix(0, length(f0), p)
  for (j in seq_len(p)) {
    e <- numeric(p)
    e[j] <- h
    out[, j] <- (f(x + e) - f(x - e)) / (2 * h)
  }
  out
}
ref_tensor <- function(f, b, h = 2e-4, hin = 1e-4) {
  p <- length(b)
  out <- array(0, c(p, p, p))
  for (s in seq_len(p)) {
    e <- numeric(p)
    e[s] <- h
    out[, s, ] <- (ref_jac(f, b + e, hin) - ref_jac(f, b - e, hin)) / (2 * h)
  }
  out
}

ref_corr <- function(struct, alpha, n) {
  R <- switch(struct,
    independence = diag(n),
    exchangeable = {
      R <- matrix(alpha[1L], n, n)
      diag(R) <- 1
      R
    },
    ar1 = alpha[1L]^abs(outer(seq_len(n), seq_len(n), "-")),
    fixed = ,
    unstructured = {
      R <- diag(n)
      R[lower.tri(R)] <- alpha
      R + t(R) - diag(n)
    },
    stop("unsupported structure in reference")
  )
  R
}
ref_pair_pos <- function(i, j, n) (i - 1) * (2 * n - i) / 2 + (j - i)

ref_plackett <- function(p, q, psi) {
  if (abs(psi - 1) < 1e-300) {
    return(p * q)
  }
  f <- 1 + (psi - 1) * (p + q)
  2 * psi * p * q / (f + sqrt(f^2 - 4 * psi * (psi - 1) * p * q))
}

ref_sim <- function(kind = c("cc", "or"), times, family = binomial(), link = "logit",
                    struct = "exchangeable", alpha = 0.3, phi = 1, weights = NULL,
                    beta = c(0.2, 0.3, -0.4), seed = 1L) {
  kind <- match.arg(kind)
  family <- get(family$family, mode = "function")(link = link)
  withr::local_seed(seed)
  N <- length(times)
  n_i <- lengths(times)
  id <- rep(seq_len(N), n_i)
  tt <- unlist(times)
  X <- cbind(1, rep(rnorm(N), n_i), rnorm(length(id)))
  lk <- make.link(link)
  mu <- lk$linkinv(drop(X %*% beta))
  y <- switch(family$family,
    binomial = rbinom(length(mu), 1, mu),
    poisson = rpois(length(mu), mu),
    gaussian = mu + rnorm(length(mu)),
    Gamma = rgamma(length(mu), 3, 3 / mu),
    inverse.gaussian = mu * (1 + 0.3 * rexp(length(mu)))
  )
  w <- if (is.null(weights)) rep(1, length(y)) else rep_len(weights, length(y))
  df <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3], id = id, visit = tt, w = w)
  list(
    kind = kind, y = y, X = X, id = id, tt = tt, w = w, family = family, link = link,
    lk = lk, struct = struct, alpha = alpha, phi = phi, beta = beta, N = N,
    tmax = max(tt), df = df
  )
}

ref_V <- function(d, mu, idx) {
  if (d$kind == "cc") {
    R <- ref_corr(d$struct, d$alpha, d$tmax)
    s <- sqrt(d$family$variance(mu) / d$w[idx])
    return(d$phi * outer(s, s) * R[d$tt[idx], d$tt[idx], drop = FALSE])
  }
  m <- length(mu)
  V <- diag(mu * (1 - mu), m)
  tt <- d$tt[idx]
  if (m > 1L) {
    for (i in 1:(m - 1L)) {
      for (j in (i + 1L):m) {
        psi <- d$alpha[ref_pair_pos(tt[i], tt[j], d$tmax)]
        V[i, j] <- V[j, i] <- ref_plackett(mu[i], mu[j], psi) - mu[i] * mu[j]
      }
    }
  }
  V / sqrt(outer(d$w[idx], d$w[idx]))
}
ref_cluster <- function(d, b, idx) {
  X <- d$X[idx, , drop = FALSE]
  eta <- drop(X %*% b)
  mu <- d$lk$linkinv(eta)
  D <- d$lk$mu.eta(eta) * X
  V <- ref_V(d, mu, idx)
  list(mu = mu, D = D, V = V, G = t(D) %*% solve(V))
}
ref_clusters <- function(d) split(seq_along(d$id), d$id)

ref_U <- function(d, b = d$beta) {
  U <- 0
  for (idx in ref_clusters(d)) {
    g <- ref_cluster(d, b, idx)
    U <- U + drop(g$G %*% (d$y[idx] - g$mu))
  }
  U
}
ref_Sigma0 <- function(d, b = d$beta) {
  S <- 0
  for (idx in ref_clusters(d)) {
    g <- ref_cluster(d, b, idx)
    S <- S + g$G %*% g$D
  }
  S
}
ref_covariances <- function(d, b = d$beta) {
  p <- length(b)
  S0 <- 0
  M <- 0
  for (idx in ref_clusters(d)) {
    g <- ref_cluster(d, b, idx)
    u <- drop(g$G %*% (d$y[idx] - g$mu))
    S0 <- S0 + g$G %*% g$D
    M <- M + tcrossprod(u)
  }
  naive <- solve(S0)
  robust <- naive %*% M %*% naive
  K <- length(ref_clusters(d))
  N <- length(d$id)
  kappa <- ((N - 1) / (N - p)) * (K / (K - 1))
  lambda <- min(0.5, p / (K - p))
  xi <- max(1, sum(diag(naive %*% M)) / p)
  list(naive = naive, robust = robust, bc = kappa * robust + lambda * xi * naive)
}

ref_lambda <- function(d, method, b = d$beta, a = 0.5) {
  p <- length(b)
  cl <- ref_clusters(d)
  if (method == "pgee-jeffreys") {
    ld <- function(bb) {
      S <- 0
      for (idx in cl) {
        g <- ref_cluster(d, bb, idx)
        S <- S + g$G %*% g$D
      }
      as.numeric(determinant(S)$modulus)
    }
    return(a * ref_jac(function(bb) ld(bb), b, 1e-5)[1, ])
  }
  S0 <- 0
  M <- 0
  Jsum <- 0
  dN <- dR <- dE <- uA <- uO <- array(0, c(p, p, p))
  for (idx in cl) {
    g <- ref_cluster(d, b, idx)
    S <- d$y[idx] - g$mu
    Ui <- drop(g$G %*% S)
    S0 <- S0 + g$G %*% g$D
    M <- M + tcrossprod(Ui)
    Gfun <- function(bb) ref_cluster(d, bb, idx)$G
    A <- ref_jac(function(bb) drop(Gfun(bb) %*% S), b)
    Uf <- function(bb) {
      gg <- ref_cluster(d, bb, idx)
      drop(gg$G %*% (d$y[idx] - gg$mu))
    }
    Ji <- ref_jac(Uf, b)
    Jsum <- Jsum + Ji
    mu0 <- g$mu
    Ef <- function(bb) {
      gg <- ref_cluster(d, bb, idx)
      drop(gg$G %*% (mu0 - gg$mu))
    }
    uA <- uA + ref_tensor(Ef, b)
    if (method == "brgee-empirical") uO <- uO + ref_tensor(Uf, b)
    for (r in 1:p) {
      for (s in 1:p) {
        for (t in 1:p) {
          dR[r, s, t] <- dR[r, s, t] + A[r, s] * Ui[t]
          dE[r, s, t] <- dE[r, s, t] + Ji[r, s] * Ui[t]
        }
      }
    }
    for (s in 1:p) {
      e <- numeric(p)
      e[s] <- 1e-5
      dGs <- (Gfun(b + e) - Gfun(b - e)) / 2e-5
      dN[, s, ] <- dN[, s, ] + dGs %*% g$V %*% t(g$G)
    }
  }
  Si <- solve(S0)
  switch(method,
    "brgee-robust" = {
      SG <- Si %*% M %*% Si
      sapply(1:p, function(r) -(sum(t(Si) * dR[r, , ]) + 0.5 * sum(SG * uA[r, , ])))
    },
    "brgee-naive" =
      sapply(1:p, function(r) -(sum(t(Si) * dN[r, , ]) + 0.5 * sum(Si * uA[r, , ]))),
    "brgee-empirical" = {
      Fi <- solve(-Jsum)
      SG <- Fi %*% M %*% t(Fi)
      sapply(1:p, function(r) -(sum(t(Fi) * dE[r, , ]) + 0.5 * sum(SG * uO[r, , ])))
    },
    stop("unknown method")
  )
}

ref_times_unbalanced <- c(
  rep(list(1:4), 4), rep(list(c(1, 2, 4)), 3), rep(list(c(2, 3, 4)), 2),
  rep(list(c(1, 3)), 2), rep(list(1:3), 2)
)

ref_fit_unquiet <- function(
    d, method = "gee",
    control = geer_control(tolerance = 1e-10, maxiter = 200),
    beta_start = NULL) {
  f <- y ~ x1 + x2
  if (d$kind == "cc") {
    n <- d$tmax
    geewa(f,
      data = d$df, id = id, repeated = visit, family = d$family,
      corstr = "fixed", alpha_vector = ref_corr(d$struct, d$alpha, n)[lower.tri(diag(n))],
      phi_fixed = TRUE, phi_value = d$phi, weights = w, method = method,
      control = control, beta_start = beta_start
    )
  } else {
    geewa_binary(f,
      data = d$df, id = id, repeated = visit, link = d$link,
      orstr = "fixed", alpha_vector = d$alpha, weights = w, method = method,
      control = control, beta_start = beta_start
    )
  }
}

ref_fit <- function(...) {
  withCallingHandlers(ref_fit_unquiet(...), warning = function(w) {
    if (grepl("non-integer #successes", conditionMessage(w),
      fixed = TRUE
    )) {
      invokeRestart("muffleWarning")
    }
  })
}

expect_rel <- function(object, expected, tolerance, label = NULL) {
  err <- max(abs(object - expected) / abs(expected))
  expect_true(is.finite(err) && err <= tolerance,
    label = sprintf(
      "%s: max relative error %.3g <= %.3g",
      if (is.null(label)) "expect_rel" else label, err, tolerance
    )
  )
}
