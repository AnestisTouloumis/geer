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
ref_fit_covariance <- function(object, indices) {
  mu <- object$fitted.values[indices]
  w <- object$prior.weights[indices]
  tt <- as.integer(object$repeated[indices])
  tmax <- max(as.integer(object$repeated))
  structure <- object$association_structure
  if (identical(object$fit_function, "geewa_binary")) {
    m <- length(mu)
    v_matrix <- diag(mu * (1 - mu), m)
    if (m > 1L && !identical(structure, "independence")) {
      for (i in 1:(m - 1L)) {
        for (j in (i + 1L):m) {
          psi <- object$alpha[ref_pair_pos(tt[i], tt[j], tmax)]
          v_matrix[i, j] <- v_matrix[j, i] <-
            ref_plackett(mu[i], mu[j], psi) - mu[i] * mu[j]
        }
      }
    }
    return(v_matrix / sqrt(outer(w, w)))
  }
  r_matrix <- ref_corr(structure, object$alpha, tmax)
  s <- sqrt(object$phi * object$family$variance(mu) / w)
  outer(s, s) * r_matrix[tt, tt, drop = FALSE]
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

times8 <- list(1:3, c(1, 2, 4), 2:4, 1:4, c(1, 4), 3:4, 1:4, 1:3)

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

ref_small_sample_covariances <- function(object) {
  model_matrix <- object$x
  d_all <- model_matrix * object$family$mu.eta(object$linear.predictors)
  e_all <- object$y - object$fitted.values
  clusters <- split(seq_along(object$id), object$id)
  clusters_no <- length(clusters)
  p <- ncol(model_matrix)
  repeated_max <- max(as.integer(object$repeated))

  pieces <- lapply(clusters, function(indices) {
    v_matrix <- ref_fit_covariance(object, indices)
    d_matrix <- d_all[indices, , drop = FALSE]
    list(
      d = d_matrix,
      vd = solve(v_matrix, d_matrix),
      e = e_all[indices],
      sd = sqrt(diag(v_matrix)),
      positions = as.integer(object$repeated[indices])
    )
  })
  bread <- solve(Reduce(`+`, lapply(pieces, function(z) crossprod(z$vd, z$d))))
  u_list <- lapply(pieces, function(z) drop(crossprod(z$vd, z$e)))
  meat_all <- Reduce(`+`, lapply(u_list, tcrossprod))
  principal_sqrt <- function(a) {
    decomposition <- eigen(a)
    values <- sqrt(as.complex(decomposition$values))
    Re(decomposition$vectors %*% diag(values, nrow(a)) %*%
         solve(decomposition$vectors))
  }
  leverage_inverses <- lapply(pieces, function(z) {
    h_matrix <- z$d %*% bread %*% t(z$vd)
    md <- solve(diag(length(z$e)) - h_matrix)
    list(md = md, kc = principal_sqrt(md))
  })
  sandwich <- function(meat) bread %*% meat %*% bread
  corrected_meat <- function(type) {
    Reduce(`+`, Map(function(z, l) {
      tcrossprod(crossprod(z$vd, l[[type]] %*% z$e))
    }, pieces, leverage_inverses))
  }
  meat_kc <- corrected_meat("kc")
  meat_md <- corrected_meat("md")
  meat_fg <- Reduce(`+`, Map(function(z, u) {
    leverage <- diag(crossprod(z$vd, z$d) %*% bread)
    tcrossprod(u / sqrt(1 - pmin(0.75, leverage)))
  }, pieces, u_list))
  meat_fz <- Reduce(`+`, Map(function(z, l, u) {
    inner <- tcrossprod(z$e) -
      z$d %*% bread %*% (meat_all - tcrossprod(u)) %*% bread %*% t(z$d)
    g_matrix <- crossprod(z$vd, l$md)
    g_matrix %*% inner %*% t(g_matrix)
  }, pieces, leverage_inverses, u_list))

  complete <- all(vapply(pieces, function(z) {
    identical(sort(z$positions), seq_len(repeated_max))
  }, logical(1)))
  pooled_meat <- function(transform) {
    pooled <- matrix(0, repeated_max, repeated_max)
    for (i in seq_along(pieces)) {
      z <- pieces[[i]]
      r <- transform(z, leverage_inverses[[i]]) / z$sd
      pooled[z$positions, z$positions] <-
        pooled[z$positions, z$positions] + tcrossprod(r)
    }
    pooled <- pooled / clusters_no
    Reduce(`+`, lapply(pieces, function(z) {
      w_matrix <- z$sd * z$vd
      t(w_matrix) %*% pooled[z$positions, z$positions] %*% w_matrix
    }))
  }
  na_matrix <- matrix(NA_real_, p, p)
  pan <- if (complete) sandwich(pooled_meat(function(z, l) z$e)) else na_matrix
  list(
    `kauermann-carroll` = sandwich(meat_kc),
    `mancl-derouen` = sandwich(meat_md),
    `ford-westgate` = sandwich((meat_kc + meat_md) / 2),
    `fay-graubard` = sandwich(meat_fg),
    `fan-zhang-zhang` = sandwich(meat_fz),
    pan = pan,
    `gosho-sato-takeuchi` = if (complete && clusters_no > p) {
      clusters_no / (clusters_no - p) * pan
    } else {
      na_matrix
    },
    `wang-long` = if (complete) {
      sandwich(pooled_meat(function(z, l) l$md %*% z$e))
    } else {
      na_matrix
    },
    `westgate-burchett` = if (complete) {
      sandwich(pooled_meat(function(z, l) l$kc %*% z$e))
    } else {
      na_matrix
    }
  )
}

ref_influence <- function(object) {
  d_all <- object$x * object$family$mu.eta(object$linear.predictors)
  e_all <- object$y - object$fitted.values
  clusters <- split(seq_along(object$id), object$id)
  pieces <- lapply(clusters, function(indices) {
    v_matrix <- ref_fit_covariance(object, indices)
    d_matrix <- d_all[indices, , drop = FALSE]
    list(v = v_matrix, d = d_matrix, vd = solve(v_matrix, d_matrix),
         e = e_all[indices])
  })
  bread <- solve(Reduce(`+`, lapply(pieces, function(z) crossprod(z$vd, z$d))))
  observation_dfbeta <- function(z, j) {
    others <- seq_along(z$e)[-j]
    a <- if (length(others)) {
      z$v[j, others] %*% solve(z$v[others, others])
    } else {
      matrix(0, 1, 0)
    }
    d_tilde <- z$d[j, ] - drop(a %*% z$d[others, , drop = FALSE])
    e_tilde <- z$e[j] - sum(a * z$e[others])
    v_tilde <- z$v[j, j] - sum(a * z$v[others, j])
    h_tilde <- drop(d_tilde %*% bread %*% d_tilde) / v_tilde
    drop(bread %*% d_tilde) * e_tilde / (v_tilde * (1 - h_tilde))
  }
  rows <- lapply(pieces, function(z) {
    h_matrix <- z$d %*% bread %*% t(z$vd)
    z_vector <- solve(z$v, z$e)
    list(
      hat = diag(h_matrix),
      cluster = drop(bread %*% t(z$vd) %*% solve(diag(length(z$e)) - h_matrix,
                                                   z$e)),
      observation = t(vapply(seq_along(z$e), observation_dfbeta,
                             numeric(ncol(z$d)), z = z)),
      cluster_delta = drop(crossprod(z$d, z_vector)),
      observation_delta = z$d * z_vector,
      response_delta = sqrt(diag(z$v)) * z$vd
    )
  })
  list(
    hat = unlist(lapply(rows, `[[`, "hat"), use.names = FALSE),
    cluster_dfbeta = do.call(rbind, lapply(rows, `[[`, "cluster")),
    observation_dfbeta = do.call(rbind, lapply(rows, `[[`, "observation")),
    cluster_delta = do.call(rbind, lapply(rows, `[[`, "cluster_delta")),
    observation_delta = do.call(rbind, lapply(rows, `[[`, "observation_delta")),
    response_delta = do.call(rbind, lapply(rows, `[[`, "response_delta")),
    bread = bread
  )
}

ref_estimating_equations <- function(object, beta) {
  eta <- drop(object$x %*% beta) + as.numeric(object$offset)
  mu <- object$family$linkinv(eta)
  shifted <- object
  shifted$fitted.values <- mu
  d_all <- object$x * object$family$mu.eta(eta)
  Reduce(`+`, lapply(split(seq_along(object$id), object$id), function(i) {
    v_matrix <- ref_fit_covariance(shifted, i)
    drop(crossprod(d_all[i, , drop = FALSE],
                   solve(v_matrix, object$y[i] - mu[i])))
  }))
}

ref_gof_residual <- function(fit, statistic, covariance) {
  mu <- fit$fitted.values
  r <- fit$y - mu
  v <- mu * (1 - mu)
  d_matrix <- fit$x * fit$family$mu.eta(fit$linear.predictors)
  clusters <- split(seq_along(fit$id), fit$id)
  n <- length(mu)
  positions_no <- max(fit$repeated)
  standardized <- matrix(NA_real_, length(clusters), positions_no)
  for (k in seq_along(clusters)) {
    i <- clusters[[k]]
    standardized[k, fit$repeated[i]] <- r[i] / sqrt(v[i])
  }
  pooled <- matrix(0, positions_no, positions_no)
  for (j in seq_len(positions_no)) {
    for (l in seq_len(positions_no)) {
      pooled[j, l] <- mean(standardized[, j] * standardized[, l], na.rm = TRUE)
    }
  }
  v_inverse <- matrix(0, n, n)
  cov_y <- matrix(0, n, n)
  for (i in clusters) {
    v_inverse[i, i] <- solve(ref_fit_covariance(fit, i))
    positions <- fit$repeated[i]
    cov_y[i, i] <- if (identical(covariance, "empirical")) {
      tcrossprod(r[i])
    } else {
      diag(sqrt(v[i])) %*% pooled[positions, positions] %*% diag(sqrt(v[i]))
    }
  }
  sensitivity <- crossprod(d_matrix, v_inverse %*% d_matrix)
  h_matrix <- d_matrix %*% solve(sensitivity, t(d_matrix) %*% v_inverse)
  if (identical(statistic, "pearson")) {
    value <- sum(r^2 / v)
    mean_value <- n
    weights <- (1 - 2 * mu) / v
  } else {
    value <- sum(r^2)
    mean_value <- sum(v)
    weights <- 1 - 2 * mu
  }
  adjusted <- crossprod(diag(n) - h_matrix, weights)
  variance <- drop(crossprod(adjusted, cov_y %*% adjusted))
  z <- (value - mean_value) / sqrt(variance)
  c(z = z, p = 2 * pnorm(-abs(z)))
}

ref_pearson <- function(d, b) {
  mu <- d$lk$linkinv(drop(d$X %*% b))
  (d$y - mu) * sqrt(d$w / d$family$variance(mu))
}

ref_alpha <- function(d, r, phi, p, struct, m = 1L) {
  cl <- ref_clusters(d)
  n_occ <- d$tmax
  pairs <- do.call(rbind, lapply(cl, function(idx) {
    n <- length(idx)
    if (n < 2) {
      return(NULL)
    }
    jk <- t(combn(n, 2))
    data.frame(
      s = d$tt[idx][jk[, 1]], t = d$tt[idx][jk[, 2]],
      prod = r[idx][jk[, 1]] * r[idx][jk[, 2]]
    )
  }))
  pairs$lag <- pairs$t - pairs$s
  den <- function(n) (n - p) * phi
  switch(struct,
    exchangeable = sum(pairs$prod) / den(nrow(pairs)),
    ar1 = sum(pairs$prod[pairs$lag == 1]) / den(sum(pairs$lag == 1)),
    toeplitz = sapply(
      seq_len(n_occ - 1),
      function(l) sum(pairs$prod[pairs$lag == l]) / den(sum(pairs$lag == l))
    ),
    "m-dependent" = sapply(
      seq_len(m),
      function(l) sum(pairs$prod[pairs$lag == l]) / den(sum(pairs$lag == l))
    ),
    unstructured = unlist(lapply(1:(n_occ - 1), function(s) {
      sapply((s + 1):n_occ, function(t) {
        sel <- pairs$s == s & pairs$t == t
        sum(pairs$prod[sel]) / den(sum(sel))
      })
    }))
  )
}

ref_corr_alpha <- function(d) {
  if (d$struct == "independence") numeric(0) else d$alpha
}

ref_frechet_cov <- function(p, q) {
  c(max(0, p + q - 1) - p * q, min(p, q) - p * q)
}

ref_or_counts <- function(y, id, tt, w, add, n_occ, exchangeable = FALSE) {
  out <- numeric(0)
  for (s in 1:(n_occ - 1)) {
    for (t in (s + 1):n_occ) {
      n <- matrix(add, 2, 2)
      for (i in unique(id)) {
        a <- which(id == i & tt == s)
        b <- which(id == i & tt == t)
        if (length(a) == 1 && length(b) == 1) {
          n[y[a] + 1, y[b] + 1] <- n[
            y[a] + 1,
            y[b] + 1
          ] + min(w[a], w[b])
        }
      }
      out <- c(out, n[1, 1] * n[2, 2] / (n[1, 2] * n[2, 1]))
    }
  }
  if (exchangeable) rep(exp(mean(log(out))), length(out)) else out
}

compute_reference_ad <- function(x, group) {
  n <- length(x)
  support <- sort(unique(x))
  support <- support[-length(support)]
  counts <- vapply(support, function(value) sum(x == value), numeric(1))
  total <- vapply(support, function(value) sum(x <= value), numeric(1))
  sum(vapply(unique(group), function(index) {
    size <- sum(group == index)
    cumulative <- vapply(support, function(value) {
      sum(x[group == index] <= value)
    }, numeric(1))
    sum(counts * (n * cumulative - size * total)^2 /
          (total * (n - total))) / size
  }, numeric(1))) / n
}

compute_reference_monotone <- function(x) {
  observed <- !is.na(x[, 2L])
  design <- cbind(1, x[observed, 1L])
  coefficients <- solve(crossprod(design),
                        crossprod(design, x[observed, 2L]))
  residuals <- x[observed, 2L] - drop(design %*% coefficients)
  variance <- mean((x[, 1L] - mean(x[, 1L]))^2)
  slope <- coefficients[2L]
  list(
    mu = c(mean(x[, 1L]), coefficients[1L] + slope * mean(x[, 1L])),
    sigma = matrix(c(variance, slope * variance, slope * variance,
                     mean(residuals^2) + slope^2 * variance), 2L)
  )
}

compute_reference_or_score <- function(beta, response, odds_ratio) {
  probabilities <- plogis(beta)
  covariance <- diag(probabilities * (1 - probabilities))
  joint <- ref_plackett(probabilities[1L], probabilities[2L], odds_ratio)
  covariance[1L, 2L] <- covariance[2L, 1L] <-
    joint - prod(probabilities)
  probabilities * (1 - probabilities) *
    drop(solve(covariance, response - probabilities))
}

compute_reference_empirical <- function(response, offset) {
  beta <- uniroot(function(value) {
    sum(response - plogis(value + offset))
  }, c(-20, 20), tol = 1e-12)$root
  probabilities <- plogis(beta + offset)
  residuals <- response - probabilities
  variances <- probabilities * (1 - probabilities)
  information <- sum(variances)
  first <- -sum(variances * residuals)
  second <- -sum(variances * (1 - 2 * probabilities))
  meat <- sum(residuals^2)
  adjustment <- -(first / information +
                    0.5 * meat * second / information^2)
  printed <- first / information -
    0.5 * meat * second / information^2
  list(beta = beta, corrected = beta + adjustment / information,
       printed = beta + printed / information)
}
