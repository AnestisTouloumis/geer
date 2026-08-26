compute_n_association_parameters <- function(object) {
  if (identical(object$association_structure, "independence")) {
    0L
  } else if (object$association_structure %in% c("exchangeable", "ar1")) {
    1L
  } else {
    length(object$alpha)
  }
}


compute_n_estimated_association_parameters <- function(object) {
  if (object$association_structure %in% c("independence", "fixed")) {
    0L
  } else if (object$association_structure %in% c("exchangeable", "ar1")) {
    1L
  } else {
    length(object$alpha)
  }
}


compute_qicc <- function(qic, p, association_params_no, clusters_no) {
  denominator <- clusters_no - p - association_params_no - 1
  if (!is.finite(denominator) || denominator <= 0) {
    return(NA_real_)
  }
  # Hardin and Hilbe (2013), Equation 4.30, subtract this correction.
  correction <- 2 * (p + association_params_no) *
    (p + association_params_no + 1) / denominator
  qic - correction
}


compute_gaussian_pseudolikelihood_criteria <- function(object,
                                                       gpc,
                                                       p,
                                                       association_params_no) {
  parameter_count <- p + association_params_no
  gaussian_deviance <- object$obs_no * log(2 * pi) - 2 * gpc
  list(
    AGPC = gaussian_deviance + 2 * parameter_count,
    SGPC = gaussian_deviance + log(object$clusters_no) * parameter_count
  )
}


compute_bivariate_probability_or <- function(row_prob, col_prob, odds_ratio) {
  ans_independence <- row_prob * col_prob
  tol <- 1e-8
  if (row_prob > 1 - tol || col_prob > 1 - tol ||
      row_prob < tol || col_prob < tol ||
      abs(odds_ratio - 1) < tol) {
    return(ans_independence)
  }
  f_value <- 1 - (1 - odds_ratio) * (row_prob + col_prob)
  root_value <- max(
    0,
    f_value^2 - 4 * odds_ratio * (odds_ratio - 1) * ans_independence
  )
  (f_value - sqrt(root_value)) / (2 * (odds_ratio - 1))
}


compute_upper_triangular_pair_index <- function(row, col, dimension) {
  as.integer((row - 1L) * (2L * dimension - row) / 2L + (col - row))
}


compute_or_working_covariance <- function(object, indices) {
  mu <- object$fitted.values[indices]
  weights <- object$prior.weights[indices]
  repeated <- as.integer(object$repeated[indices])
  cluster_size <- length(indices)
  ans <- diag(mu * (1 - mu) / weights, nrow = cluster_size)
  if (cluster_size <= 1L) {
    return(ans)
  }
  repeated_max <- max(as.integer(object$repeated))
  alpha <- get_or_alpha(object)
  for (j in seq_len(cluster_size - 1L)) {
    for (k in seq.int(j + 1L, cluster_size)) {
      row_time <- repeated[[j]]
      col_time <- repeated[[k]]
      if (row_time > col_time) {
        tmp <- row_time
        row_time <- col_time
        col_time <- tmp
      }
      pair_index <- compute_upper_triangular_pair_index(
        row_time,
        col_time,
        repeated_max
      )
      odds_ratio <- alpha[[pair_index]]
      joint <- compute_bivariate_probability_or(mu[[j]], mu[[k]], odds_ratio)
      covariance <- (joint - mu[[j]] * mu[[k]]) /
        sqrt(weights[[j]] * weights[[k]])
      ans[j, k] <- covariance
      ans[k, j] <- covariance
    }
  }
  ans
}


compute_cc_working_covariance <- function(object, indices) {
  mu <- object$fitted.values[indices]
  weights <- object$prior.weights[indices]
  repeated <- as.integer(object$repeated[indices])
  repeated_max <- max(as.integer(object$repeated))
  correlation <- if (identical(object$association_structure, "independence")) {
    diag(repeated_max)
  } else {
    get_correlation_matrix(
      object$association_structure,
      object$alpha,
      repeated_max
    )
  }
  marginal_sd <- sqrt(object$phi * object$family$variance(mu) / weights)
  correlation[repeated, repeated, drop = FALSE] * tcrossprod(marginal_sd)
}


compute_working_covariance_for_criteria <- function(object, indices) {
  if (is_geewa_fit(object)) {
    compute_cc_working_covariance(object, indices)
  } else {
    compute_or_working_covariance(object, indices)
  }
}


compute_ghyc_pac <- function(object) {
  repeated_max <- max(as.integer(object$repeated))
  empirical_sum <- matrix(0, repeated_max, repeated_max)
  working_sum <- matrix(0, repeated_max, repeated_max)
  pair_counts <- matrix(0, repeated_max, repeated_max)
  residuals <- object$y - object$fitted.values
  cluster_indices <- split(seq_along(object$id), object$id)

  for (indices in cluster_indices) {
    repeated <- as.integer(object$repeated[indices])
    empirical_sum[repeated, repeated] <-
      empirical_sum[repeated, repeated, drop = FALSE] +
      tcrossprod(residuals[indices])
    working_sum[repeated, repeated] <-
      working_sum[repeated, repeated, drop = FALSE] +
      compute_working_covariance_for_criteria(object, indices)
    pair_counts[repeated, repeated] <-
      pair_counts[repeated, repeated, drop = FALSE] + 1
  }

  if (any(pair_counts <= 0)) {
    return(list(GHYC = NA_real_, PAC = NA_real_))
  }

  empirical_mean <- empirical_sum / pair_counts
  working_mean <- working_sum / pair_counts

  working_inverse <- tryCatch(
    solve(working_mean),
    error = function(e) NULL
  )
  ghyc <- if (is.null(working_inverse)) {
    NA_real_
  } else {
    discrepancy <- empirical_mean %*% working_inverse - diag(repeated_max)
    value <- sum(diag(discrepancy %*% discrepancy))
    if (is.finite(value)) value else NA_real_
  }

  working_det <- tryCatch(det(working_mean), error = function(e) NA_real_)
  empirical_det <- tryCatch(det(empirical_mean), error = function(e) NA_real_)
  pac <- if (!is.finite(working_det) || abs(working_det) <= .Machine$double.eps ||
             !is.finite(empirical_det)) {
    NA_real_
  } else {
    value <- abs(empirical_det / working_det - 1)
    if (is.finite(value)) value else NA_real_
  }

  list(GHYC = ghyc, PAC = pac)
}


compute_independence_naive_inverse <- function(object,
                                               mu = object$fitted.values,
                                               eta = object$linear.predictors,
                                               phi = object$phi) {
  get_naive_matrix_inverse_independence(
    object$x,
    object$id,
    object$family$link,
    object$family$family,
    mu,
    eta,
    phi,
    object$prior.weights
  )
}


compute_quasi_loglikelihood_values <- function(y,
                                               mu,
                                               weights,
                                               family_name,
                                               phi) {
  eps <- sqrt(.Machine$double.eps)
  ans <- switch(
    family_name,
    gaussian = -sum(weights * (y - mu)^2) / 2,
    binomial = {
      mu_safe <- pmin(pmax(mu, eps), 1 - eps)
      sum(weights * (y * stats::qlogis(mu_safe) + log1p(-mu_safe)))
    },
    poisson = {
      mu_safe <- pmax(mu, eps)
      sum(weights * (y * log(mu_safe) - mu_safe))
    },
    Gamma = {
      mu_safe <- pmax(mu, eps)
      -sum(weights * (y / mu_safe + log(mu_safe)))
    },
    inverse.gaussian = {
      mu_safe <- pmax(mu, eps)
      -sum(weights * (mu_safe - 0.5 * y) / mu_safe^2)
    },
    stop("'family' is not a recognized distribution", call. = FALSE)
  )
  ans / phi
}


compute_quasi_loglikelihood <- function(object) {
  compute_quasi_loglikelihood_values(
    y = object$y,
    mu = object$fitted.values,
    weights = object$prior.weights,
    family_name = object$family$family,
    phi = object$phi
  )
}


geer_phi_is_fixed <- function(object) {
  if (!is_geewa_fit(object)) {
    return(TRUE)
  }
  phi_fixed <- object$call$phi_fixed
  if (is.null(phi_fixed)) {
    return(FALSE)
  }
  if (is.logical(phi_fixed) && length(phi_fixed) == 1L && !is.na(phi_fixed)) {
    return(isTRUE(phi_fixed))
  }
  FALSE
}


compute_independence_gee_quantities <- function(object) {
  glm_control <- stats::glm.control()
  if (!is.null(object$control$maxiter) && is.finite(object$control$maxiter)) {
    glm_control$maxit <- max(glm_control$maxit, as.integer(object$control$maxiter))
  }
  if (!is.null(object$control$tolerance) &&
      is.finite(object$control$tolerance) && object$control$tolerance > 0) {
    glm_control$epsilon <- min(glm_control$epsilon, object$control$tolerance)
  }
  fit_independence <- suppressWarnings(
    stats::glm.fit(
      x = object$x,
      y = object$y,
      weights = object$prior.weights,
      offset = object$offset,
      family = object$family,
      start = as.numeric(object$coefficients),
      control = glm_control,
      intercept = FALSE
    )
  )
  beta <- as.numeric(fit_independence$coefficients)
  if (!isTRUE(fit_independence$converged) || any(!is.finite(beta))) {
    stop(
      "QICHH could not be computed because the independence GEE fit did not converge to finite coefficients",
      call. = FALSE
    )
  }
  eta <- as.numeric(object$x %*% beta + object$offset)
  mu <- as.numeric(object$family$linkinv(eta))
  if (any(!is.finite(mu))) {
    stop(
      "QICHH could not be computed because the independence fitted means are non-finite",
      call. = FALSE
    )
  }
  family_name <- object$family$family
  if (identical(family_name, "binomial")) {
    phi <- 1
  } else if (geer_phi_is_fixed(object)) {
    phi <- object$phi
  } else {
    variance <- object$family$variance(mu)
    if (any(!is.finite(variance)) || any(variance <= 0)) {
      stop(
        "QICHH could not be computed because the independence variance function is non-positive or non-finite",
        call. = FALSE
      )
    }
    denominator <- object$obs_no - length(beta)
    if (denominator <= 0) {
      stop(
        "QICHH could not be computed because the residual degrees of freedom are non-positive",
        call. = FALSE
      )
    }
    phi <- sum(object$prior.weights * (object$y - mu)^2 / variance) / denominator
    phi <- max(phi, 10 * .Machine$double.eps)
  }
  naive_inverse <- compute_independence_naive_inverse(
    object,
    mu = mu,
    eta = eta,
    phi = phi
  )
  quasi_loglikelihood <- compute_quasi_loglikelihood_values(
    y = object$y,
    mu = mu,
    weights = object$prior.weights,
    family_name = family_name,
    phi = phi
  )
  list(
    coefficients = beta,
    linear.predictors = eta,
    fitted.values = mu,
    phi = phi,
    naive_inverse = naive_inverse,
    quasi_loglikelihood = quasi_loglikelihood
  )
}


compute_eqic_adjusted_variance <- function(mu, family_name, k = 1 / 6) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("'k' must be a single positive finite number", call. = FALSE)
  }
  variance <- switch(
    family_name,
    gaussian = rep.int(1, length(mu)),
    poisson = mu + k,
    Gamma = (mu + k)^2,
    inverse.gaussian = (mu + k)^3,
    binomial = (mu + k) * (1 - mu + k),
    stop("EQIC is not available for this family", call. = FALSE)
  )
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop(
      "EQIC could not be computed because the adjusted variance function is non-positive or non-finite",
      call. = FALSE
    )
  }
  variance
}


compute_eqic_adjusted_deviance <- function(y, mu, family_name, k = 1 / 6) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k <= 0) {
    stop("'k' must be a single positive finite number", call. = FALSE)
  }
  eps <- sqrt(.Machine$double.eps)
  switch(
    family_name,
    gaussian = (y - mu)^2,
    poisson = {
      z <- y + k
      mu_adjusted <- pmax(mu + k, eps)
      2 * (z * log(z / mu_adjusted) - (z - mu_adjusted))
    },
    Gamma = {
      z <- y + k
      mu_adjusted <- pmax(mu + k, eps)
      ratio <- z / mu_adjusted
      2 * (ratio - log(ratio) - 1)
    },
    inverse.gaussian = {
      z <- y + k
      mu_adjusted <- pmax(mu + k, eps)
      1 / z + z / mu_adjusted^2 - 2 / mu_adjusted
    },
    binomial = {
      mu_safe <- pmin(pmax(mu, eps), 1 - eps)
      denominator <- 1 + 2 * k
      2 / denominator * (
        (y + k) * log((y + k) / (mu_safe + k)) +
          (1 - y + k) * log((1 - y + k) / (1 - mu_safe + k))
      )
    },
    stop("EQIC is not available for this family", call. = FALSE)
  )
}


compute_eqic <- function(object, beta_covariance, k = 1 / 6) {
  family_name <- object$family$family
  deviance_contributions <- compute_eqic_adjusted_deviance(
    y = object$y,
    mu = object$fitted.values,
    family_name = family_name,
    k = k
  )
  deviance <- sum(object$prior.weights * deviance_contributions)
  if (!is.finite(deviance) || deviance < 0) {
    stop("EQIC could not be computed because the adjusted deviance is invalid", call. = FALSE)
  }
  if (identical(family_name, "binomial")) {
    phi <- 1
  } else if (geer_phi_is_fixed(object)) {
    phi <- object$phi
  } else {
    phi <- deviance / object$obs_no
    phi <- max(phi, 10 * .Machine$double.eps)
  }
  adjusted_variance <- compute_eqic_adjusted_variance(
    mu = object$fitted.values,
    family_name = family_name,
    k = k
  )
  log_variance_term <- sum(
    object$prior.weights * log(2 * pi * phi * adjusted_variance)
  )
  eqic_independence_inverse <- compute_independence_naive_inverse(
    object,
    phi = phi
  )
  eqic_cic <- sum(eqic_independence_inverse * beta_covariance)
  deviance / phi + log_variance_term + 2 * eqic_cic
}


compute_gee_criteria <- function(object,
                                 cov_type,
                                 digits = NULL,
                                 include_extended = TRUE) {
  quasi_loglikelihood <- compute_quasi_loglikelihood(object)
  sc_wc_stats <- if (is_geewa_fit(object)) {
    get_gee_criteria_sc_cw(
      object$y,
      object$id,
      object$repeated,
      object$family$family,
      object$fitted.values,
      object$association_structure,
      object$alpha,
      object$phi,
      object$prior.weights
    )
  } else {
    get_gee_criteria_sc_cw_or(
      object$y,
      object$id,
      object$repeated,
      object$fitted.values,
      get_or_alpha(object),
      object$prior.weights
    )
  }
  sc_wc_stats <- unlist(sc_wc_stats, use.names = FALSE)
  naive_covariance <- vcov(object, cov_type = "naive")
  beta_covariance <- vcov(object, cov_type = cov_type)
  independence_inverse <- compute_independence_naive_inverse(object)
  p <- length(object$coefficients)
  association_params_no <- compute_n_association_parameters(object)
  estimated_association_params_no <- compute_n_estimated_association_parameters(object)
  gessc <- sc_wc_stats[[1L]] / (object$obs_no - p - association_params_no)
  gpc <- sc_wc_stats[[2L]]
  penalized_gpc <- compute_gaussian_pseudolikelihood_criteria(
    object = object,
    gpc = gpc,
    p = p,
    association_params_no = estimated_association_params_no
  )
  qic_u <- 2 * (p - quasi_loglikelihood)
  cic <- sum(independence_inverse * beta_covariance)
  qic <- 2 * (cic - quasi_loglikelihood)
  qicc <- compute_qicc(
    qic = qic,
    p = p,
    association_params_no = estimated_association_params_no,
    clusters_no = object$clusters_no
  )
  q_matrix <- tryCatch(
    solve(naive_covariance, beta_covariance),
    error = function(e) {
      stop(
        "failed to compute RJC because the naive covariance matrix is singular or invalid",
        call. = FALSE
      )
    }
  )
  rjc_trace <- sum(diag(q_matrix)) / p
  rjc_frobenius <- sum(q_matrix * t(q_matrix)) / p
  rjc <- sqrt((1 - rjc_trace)^2 + (1 - rjc_frobenius)^2)
  if (isTRUE(include_extended)) {
    independence_quantities <- compute_independence_gee_quantities(object)
    qichh_penalty <- sum(independence_quantities$naive_inverse * beta_covariance)
    qichh <- 2 * (qichh_penalty - independence_quantities$quasi_loglikelihood)
    eqic <- compute_eqic(object, beta_covariance = beta_covariance)
    covariance_match <- compute_ghyc_pac(object)
    ans <- data.frame(
      QIC = qic,
      QICHH = qichh,
      QICC = qicc,
      CIC = cic,
      RJC = rjc,
      QICu = qic_u,
      EQIC = eqic,
      GESSC = gessc,
      GPC = gpc,
      AGPC = penalized_gpc$AGPC,
      SGPC = penalized_gpc$SGPC,
      GHYC = covariance_match$GHYC,
      PAC = covariance_match$PAC,
      Parameters = p
    )
    numeric_cols <- c(
      "QIC", "QICHH", "QICC", "CIC", "RJC", "QICu", "EQIC",
      "GESSC", "GPC", "AGPC", "SGPC", "GHYC", "PAC"
    )
  } else {
    ans <- data.frame(
      QIC = qic,
      CIC = cic,
      RJC = rjc,
      QICu = qic_u,
      GESSC = gessc,
      GPC = gpc,
      Parameters = p
    )
    numeric_cols <- c("QIC", "CIC", "RJC", "QICu", "GESSC", "GPC")
  }
  if (!is.null(digits)) {
    digits <- check_nonnegative_integerish(digits, "digits")
    ans[, numeric_cols] <- lapply(ans[, numeric_cols, drop = FALSE], round, digits = digits)
    ans[, "Parameters"] <- as.integer(ans[, "Parameters"])
  }
  ans
}


compute_gee_cic <- function(object, cov_type) {
  independence_inverse <- compute_independence_naive_inverse(object)
  beta_covariance <- vcov(object, cov_type = cov_type)
  sum(independence_inverse * beta_covariance)
}
