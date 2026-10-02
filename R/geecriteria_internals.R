## A criterion that cannot be evaluated for a particular fit is reported as NA
## rather than aborting the whole table: geecriteria() compares several models
## at once, and a single degenerate candidate should not hide the remaining
## criteria or the remaining models. Failures of argument validation remain
## ordinary errors, so only this condition class is converted.
geer_criterion_unavailable <- function(message) {
  stop(structure(
    list(message = message, call = NULL),
    class = c("geer_criterion_unavailable", "error", "condition")
  ))
}


compute_criterion_or_na <- function(expr) {
  tryCatch(
    expr,
    geer_criterion_unavailable = function(cnd) NA_real_
  )
}


geer_criteria_columns <- c(
  "QIC", "QICHH", "QICC", "CIC", "RJC", "QICu", "EQIC", "GESSC",
  "GPC", "AGPC", "SGPC", "GHYC", "PAC", "PT", "WR", "RMR"
)


geer_criteria_columns_basic <- c(
  "QIC", "CIC", "RJC", "QICu", "GESSC", "GPC"
)


normalize_geer_criteria <- function(criteria) {
  if (!is.character(criteria) || length(criteria) < 1L || anyNA(criteria)) {
    stop(
      "'criteria' must be a non-missing character vector",
      call. = FALSE
    )
  }
  if (any(toupper(criteria) == "ALL")) {
    if (length(criteria) > 1L) {
      stop(
        "'criteria' must be either \"all\" or a selection of criterion names, not both",
        call. = FALSE
      )
    }
    return(geer_criteria_columns)
  }
  ## Matching ignores case so that the usual lower-case spellings are accepted,
  ## but the canonical names are returned and the requested order is kept.
  matched <- match(toupper(criteria), toupper(geer_criteria_columns))
  if (anyNA(matched)) {
    stop(
      sprintf(
        "unknown entries in 'criteria': %s. Available criteria are %s, or \"all\"",
        paste(criteria[is.na(matched)], collapse = ", "),
        paste(geer_criteria_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  geer_criteria_columns[unique(matched)]
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
  cluster_indices <- split(seq_along(object$id), object$id)

  ## Gosho, Hamada and Yoshimura (2011) and Pardo and Alonso (2019) both define
  ## their criteria for balanced designs: the empirical and working covariance
  ## matrices are sums of cluster-level matrices, which are only conformable
  ## when every cluster contributes the same repeated positions. Both criteria
  ## are therefore reported only when each cluster observes every position
  ## exactly once.
  positions <- seq_len(repeated_max)
  balanced <- all(vapply(
    cluster_indices,
    function(indices) {
      identical(sort(as.integer(object$repeated[indices])), positions)
    },
    logical(1L)
  ))
  if (!balanced) {
    return(list(GHYC = NA_real_, PAC = NA_real_))
  }

  clusters_no <- length(cluster_indices)
  empirical_sum <- matrix(0, repeated_max, repeated_max)
  working_sum <- matrix(0, repeated_max, repeated_max)
  residuals <- object$y - object$fitted.values

  for (indices in cluster_indices) {
    repeated <- as.integer(object$repeated[indices])
    empirical_sum[repeated, repeated] <-
      empirical_sum[repeated, repeated, drop = FALSE] +
      tcrossprod(residuals[indices])
    working_sum[repeated, repeated] <-
      working_sum[repeated, repeated, drop = FALSE] +
      compute_working_covariance_for_criteria(object, indices)
  }

  empirical_mean <- empirical_sum / clusters_no
  working_mean <- working_sum / clusters_no

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


## Fitted means are clamped away from the boundaries of their support before
## the logarithms are taken. The quasi-likelihood diverges there, and a
## non-finite first term would propagate to QIC, QICu, QICC and QICHH at once.
## The clamp keeps those criteria finite and comparable across working
## association structures, for which the first term is very nearly constant.
## See the documentation of geecriteria() for the caveat this carries when
## candidate models differ in whether they produce boundary fitted values.
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
    geer_criterion_unavailable(
      "QICHH could not be computed because the independence GEE fit did not converge to finite coefficients"
    )
  }
  eta <- as.numeric(object$x %*% beta + object$offset)
  mu <- as.numeric(object$family$linkinv(eta))
  if (any(!is.finite(mu))) {
    geer_criterion_unavailable(
      "QICHH could not be computed because the independence fitted means are non-finite"
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
      geer_criterion_unavailable(
        "QICHH could not be computed because the independence variance function is non-positive or non-finite"
      )
    }
    denominator <- object$obs_no - length(beta)
    if (denominator <= 0) {
      geer_criterion_unavailable(
        "QICHH could not be computed because the residual degrees of freedom are non-positive"
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
    geer_criterion_unavailable("EQIC is not available for this family")
  )
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    geer_criterion_unavailable(
      "EQIC could not be computed because the adjusted variance function is non-positive or non-finite"
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
    geer_criterion_unavailable("EQIC is not available for this family")
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
    geer_criterion_unavailable(
      "EQIC could not be computed because the adjusted deviance is invalid"
    )
  }
  if (identical(family_name, "binomial")) {
    phi <- 1
  } else if (geer_phi_is_fixed(object)) {
    phi <- object$phi
  } else {
    phi <- deviance / object$obs_no
  }
  if (!is.finite(phi) || phi <= 0) {
    geer_criterion_unavailable(
      "EQIC could not be computed because the dispersion estimate is not a usable positive value"
    )
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


compute_jang_criteria <- function(beta_covariance, independence_inverse) {
  failed <- list(PT = NA_real_, WR = NA_real_, RMR = NA_real_)
  ## The generalized eigenvalues of the sandwich covariance with respect to the
  ## model-based independence covariance are the eigenvalues of
  ## beta_covariance %*% independence_inverse. That product of two symmetric
  ## matrices is not itself symmetric, so it is rewritten as the similar
  ## symmetric matrix L %*% beta_covariance %*% t(L), where t(L) %*% L is the
  ## Cholesky factorization of independence_inverse. The eigenvalues are then
  ## real by construction, and the factorization doubles as a check that the
  ## independence information matrix is positive definite.
  chol_factor <- tryCatch(chol(independence_inverse), error = function(e) NULL)
  if (is.null(chol_factor)) {
    return(failed)
  }
  symmetric_form <- chol_factor %*% beta_covariance %*% t(chol_factor)
  symmetric_form <- (symmetric_form + t(symmetric_form)) / 2
  eigenvalues <- tryCatch(
    eigen(symmetric_form, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eigenvalues) || any(!is.finite(eigenvalues)) ||
      any(eigenvalues <= 0)) {
    return(failed)
  }
  ratios <- eigenvalues / (1 + eigenvalues)
  list(
    PT = sum(ratios),
    WR = prod(ratios),
    RMR = max(ratios)
  )
}


compute_gee_criteria <- function(object,
                                 cov_type,
                                 digits = NULL,
                                 include_extended = TRUE,
                                 criteria = geer_criteria_columns) {
  ## Only the requested criteria are evaluated; the remaining columns are
  ## filled with NA so that the shape of the returned table never depends on
  ## the selection. This matters because several of the blocks below are
  ## expensive: QICHH refits the model under working independence, GHYC and PAC
  ## loop over clusters, the eigenvalue criteria factorize a p by p matrix, and
  ## the sandwich covariance itself costs one deletion refit per cluster when
  ## cov_type is "jackknife".
  wanted <- function(...) any(c(...) %in% criteria)
  needs_quasi_loglikelihood <- wanted("QIC", "QICu", "QICC")
  needs_covariance <- wanted(
    "QIC", "QICHH", "QICC", "CIC", "RJC", "EQIC", "PT", "WR", "RMR"
  )
  needs_sc_wc <- wanted("GESSC", "GPC", "AGPC", "SGPC")

  quasi_loglikelihood <- if (needs_quasi_loglikelihood) {
    compute_quasi_loglikelihood(object)
  } else {
    NA_real_
  }
  sc_wc_stats <- if (!needs_sc_wc) {
    c(sc = NA_real_, gp = NA_real_)
  } else if (is_geewa_fit(object)) {
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
  ## The C++ helpers return a named list; the names are kept so that the two
  ## quantities are addressed by name rather than by position.
  sc_wc_stats <- unlist(sc_wc_stats)
  naive_covariance <- if (wanted("RJC")) {
    stats::vcov(object, cov_type = "naive")
  } else {
    NULL
  }
  beta_covariance <- if (needs_covariance) {
    stats::vcov(object, cov_type = cov_type)
  } else {
    NULL
  }
  independence_inverse <- if (needs_covariance) {
    compute_independence_naive_inverse(object)
  } else {
    NULL
  }
  p <- length(object$coefficients)
  ## Every criterion that discounts model complexity counts only the
  ## association parameters that were actually estimated, so a supplied
  ## correlation or odds-ratio structure contributes none.
  association_params_no <- compute_n_estimated_association_parameters(object)
  gessc <- if (wanted("GESSC")) {
    sc_wc_stats[["sc"]] / (object$obs_no - p - association_params_no)
  } else {
    NA_real_
  }
  gpc <- sc_wc_stats[["gp"]]
  penalized_gpc <- if (wanted("AGPC", "SGPC")) {
    compute_gaussian_pseudolikelihood_criteria(
      object = object,
      gpc = gpc,
      p = p,
      association_params_no = association_params_no
    )
  } else {
    list(AGPC = NA_real_, SGPC = NA_real_)
  }
  qic_u <- if (wanted("QICu")) 2 * (p - quasi_loglikelihood) else NA_real_
  ## Same quantity as compute_gee_cic(), which exists so that step_p(), add1()
  ## and drop1() can obtain CIC on its own without building the whole table.
  ## It is computed inline here because both matrices are already available:
  ## calling the standalone helper would re-derive the independence
  ## information matrix and refit the sandwich covariance, which for
  ## cov_type = "jackknife" means one deletion refit per cluster. Keep the two
  ## expressions in step.
  cic <- if (wanted("QIC", "CIC", "QICC")) {
    sum(independence_inverse * beta_covariance)
  } else {
    NA_real_
  }
  qic <- if (wanted("QIC", "QICC")) 2 * (cic - quasi_loglikelihood) else NA_real_
  qicc <- if (wanted("QICC")) {
    compute_qicc(
      qic = qic,
      p = p,
      association_params_no = association_params_no,
      clusters_no = object$clusters_no
    )
  } else {
    NA_real_
  }
  if (!wanted("CIC")) {
    cic <- NA_real_
  }
  if (!wanted("QIC")) {
    qic <- NA_real_
  }
  rjc <- if (!wanted("RJC")) NA_real_ else compute_criterion_or_na({
    q_matrix <- tryCatch(
      solve(naive_covariance, beta_covariance),
      error = function(e) {
        geer_criterion_unavailable(
          "RJC could not be computed because the naive covariance matrix is singular or invalid"
        )
      }
    )
    ## sum(Q * t(Q)) is the trace of Q squared, as in Hin, Carey and Wang
    ## (2007), and not the Frobenius norm, which would be sum(Q^2). The two
    ## coincide only for symmetric Q, and Q is not symmetric here.
    rjc_trace <- sum(diag(q_matrix)) / p
    rjc_squared_trace <- sum(q_matrix * t(q_matrix)) / p
    sqrt((1 - rjc_trace)^2 + (1 - rjc_squared_trace)^2)
  })
  if (isTRUE(include_extended)) {
    eigenvalue_criteria <- if (wanted("PT", "WR", "RMR")) {
      compute_jang_criteria(
        beta_covariance = beta_covariance,
        independence_inverse = independence_inverse
      )
    } else {
      list(PT = NA_real_, WR = NA_real_, RMR = NA_real_)
    }
    qichh <- if (!wanted("QICHH")) NA_real_ else compute_criterion_or_na({
      independence_quantities <- compute_independence_gee_quantities(object)
      qichh_penalty <- sum(
        independence_quantities$naive_inverse * beta_covariance
      )
      2 * (qichh_penalty - independence_quantities$quasi_loglikelihood)
    })
    eqic <- if (!wanted("EQIC")) {
      NA_real_
    } else {
      compute_criterion_or_na(
        compute_eqic(object, beta_covariance = beta_covariance)
      )
    }
    covariance_match <- if (wanted("GHYC", "PAC")) {
      compute_ghyc_pac(object)
    } else {
      list(GHYC = NA_real_, PAC = NA_real_)
    }
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
      PT = eigenvalue_criteria$PT,
      WR = eigenvalue_criteria$WR,
      RMR = eigenvalue_criteria$RMR,
      Parameters = p
    )
    numeric_cols <- geer_criteria_columns
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
    numeric_cols <- geer_criteria_columns_basic
  }
  if (!is.null(digits)) {
    digits <- check_nonnegative_integerish(digits, "digits")
    ans[, numeric_cols] <- lapply(ans[, numeric_cols, drop = FALSE], round, digits = digits)
    ans[, "Parameters"] <- as.integer(ans[, "Parameters"])
  }
  ans
}


## Standalone CIC for callers that need this criterion alone. The same
## expression appears inline in compute_gee_criteria(), where both matrices
## have already been computed; keep the two in step.
compute_gee_cic <- function(object, cov_type) {
  independence_inverse <- compute_independence_naive_inverse(object)
  beta_covariance <- stats::vcov(object, cov_type = cov_type)
  sum(independence_inverse * beta_covariance)
}
