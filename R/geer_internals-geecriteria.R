compute_n_association_parameters <- function(object) {
  if (identical(object$association_structure, "independence")) {
    0L
  } else if (object$association_structure %in% c("exchangeable", "ar1")) {
    1L
  } else {
    length(object$alpha)
  }
}


compute_independence_naive_inverse <- function(object) {
  get_naive_matrix_inverse_independence(
    object$x,
    object$id,
    object$family$link,
    object$family$family,
    object$fitted.values,
    object$linear.predictors,
    object$phi,
    object$prior.weights
  )
}


compute_quasi_loglikelihood <- function(object) {
  y <- object$y
  mu <- object$fitted.values
  weights <- object$prior.weights
  family_name <- object$family$family
  phi <- object$phi
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


compute_gee_criteria <- function(object, cov_type, digits = NULL) {
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
  gessc <- sc_wc_stats[[1L]] / (object$obs_no - p - association_params_no)
  gpc <- sc_wc_stats[[2L]]
  qic_u <- 2 * (p - quasi_loglikelihood)
  cic <- sum(independence_inverse * beta_covariance)
  qic <- 2 * (cic - quasi_loglikelihood)
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
  ans <- data.frame(
    QIC = qic,
    CIC = cic,
    RJC = rjc,
    QICu = qic_u,
    GESSC = gessc,
    GPC = gpc,
    Parameters = p
  )
  if (!is.null(digits)) {
    digits <- check_nonnegative_integerish(digits, "digits")
    numeric_cols <- c("QIC", "CIC", "RJC", "QICu", "GESSC", "GPC")
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
