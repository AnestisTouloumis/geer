.compute_or_alpha <- function(object) {
  if (length(object$alpha) == 1L) {
    rep(object$alpha, choose(max(as.integer(object$repeated)), 2))
  } else {
    object$alpha
  }
}


.compute_association_params_no <- function(object) {
  if (identical(object$association_structure, "independence")) {
    0L
  } else if (object$association_structure %in% c("exchangeable", "ar1")) {
    1L
  } else {
    length(object$alpha)
  }
}


.get_independence_naive_inverse <- function(object) {
  get_naive_matrix_inverse_independence(
    object$y,
    object$x,
    object$id,
    object$family$link,
    object$family$family,
    object$fitted.values,
    object$linear.predictors,
    object$phi
  )
}


compute_quasi_loglikelihood <- function(object) {
  y <- object$y
  mu <- object$fitted.values
  wt <- object$prior.weights
  mdis <- object$family$family
  phi <- object$phi
  eps <- sqrt(.Machine$double.eps)
  ans <- switch(
    mdis,
    gaussian = -sum(wt * (y - mu)^2) / 2,
    binomial = {
      mu_safe <- pmin(pmax(mu, eps), 1 - eps)
      sum(wt * (y * stats::qlogis(mu_safe) + log1p(-mu_safe)))
    },
    poisson = {
      mu_safe <- pmax(mu, eps)
      sum(wt * (y * log(mu_safe) - mu_safe))
    },
    Gamma = {
      mu_safe <- pmax(mu, eps)
      -sum(wt * (y / mu_safe + log(mu_safe)))
    },
    inverse.gaussian = {
      mu_safe <- pmax(mu, eps)
      -sum(wt * (mu_safe - 0.5 * y) / mu_safe^2)
    },
    stop("distribution not recognized", call. = FALSE)
  )
  ans / phi
}


## compute all criteria
compute_criteria <- function(object, cov_type, digits = NULL) {
  quasi_loglikelihood <- compute_quasi_loglikelihood(object)
  sc_wc_stats <- if (identical(as.character(object$call[[1L]]), "geewa")) {
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
      .compute_or_alpha(object),
      object$prior.weights
    )
  }
  sc_wc_stats <- unlist(sc_wc_stats, use.names = FALSE)
  naive_covariance <- vcov(object, cov_type = "naive")
  beta_covariance <- vcov(object, cov_type = cov_type)
  independence_inverse <- .get_independence_naive_inverse(object)
  p <- length(object$coefficients)
  association_params_no <- .compute_association_params_no(object)
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
  c1 <- sum(diag(q_matrix)) / p
  c2 <- sum(q_matrix * t(q_matrix)) / p
  rjc <- sqrt((1 - c1)^2 + (1 - c2)^2)
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
    ans[, 1:6] <- apply(ans[, 1:6, drop = FALSE], 2, round, digits = digits)
    ans[, 7] <- round(ans[, 7], digits = 0)
  }
  ans
}


## extract only CIC criterion
extract_cic <- function(object, cov_type) {
  independence_inverse <- .get_independence_naive_inverse(object)
  beta_covariance <- vcov(object, cov_type = cov_type)
  sum(independence_inverse * beta_covariance)
}
