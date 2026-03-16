compute_quasi_loglikelihood <- function(object) {
  y  <- object$y
  mu <- object$fitted.values
  wt <- object$prior.weights
  mdis <- object$family$family
  phi <- object$phi
  ans <-
    switch(mdis,
           gaussian = -sum(wt * (y - mu)^2)/2,
           binomial = sum(wt * (y * qlogis(mu) + log(1 - mu))),
           poisson  = sum(wt * (y * log(mu) - mu)),
           Gamma = -sum(wt * (y/mu + log(mu))),
           inverse.gaussian = -sum(wt * (mu - 0.5*y)/mu^2),
           stop("Error: distribution not recognized")
    )
  ans/phi
}


## compute all criteria
compute_criteria <- function(object, cov_type, digits) {
  quasi_loglikelihood <- compute_quasi_loglikelihood(object)
  if (as.character(object$call[1]) == "geewa") {
    sc_wc_stats <-
      get_gee_criteria_sc_cw(object$y,
                             object$id,
                             object$repeated,
                             object$family$family,
                             object$fitted.values,
                             object$association_structure,
                             object$alpha,
                             object$phi,
                             object$prior.weights)
  } else {
    if (length(object$alpha) == 1) {
      object$alpha <- rep(object$alpha, choose(max(object$repeated), 2))
    }
    sc_wc_stats <-
      get_gee_criteria_sc_cw_or(object$y,
                                object$id,
                                object$repeated,
                                object$fitted.values,
                                object$alpha,
                                object$prior.weights)
  }
  if (object$association == "independence") {
    association_params_no <- 0
  } else if (object$association == "exchangeable" || object$association == "ar1") {
    association_params_no <- 1
  } else {
    association_params_no <- length(object$alpha)
  }
  naive_covariance <- vcov(object,
                           cov_type = "naive")
  beta_covariance <- vcov(object,
                          cov_type = cov_type)
  sc_wc_stats <- unlist(sc_wc_stats)
  independence_naive_covariance_inverse <-
    get_naive_matrix_inverse_independence(object$y,
                                          object$x,
                                          object$id,
                                          object$family$link,
                                          object$family$family,
                                          object$fitted.values,
                                          object$linear.predictors,
                                          object$phi)
  p <- length(object$coeff)
  gessc <- sc_wc_stats[[1]]/(object$obs_no - p - association_params_no)
  gpc <- sc_wc_stats[[2]]
  qic_u <- 2 * (p - quasi_loglikelihood)
  cic <- sum(independence_naive_covariance_inverse * beta_covariance)
  qic <- 2 * (cic - quasi_loglikelihood)
  q_matrix <- solve(naive_covariance, beta_covariance)
  c1 <- sum(diag(q_matrix)) / p
  c2 <- sum(q_matrix * t(q_matrix)) / p
  rjc <- sqrt(((1 - c1) ^ 2) + ((1 - c2) ^ 2))
  ans <- data.frame(QIC = qic,
                    CIC = cic,
                    RJC = rjc,
                    QICu = qic_u,
                    GESSC = gessc,
                    GPC = gpc,
                    Parameters = p)
  ans[, 1:6] <- apply(ans[, 1:6], 2, round, digits = digits)
  ans[, 7] <- round(ans[, 7], digits = 0)
  ans
}


## extract only cic criterion
extract_cic <- function(object, cov_type) {
  cov_mat_independence_inverse <-
    get_naive_matrix_inverse_independence(object$y,
                                          object$x,
                                          object$id,
                                          object$family$link,
                                          object$family$family,
                                          object$fitted.values,
                                          object$linear.predictors,
                                          object$phi)
  cov_mat <- vcov(object, cov_type = cov_type)
  ans <- sum(cov_mat_independence_inverse * cov_mat)
  ans
}
