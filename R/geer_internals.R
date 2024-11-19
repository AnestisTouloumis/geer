format_perc <- function(probs, digits) {
  paste(
    format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
    "%"
  )
}


compute_quasi_loglikelihood <- function(object) {
  mu <- object$fitted.values
  y  <- object$y
  marginal_distribution <- object$family$family
  ans <-
    switch(marginal_distribution,
           gaussian = sum(((y - mu)^2)/-2),
           binomial = sum(y * log(mu) + (1 - y) * log(1 - mu)),
           poisson  = sum((y * log(mu)) - mu),
           Gamma    = sum(-y/(mu - log(mu))),
           stop("Error: distribution not recognized")
    )/object$phi
  ans
}

compute_criteria <- function(object, cov_type, digits) {
  quasi_loglikelihood <- compute_quasi_loglikelihood(object)
  if (as.character(object$call[1]) == "geewa") {
    quasi_loglikelihood_independence <-
      compute_quasi_loglikelihood(
        update(object,
               correlation_structure = "independence")
      )
    sc_wc_stats <-
      get_gee_criteria_sc_cw(object$y,
                             object$id,
                             object$repeated,
                             object$family$family,
                             object$fitted.values,
                             object$association_structure,
                             object$alpha,
                             object$phi)
  } else {
    quasi_loglikelihood_independence <-
      compute_quasi_loglikelihood(
        update(object,
               or_structure = "independence")
      )
    sc_wc_stats <-
      get_gee_criteria_sc_cw_or(object$y,
                                object$id,
                                object$repeated,
                                object$fitted.values,
                                object$alpha)
  }
  if (object$association == "indepedence") {
    association_params_no <- 0
  } else if (object$association == "exchangeable") {
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
                                          object$model_matrix,
                                          object$id,
                                          object$family$link,
                                          object$family$family,
                                          object$fitted.values,
                                          object$linear.predictors,
                                          object$phi)
  p <- length(object$coeff)
  gess <- sc_wc_stats[[1]]/(object$obs_no - p - association_params_no)
  cwc <- sc_wc_stats[[2]]
  qic_u <- 2 * (p - quasi_loglikelihood)
  cic <- sum(diag(independence_naive_covariance_inverse %*% beta_covariance))
  qic <- 2 * (cic - quasi_loglikelihood)
  qic_hh <- 2 * (cic - quasi_loglikelihood_independence)
  q_matrix <- solve(naive_covariance) %*% beta_covariance
  c1 <- sum(diag(q_matrix)) / p
  c2 <- sum(q_matrix * t(q_matrix)) / p
  rjc <- sqrt(((1 - c1) ^ 2) + ((1 - c2) ^ 2))
  ans <- data.frame(QIC = qic,
                    QICHH = qic_hh,
                    CIC = cic,
                    RJC = rjc,
                    QICU = qic_u,
                    GESS = gess,
                    CWC = cwc,
                    Parameters = p)
  ans[, 1:7] <- apply(ans[, 1:7], 2, round, digits = digits)
  ans[, 8] <- round(ans[, 8], digits = 0)
  ans
}
