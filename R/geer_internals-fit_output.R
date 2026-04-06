build_geer_output <- function(geesolver_fit,
                              xnames,
                              qr_model_matrix,
                              family,
                              weights,
                              y,
                              model_matrix,
                              model_frame,
                              id,
                              repeated,
                              call,
                              data,
                              model_terms,
                              control,
                              method,
                              association_structure) {
  fit <- list()
  fit$coefficients <- as.numeric(geesolver_fit$beta_hat)
  names(fit$coefficients) <- xnames
  fit$residuals <- as.numeric(geesolver_fit$residuals)
  fit$fitted.values <- as.numeric(geesolver_fit$fitted)
  fit$qr <- qr_model_matrix
  fit$rank <- qr_model_matrix$rank
  fit$family <- family
  fit$linear.predictors <- as.numeric(geesolver_fit$eta)
  fit$iter <- ncol(geesolver_fit$beta_mat) - 1L
  fit$prior.weights <- weights
  fit$df.residual <- nrow(model_matrix) - ncol(model_matrix)
  fit$y <- y
  fit$x <- model_matrix
  fit$na.action <- attr(model_frame, "na.action")
  fit$id <- as.numeric(id)
  fit$repeated <- as.numeric(repeated)
  fit$call <- call
  fit$formula <- fit$call$formula
  fit$terms <- model_terms
  fit$data <- data
  fit$offset <- geesolver_fit$offset
  fit$control <- control
  fit$method <- method
  fit$contrasts <- attr(model_matrix, "contrasts")
  fit$xlevels <- .getXlevels(attr(model_frame, "terms"), model_frame)
  fit$naive_covariance <- geesolver_fit$naive_covariance
  dimnames(fit$naive_covariance) <- list(xnames, xnames)
  fit$robust_covariance <- geesolver_fit$robust_covariance
  dimnames(fit$robust_covariance) <- list(xnames, xnames)
  fit$bias_corrected_covariance <- geesolver_fit$bc_covariance
  dimnames(fit$bias_corrected_covariance) <- list(xnames, xnames)
  fit$association_structure <- association_structure
  fit$alpha <- as.numeric(geesolver_fit$alpha)
  fit$phi <- geesolver_fit$phi
  fit$obs_no <- nrow(model_matrix)
  fit$clusters_no <- length(unique(id))
  cluster_sizes <- vapply(split(repeated, id), length, integer(1))
  fit$min_cluster_size <- min(cluster_sizes)
  fit$max_cluster_size <- max(cluster_sizes)
  fit
}
