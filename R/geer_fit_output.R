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
                              formula,
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
  fit$df.residual <- nrow(model_matrix) - fit$rank
  fit$y <- y
  fit$x <- model_matrix
  fit$na.action <- attr(model_frame, "na.action")
  fit$id <- as.numeric(id)
  fit$repeated <- as.numeric(repeated)
  fit$call <- call
  fit$formula <- formula
  fit$terms <- model_terms
  fit$data <- data
  fit$offset <- geesolver_fit$offset
  fit$control <- control
  fit$method <- method
  fit$contrasts <- attr(model_matrix, "contrasts")
  fit$xlevels <- stats::.getXlevels(attr(model_frame, "terms"), model_frame)
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


finalize_geer_fit <- function(fit,
                              geesolver_fit,
                              tolerance,
                              method,
                              family,
                              association_structure,
                              repeated,
                              fit_function) {
  fit$fit_function <- fit_function
  fit$converged <- (geesolver_fit$criterion[fit$iter] <= tolerance)
  if (method %in% geer_onestep_methods) {
    fit$converged <- TRUE
  }
  independence_value <- if (identical(fit_function, "geewa_binary")) 1 else 0
  fit$alpha <- if (identical(association_structure, "independence")) {
    independence_value
  } else {
    as.numeric(geesolver_fit$alpha)
  }
  if (association_structure %in% c("unstructured", "fixed")) {
    pairs_matrix <- utils::combn(max(repeated), 2)
    names(fit$alpha) <-
      paste0("alpha_", pairs_matrix[1, ], ".", pairs_matrix[2, ])
  }
  if (association_structure %in% c("toeplitz", "m-dependent")) {
    names(fit$alpha) <- paste0("alpha_lag", seq_along(fit$alpha))
  }
  if (!fit$converged) {
    warning(fit_function, ": algorithm did not converge", call. = FALSE)
  }
  ## Matches the boundary threshold and wording used by stats::glm.fit().
  eps <- 10 * .Machine$double.eps
  if (identical(family$family, "binomial")) {
    if (any(fit$fitted.values > 1 - eps) || any(fit$fitted.values < eps)) {
      warning(
        fit_function, ": fitted probabilities numerically 0 or 1 occurred",
        call. = FALSE
      )
    }
  }
  if (identical(family$family, "poisson")) {
    if (any(fit$fitted.values < eps)) {
      warning(
        fit_function, ": fitted rates numerically 0 occurred",
        call. = FALSE
      )
    }
  }
  fit
}
