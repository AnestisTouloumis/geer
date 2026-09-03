jackknife_fit_origin <- function(object) {
  if (!is.null(object$fit_function)) {
    return(object$fit_function)
  }
  call_head <- as.character(object$call[[1L]])
  call_head <- call_head[[length(call_head)]]
  if (call_head %in% c("geewa", "geewa_binary")) {
    return(call_head)
  }
  stop(
    "jackknife covariance could not determine whether the model was fitted by 'geewa' or 'geewa_binary'",
    call. = FALSE
  )
}


jackknife_pair_subset <- function(alpha, full_max, subset_max) {
  if (subset_max < 2L) {
    return(numeric(0))
  }
  if (subset_max == full_max) {
    return(as.numeric(alpha))
  }
  full_pairs <- utils::combn(seq_len(full_max), 2L)
  subset_pairs <- utils::combn(seq_len(subset_max), 2L)
  full_keys <- paste(full_pairs[1L, ], full_pairs[2L, ], sep = ":")
  subset_keys <- paste(subset_pairs[1L, ], subset_pairs[2L, ], sep = ":")
  indices <- match(subset_keys, full_keys)
  if (anyNA(indices) || length(alpha) != length(full_keys)) {
    stop("jackknife covariance failed to map the fitted association parameters", call. = FALSE)
  }
  as.numeric(alpha[indices])
}


jackknife_cc_alpha <- function(object, repeated) {
  structure <- object$association_structure
  alpha <- as.numeric(object$alpha)
  subset_max <- max(repeated)
  full_max <- max(object$repeated)

  if (structure %in% c("unstructured", "fixed")) {
    return(jackknife_pair_subset(alpha, full_max, subset_max))
  }
  if (structure %in% c("toeplitz", "m-dependent")) {
    keep <- min(length(alpha), max(subset_max - 1L, 0L))
    return(if (keep > 0L) alpha[seq_len(keep)] else numeric(0))
  }
  alpha
}


jackknife_or_alpha <- function(object, repeated) {
  ## fit_bingee_or() always indexes alpha_vector by pair position, so it must
  ## have length choose(max(repeated), 2) for every odds-ratio structure. An
  ## independence fit stores the scalar 1 rather than a pair vector, so it is
  ## expanded here; passing the scalar through would index out of bounds.
  if (identical(object$association_structure, "independence")) {
    return(rep.int(1, choose(max(repeated), 2L)))
  }
  jackknife_pair_subset(
    alpha = object$alpha,
    full_max = max(object$repeated),
    subset_max = max(repeated)
  )
}


jackknife_last_criterion <- function(fit) {
  index <- ncol(fit$beta_mat) - 1L
  if (index < 1L || length(fit$criterion) < index) {
    return(Inf)
  }
  as.numeric(fit$criterion[[index]])
}


jackknife_check_convergence <- function(fit, tolerance, cluster_label, stage = NULL) {
  criterion <- jackknife_last_criterion(fit)
  if (!is.finite(criterion) || criterion > tolerance) {
    stage_text <- if (is.null(stage)) "" else paste0(" during ", stage)
    stop(
      sprintf(
        "jackknife covariance failed for cluster '%s'%s: leave-one-cluster fit did not converge",
        cluster_label,
        stage_text
      ),
      call. = FALSE
    )
  }
  invisible(fit)
}


jackknife_refit_cc <- function(object, keep, cluster_label) {
  y <- object$y[keep]
  x <- object$x[keep, , drop = FALSE]
  id <- as.numeric(factor(object$id[keep]))
  repeated <- object$repeated[keep]
  weights <- object$prior.weights[keep]
  offset <- object$offset[keep]
  control <- object$control
  tolerance <- control$tolerance
  method <- object$method
  alpha <- jackknife_cc_alpha(object, repeated)
  alpha_fixed <- 1L
  mdependence <- if (identical(object$association_structure, "m-dependent")) {
    length(alpha)
  } else {
    1L
  }
  use_p <- if (!is.null(object$use_p)) isTRUE(object$use_p) else TRUE
  use_params <- if (use_p) ncol(x) else 0L
  phi_fixed <- if (!is.null(object$phi_fixed)) isTRUE(object$phi_fixed) else FALSE
  phi_value <- if (phi_fixed) object$phi else 1
  beta_start <- as.numeric(object$coefficients)

  fit_once <- function(beta, fit_method, maxiter = control$maxiter,
                       step_maxiter = control$step_maxiter,
                       step_multiplier = control$step_multiplier,
                       alpha_value = alpha,
                       fixed_alpha = alpha_fixed,
                       corstr = object$association_structure,
                       phi = phi_value,
                       fixed_phi = as.integer(phi_fixed)) {
    fit_geesolver_cc(
      y, x, id, repeated, weights,
      object$family$link, object$family$family, beta, offset,
      maxiter, tolerance, step_maxiter, step_multiplier,
      control$jeffreys_power, fit_method, use_params,
      alpha_value, fixed_alpha, corstr, mdependence,
      phi, fixed_phi
    )
  }

  if (method %in% geer_bcgee_methods) {
    first <- fit_once(beta_start, "gee")
    jackknife_check_convergence(first, tolerance, cluster_label, "the preliminary GEE fit")
    final <- fit_once(
      as.numeric(first$beta_hat),
      sub("bcgee", "brgee", method),
      maxiter = 1L,
      step_maxiter = 1L,
      step_multiplier = 1L,
      alpha_value = alpha,
      fixed_alpha = 1L,
      phi = first$phi,
      fixed_phi = 1L
    )
  } else if (identical(method, "hpgee-jeffreys")) {
    first <- fit_once(
      beta_start,
      "pgee-jeffreys",
      alpha_value = 0,
      fixed_alpha = 1L,
      corstr = "independence"
    )
    jackknife_check_convergence(first, tolerance, cluster_label, "the preliminary independence PGEE fit")
    final <- fit_once(
      as.numeric(first$beta_hat),
      "gee",
      maxiter = 1L,
      step_maxiter = 1L,
      step_multiplier = 1L,
      alpha_value = alpha,
      fixed_alpha = 1L
    )
  } else if (identical(method, "opgee-jeffreys")) {
    first <- fit_once(
      beta_start,
      "pgee-jeffreys",
      alpha_value = 0,
      fixed_alpha = 1L,
      corstr = "independence"
    )
    jackknife_check_convergence(first, tolerance, cluster_label, "the preliminary independence PGEE fit")
    final <- fit_once(
      as.numeric(first$beta_hat),
      "pgee-jeffreys",
      maxiter = 1L,
      step_maxiter = 1L,
      step_multiplier = 1L,
      alpha_value = alpha,
      fixed_alpha = 1L
    )
  } else {
    final <- fit_once(beta_start, method)
    jackknife_check_convergence(final, tolerance, cluster_label)
  }

  beta <- as.numeric(final$beta_hat)
  if (length(beta) != ncol(x) || any(!is.finite(beta))) {
    stop(
      sprintf(
        "jackknife covariance failed for cluster '%s': invalid leave-one-cluster coefficient estimate",
        cluster_label
      ),
      call. = FALSE
    )
  }
  beta
}


jackknife_refit_or <- function(object, keep, cluster_label) {
  y <- object$y[keep]
  x <- object$x[keep, , drop = FALSE]
  id <- as.numeric(factor(object$id[keep]))
  repeated <- object$repeated[keep]
  weights <- object$prior.weights[keep]
  offset <- object$offset[keep]
  control <- object$control
  tolerance <- control$tolerance
  method <- object$method
  alpha <- jackknife_or_alpha(object, repeated)
  alpha_independence <- rep.int(1, choose(max(repeated), 2L))
  beta_start <- as.numeric(object$coefficients)

  fit_once <- function(beta, fit_method, alpha_value,
                       maxiter = control$maxiter,
                       step_maxiter = control$step_maxiter,
                       step_multiplier = control$step_multiplier) {
    fit_bingee_or(
      y, x, id, repeated, weights, object$family$link,
      beta, offset, maxiter, tolerance,
      step_maxiter, step_multiplier,
      control$jeffreys_power, fit_method, alpha_value
    )
  }

  if (method %in% geer_bcgee_methods) {
    first <- fit_once(beta_start, "gee", alpha)
    jackknife_check_convergence(first, tolerance, cluster_label, "the preliminary GEE fit")
    final <- fit_once(
      as.numeric(first$beta_hat),
      sub("bcgee", "brgee", method),
      alpha,
      maxiter = 1L,
      step_maxiter = 1L,
      step_multiplier = 1L
    )
  } else if (identical(method, "hpgee-jeffreys")) {
    first <- fit_once(beta_start, "pgee-jeffreys", alpha_independence)
    jackknife_check_convergence(first, tolerance, cluster_label, "the preliminary independence PGEE fit")
    final <- fit_once(
      as.numeric(first$beta_hat),
      "gee",
      alpha,
      maxiter = 1L,
      step_maxiter = 1L,
      step_multiplier = 1L
    )
  } else if (identical(method, "opgee-jeffreys")) {
    first <- fit_once(beta_start, "pgee-jeffreys", alpha_independence)
    jackknife_check_convergence(first, tolerance, cluster_label, "the preliminary independence PGEE fit")
    final <- fit_once(
      as.numeric(first$beta_hat),
      "pgee-jeffreys",
      alpha,
      maxiter = 1L,
      step_maxiter = 1L,
      step_multiplier = 1L
    )
  } else {
    final <- fit_once(beta_start, method, alpha)
    jackknife_check_convergence(final, tolerance, cluster_label)
  }

  beta <- as.numeric(final$beta_hat)
  if (length(beta) != ncol(x) || any(!is.finite(beta))) {
    stop(
      sprintf(
        "jackknife covariance failed for cluster '%s': invalid leave-one-cluster coefficient estimate",
        cluster_label
      ),
      call. = FALSE
    )
  }
  beta
}


compute_jackknife_delete_estimates <- function(object) {
  object <- check_geer_object(object)
  cluster_indices <- split(seq_along(object$id), object$id)
  cluster_labels <- names(cluster_indices)
  n_clusters <- length(cluster_indices)
  if (n_clusters < 2L) {
    stop("jackknife covariance requires at least two clusters", call. = FALSE)
  }
  p <- length(object$coefficients)
  estimates <- matrix(
    NA_real_,
    nrow = n_clusters,
    ncol = p,
    dimnames = list(cluster_labels, names(object$coefficients))
  )
  fit_origin <- jackknife_fit_origin(object)

  for (i in seq_along(cluster_indices)) {
    keep <- rep.int(TRUE, length(object$id))
    keep[cluster_indices[[i]]] <- FALSE
    estimates[i, ] <- if (identical(fit_origin, "geewa_binary")) {
      jackknife_refit_or(object, keep, cluster_labels[[i]])
    } else {
      jackknife_refit_cc(object, keep, cluster_labels[[i]])
    }
  }
  estimates
}


compute_jackknife_covariance <- function(object) {
  object <- check_geer_object(object)
  delete_estimates <- compute_jackknife_delete_estimates(object)
  clusters_no <- nrow(delete_estimates)
  centered <- sweep(delete_estimates, 2L, colMeans(delete_estimates), `-`)
  ## The (K - 1) / K factor is the usual finite-sample correction of the
  ## Quenouille-Tukey jackknife, obtained from the sample variance of the
  ## pseudo-values K * beta_hat - (K - 1) * beta_hat_(i).
  ans <- ((clusters_no - 1L) / clusters_no) * crossprod(centered)
  ans <- 0.5 * (ans + t(ans))
  coefficient_names <- names(object$coefficients)
  dimnames(ans) <- list(coefficient_names, coefficient_names)
  ans
}
