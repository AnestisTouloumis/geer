extract_jackknife_fit_function <- function(object) {
  fit_function <- get_geer_fit_function(object)
  if (!is.na(fit_function)) {
    return(fit_function)
  }
  stop(
    "jackknife covariance could not determine whether the model was fitted by 'geewa' or 'geewa_binary'",
    call. = FALSE
  )
}


select_jackknife_pair_subset <- function(alpha, full_max, subset_max) {
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


compute_jackknife_alpha_cc <- function(object, repeated) {
  structure <- object$association_structure
  alpha <- as.numeric(object$alpha)
  subset_max <- max(repeated)
  full_max <- max(object$repeated)

  if (structure %in% c("unstructured", "fixed")) {
    return(select_jackknife_pair_subset(alpha, full_max, subset_max))
  }
  if (structure %in% c("toeplitz", "m-dependent")) {
    keep <- min(length(alpha), max(subset_max - 1L, 0L))
    return(if (keep > 0L) alpha[seq_len(keep)] else numeric(0))
  }
  alpha
}


compute_jackknife_alpha_or <- function(object, repeated) {
  ## fit_geesolver_or() always indexes alpha_vector by pair position, so it must
  ## have length choose(max(repeated), 2) for every odds-ratio structure. An
  ## independence fit stores the scalar 1 rather than a pair vector, so it is
  ## expanded here; passing the scalar through would index out of bounds.
  if (identical(object$association_structure, "independence")) {
    return(rep.int(1, choose(max(repeated), 2L)))
  }
  select_jackknife_pair_subset(
    alpha = object$alpha,
    full_max = max(object$repeated),
    subset_max = max(repeated)
  )
}


extract_jackknife_last_criterion <- function(fit) {
  index <- ncol(fit$beta_mat) - 1L
  if (index < 1L || length(fit$criterion) < index) {
    return(Inf)
  }
  as.numeric(fit$criterion[[index]])
}


check_jackknife_convergence <- function(fit, tolerance, cluster_label, stage = NULL) {
  criterion <- extract_jackknife_last_criterion(fit)
  if (!is.finite(criterion) || criterion > tolerance) {
    stage_text <- if (is.null(stage)) "" else paste0(" during ", stage)
    reason <- geer_solver_failure(fit)
    reason_text <- if (nzchar(reason)) paste0(" (", reason, ")") else ""
    stop(
      sprintf(
        "jackknife covariance failed for cluster '%s'%s: leave-one-cluster fit did not converge%s",
        cluster_label,
        stage_text,
        reason_text
      ),
      call. = FALSE
    )
  }
  invisible(fit)
}


## One-step methods skip the convergence check of the final pass, so a
## numerical failure reported by the solver is checked separately.
check_jackknife_solver_failure <- function(fit, cluster_label) {
  reason <- geer_solver_failure(fit)
  if (nzchar(reason)) {
    stop(
      sprintf(
        "jackknife covariance failed for cluster '%s': %s",
        cluster_label,
        reason
      ),
      call. = FALSE
    )
  }
  invisible(fit)
}


refit_jackknife_cc <- function(object, keep, cluster_label) {
  y <- object$y[keep]
  x <- object$x[keep, , drop = FALSE]
  id <- as.numeric(factor(object$id[keep]))
  repeated <- object$repeated[keep]
  weights <- object$prior.weights[keep]
  offset <- object$offset[keep]
  control <- object$control
  tolerance <- control$tolerance
  method <- object$method
  alpha <- compute_jackknife_alpha_cc(object, repeated)
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

  fit_pass <- function(beta, pass, previous) {
    iterations <- if (pass$one_step) {
      list(1L, 1L, 1)
    } else {
      list(control$maxiter, control$step_maxiter, control$step_multiplier)
    }
    if (pass$carry_nuisance) {
      phi_pass <- previous$phi
      phi_fixed_pass <- 1L
    } else {
      phi_pass <- phi_value
      phi_fixed_pass <- as.integer(phi_fixed)
    }
    fit_geesolver_cc(
      y, x, id, repeated, weights,
      object$family$link, object$family$family, as.numeric(beta), offset,
      iterations[[1L]], tolerance, iterations[[2L]], iterations[[3L]],
      control$jeffreys_power, pass$method, use_params,
      if (pass$independence) 0 else alpha, alpha_fixed,
      if (pass$independence) "independence" else object$association_structure,
      mdependence, phi_pass, phi_fixed_pass,
      as.integer(isTRUE(pass$hold_nuisance))
    )
  }
  final <- run_geer_estimation_passes(
    fit_pass = fit_pass,
    method = method,
    beta_start = beta_start,
    check_first = function(fit, method) {
      stage <- if (method %in% geer_bcgee_methods) {
        "the preliminary GEE fit"
      } else {
        "the preliminary independence PGEE fit"
      }
      check_jackknife_convergence(fit, tolerance, cluster_label, stage)
    }
  )
  if (!(method %in% geer_onestep_methods)) {
    check_jackknife_convergence(final, tolerance, cluster_label)
  }
  check_jackknife_solver_failure(final, cluster_label)

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


refit_jackknife_or <- function(object, keep, cluster_label) {
  y <- object$y[keep]
  x <- object$x[keep, , drop = FALSE]
  id <- as.numeric(factor(object$id[keep]))
  repeated <- object$repeated[keep]
  weights <- object$prior.weights[keep]
  offset <- object$offset[keep]
  control <- object$control
  tolerance <- control$tolerance
  method <- object$method
  alpha <- compute_jackknife_alpha_or(object, repeated)
  alpha_independence <- rep.int(1, choose(max(repeated), 2L))
  beta_start <- as.numeric(object$coefficients)

  fit_pass <- function(beta, pass, previous) {
    iterations <- if (pass$one_step) {
      list(1L, 1L, 1)
    } else {
      list(control$maxiter, control$step_maxiter, control$step_multiplier)
    }
    fit_geesolver_or(
      y, x, id, repeated, weights, object$family$link,
      as.numeric(beta), offset,
      iterations[[1L]], tolerance, iterations[[2L]], iterations[[3L]],
      control$jeffreys_power, pass$method,
      if (pass$independence) alpha_independence else alpha
    )
  }
  final <- run_geer_estimation_passes(
    fit_pass = fit_pass,
    method = method,
    beta_start = beta_start,
    check_first = function(fit, method) {
      stage <- if (method %in% geer_bcgee_methods) {
        "the preliminary GEE fit"
      } else {
        "the preliminary independence PGEE fit"
      }
      check_jackknife_convergence(fit, tolerance, cluster_label, stage)
    }
  )
  if (!(method %in% geer_onestep_methods)) {
    check_jackknife_convergence(final, tolerance, cluster_label)
  }
  check_jackknife_solver_failure(final, cluster_label)

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
  fit_origin <- extract_jackknife_fit_function(object)

  for (i in seq_along(cluster_indices)) {
    keep <- rep.int(TRUE, length(object$id))
    keep[cluster_indices[[i]]] <- FALSE
    estimates[i, ] <- if (identical(fit_origin, "geewa_binary")) {
      refit_jackknife_or(object, keep, cluster_labels[[i]])
    } else {
      refit_jackknife_cc(object, keep, cluster_labels[[i]])
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


## The jackknife covariance needs one full refit per cluster, and summary(),
## confint(), tidy(), predict(se.fit = TRUE) and every test that uses it would
## otherwise repeat that work. The result is kept in the 'cache' environment
## created with the fit. It is reused only while the stored coefficients are
## identical to the current ones, so a modified copy of the fit (for example
## through set_coef()) is recomputed rather than served stale.
get_cached_jackknife_covariance <- function(object) {
  cache <- object$cache
  if (is.environment(cache)) {
    stored <- cache$jackknife
    if (!is.null(stored) && identical(stored$coefficients, object$coefficients)) {
      return(stored$covariance)
    }
  }
  covariance <- compute_jackknife_covariance(object)
  if (is.environment(cache)) {
    assign(
      "jackknife",
      list(coefficients = object$coefficients, covariance = covariance),
      envir = cache
    )
  }
  covariance
}
