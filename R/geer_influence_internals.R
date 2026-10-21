compute_influence_quantities <- function(object) {
  if (is_geewa_fit(object)) {
    get_influence_quantities_cc(
      object$y, object$x, object$id, object$repeated, object$prior.weights,
      object$family$link, object$family$family,
      object$fitted.values, object$linear.predictors,
      object$association_structure, object$alpha, object$phi
    )
  } else {
    get_influence_quantities_or(
      object$y, object$x, object$id, object$repeated, object$prior.weights,
      object$family$link,
      object$fitted.values, object$linear.predictors,
      get_alpha_or(object)
    )
  }
}


get_cached_influence_values <- function(object) {
  get_cached_value(
    object,
    "influence",
    function() compute_influence_quantities(object)
  )
}


get_influence_cluster_order <- function(object) {
  codes <- names(split(seq_along(object$id), object$id))
  list(
    labels = get_cluster_labels(object),
    index = match(codes, as.character(unique(object$id)))
  )
}


warn_influence_na <- function(values, unit) {
  missing <- rownames(values)[!stats::complete.cases(values)]
  if (!length(missing)) {
    return(invisible(NULL))
  }
  label <- if (identical(unit, "clusters")) "cluster" else "observation"
  warning(
    sprintf(
      paste0(
        "deleting %s %s leaves coefficients that cannot be estimated, so ",
        "the corresponding values are NA"
      ),
      if (length(missing) > 1L) paste0(label, "s") else label,
      paste(missing, collapse = ", ")
    ),
    call. = FALSE
  )
}


check_influence_method <- function(object, unit, method, context) {
  if (identical(unit, "observations") && identical(method, "exact")) {
    stop(
      sprintf(
        paste0(
          "%s failed: method = \"exact\" is available only for ",
          "unit = \"clusters\""
        ),
        context
      ),
      call. = FALSE
    )
  }
  if (identical(method, "one-step") && !identical(object$method, "gee")) {
    stop(
      sprintf(
        paste0(
          "%s failed: method = \"one-step\" requires a fit with ",
          "method = \"gee\"; use method = \"exact\" with unit = \"clusters\""
        ),
        context
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}


compute_influence_hatvalues <- function(object, unit) {
  hat_vector <- as.numeric(get_cached_influence_values(object)$hat)
  if (identical(unit, "observations")) {
    names(hat_vector) <- get_observation_labels(object)
    return(hat_vector)
  }
  ans <- vapply(
    split(hat_vector, object$id),
    sum,
    numeric(1)
  )
  names(ans) <- get_cluster_labels(object)
  ans
}


compute_influence_dfbeta <- function(object, unit, method, context) {
  check_influence_method(object, unit, method, context)
  coef_names <- names(object$coefficients)
  quantities <- get_cached_influence_values(object)
  if (identical(unit, "observations")) {
    ans <- quantities$observation_dfbeta
    dimnames(ans) <- list(get_observation_labels(object), coef_names)
    warn_influence_na(ans, unit)
    return(ans)
  }
  cluster_order <- get_influence_cluster_order(object)
  ans <- quantities$cluster_dfbeta[cluster_order$index, , drop = FALSE]
  dimnames(ans) <- list(cluster_order$labels, coef_names)
  if (identical(method, "exact")) {
    singular <- !stats::complete.cases(ans)
    delete_estimates <- if (any(singular)) {
      compute_jackknife_deletes(object, context, skip = singular)
    } else {
      get_cached_jackknife_estimates(object, context)
    }
    ans <- sweep(-delete_estimates, 2L, object$coefficients, `+`)
  }
  warn_influence_na(ans, unit)
  ans
}


compute_cooks_distance <- function(object,
                                   unit,
                                   method,
                                   cov_type,
                                   parm,
                                   context) {
  dfbeta_matrix <- compute_influence_dfbeta(object, unit, method, context)
  vcov_matrix <- stats::vcov(object, cov_type = cov_type)[
    parm, parm, drop = FALSE
  ]
  inverse_matrix <- tryCatch(
    solve(vcov_matrix),
    error = function(e) {
      stop(
        sprintf(
          paste0(
            "%s failed: the %s covariance matrix of the selected ",
            "coefficients is singular"
          ),
          context,
          cov_type
        ),
        call. = FALSE
      )
    }
  )
  selected <- dfbeta_matrix[, parm, drop = FALSE]
  ans <- rowSums((selected %*% inverse_matrix) * selected) / length(parm)
  names(ans) <- rownames(dfbeta_matrix)
  ans
}


compute_influence_dfbetas <- function(object,
                                      unit,
                                      method,
                                      cov_type,
                                      context) {
  dfbeta_matrix <- compute_influence_dfbeta(object, unit, method, context)
  standard_errors <- compute_se_or_stop(
    stats::vcov(object, cov_type = cov_type),
    parm = colnames(dfbeta_matrix),
    context = context
  )
  sweep(dfbeta_matrix, 2L, standard_errors, `/`)
}


compute_local_influence_matrix <- function(sensitivity, parm) {
  inverse_matrix <- tryCatch(
    solve(sensitivity),
    error = function(e) NULL
  )
  nuisance <- setdiff(colnames(sensitivity), parm)
  nuisance_inverse <- if (length(nuisance) && !is.null(inverse_matrix)) {
    tryCatch(
      solve(sensitivity[nuisance, nuisance, drop = FALSE]),
      error = function(e) NULL
    )
  }
  if (is.null(inverse_matrix) ||
      (length(nuisance) && is.null(nuisance_inverse))) {
    stop(
      paste0(
        "local_influence failed: the derivative of the estimating ",
        "equations is singular"
      ),
      call. = FALSE
    )
  }
  if (length(nuisance)) {
    inverse_matrix[nuisance, nuisance] <-
      inverse_matrix[nuisance, nuisance] - nuisance_inverse
  }
  inverse_matrix
}


compute_local_direction <- function(delta_matrix, inverse_matrix) {
  decomposition <- qr(delta_matrix)
  rank <- decomposition$rank
  pivot <- decomposition$pivot[seq_len(rank)]
  r_matrix <- qr.R(decomposition)[seq_len(rank), seq_len(rank), drop = FALSE]
  curvature <- r_matrix %*% inverse_matrix[pivot, pivot, drop = FALSE] %*%
    t(r_matrix)
  eigen_decomposition <- eigen(
    0.5 * (curvature + t(curvature)),
    symmetric = TRUE
  )
  leading <- which.max(abs(eigen_decomposition$values))
  q_matrix <- qr.Q(decomposition)[, seq_len(rank), drop = FALSE]
  abs(drop(q_matrix %*% eigen_decomposition$vectors[, leading]))
}


compute_local_influence <- function(object, perturbation, type, parm) {
  quantities <- get_cached_influence_values(object)
  coef_names <- names(object$coefficients)
  sensitivity <- quantities$sensitivity
  dimnames(sensitivity) <- list(coef_names, coef_names)
  inverse_matrix <- compute_local_influence_matrix(sensitivity, parm)
  delta_matrix <- switch(
    perturbation,
    "cluster-weights" = {
      cluster_order <- get_influence_cluster_order(object)
      ans <- quantities$cluster_delta[cluster_order$index, , drop = FALSE]
      rownames(ans) <- cluster_order$labels
      ans
    },
    "observation-weights" = quantities$observation_delta,
    response = quantities$response_delta
  )
  if (!identical(perturbation, "cluster-weights")) {
    rownames(delta_matrix) <- get_observation_labels(object)
  }
  colnames(delta_matrix) <- coef_names
  ans <- if (identical(type, "total")) {
    2 * abs(rowSums((delta_matrix %*% inverse_matrix) * delta_matrix))
  } else {
    compute_local_direction(delta_matrix, inverse_matrix)
  }
  names(ans) <- rownames(delta_matrix)
  ans
}
