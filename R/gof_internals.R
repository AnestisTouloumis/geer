check_gof_fit <- function(object, context) {
  binary_family <- identical(get_geer_fit_function(object), "geewa_binary") ||
    object$family$family %in% c("binomial", "quasibinomial")
  y <- object$y
  if (!binary_family || anyNA(y) || !all(y == 0 | y == 1) ||
      !all(object$prior.weights == 1)) {
    stop(
      sprintf(
        paste0(
          "%s failed: goodness-of-fit tests require a binary (0/1) response ",
          "with unit prior weights"
        ),
        context
      ),
      call. = FALSE
    )
  }
  mu <- object$fitted.values
  if (any(!is.finite(mu) | mu <= 0 | mu >= 1)) {
    stop(
      sprintf(
        paste0(
          "%s failed: some fitted probabilities are 0 or 1, which suggests ",
          "separation"
        ),
        context
      ),
      call. = FALSE
    )
  }
  if (!identical(object$method, "gee")) {
    warning(
      sprintf(
        paste0(
          "the goodness-of-fit tests were developed for method = \"gee\"; ",
          "for method = \"%s\" they are evaluated at the fitted values, and ",
          "their null distribution has not been studied"
        ),
        object$method
      ),
      call. = FALSE
    )
  }
  invisible(object)
}


check_gof_cov_type <- function(cov_type, context) {
  if (identical(cov_type, "jackknife")) {
    stop(
      sprintf(
        paste0(
          "%s failed: cov_type = \"jackknife\" is not available for ",
          "goodness-of-fit tests"
        ),
        context
      ),
      call. = FALSE
    )
  }
  cov_type
}


check_gof_by_repeated <- function(by_repeated) {
  if (!is.logical(by_repeated) || length(by_repeated) != 1L ||
      is.na(by_repeated)) {
    stop("'by_repeated' must be TRUE or FALSE", call. = FALSE)
  }
  invisible(by_repeated)
}


build_gof_risk_groups <- function(object, groups, context) {
  mu <- object$fitted.values
  breaks <- unique(stats::quantile(
    mu,
    probs = seq(0, 1, length.out = groups + 1L),
    names = FALSE
  ))
  if (length(breaks) < 3L) {
    stop(
      sprintf(
        paste0(
          "%s failed: the fitted probabilities take too few distinct values ",
          "to form at least two groups"
        ),
        context
      ),
      call. = FALSE
    )
  }
  risk <- cut(mu, breaks = breaks, include.lowest = TRUE)
  if (nlevels(risk) < groups) {
    warning(
      sprintf(
        paste0(
          "only %d groups of fitted probabilities could be formed because ",
          "of tied fitted probabilities, instead of the %d requested"
        ),
        nlevels(risk),
        groups
      ),
      call. = FALSE
    )
  }
  risk
}


extract_gof_region_variables <- function(object, regions, env, context) {
  variables <- as.list(attr(stats::terms(regions), "variables"))[-1L]
  if (!length(variables)) {
    stop(
      sprintf(
        "%s failed: 'regions' must contain at least one variable",
        context
      ),
      call. = FALSE
    )
  }
  extracted <- extract_fit_data_variables(
    object,
    variables,
    enclos = environment(regions),
    env = env
  )
  if (identical(extracted$status, "incomplete")) {
    stop(
      sprintf(
        paste0(
          "%s failed: the variables in 'regions' have missing values for ",
          "observations used to fit the model; supply a vector with one ",
          "value per observation instead"
        ),
        context
      ),
      call. = FALSE
    )
  }
  if (!identical(extracted$status, "ok")) {
    stop(
      sprintf(
        paste0(
          "%s failed: 'regions' cannot be evaluated in the data used to fit ",
          "the model; supply a vector with one value per observation instead"
        ),
        context
      ),
      call. = FALSE
    )
  }
  extracted$values
}


build_gof_regions <- function(object, regions, env, context) {
  values <- if (inherits(regions, "formula") && length(regions) == 2L) {
    extract_gof_region_variables(object, regions, env, context)
  } else if (is.atomic(regions) && length(regions) == object$obs_no &&
             !anyNA(regions)) {
    list(regions)
  } else {
    stop(
      sprintf(
        paste0(
          "%s failed: 'regions' must be a one-sided formula or a vector with ",
          "one non-missing value per observation, in the order stored in ",
          "the fit"
        ),
        context
      ),
      call. = FALSE
    )
  }
  region <- interaction(lapply(values, factor), drop = TRUE, lex.order = TRUE)
  if (nlevels(region) < 2L) {
    stop(
      sprintf("%s failed: 'regions' must define at least two regions", context),
      call. = FALSE
    )
  }
  region
}


check_gof_regions <- function(object, region, by_repeated) {
  count_clusters <- function(x) length(unique(x))
  clusters <- vapply(split(object$id, region), count_clusters, numeric(1))
  cells <- if (by_repeated) {
    interaction(region, object$repeated, drop = TRUE)
  } else {
    region
  }
  events <- vapply(split(object$y, cells), sum, numeric(1))
  sizes <- tabulate(cells, nlevels(cells))
  problems <- c(
    if (any(clusters < 10)) "some regions have fewer than 10 clusters",
    if (mean(clusters < 25) > 0.25) {
      "more than 25% of the regions have fewer than 25 clusters"
    },
    if (any(events == 0 | events == sizes)) {
      if (by_repeated) {
        "some region by position cells have no events or no non-events"
      } else {
        "some regions have no events or no non-events"
      }
    }
  )
  if (length(problems)) {
    warning(
      sprintf(
        paste0(
          "the chi-squared approximation may be unreliable (Barnhart and ",
          "Williamson, 1998): %s"
        ),
        paste(problems, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}


build_gof_alternative <- function(object, partition, by_repeated) {
  extra <- stats::model.matrix(~ 0 + partition)
  if (by_repeated) {
    occasion <- factor(object$repeated)
    extra <- cbind(
      extra,
      stats::model.matrix(~ 0 + occasion),
      stats::model.matrix(~ 0 + partition:occasion)
    )
  }
  coef_no <- ncol(object$x)
  decomposition <- qr(cbind(object$x, extra))
  kept <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  kept <- kept[kept > coef_no] - coef_no
  if (!length(kept)) {
    return(NULL)
  }
  extra <- extra[, kept, drop = FALSE]
  colnames(extra) <- paste0(".gof", seq_len(ncol(extra)))
  alternative <- object
  alternative$x <- cbind(object$x, extra)
  alternative$coefficients <- c(
    object$coefficients,
    stats::setNames(numeric(ncol(extra)), colnames(extra))
  )
  alternative$cache <- new.env(parent = emptyenv())
  alternative
}


compute_gof_score_test <- function(object,
                                   partition,
                                   by_repeated,
                                   cov_type,
                                   context) {
  alternative <- build_gof_alternative(object, partition, by_repeated)
  if (is.null(alternative)) {
    stop(
      sprintf(
        paste0(
          "%s failed: the partition adds no column that is linearly ",
          "independent of the model matrix"
        ),
        context
      ),
      call. = FALSE
    )
  }
  compute_score_test(object, alternative, cov_type = cov_type)
}


build_gof_table <- function(object, partition) {
  data.frame(
    observations = tabulate(partition, nlevels(partition)),
    clusters = vapply(
      split(object$id, partition),
      function(x) length(unique(x)),
      numeric(1)
    ),
    observed = vapply(split(object$y, partition), sum, numeric(1)),
    expected = vapply(
      split(object$fitted.values, partition),
      sum,
      numeric(1)
    ),
    row.names = levels(partition)
  )
}


build_gof_response_covariance <- function(object, residuals, variances, type) {
  cluster_indices <- split(seq_along(object$id), object$id)
  if (identical(type, "empirical")) {
    return(lapply(cluster_indices, function(i) tcrossprod(residuals[i])))
  }
  positions <- as.integer(object$repeated)
  positions_no <- max(positions)
  standardized <- residuals / sqrt(variances)
  cross_sum <- matrix(0, positions_no, positions_no)
  pair_counts <- matrix(0, positions_no, positions_no)
  for (i in cluster_indices) {
    cross_sum[positions[i], positions[i]] <-
      cross_sum[positions[i], positions[i]] + tcrossprod(standardized[i])
    pair_counts[positions[i], positions[i]] <-
      pair_counts[positions[i], positions[i]] + 1
  }
  pooled <- cross_sum / pair_counts
  lapply(cluster_indices, function(i) {
    scale <- sqrt(variances[i])
    scale * pooled[positions[i], positions[i], drop = FALSE] *
      rep(scale, each = length(i))
  })
}


compute_gof_residual_parts <- function(object, covariance) {
  mu <- object$fitted.values
  residuals <- object$y - mu
  variances <- mu * (1 - mu)
  cluster_indices <- split(seq_along(object$id), object$id)
  d_matrix <- object$x * object$family$mu.eta(object$linear.predictors)
  list(
    mu = mu,
    residuals = residuals,
    variances = variances,
    cluster_indices = cluster_indices,
    d_matrix = d_matrix,
    naive_covariance = object$naive_covariance,
    inverse_v_d = lapply(cluster_indices, function(i) {
      solve(
        compute_working_covariance_for_criteria(object, i),
        d_matrix[i, , drop = FALSE]
      )
    }),
    covariances = build_gof_response_covariance(
      object,
      residuals,
      variances,
      covariance
    )
  )
}


compute_gof_residual_test <- function(components, statistic, context) {
  mu <- components$mu
  residuals <- components$residuals
  variances <- components$variances
  if (identical(statistic, "pearson")) {
    value <- sum(residuals^2 / variances)
    mean_value <- length(mu)
    weights <- (1 - 2 * mu) / variances
  } else {
    value <- sum(residuals^2)
    mean_value <- sum(variances)
    weights <- 1 - 2 * mu
  }
  projected <- components$naive_covariance %*%
    crossprod(components$d_matrix, weights)
  cluster_indices <- components$cluster_indices
  variance <- sum(vapply(
    seq_along(cluster_indices),
    function(k) {
      i <- cluster_indices[[k]]
      adjusted <- weights[i] - drop(components$inverse_v_d[[k]] %*% projected)
      drop(crossprod(adjusted, components$covariances[[k]] %*% adjusted))
    },
    numeric(1)
  ))
  if (!is.finite(variance) || variance <= 0) {
    stop(
      sprintf(
        paste0(
          "%s failed: the approximate variance of the statistic is not ",
          "positive"
        ),
        context
      ),
      call. = FALSE
    )
  }
  z_value <- (value - mean_value) / sqrt(variance)
  list(
    value = value,
    mean = mean_value,
    variance = variance,
    z = z_value,
    p_value = 2 * stats::pnorm(-abs(z_value))
  )
}
