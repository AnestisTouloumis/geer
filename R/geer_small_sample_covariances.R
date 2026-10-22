compute_ss_covariances <- function(
    object,
    fitted_values = object$fitted.values,
    linear_predictors = object$linear.predictors) {
  covariances <- if (is_geewa_fit(object)) {
    get_small_sample_covariances_cc(
      object$y, object$x, object$id, object$repeated, object$prior.weights,
      object$family$link, object$family$family,
      fitted_values, linear_predictors,
      object$association_structure, object$alpha, object$phi
    )
  } else {
    get_small_sample_covariances_or(
      object$y, object$x, object$id, object$repeated, object$prior.weights,
      object$family$link,
      fitted_values, linear_predictors,
      get_alpha_or(object)
    )
  }
  singular_ids <- as.numeric(covariances$singular_clusters)
  covariances$singular_clusters <- NULL
  coef_names <- names(object$coefficients)
  covariances <- lapply(covariances, function(x) {
    dimnames(x) <- list(coef_names, coef_names)
    x
  })
  attr(covariances, "singular_clusters") <-
    get_cluster_labels(object)[match(singular_ids, unique(object$id))]
  covariances
}


get_cached_ss_covariances <- function(object) {
  get_cached_value(
    object,
    "small_sample",
    function() compute_ss_covariances(object)
  )
}


select_small_sample_covariance <- function(covariances,
                                           cov_type,
                                           object,
                                           coef_no,
                                           context) {
  if (identical(cov_type, "gosho-sato-takeuchi") &&
      object$clusters_no <= coef_no) {
    stop(
      sprintf(
        paste0(
          "%s failed: the gosho-sato-takeuchi covariance requires more ",
          "clusters than coefficients"
        ),
        context
      ),
      call. = FALSE
    )
  }
  if (cov_type %in% geer_pooled_cov_types && !is_balanced_design(object)) {
    stop(
      sprintf(
        paste0(
          "%s failed: the %s covariance requires a balanced design, in which ",
          "every cluster observes each repeated position once"
        ),
        context,
        cov_type
      ),
      call. = FALSE
    )
  }
  singular_clusters <- attr(covariances, "singular_clusters")
  if (cov_type %in% geer_leverage_cov_types && length(singular_clusters)) {
    stop(
      sprintf(
        paste0(
          "%s failed: the %s covariance is undefined because I - H is ",
          "singular for cluster(s) %s, whose leverage is 1"
        ),
        context,
        cov_type,
        paste(singular_clusters, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  covariances[[cov_type]]
}
