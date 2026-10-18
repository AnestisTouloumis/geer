compute_small_sample_covariances <- function(
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
  coef_names <- names(object$coefficients)
  lapply(covariances, function(x) {
    dimnames(x) <- list(coef_names, coef_names)
    x
  })
}


get_cached_small_sample_covariances <- function(object) {
  cache <- object$cache
  if (is.environment(cache)) {
    stored <- cache$small_sample
    if (!is.null(stored) &&
        identical(stored$coefficients, object$coefficients)) {
      return(stored$covariances)
    }
  }
  covariances <- compute_small_sample_covariances(object)
  if (is.environment(cache)) {
    assign(
      "small_sample",
      list(coefficients = object$coefficients, covariances = covariances),
      envir = cache
    )
  }
  covariances
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
  covariances[[cov_type]]
}
