new_geer <- function(x) {
  if (!is.list(x)) {
    stop("'x' must be a list", call. = FALSE)
  }

  required <- c(
    "coefficients", "residuals", "fitted.values", "qr", "rank", "family",
    "linear.predictors", "iter", "prior.weights", "df.residual", "y", "x",
    "id", "repeated", "converged", "call", "formula", "terms",
    "data", "offset", "control", "method",
    "naive_covariance", "robust_covariance", "bias_corrected_covariance",
    "association_structure", "alpha", "phi", "obs_no", "clusters_no",
    "min_cluster_size", "max_cluster_size"
  )

  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "Input is missing required components: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!("na.action" %in% names(x))) {
    x$na.action <- NULL
  }

  if (!("contrasts" %in% names(x))) {
    x$contrasts <- NULL
  }

  if (!("xlevels" %in% names(x))) {
    x$xlevels <- list()
  }

  class(x) <- unique(c("geer", class(x)))
  x
}


validate_geer <- function(x) {
  if (!is.list(x)) {
    stop("a 'geer' object must be a list", call. = FALSE)
  }
  if (!inherits(x, "geer")) {
    stop("object must inherit from class 'geer'", call. = FALSE)
  }

  required <- c(
    "coefficients",
    "fitted.values",
    "residuals",
    "family",
    "formula",
    "terms",
    "call",
    "x",
    "y",
    "id",
    "alpha",
    "association_structure",
    "method",
    "df.residual"
  )
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields)) {
    stop(
      sprintf(
        "invalid 'geer' object: missing required component(s): %s",
        paste(missing_fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(x$coefficients)) {
    stop("'coefficients' must be numeric", call. = FALSE)
  }
  if (is.null(names(x$coefficients)) || !length(names(x$coefficients))) {
    stop("'coefficients' must be a named numeric vector", call. = FALSE)
  }

  if (!is.matrix(x$x)) {
    stop("'x' must be a matrix", call. = FALSE)
  }
  if (ncol(x$x) != length(x$coefficients)) {
    stop("number of columns of 'x' must match length of 'coefficients'", call. = FALSE)
  }

  if (!is.numeric(x$fitted.values)) {
    stop("'fitted.values' must be numeric", call. = FALSE)
  }
  if (!is.numeric(x$residuals)) {
    stop("'residuals' must be numeric", call. = FALSE)
  }

  n_obs <- NROW(x$fitted.values)

  if (NROW(x$x) != n_obs) {
    stop("number of rows of 'x' must match length of 'fitted.values'", call. = FALSE)
  }
  if (NROW(x$id) != n_obs) {
    stop("length of 'id' must match number of observations", call. = FALSE)
  }
  if (NROW(x$y) != n_obs) {
    stop("response 'y' must match number of observations", call. = FALSE)
  }

  if ("obs_no" %in% names(x)) {
    if (!is.numeric(x$obs_no) || length(x$obs_no) != 1L ||
        !is.finite(x$obs_no) || x$obs_no < 0 ||
        abs(x$obs_no - round(x$obs_no)) > sqrt(.Machine$double.eps)) {
      stop("'obs_no' must be a single non-negative integer", call. = FALSE)
    }

    if (x$obs_no != n_obs) {
      stop("'obs_no' is inconsistent with observation-level components", call. = FALSE)
    }
  }

  if ("repeated" %in% names(x) && !is.null(x$repeated)) {
    if (NROW(x$repeated) != n_obs) {
      stop("'repeated' must match number of observations", call. = FALSE)
    }
  }

  if (!inherits(x$terms, "terms")) {
    stop("'terms' must be a terms object", call. = FALSE)
  }

  if (!is.call(x$call)) {
    stop("'call' must be a call", call. = FALSE)
  }

  if (!("formula" %in% names(x$call))) {
    stop("'call' must contain a formula component", call. = FALSE)
  }

  if (!identical(
    paste(deparse(x$formula), collapse = ""),
    paste(deparse(x$call$formula), collapse = "")
  )) {
    stop("'formula' and 'call$formula' must be identical", call. = FALSE)
  }

  if (!is.list(x$family)) {
    stop("'family' must be a family object", call. = FALSE)
  }
  if (is.null(x$family$family) || !is.character(x$family$family) || length(x$family$family) != 1L) {
    stop("'family$family' must be a length-1 character value", call. = FALSE)
  }
  if (is.null(x$family$link) || !is.character(x$family$link) || length(x$family$link) != 1L) {
    stop("'family$link' must be a length-1 character value", call. = FALSE)
  }
  if (is.null(x$family$linkinv) || !is.function(x$family$linkinv)) {
    stop("'family$linkinv' must be a function", call. = FALSE)
  }

  if (!is.numeric(x$alpha)) {
    stop("'alpha' must be numeric", call. = FALSE)
  }
  if (!is.character(x$association_structure) || length(x$association_structure) != 1L) {
    stop("'association_structure' must be a length-1 character value", call. = FALSE)
  }
  if (!is.character(x$method) || length(x$method) != 1L) {
    stop("'method' must be a length-1 character value", call. = FALSE)
  }

  if (!is.numeric(x$df.residual) || length(x$df.residual) != 1L || !is.finite(x$df.residual)) {
    stop("'df.residual' must be a single finite numeric value", call. = FALSE)
  }

  if ("rank" %in% names(x)) {
    if (!is.numeric(x$rank) || length(x$rank) != 1L ||
        !is.finite(x$rank) || x$rank < 0 ||
        abs(x$rank - round(x$rank)) > sqrt(.Machine$double.eps)) {
      stop("'rank' must be a single non-negative integer", call. = FALSE)
    }
  }

  if ("converged" %in% names(x)) {
    if (!is.logical(x$converged) || length(x$converged) != 1L || is.na(x$converged)) {
      stop("'converged' must be a single non-missing logical value", call. = FALSE)
    }
  }

  if ("phi" %in% names(x) && !is.null(x$phi)) {
    if (!is.numeric(x$phi) || length(x$phi) != 1L || !is.finite(x$phi)) {
      stop("'phi' must be a single finite numeric value", call. = FALSE)
    }
  }

  if (!is.null(colnames(x$x)) && !identical(colnames(x$x), names(x$coefficients))) {
    stop("column names of 'x' must match names of 'coefficients'", call. = FALSE)
  }

  if (!is.null(x$xlevels) && !is.list(x$xlevels)) {
    stop("'xlevels' must be NULL or a named list", call. = FALSE)
  }
  if (is.list(x$xlevels) && length(x$xlevels) > 0L && is.null(names(x$xlevels))) {
    stop("'xlevels' must be a named list", call. = FALSE)
  }

  if (!is.null(x$contrasts) && !is.list(x$contrasts)) {
    stop("'contrasts' must be NULL or a named list", call. = FALSE)
  }
  if (is.list(x$contrasts) && length(x$contrasts) > 0L && is.null(names(x$contrasts))) {
    stop("'contrasts' must be a named list", call. = FALSE)
  }

  if ("clusters_no" %in% names(x)) {
    if (!is.numeric(x$clusters_no) || length(x$clusters_no) != 1L ||
        !is.finite(x$clusters_no) || x$clusters_no < 1 ||
        abs(x$clusters_no - round(x$clusters_no)) > sqrt(.Machine$double.eps)) {
      stop("'clusters_no' must be a single positive integer", call. = FALSE)
    }

    if (length(unique(x$id)) != x$clusters_no) {
      stop("'clusters_no' is inconsistent with 'id'", call. = FALSE)
    }
  }

  if ("min_cluster_size" %in% names(x)) {
    if (!is.numeric(x$min_cluster_size) || length(x$min_cluster_size) != 1L ||
        !is.finite(x$min_cluster_size) || x$min_cluster_size < 1 ||
        abs(x$min_cluster_size - round(x$min_cluster_size)) > sqrt(.Machine$double.eps)) {
      stop("'min_cluster_size' must be a single positive integer", call. = FALSE)
    }
  }

  if ("max_cluster_size" %in% names(x)) {
    if (!is.numeric(x$max_cluster_size) || length(x$max_cluster_size) != 1L ||
        !is.finite(x$max_cluster_size) || x$max_cluster_size < 1 ||
        abs(x$max_cluster_size - round(x$max_cluster_size)) > sqrt(.Machine$double.eps)) {
      stop("'max_cluster_size' must be a single positive integer", call. = FALSE)
    }
  }

  if (all(c("min_cluster_size", "max_cluster_size") %in% names(x))) {
    if (x$min_cluster_size > x$max_cluster_size) {
      stop("'min_cluster_size' cannot exceed 'max_cluster_size'", call. = FALSE)
    }
  }

  cov_fields <- c(
    "naive_covariance",
    "robust_covariance",
    "bias_corrected_covariance",
    "df_adjusted_covariance"
  )

  for (nm in intersect(cov_fields, names(x))) {
    V <- x[[nm]]
    if (is.null(V)) {
      next
    }
    if (!is.matrix(V) || !is.numeric(V)) {
      stop(sprintf("'%s' must be a numeric matrix", nm), call. = FALSE)
    }
    if (nrow(V) != length(x$coefficients) || ncol(V) != length(x$coefficients)) {
      stop(
        sprintf("'%s' must have dimensions compatible with 'coefficients'", nm),
        call. = FALSE
      )
    }
    if (is.null(rownames(V)) || is.null(colnames(V))) {
      stop(sprintf("'%s' must have row and column names", nm), call. = FALSE)
    }
    if (!identical(rownames(V), names(x$coefficients)) ||
        !identical(colnames(V), names(x$coefficients))) {
      stop(
        sprintf("dimnames of '%s' must match names of 'coefficients'", nm),
        call. = FALSE
      )
    }
  }

  x
}
