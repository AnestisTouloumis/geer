#' @title marginaleffects Support for geer Objects
#'
#' @description
#' Methods that allow fitted \code{geer} objects to work with the
#' \pkg{marginaleffects} package.
#'
#' These functions are used internally by \pkg{marginaleffects} and are not
#' normally called directly by users. Once both packages are installed,
#' functions such as \code{marginaleffects::predictions()},
#' \code{marginaleffects::comparisons()}, and
#' \code{marginaleffects::slopes()} can be applied to fitted \code{geer}
#' objects.
#'
#' By default, \pkg{marginaleffects} inference for \code{geer} objects uses
#' the bias-corrected covariance matrix, matching the default inference in
#' \pkg{geer}. Supplying \code{vcov = "robust"} to a \pkg{marginaleffects}
#' function uses the sandwich covariance matrix. Other covariance estimators
#' can be supplied through a function, for example
#' \code{vcov = function(x) vcov(x, cov_type = "naive")}.
#'
#' @param model a fitted model object of class \code{"geer"}.
#' @param x a fitted model object of class \code{"geer"}.
#' @param coefs a named numeric vector of regression coefficients.
#' @param vcov covariance specification supplied by \pkg{marginaleffects}.
#'   \code{NULL} or \code{TRUE} uses the bias-corrected covariance matrix,
#'   \code{"robust"} uses the sandwich covariance matrix, and a square
#'   covariance matrix is returned unchanged after validation.
#' @param newdata a data frame containing predictor values at which to compute
#'   predictions.
#' @param type prediction scale. Supported values are \code{"response"} and
#'   \code{"link"}. \code{NULL} defaults to \code{"response"}.
#' @param ... additional arguments passed through by \pkg{marginaleffects}.
#'
#' @return
#' These methods return objects in the formats required internally by
#' \pkg{marginaleffects}: a named coefficient vector, a modified \code{geer}
#' object, a named covariance matrix, a model data frame, or a data frame
#' with columns \code{rowid} and \code{estimate}.
#'
#' @seealso \code{\link{predict.geer}}, \code{\link{vcov.geer}},
#'   \code{\link{coef.geer}}.
#'
#' @examples
#' if (requireNamespace("marginaleffects", quietly = TRUE)) {
#'   data("cerebrovascular", package = "geer")
#'
#'   fit <- geewa_binary(
#'     formula = ecg ~ treatment + factor(period),
#'     link = "logit",
#'     data = cerebrovascular,
#'     id = id,
#'     orstr = "exchangeable"
#'   )
#'
#'   marginaleffects::avg_predictions(fit)
#'   marginaleffects::avg_comparisons(fit, variables = "treatment")
#'   marginaleffects::avg_predictions(fit, vcov = "robust")
#' }
#'
#' @name marginaleffects-support
NULL


validate_marginaleffects_vcov <- function(vcov_matrix, coef_names) {
  if (!is.matrix(vcov_matrix) || !is.numeric(vcov_matrix)) {
    stop("'vcov' must be a numeric matrix", call. = FALSE)
  }
  if (nrow(vcov_matrix) != ncol(vcov_matrix)) {
    stop("'vcov' must be a square matrix", call. = FALSE)
  }
  if (nrow(vcov_matrix) != length(coef_names)) {
    stop(
      "'vcov' has dimensions incompatible with the model coefficients",
      call. = FALSE
    )
  }
  if (is.null(rownames(vcov_matrix)) || is.null(colnames(vcov_matrix))) {
    dimnames(vcov_matrix) <- list(coef_names, coef_names)
  }
  missing_rows <- setdiff(coef_names, rownames(vcov_matrix))
  missing_cols <- setdiff(coef_names, colnames(vcov_matrix))
  if (length(missing_rows) || length(missing_cols)) {
    stop("names of 'vcov' must match the model coefficient names", call. = FALSE)
  }
  vcov_matrix[coef_names, coef_names, drop = FALSE]
}


#' @rdname marginaleffects-support
#' @exportS3Method insight::get_data
get_data.geer <- function(x, ...) {
  x <- check_geer_object(x)
  if (!is.data.frame(x$data)) {
    stop(
      "the fitted 'geer' object does not contain a recoverable data frame",
      call. = FALSE
    )
  }
  x$data
}


#' @rdname marginaleffects-support
#' @exportS3Method marginaleffects::get_coef
get_coef.geer <- function(model, ...) {
  model <- check_geer_object(model)
  coef(model)
}


#' @rdname marginaleffects-support
#' @exportS3Method marginaleffects::set_coef
set_coef.geer <- function(model, coefs, ...) {
  model <- check_geer_object(model)
  if (!is.numeric(coefs) || is.null(names(coefs))) {
    stop("'coefs' must be a named numeric vector", call. = FALSE)
  }
  coef_names <- names(model$coefficients)
  if (length(coefs) != length(coef_names) ||
      !setequal(names(coefs), coef_names)) {
    stop("names of 'coefs' must match the model coefficient names", call. = FALSE)
  }
  model$coefficients <- as.numeric(coefs[coef_names])
  names(model$coefficients) <- coef_names
  model
}


#' @rdname marginaleffects-support
#' @exportS3Method marginaleffects::get_vcov
get_vcov.geer <- function(model, vcov = NULL, ...) {
  model <- check_geer_object(model)
  coef_names <- names(model$coefficients)

  if (is.null(vcov) || isTRUE(vcov)) {
    out <- stats::vcov(model, cov_type = "bias-corrected")
  } else if (isFALSE(vcov)) {
    return(NULL)
  } else if (is.character(vcov) && length(vcov) == 1L && !is.na(vcov)) {
    if (!identical(tolower(vcov), "robust")) {
      stop(
        paste0(
          "unsupported 'vcov' specification for a 'geer' model; use TRUE ",
          "for the bias-corrected covariance, 'robust' for the sandwich ",
          "covariance, or supply a covariance matrix/function"
        ),
        call. = FALSE
      )
    }
    out <- stats::vcov(model, cov_type = "robust")
  } else if (is.function(vcov)) {
    out <- vcov(model)
  } else {
    out <- vcov
  }

  validate_marginaleffects_vcov(out, coef_names)
}


#' @rdname marginaleffects-support
#' @exportS3Method marginaleffects::get_predict
get_predict.geer <- function(model,
                             newdata = NULL,
                             type = NULL,
                             ...) {
  model <- check_geer_object(model)
  if (is.null(type)) {
    type <- "response"
  } else {
    type <- match.arg(type, c("response", "link"))
  }

  if (is.null(newdata)) {
    if (!is.data.frame(model$data)) {
      stop(
        paste0(
          "automatic prediction data recovery requires a model fitted with ",
          "a data frame; supply 'newdata' explicitly"
        ),
        call. = FALSE
      )
    }
    newdata <- model$data
  }
  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  estimate <- stats::predict(
    model,
    newdata = newdata,
    type = type,
    ...
  )
  if (is.list(estimate) && !is.null(estimate$fit)) {
    estimate <- estimate$fit
  }
  if (!is.numeric(estimate) || length(estimate) != nrow(newdata)) {
    stop("'predict.geer' returned an invalid prediction vector", call. = FALSE)
  }

  rowid <- if ("rowid" %in% names(newdata)) {
    newdata$rowid
  } else {
    seq_len(nrow(newdata))
  }

  data.frame(
    rowid = rowid,
    estimate = as.numeric(estimate)
  )
}
