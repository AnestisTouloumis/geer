check_prediction_offset_variables <- function(expression, newdata, source) {
  needed <- all.vars(expression)
  absent <- setdiff(needed, names(newdata))
  if (length(absent)) {
    stop(
      sprintf(
        paste0(
          "prediction offset variable(s) %s supplied via the %s are not ",
          "present in 'newdata'; they must be supplied in 'newdata' rather ",
          "than taken from the environment in which the model was fitted"
        ),
        paste(sprintf("'%s'", absent), collapse = ", "),
        source
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}


evaluate_prediction_offset <- function(object, newdata, n_obs) {
  offset_formula <- rep.int(0, n_obs)
  pred_terms <- stats::delete.response(object$terms)
  offset_position <- attr(pred_terms, "offset")
  if (length(offset_position)) {
    variable_list <- as.list(attr(pred_terms, "variables"))[-1L]
    for (expression in variable_list[offset_position]) {
      check_prediction_offset_variables(expression, newdata, "model formula")
    }
  }
  mf <- stats::model.frame(
    pred_terms,
    newdata,
    xlev = object$xlevels,
    na.action = stats::na.pass
  )
  offset_from_formula <- stats::model.offset(mf)
  if (!is.null(offset_from_formula)) {
    if (length(offset_from_formula) == 1L) {
      offset_from_formula <- rep.int(as.numeric(offset_from_formula), n_obs)
    }
    if (!is.numeric(offset_from_formula) ||
        length(offset_from_formula) != n_obs ||
        any(!is.finite(offset_from_formula))) {
      stop(
        paste0(
          "formula-based prediction offset must be a finite numeric vector ",
          "of the correct length"
        ),
        call. = FALSE
      )
    }
    offset_formula <- as.numeric(offset_from_formula)
  }
  offset_argument <- rep.int(0, n_obs)
  call_offset <- object$call$offset
  if (!is.null(call_offset)) {
    check_prediction_offset_variables(call_offset, newdata, "'offset' argument")
    offset_env <- environment(stats::formula(object))
    offset_from_argument <- tryCatch(
      eval(call_offset, envir = newdata, enclos = offset_env),
      error = function(e) {
        stop(
          paste0(
            "could not evaluate the offset supplied via the 'offset' ",
            "argument in 'newdata'"
          ),
          call. = FALSE
        )
      }
    )
    if (length(offset_from_argument) == 1L) {
      offset_from_argument <- rep.int(as.numeric(offset_from_argument), n_obs)
    }
    if (!is.numeric(offset_from_argument) ||
        length(offset_from_argument) != n_obs ||
        any(!is.finite(offset_from_argument))) {
      stop(
        paste0(
          "prediction offset supplied via the 'offset' argument must be a ",
          "finite numeric vector of the correct length"
        ),
        call. = FALSE
      )
    }
    offset_argument <- as.numeric(offset_from_argument)
  }
  offset_formula + offset_argument
}


#' @title
#' Extract Model Fitted Values from a geer Object
#'
#' @rdname fitted.geer
#' @method fitted geer
#'
#' @description
#' Extracts the fitted mean values on the response scale from a fitted
#' \code{geer} object. \code{fitted.values} is an alias for \code{fitted}.
#'
#' @inheritParams coef.geer
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @return
#' A numeric vector of fitted means, one per observation used in the fit, in
#' the order given under Observation order.
#'
#' @section Observation order:
#' The values follow the internal order of the fitted model: the complete
#' cases, sorted by cluster and, within cluster, by \code{repeated}. Unless the
#' data were already in this order and had no missing values, assigning the
#' values to the data, as in \code{data$r <- residuals(fit)}, misaligns the
#' rows without warning. \code{fit$id} and \code{fit$repeated} are stored in
#' the same order and can be used to align them.
#'
#' @seealso \code{\link{residuals.geer}}, \code{\link{predict.geer}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' head(fitted(fit))
#'
#' @export
fitted.geer <- function(object, ...) {
  check_unused_dots(list(...), "fitted.geer")
  object <- check_geer_object(object)
  object$fitted.values
}


#' @title
#' Predictions from a geer Object
#'
#' @method predict geer
#'
#' @description
#' Generates predictions from a fitted \code{geer} object, optionally with
#' approximate standard errors.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param newdata optional data frame containing the variables used for
#'   prediction. Defaults to \code{NULL}, which uses the data used to fit the
#'   model.
#' @param type character string specifying the scale of the predictions: linear
#'   predictor (\code{"link"}) or response (\code{"response"}). Defaults to
#'   \code{"link"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for the standard errors when \code{se.fit = TRUE}: one of the
#'   estimators described in \code{\link{vcov.geer}}. Defaults to
#'   \code{"bias-corrected"}.
#' @param se.fit logical indicating whether standard errors are returned.
#'   Defaults to \code{FALSE}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @details
#' The linear predictor is computed from the model matrix of \code{newdata},
#' or of the fitted data, and the offset, and transformed by the inverse link
#' for \code{type = "response"}. Standard errors are obtained by the delta
#' method. When \code{newdata} is supplied and the model has an offset, the
#' variables of the offset must be in \code{newdata}.
#'
#' @return
#' For \code{se.fit = FALSE}, a numeric vector of predictions. For
#' \code{se.fit = TRUE}, a list with components:
#' \item{fit}{the predictions.}
#' \item{se.fit}{their standard errors.}
#' With \code{newdata}, the results are named by its row names; otherwise they
#' are unnamed.
#'
#' @section Observation order:
#' With \code{newdata}, the predictions follow its rows, and rows with missing
#' values give \code{NA}. Without it, they follow the internal order of the
#' fitted model: the complete cases, sorted by cluster and, within cluster, by
#' \code{repeated}. To obtain predictions in the order of the data, supply the
#' data as \code{newdata}.
#'
#' @seealso \code{\link{fitted.geer}}, \code{\link{residuals.geer}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}},
#'   \code{\link[stats]{predict}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' head(predict(fit))
#' head(predict(fit, type = "response"))
#'
#' ## Predictions for new data, with standard errors
#' pred <- predict(fit, newdata = cerebrovascular[1:5, ], type = "response",
#'                 se.fit = TRUE)
#' pred
#'
#' @export
predict.geer <- function(object,
                         newdata = NULL,
                         type = c("link", "response"),
                         cov_type = "bias-corrected",
                         se.fit = FALSE,
                         ...) {
  check_unused_dots(list(...), "predict.geer")
  object <- check_geer_object(object)
  type <- match.arg(type)
  cov_type <- match_cov_type(cov_type)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("'se.fit' must be TRUE or FALSE", call. = FALSE)
  }
  se.fit <- isTRUE(se.fit)
  coef_names <- names(object$coefficients)
  if (is.null(coef_names) || !length(coef_names)) {
    stop("'coefficients' must be named", call. = FALSE)
  }
  if (is.null(newdata)) {
    eta_vector <- object$linear.predictors
    mu_vector <- object$fitted.values
    out <- if (type == "link") eta_vector else mu_vector
    if (!se.fit) {
      return(out)
    }
    vcov_matrix <- stats::vcov(object, cov_type = cov_type)
    vcov_matrix <- vcov_matrix[coef_names, coef_names, drop = FALSE]
    design_matrix <- object$x
    if (is.null(colnames(design_matrix))) {
      stop("'object$x' must have column names", call. = FALSE)
    }
    design_matrix <- design_matrix[, coef_names, drop = FALSE]
    se <- compute_standard_errors_or_na(
      rowSums((design_matrix %*% vcov_matrix) * design_matrix),
      "predict.geer",
      what = "prediction"
    )
    if (type == "response") {
      se <- se * abs(object$family$mu.eta(eta_vector))
    }
    ans <- list(fit = out, se.fit = as.numeric(se))
    return(ans)
  }
  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  pred_terms <- stats::delete.response(object$terms)
  mf <- stats::model.frame(
    pred_terms,
    newdata,
    xlev = object$xlevels,
    na.action = stats::na.pass
  )
  design_matrix <- stats::model.matrix(
    pred_terms,
    mf,
    contrasts.arg = object$contrasts
  )
  if (is.null(colnames(design_matrix))) {
    stop("prediction model matrix must have column names", call. = FALSE)
  }
  missing_cols <- setdiff(coef_names, colnames(design_matrix))
  if (length(missing_cols)) {
    stop(
      sprintf(
        "prediction data are missing required model columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  design_matrix <- design_matrix[, coef_names, drop = FALSE]
  eta_vector <- drop(design_matrix %*% object$coefficients)
  offset_vector <- evaluate_prediction_offset(
    object = object,
    newdata = newdata,
    n_obs = nrow(design_matrix)
  )
  eta_vector <- eta_vector + offset_vector
  if (!se.fit) {
    if (type == "link") {
      return(eta_vector)
    }
    return(object$family$linkinv(eta_vector))
  }
  vcov_matrix <- stats::vcov(object, cov_type = cov_type)
  vcov_matrix <- vcov_matrix[coef_names, coef_names, drop = FALSE]
  se <- compute_standard_errors_or_na(
    rowSums((design_matrix %*% vcov_matrix) * design_matrix),
    "predict.geer", what = "prediction"
  )
  if (type == "response") {
    mu_vector <- object$family$linkinv(eta_vector)
    se <- se * abs(object$family$mu.eta(eta_vector))
    ans <- list(fit = mu_vector, se.fit = se)
    return(ans)
  }
  ans <- list(fit = eta_vector, se.fit = se)
  ans
}


compute_mahalanobis_residuals <- function(object) {
  cluster_index <- split(seq_along(object$id), object$id)
  cluster_ids <- names(cluster_index)
  ans <- vapply(
    cluster_ids,
    FUN.VALUE = numeric(1),
    FUN = function(cluster_id) {
      indices <- cluster_index[[cluster_id]]
      residual_vector <- object$y[indices] - object$fitted.values[indices]
      working_covariance <- compute_working_covariance_for_criteria(
        object,
        indices
      )
      chol_factor <- tryCatch(
        chol(working_covariance),
        error = function(e) NULL
      )
      if (is.null(chol_factor)) {
        stop(
          sprintf(
            paste0(
              "Mahalanobis residual could not be computed for cluster '%s' ",
              "because its working covariance matrix is not positive definite"
            ),
            as.character(cluster_id)
          ),
          call. = FALSE
        )
      }
      solved <- forwardsolve(t(chol_factor), residual_vector)
      value <- sum(solved^2) / length(indices)
      if (!is.finite(value)) {
        stop(
          sprintf(
            paste0(
              "Mahalanobis residual could not be computed for cluster '%s' ",
              "because the result is non-finite"
            ),
            as.character(cluster_id)
          ),
          call. = FALSE
        )
      }
      value
    }
  )
  ans
}


#' @title
#' Residuals from a geer Object
#'
#' @method residuals geer
#'
#' @description
#' Extracts residuals of different types from a fitted \code{geer} object.
#'
#' @inheritParams coef.geer
#' @param type character string specifying the type of residuals: raw
#'   (\code{"working"}, or its alias \code{"response"}), Pearson
#'   (\code{"pearson"}), deviance (\code{"deviance"}), or cluster-level
#'   Mahalanobis (\code{"mahalanobis"}) residuals. See Details. Defaults to
#'   \code{"working"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @details
#' With prior weights \eqn{w_{ij}}, the Pearson and deviance residuals are
#' \deqn{r^P_{ij} = \frac{y_{ij} - \hat\mu_{ij}}{\sqrt{\hat\phi
#' V(\hat\mu_{ij})/w_{ij}}}, \qquad
#' r^D_{ij} = \mathrm{sign}(y_{ij} - \hat\mu_{ij})
#' \sqrt{d_{ij}/\hat\phi},}{r^P_ij = (y_ij - mu_ij)/sqrt(phi V(mu_ij)/w_ij),   r^D_ij = sign(y_ij - mu_ij) sqrt(d_ij/phi),}
#' where \eqn{d_{ij}} is the weighted deviance contribution of the family, as
#' in \pkg{glmtoolbox} (Vanegas, Rondon and Paula, 2023). The Mahalanobis
#' residual of cluster \eqn{i} is
#' \deqn{r^M_i = \frac{1}{n_i}(y_i - \hat\mu_i)^\top \hat V_i^{-1}
#' (y_i - \hat\mu_i),}{r^M_i = (1/n_i)(y_i - mu_i)' V_i^(-1) (y_i - mu_i),}
#' where \eqn{\hat V_i}{V_i} is the fitted working covariance matrix, including the
#' dispersion. It is a nonnegative cluster-level diagnostic.
#'
#' @section Observation order:
#' The values follow the internal order of the fitted model: the complete
#' cases, sorted by cluster and, within cluster, by \code{repeated}. Unless the
#' data were already in this order and had no missing values, assigning the
#' values to the data, as in \code{data$r <- residuals(fit)}, misaligns the
#' rows without warning. \code{fit$id} and \code{fit$repeated} are stored in
#' the same order and can be used to align them.
#'
#' @return
#' A numeric vector with one value per observation, in the order given under
#' Observation order, or, for \code{type = "mahalanobis"}, one value per
#' cluster, named by the cluster identifier.
#'
#' @references
#' Vanegas, L.H., Rondon, L.M. and Paula, G.A. (2023) Generalized Estimating
#' Equations using the new R package glmtoolbox. \emph{The R Journal},
#' \bold{15}, 105--133.
#'
#' @seealso \code{\link{fitted.geer}}, \code{\link{predict.geer}},
#'   \code{\link{runs_test}}, \code{\link{geewa}}, \code{\link{geewa_binary}},
#'   \code{\link[stats]{residuals}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' head(residuals(fit))
#' head(residuals(fit, type = "pearson"))
#' head(residuals(fit, type = "deviance"))
#' head(residuals(fit, type = "mahalanobis"))
#'
#' @export
residuals.geer <- function(object,
                           type = c(
                             "working", "pearson", "deviance", "mahalanobis",
                             "response"
                           ),
                           ...) {
  check_unused_dots(list(...), "residuals.geer")
  object <- check_geer_object(object)
  type <- match.arg(type)
  if (type %in% c("pearson", "deviance")) {
    y <- as.numeric(object$y)
    mu <- object$fitted.values
    weights <- object$prior.weights
  }
  if (!identical(type, "working")) {
    if (!is.numeric(object$phi) || length(object$phi) != 1L ||
        !is.finite(object$phi) || object$phi <= 0) {
      stop(
        sprintf(
          "'%s' residuals require a finite positive dispersion estimate",
          type
        ),
        call. = FALSE
      )
    }
  }
  ans <- switch(
    type,
    working = object$residuals,
    response = object$residuals,
    pearson = {
      pr <- get_pearson_residuals(object$family$family, y, mu, weights)
      as.numeric(pr) / sqrt(object$phi)
    },
    deviance = {
      dev_contribution <- object$family$dev.resids(y, mu, weights)
      if (any(!is.finite(dev_contribution)) || any(dev_contribution < 0)) {
        stop(
          paste0(
            "'deviance' residuals require finite non-negative deviance ",
            "contributions from the fitted family"
          ),
          call. = FALSE
        )
      }
      sign(y - mu) * sqrt(dev_contribution / object$phi)
    },
    mahalanobis = compute_mahalanobis_residuals(object)
  )
  ans
}
