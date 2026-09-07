check_prediction_offset_variables <- function(expression, newdata, source) {
  ## Offset variables are resolved in 'newdata' only. Without this check R
  ## falls back to the environment in which the model was fitted, where the
  ## original offset still lives, and an offset of the right length is then
  ## silently accepted for the wrong observations.
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
#' @param ... not used. Supplying any argument here is an error rather than
#'   being silently ignored.
#'
#' @return
#' A numeric vector of fitted mean values on the response scale, of the same
#' length as the number of observations used in fitting, in the internal
#' observation order described below rather than in the row order of the
#' supplied data.
#'
#' @section Observation order:
#' The returned vector follows the internal observation order of the fitted
#' model, not the row order of the data supplied to \code{geewa} or
#' \code{geewa_binary}. Rows are restricted to the complete cases retained by
#' the \code{na.action} in force and are then sorted by cluster identifier
#' and, within cluster, by the repeated index. Unless the supplied data were
#' already sorted that way and free of missing values, assigning the result
#' back with an expression such as \code{data$r <- residuals(fit)} will
#' attach values to the wrong rows without any warning, because the lengths
#' can still agree. Order the data by cluster and repeated index before
#' fitting, or align the result yourself using \code{fit$id} and
#' \code{fit$repeated}, which are stored in the same order as the returned
#' vector.
#'
#' @seealso \code{\link{residuals.geer}}, \code{\link{predict.geer}}.
#'
#' @examples
#' data("leprosy", package = "geer")
#' fit <- geewa(
#'   formula = bacilli ~ factor(period) + factor(period):treatment,
#'   family = poisson(link = "log"),
#'   data = leprosy,
#'   id = id
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
#' @param newdata optional data frame in which to look for variables used for
#'   prediction. If omitted, predictions are made on the data used for fitting.
#' @param type character string specifying the scale of the predictions.
#'   Options are \code{"link"} for the linear predictor scale and
#'   \code{"response"} for the response scale. Defaults to \code{"link"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used to compute approximate standard errors when \code{se.fit = TRUE}.
#'   Options are the bias-corrected estimator (\code{"bias-corrected"}),
#'   the sandwich or robust estimator (\code{"robust"}), the
#'   degrees-of-freedom adjusted estimator (\code{"df-adjusted"}), the
#'   leave-one-cluster jackknife estimator (\code{"jackknife"}), and the
#'   model-based or naive estimator (\code{"naive"}). Defaults to
#'   \code{"bias-corrected"}.
#' @param se.fit logical indicating whether approximate standard errors are to
#'   be returned. Defaults to \code{FALSE}.
#' @param ... not used. Supplying any argument here is an error, so that a
#'   misspelt \code{type} or \code{cov_type} is reported rather than
#'   silently ignored.
#'
#' @details
#' Predictions are obtained by computing the model matrix for \code{newdata}
#' (or using the original fit when \code{newdata} is omitted), multiplying by
#' the estimated coefficients, and adding any offset. If
#' \code{type = "response"}, the linear predictor is transformed via the
#' inverse link function.
#'
#' When \code{se.fit = TRUE}, approximate standard errors are computed by the
#' delta method using the covariance matrix specified by \code{cov_type}. On
#' the response scale, standard errors are additionally scaled by the absolute
#' derivative of the inverse link function.
#'
#' When \code{newdata} is supplied and the fitted model has an offset, whether
#' through \code{offset()} in the formula or through the \code{offset}
#' argument, every variable the offset depends on must be present in
#' \code{newdata}. An error is signalled otherwise. Without that requirement
#' the offset would be resolved in the environment in which the model was
#' fitted, so an offset of the right length taken from the original data would
#' be accepted silently for a different set of observations.
#'
#' @return
#' If \code{se.fit = FALSE}, a numeric vector of predictions on the scale
#' specified by \code{type}. If \code{se.fit = TRUE}, a list with components:
#' \item{fit}{numeric vector of predictions.}
#' \item{se.fit}{numeric vector of approximate standard errors of the
#'   predictions.}
#'
#' When \code{newdata} is supplied both the predictions and their standard
#' errors carry the row names of \code{newdata}. When it is omitted both are
#' unnamed, because the row names stored with the fitted model refer to the
#' data before they were sorted into the internal observation order and would
#' therefore identify the wrong rows.
#'
#' @section Observation order:
#' When \code{newdata} is supplied the predictions follow the row order of
#' \code{newdata}. When it is omitted they follow the internal observation
#' order of the fitted model instead: the complete cases retained by the
#' \code{na.action} in force, sorted by cluster identifier and, within
#' cluster, by the repeated index. Unless the data supplied to \code{geewa}
#' or \code{geewa_binary} were already sorted that way and free of missing
#' values, assigning the result back to those data will attach values to the
#' wrong rows without any warning, because the lengths can still agree. Pass
#' the original data as \code{newdata} to obtain predictions in the row order
#' of the data.
#'
#' Rows of \code{newdata} are not dropped when a variable used by the model is
#' missing. The model frame is built with \code{na.action = na.pass}, so such
#' a row yields \code{NA} for the prediction, and for its standard error when
#' \code{se.fit = TRUE}, and the returned vector keeps one element per row of
#' \code{newdata}.
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
#'   orstr = "exchangeable"
#' )
#' head(predict(fit, type = "link"))
#' head(predict(fit, type = "response"))
#'
#' nd <- cerebrovascular[1:5, , drop = FALSE]
#' predict(fit, newdata = nd, type = "response")
#'
#' pred <- predict(fit, type = "response", se.fit = TRUE, cov_type = "robust")
#' head(pred$fit)
#' head(pred$se.fit)
#'
#' @export
predict.geer <- function(object,
                         newdata = NULL,
                         type = c("link", "response"),
                         cov_type = geer_cov_type_choices,
                         se.fit = FALSE,
                         ...) {
  check_unused_dots(list(...), "predict.geer")
  object <- check_geer_object(object)
  type <- match.arg(type)
  cov_type <- match.arg(cov_type)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    stop("'se.fit' must be a single logical value", call. = FALSE)
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
    se <- sqrt(rowSums((design_matrix %*% vcov_matrix) * design_matrix))
    if (type == "response") {
      se <- se * abs(object$family$mu.eta(eta_vector))
    }
    ## 'object$x' keeps the model-matrix row names of the data before sorting,
    ## which no longer identify the rows of the stored fit, so they are
    ## dropped to match the unnamed 'fit' component.
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
  se <- sqrt(rowSums((design_matrix %*% vcov_matrix) * design_matrix))
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
      ## A Cholesky factor both solves the system and certifies that the
      ## working covariance is positive definite. The quadratic form is then
      ## the squared norm of the forward-solved residual, so it is nonnegative
      ## by construction and needs no clamping.
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
#' @param type character string specifying the type of residuals to return.
#'   Options are \code{"working"} for raw residuals, \code{"pearson"} for
#'   residuals standardized by the marginal variance, \code{"deviance"}
#'   for dispersion-scaled signed square roots of the deviance contributions,
#'   and \code{"mahalanobis"} for cluster-level Mahalanobis residuals.
#'   Defaults to \code{"working"}.
#' @param ... not used. Supplying any argument here is an error, so that a
#'   misspelt \code{type} is reported rather than silently ignored.
#'
#' @details
#' Pearson residuals are computed as
#' \deqn{r^P_{ij} = \frac{y_{ij} - \hat\mu_{ij}}{\sqrt{\hat\phi
#' V(\hat\mu_{ij}) / w_{ij}}}.}
#' Deviance residuals are computed as
#' \deqn{r^D_{ij} = \mathrm{sign}(y_{ij} - \hat\mu_{ij})
#' \sqrt{d(y_{ij}, \hat\mu_{ij}, w_{ij}) / \hat\phi},}
#' where \eqn{d(y_{ij}, \hat\mu_{ij}, w_{ij})} is the non-scaled deviance
#' contribution returned by the fitted family. As in
#' \code{\link[stats]{glm}}, that contribution already includes the prior
#' weight \eqn{w_{ij}}, which is why the weight appears inside \eqn{d} here
#' but explicitly in the denominator of the Pearson residual.
#'
#' Mahalanobis residuals are computed once per cluster as
#' \deqn{r^M_i = \frac{1}{n_i}(y_i - \hat\mu_i)^T
#' \widehat{\operatorname{Var}}(Y_i)^{-1}(y_i - \hat\mu_i),}
#' using the fitted working covariance matrix for that cluster, which already
#' includes the estimated dispersion. They are therefore nonnegative
#' cluster-level diagnostics, rather than observation-level signed residuals.
#' The quadratic form is evaluated through a Cholesky factorization of the
#' working covariance matrix, and an error is signalled for any cluster whose
#' working covariance matrix is not positive definite.
#'
#' @section Observation order:
#' The returned vector follows the internal observation order of the fitted
#' model, not the row order of the data supplied to \code{geewa} or
#' \code{geewa_binary}. Rows are restricted to the complete cases retained by
#' the \code{na.action} in force and are then sorted by cluster identifier
#' and, within cluster, by the repeated index. Unless the supplied data were
#' already sorted that way and free of missing values, assigning the result
#' back with an expression such as \code{data$r <- residuals(fit)} will
#' attach values to the wrong rows without any warning, because the lengths
#' can still agree. Order the data by cluster and repeated index before
#' fitting, or align the result yourself using \code{fit$id} and
#' \code{fit$repeated}, which are stored in the same order as the returned
#' vector.
#'
#' @return
#' A numeric vector. For \code{type = "working"}, \code{"pearson"}, or
#' \code{"deviance"}, the vector has one value per fitted observation, in the
#' internal observation order described below rather than in the row order of
#' the supplied data. For \code{type = "mahalanobis"}, it has one value per
#' cluster and is named by the cluster identifier.
#'
#' @references
#' Vanegas, L.H., Rondon, L.M. and Paula, G.A. (2023) Generalized Estimating
#' Equations using the new R package glmtoolbox. \emph{The R Journal},
#' \bold{15}, 105--133.
#'
#' @seealso \code{\link{fitted.geer}}, \code{\link{predict.geer}},
#'   \code{\link{runs_test}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}, \code{\link[stats]{residuals}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' head(residuals(fit, type = "working"))
#' head(residuals(fit, type = "pearson"))
#' head(residuals(fit, type = "deviance"))
#' head(residuals(fit, type = "mahalanobis"))
#'
#' @export
residuals.geer <- function(object,
                           type = c(
                             "working", "pearson", "deviance", "mahalanobis"
                           ),
                           ...) {
  check_unused_dots(list(...), "residuals.geer")
  object <- check_geer_object(object)
  type <- match.arg(type)
  if (type %in% c("pearson", "deviance")) {
    ## object$y carries the model-frame row names while $fitted.values and
    ## $residuals do not, so the names are dropped to keep every
    ## observation-level residual type unnamed.
    y <- as.numeric(object$y)
    mu <- object$fitted.values
    weights <- object$prior.weights
  }
  if (!identical(type, "working")) {
    ## Every scaled residual divides by the estimated dispersion, so an
    ## unusable estimate is reported rather than propagated.
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
    pearson = {
      pr <- get_pearson_residuals(object$family$family, y, mu, weights)
      as.numeric(pr) / sqrt(object$phi)
    },
    deviance = {
      ## A negative deviance contribution is not a rounding artifact: it means
      ## the family has returned something the deviance cannot produce, so it
      ## is reported rather than clamped to zero.
      dev_contribution <- object$family$dev.resids(y, mu, weights)
      if (any(!is.finite(dev_contribution)) || any(dev_contribution < 0)) {
        stop(
          paste0(
            "'deviance' residuals require finite nonnegative deviance ",
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


