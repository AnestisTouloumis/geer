evaluate_prediction_offset <- function(object, newdata, n_obs) {
  offset_formula <- rep.int(0, n_obs)
  pred_terms <- delete.response(object$terms)
  mf <- model.frame(pred_terms, newdata, xlev = object$xlevels, na.action = na.pass)
  offset_from_formula <- model.offset(mf)
  if (!is.null(offset_from_formula)) {
    if (length(offset_from_formula) == 1L) {
      offset_from_formula <- rep.int(as.numeric(offset_from_formula), n_obs)
    }
    if (!is.numeric(offset_from_formula) || length(offset_from_formula) != n_obs ||
        anyNA(offset_from_formula) || any(!is.finite(offset_from_formula))) {
      stop("formula-based prediction offset must be a finite numeric vector of the correct length",
           call. = FALSE)
    }
    offset_formula <- as.numeric(offset_from_formula)
  }
  offset_argument <- rep.int(0, n_obs)
  call_offset <- object$call$offset
  if (!is.null(call_offset)) {
    offset_env <- environment(formula(object))
    offset_from_argument <- tryCatch(
      eval(call_offset, envir = newdata, enclos = offset_env),
      error = function(e) {
        stop(
          "could not evaluate the offset supplied via the 'offset' argument in 'newdata'",
          call. = FALSE
        )
      }
    )
    if (length(offset_from_argument) == 1L) {
      offset_from_argument <- rep.int(as.numeric(offset_from_argument), n_obs)
    }
    if (!is.numeric(offset_from_argument) || length(offset_from_argument) != n_obs ||
        anyNA(offset_from_argument) || any(!is.finite(offset_from_argument))) {
      stop(
        "prediction offset supplied via the 'offset' argument must be a finite numeric vector of the correct length",
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
#' @aliases fitted fitted.values
#' @method fitted geer
#'
#' @description
#' Extracts the fitted mean values on the response scale from a fitted
#' \code{geer} object. \code{fitted.values} is an alias for \code{fitted}.
#'
#' @inheritParams coef.geer
#'
#' @return
#' A numeric vector of fitted mean values on the response scale, of the same
#' length as the number of observations used in fitting.
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
  object <- check_geer_object(object)
  object$fitted.values
}


#' @title
#' Predictions from a geer Object
#'
#' @aliases predict predict.geer
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
#'   degrees-of-freedom adjusted estimator (\code{"df-adjusted"}), and the
#'   model-based or naive estimator (\code{"naive"}). Defaults to
#'   \code{"bias-corrected"}.
#' @param se.fit logical indicating whether approximate standard errors are to
#'   be returned. Defaults to \code{FALSE}.
#' @param ... additional arguments passed to or from other methods.
#'
#' @details
#' Predictions are obtained by computing the model matrix for \code{newdata}
#' (or using the original fit when \code{newdata} is omitted) and multiplying
#' by the estimated coefficients. If \code{type = "response"}, the linear
#' predictor is transformed via the inverse link function.
#'
#' When \code{se.fit = TRUE}, approximate standard errors are computed by the
#' delta method using the covariance matrix specified by \code{cov_type}. On
#' the response scale, standard errors are additionally scaled by the absolute
#' derivative of the inverse link function.
#'
#' @return
#' If \code{se.fit = FALSE}, a numeric vector of predictions on the scale
#' specified by \code{type}. If \code{se.fit = TRUE}, a list with components:
#' \item{fit}{numeric vector of predictions.}
#' \item{se.fit}{numeric vector of approximate standard errors of the
#'   predictions.}
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
                         cov_type = c("bias-corrected", "robust", "df-adjusted", "naive"),
                         se.fit = FALSE,
                         ...) {
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
    vcov_matrix <- vcov(object, cov_type = cov_type)
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
    ans <- list(fit = out, se.fit = se)
    return(ans)
  }
  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }
  pred_terms <- delete.response(object$terms)
  mf <- model.frame(pred_terms, newdata, xlev = object$xlevels, na.action = na.pass)
  design_matrix <- model.matrix(pred_terms, mf, contrasts.arg = object$contrasts)
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
    return(if (type == "link") eta_vector else object$family$linkinv(eta_vector))
  }
  vcov_matrix <- vcov(object, cov_type = cov_type)
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


#' @title
#' Residuals from a geer Object
#'
#' @aliases resid residuals residuals.geer
#' @method residuals geer
#'
#' @description
#' Extracts residuals of different types from a fitted \code{geer} object.
#'
#' @inheritParams coef.geer
#' @param type character string specifying the type of residuals to return.
#'   Options are \code{"working"} for raw residuals, \code{"pearson"} for
#'   residuals standardized by the variance function, and \code{"deviance"}
#'   for signed square roots of the deviance contributions. Defaults to
#'   \code{"working"}.
#'
#' @details
#' Residuals are computed using the marginal variance and deviance functions
#' of the family specified in the fitted model.
#'
#' @return
#' A numeric vector of residuals of the requested \code{type}, of the same
#' length as the number of observations used in fitting.
#'
#' @seealso \code{\link{fitted.geer}}, \code{\link{predict.geer}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}},
#'   \code{\link[stats]{residuals}}.
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
#'
#' @export
residuals.geer <- function(object,
                           type = c("working", "pearson", "deviance"),
                           ...) {
  object <- check_geer_object(object)
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  weights <- object$prior.weights
  ans <- switch(
    type,
    working = object$residuals,
    pearson = as.numeric(get_pearson_residuals(object$family$family, y, mu, weights)),
    deviance = {
      if (object$df.residual > 0) {
        dr <- sqrt(pmax(object$family$dev.resids(y, mu, weights), 0))
        sign_term <- ifelse(y > mu, 1, ifelse(y < mu, -1, 0))
        dr * sign_term
      } else {
        rep.int(0, length(mu))
      }
    }
  )
  ans
}


