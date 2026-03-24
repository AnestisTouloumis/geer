#' @title
#' Extract Model Fitted Values from a \code{geer} Object
#'
#' @rdname fitted
#' @aliases fitted fitted.values
#' @method fitted geer
#'
#' @description
#' Extract fitted values (mean response on the response scale) from a fitted
#' \code{geer} model. \code{fitted.values()} is an alias for \code{fitted()}.
#'
#' @inheritParams coef.geer
#'
#' @return
#' A numeric vector of fitted values extracted from \code{object}.
#'
#' @examples
#' data("leprosy")
#' fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
#'              family = poisson(link = "log"), id = id, data = leprosy)
#' head(fitted(fit))
#'
#' @export
fitted.geer <- function(object, ...){
  object <- check_geer_object(object)
  object$fitted.values
}



#' @title
#' Predictions from a \code{geer} Object
#'
#' @aliases predict predict.geer
#' @method predict geer
#'
#' @description
#' Generate predictions, optionally with standard errors, from a fitted
#' \code{geer} object.
#'
#' @inheritParams add1.geer
#' @param newdata optional data frame in which to look for variables used
#'        for prediction. If omitted, predictions are made on the data used for
#'        fitting.
#' @param type type of prediction required. Options include the scale of the
#'        linear predictors (\code{"link"}) and the scale of the response
#'        variable (\code{"response"}). By default, the scale of the linear
#'        predictors is used.
#' @param se.fit logical indicating if standard errors are required.
#'        Defaults to \code{FALSE}.
#'
#' @details
#' Predictions are obtained by computing the model matrix for \code{newdata}
#' (or the original data if \code{newdata} is missing) and multiplying by the
#' estimated coefficients. If \code{type = "response"}, predictions are
#' transformed via the inverse link function.
#'
#' When \code{se.fit = TRUE}, approximate standard errors of predictions are
#' computed using the model-based covariance matrix of the estimated
#' coefficients.
#'
#' If \code{newdata} is omitted the predictions are based on the data used for
#' the fit.
#'
#' @return
#' If \code{se.fit = FALSE}, a numeric vector of predictions.
#' If \code{se.fit = TRUE}, a list with components \code{fit} and \code{se.fit}.
#'
#' @seealso
#' \code{\link{geewa}}, \code{\link{geewa_binary}}, \code{\link[stats]{predict}}.
#'
#' @examples
#' data("cerebrovascular")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   id = id,
#'   data = cerebrovascular,
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
                         cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                         se.fit = FALSE,
                         ...) {
  object <- check_geer_object(object)
  type <- match.arg(type)
  cov_type <- match.arg(cov_type)
  se.fit <- isTRUE(se.fit)
  coef_names <- names(object$coefficients)
  if (is.null(coef_names) || !length(coef_names)) {
    stop("'object$coefficients' must be named", call. = FALSE)
  }
  if (is.null(newdata)) {
    eta <- object$linear.predictors
    mu <- object$fitted.values
    out <- if (type == "link") eta else mu
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
      se <- se * abs(object$family$mu.eta(eta))
    }

    return(list(fit = out, se.fit = se))
  }

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  tt <- delete.response(object$terms)
  mf <- model.frame(tt, newdata, xlev = object$xlevels, na.action = na.pass)
  design_matrix <- model.matrix(tt, mf, contrasts.arg = object$contrasts)

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

  eta <- drop(design_matrix %*% object$coefficients)

  off <- model.offset(mf)
  if (!is.null(off)) {
    eta <- eta + drop(off)
  }

  if (!se.fit) {
    return(if (type == "link") eta else object$family$linkinv(eta))
  }

  vcov_matrix <- vcov(object, cov_type = cov_type)
  vcov_matrix <- vcov_matrix[coef_names, coef_names, drop = FALSE]

  se <- sqrt(rowSums((design_matrix %*% vcov_matrix) * design_matrix))

  if (type == "response") {
    mu <- object$family$linkinv(eta)
    se <- se * abs(object$family$mu.eta(eta))
    return(list(fit = mu, se.fit = se))
  }

  list(fit = eta, se.fit = se)
}



#' @title
#' Residuals from a \code{geer} Object
#'
#' @aliases resid residuals residuals.geer
#' @method residuals geer
#'
#' @description
#' Extract residuals of different types from a fitted \code{geer} object.
#'
#' @inheritParams coef.geer
#' @param type character indicating whether the type of the residuals to return.
#'        Options include the working residuals (\code{"working"}), the pearson
#'        residuals (\code{"pearson"}) and the deviance residuals
#'        (\code{"deviance"}). By default, the working residuals are returned.
#'
#' @details
#' The residuals are computed according to the marginal distribution specified
#' by \code{object$family$family}.
#'
#' \itemize{
#'   \item Working residuals: raw differences between observed and fitted values.
#'   \item Pearson residuals: standardized by the variance function of the
#'         specified family.
#'   \item Deviance residuals: signed square roots of the contribution of each
#'         observation to the model deviance.
#' }
#'
#' @return
#' A numeric vector of residuals of the requested \code{type}.
#'
#' @seealso
#' \code{\link{geewa}}, \code{\link{geewa_binary}},
#' \code{\link[stats]{residuals}}
#'
#' @examples
#' data("cerebrovascular")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   id = id,
#'   data = cerebrovascular,
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
  w <- object$prior.weights
  ans <- switch(
    type,
    working = object$residuals,
    pearson = as.numeric(get_pearson_residuals(object$family$family, y, mu, w)),
    deviance = {
      if (object$df.residual > 0) {
        dr <- sqrt(pmax(object$family$dev.resids(y, mu, w), 0))
        ifelse(y > mu, dr, -dr)
      } else {
        rep.int(0, length(mu))
      }
    }
  )
  ans
}




