geer <- function(x, ...) { # nolint
  UseMethod("geer")
}

#' @export
geer.default <- function(x, ...) {
  object <- list()
  object$call <- x$call
  object$coefficients <- x$coefficients
  object$phi <- x$phi
  object$association_structure <- x$association_structure
  object$alpha <- x$alpha
  object$naive_covariance <- x$naive_covariance
  object$robust_covariance <- x$robust_covariance
  object$bias_corrected_covariance <- x$bias_corrected_covariance
  object$ee_value <- x$ee_value
  object$converged <- x$converged
  object$niter <- x$niter
  object$criterion <- x$criterion
  object$terms <- x$terms
  object$family <- x$family
  object$y <- x$y
  object$model_matrix <- x$model_matrix
  object$id <- x$id
  object$repeated <- x$repeated
  object$residuals <- x$residuals
  object$fitted.values <- x$fitted.values
  object$linear.predictors <- x$linear.predictors
  object$obs_no <- x$obs_no
  object$clusters_no <- x$clusters_no
  object$min_cluster_size <- x$min_cluster_size
  object$max_cluster_size <- x$max_cluster_size
  object$method <- x$method
  object$levels <- x$levels
  object$contrasts <- x$contrasts
  class(object) <- "geer"
  object
}


#' @export
print.geer <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nNumber of iterations :", x$niter, "\n")
  cat("Algorithm converged  :", x$converged, "\n")
}


#' @method summary geer
#' @export
summary.geer <- function(object, cov_type = "robust", ...) {
  coefficients_hat <- coef(object)
  standard_errors <- sqrt(diag(vcov(object, cov_type)))
  z_statistics  <- coefficients_hat / standard_errors
  pvalue <- 2 * (1 - pnorm(abs(z_statistics)))
  TAB <- cbind(
    Estimate = coefficients_hat,
    `Std.Error` = standard_errors,
    `z value` = z_statistics,
    `Pr(>|z|)` = pvalue
  )
  res <- list(
    coefficients = TAB,
    family = object$family,
    alpha = object$alpha,
    call = object$call,
    residuals = object$residuals,
    niter = object$niter,
    converged = object$converged,
    phi = object$phi,
    association_structure = object$association_structure,
    method = object$method,
    cov_type = cov_type
  )
  class(res) <- "summary.geer"
  res
}


#' @export
print.summary.geer <- function(x, ...) {
  cat("\ncall:\n")
  print(x$call)
  cat("\nEstimating Method   :", x$method, "\n")
  cat("Number of iterations:", x$niter, "\n")
  cat("Algorithm converged :", x$converged, "\n")
  cat("\nMarginal Model\n")
  cat("Family       :", x$family$family, "\n")
  cat("Link Function:", x$family$link, "\n")
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  if (x$cov_type == "bias_corrected") x$cov_type <- "bias-corrected"
  cat("Std.Errors are taken from the", x$cov_type, "covariance matrix.", "\n")
  cat("\nDispersion Parameter:", round(x$phi, digits = 4), "\n")
  cat("\nAssociation Structure:", x$association_structure, "\n")
  if (length(x$alpha) == 1) {
    cat("Association Parameter:", round(x$alpha, digits = 4), "\n")
  } else {
    cat("Association Parameters:\n")
    print(round(x$alpha, digits = 4))
    cat("\n")
  }
}


#' @aliases coefficients
#' @method coef geer
#' @export
coef.geer <- function(object, ...){
  coeffs <- object$coefficients
  coeffs
 }

#' @aliases fitted.values
#' @method fitted geer
#' @export
fitted.geer <- function(object, ...){
  object$fitted.values
}

#' @method confint geer
#' @export
confint.geer <- function(object, parm, level = 0.95, type = "robust", ...) {
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- format_perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov(object, type = type)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}


#' @method residuals geer
#' @export
residuals.geer <- function(object, type = c("working", "pearson"), ...) {
  type <- match.arg(type)
  r   <- object$residuals
  mu  <- object$fitted.values
  res <- switch(type,
                pearson = r/sqrt(object$family$variance(mu)),
                working = r)
  res
}


#' @title Predict Method for Generalized Estimating Equations Fits
#' @description Obtains predictions and optionally estimates standard errors of those predictions from a fitted generalized estimating equations object.
#' @param object a fitted object of class inheriting from \code{geewa} or \code{geewa_binary}.
#' @param newdata	optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param type the type of prediction required. The default \code{"link"} is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable.
#' @param se.fit logical switch indicating if standard errors are required. As default, \code{se.fit = FALSE}.
#' @param covtype an optional character string indicating the type of estimator which should be used to the variance-covariance matrix of the interest parameters. The available options are: robust or sandwich estimator ("robust"), bias-corrected estimator ("bias-corrected"), and the model-based or naive estimator ("model"). As default, \code{covtype = "robust"}.
#' @param ... further arguments passed to or from other methods.
#' @return If \code{se.fit = FALSE}, a vector or matrix of predictions.
#'
#' If \code{se.fit = TRUE}, a list with components
#' \item{fit}{Predictions, as for \code{se.fit = FALSE}.}
#' \item{se.fit}{Estimated standard errors.}
#' @method predict geer
#' @export
predict.geer <- function(object, newdata = NULL, type = c("link", "response"),
                         covtype = c("robust", "bias-corrected", "naive"),
                         se.fit = FALSE, ...){
  type <- match.arg(type)
  covtype <- match.arg(covtype)
  if (missing(newdata)) {
    if (type == "link") {
      predicts <- object$linear.predictors
    } else {
      predicts <- object$fitted.values
    }
    if (se.fit) {
      covariance_matrix <- vcov(object, type = covtype)
      model_matrix <- object$model_matrix
      se.fit <-
        sqrt(apply(tcrossprod(model_matrix, covariance_matrix) * model_matrix, 1, sum))
      if (type == "response")
        se.fit <- se.fit * abs(object$family$mu.eta(object$linear.predictors))
      names(predicts) <- names(se.fit)
      predicts <- list(fit = predicts, se.fit = se.fit)
    }
  } else {
    newdata <- data.frame(newdata)
    model_frame <-
      model.frame(delete.response(object$terms), newdata, xlev = object$levels)
    model_matrix <-
      model.matrix(delete.response(object$terms), model_frame, contrasts = object$contrasts)
    predicts <- c(model_matrix %*% object$coefficients)
    offset_term <- model.offset(model_frame)
    if (!is.null(offset_term)) predicts <- predicts + c(offset_term)
    if (se.fit) {
      covariance_matrix <- vcov(object, covtype)
      se.fit <-
        sqrt(apply(tcrossprod(model_matrix, covariance_matrix) * model_matrix, 1, sum))
      if (type == "response") {
        se.fit <- se.fit * abs(object$family$mu.eta(predicts))
        predicts <- object$family$linkinv(predicts)
      }
      names(predicts) <- names(se.fit)
      predicts <- list(fit = predicts, se.fit = se.fit)
    } else if (type == "response") {
      predicts <- object$family$linkinv(predicts)
    }
  }
  predicts
}
