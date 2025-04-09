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
  object$df.residuals <- x$df.residuals
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

#' @method model.matrix geer
#' @export
model.matrix.geer <-	function(object,...){
  model.matrix(object = object$terms, data = object$data, contrasts = object$contrasts)
  }


#' @title
#' Extract Model Residuals from a \code{geer} Object
#'
#' @aliases resid residuals residuals.geer
#'
#' @method residuals geer
#'
#' @description
#' Extract model residuals from \code{geer} objects.
#'
#' @param object an object representing a model of the class \code{geer}.
#' @param type character indicating whether the type of the residuals to be
#'             returned. Options include the working residuals (\code{type =
#'             "working"}), the pearson residuals (\code{type = "pearson"}) and
#'             the deviance residuals (\code{type = "deviance"}). By default, the
#'             \code{"working"} residuals are returned.
#' @param ... other arguments.
#'
#' @details
#' If \code{type = "working"}, then the raw residuals (\code{observed - fitted})
#' are returned.
#'
#' If \code{type = "pearson"}, then the pearson residuals are returned. The
#' marginal distribution of the responses is defined by \code{object$family}.
#'
#' If \code{type = "deviance"}, then the deviance residuals are returned. The
#' marginal distribution of the responses is defined by \code{object$family}.
#'
#' @return
#' A vector with the observed residuals of type \code{type}.
#'
#' @export
residuals.geer <- function(object,
                           type = "working",
                           ...) {
  icheck <- pmatch(type,
                   c("working", "pearson", "deviance"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown type for the residuals")
  response_vector <- object$y
  raw_residuals <- object$residuals
  mu_vector <- object$fitted.values
  weight_vector <- object$weights
  ans <- switch(type,
                deviance = if (object$df.residual > 0) {
                  deviance_res <-
                    sqrt(pmax((object$family$dev.resids)(response_vector, mu_vector, weight_vector), 0))
                  ifelse(response_vector > mu_vector, deviance_res, -deviance_res)
                  } else rep.int(0, length(mu_vector)),
                pearson =
                  c(get_pearson_residuals(object$family$family, response_vector,
                                        mu_vector, weight_vector)),
                working = raw_residuals
                )
  ans
}



#' @title
#' Confidence Intervals for Model Parameters from a \code{geer} Object
#'
#' @aliases confint confint.geer
#' @method confint geer
#'
#'
#' @description
#' Confidence intervals for one or more parameters in a fitted model.
#'
#'
#' @inheritParams anova.geer
#' @inheritParams stats::confint
#' @param ... additional argument(s) for methods.
#'
#'
#' @details
#' The references in \code{\link{vcov}} include the formulae for the covariance type
#' implied by \code{cov_type}.
#'
#'
#' @inherit stats::confint.default return
#'
#'
#' @export
confint.geer <- function(object, parm, level = 0.95, cov_type = "robust", ...) {
  coefficients <- coef(object)
  coefficients_names <- names(coefficients)
  if (missing(parm)) {
    parm <- coefficients_names
  } else if (is.numeric(parm)) {
    parm <- coefficients_names[parm]
  }
  alpha <- (1 - level) / 2
  alpha <- c(alpha, 1 - alpha)
  pct <- format_perc(alpha, 3)
  percentiles <- qnorm(alpha)
  ans <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  standard_errors <- sqrt(diag(vcov(object, cov_type = cov_type)))[parm]
  ans[] <- coefficients[parm] + standard_errors %o% percentiles
  ans
}



#' @title
#' Model Predictions for a \code{geer} Object
#'
#' @aliases predict predict.geer
#' @method predict geer
#'
#'
#' @description
#' Obtains predictions and optionally estimates standard errors of those
#' predictions from a fitted generalized estimating equations object of the class
#' \code{geer}.
#'
#'
#' @param object an object representing a model of the class \code{geer}.
#' @param newdata	optionally, a data frame in which to look for variables with
#'        which to predict. If omitted, the fitted linear predictors are used.
#' @param type the type of prediction required. The default
#'        \code{"link"} is on the scale of the linear predictors;
#'        the alternative \code{"response"} is on the scale of the response variable.
#' @param se.fit logical switch indicating if standard errors are required.
#'        By default, the standard errors are not required.
#' @param cov_type an optional character string indicating the type of estimator
#'        which should be used to the variance-covariance matrix of the interest
#'        parameters. Options include the sandwich or robust estimator
#'        (\code{"robust"}), the bias-corrected estimator (\code{"bias-corrected"}),
#'        the degrees of freedom adjusted estimator (\code{"bias-corrected"})
#'        and the model-based or naive estimator (\code{"naive"}). By default,
#'        the robust covariance estimator is used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' If \code{newdata} is omitted the predictions are based on the data used for
#' the fit.
#'
#' @return
#' If \code{se.fit = FALSE}, a vector or matrix of predictions.
#'
#' If \code{se.fit = TRUE}, a list with components
#' \item{fit}{Predictions, as for \code{se.fit = FALSE}.}
#' \item{se.fit}{Estimated standard errors.}
#'
#'
#' @export
predict.geer <- function(object, newdata = NULL, type = "link",
                         cov_type = "robust", se.fit = FALSE, ...){
  icheck <- pmatch(type,
                   c("link", "response"),
                   nomatch = 0)
  if (icheck == 0) stop("unknown type of prediction")
  icheck <- pmatch(cov_type,
                   c("robust", "df-adjusted", "bias-corrected", "naive"),
                   nomatch = 0)
  if (icheck == 0) stop("unknown covariance matrix")
  if (missing(newdata)) {
    if (type == "link") {
      predicts <- object$linear.predictors
    } else {
      predicts <- object$fitted.values
    }
    if (se.fit) {
      covariance_matrix <- vcov(object, type = cov_type)
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
      covariance_matrix <- vcov(object, cov_type)
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
