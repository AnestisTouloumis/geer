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
  object$residuals <- x$residuals
  object$fitted_values <- x$fitted_values
  object$linear_predictors <- x$linear_predictors
  object$obs_no <- x$obs_no
  object$clusters_no <- x$clusters_no
  object$min_cluster_size <- x$min_cluster_size
  object$max_cluster_size <- x$max_cluster_size
  object$method <- x$method
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
summary.geer <- function(object, type = "robust", ...) {
  coefficients_hat <- coef(object)
  standard_errors <- sqrt(diag(vcov(object, type)))
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
    type = type
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
  if (x$type == "bias_corrected") x$type <- "bias-corrected"
  cat("Std.Errors are taken from the", x$type, "covariance matrix.", "\n")
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
  object$fitted_values
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
  mu  <- object$fitted_values
  res <- switch(type,
                pearson = r/sqrt(object$family$variance(mu)),
                working = r)
  res
}
