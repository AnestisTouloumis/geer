#' @title
#' Print a geer Object
#'
#' @description
#' Prints the call, estimated regression coefficients, and basic fitting
#' information for a fitted \code{geer} object.
#'
#' @param x a fitted model object of class \code{"geer"}.
#' @param ... additional arguments passed to or from other methods.
#'
#' @return
#' The input object \code{x} is returned invisibly.
#'
#' @seealso \code{\link{summary.geer}}, \code{\link{coef.geer}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   corstr = "exchangeable"
#' )
#' print(fit)
#'
#' @export
print.geer <- function(x, ...) {
  x <- check_geer_object(x, "x")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  if (!is.null(x$coefficients)) {
    print(x$coefficients)
  } else {
    cat("<none>\n")
  }
  cat("\nNumber of iterations :", x$iter, "\n")
  cat("Algorithm converged  :", x$converged, "\n")
  invisible(x)
}


#' @title
#' Summarize a geer Object
#'
#' @description
#' Produces a coefficient table and basic model summary information for a
#' fitted \code{geer} object.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param cov_type character string specifying the covariance estimator used to
#'   compute standard errors, z-statistics, and p-values. Options are \code{"bias-corrected"},
#'   \code{"robust"}, \code{"df-adjusted"}, and
#'   \code{"naive"}. Defaults to \code{"bias-corrected"}.
#' @param ... additional arguments passed to or from other methods. Currently
#'   unused.
#'
#' @return
#' An object of class \code{"summary.geer"}, a list with components:
#' \item{coefficients}{a matrix with columns \code{Estimate},
#'   \code{Std. Error}, \code{z value}, and \code{Pr(>|z|)}, with one row per
#'   regression parameter.}
#' \item{family}{the \code{\link[stats]{family}} object used.}
#' \item{alpha}{the estimated or fixed association parameters.}
#' \item{call}{the matched call.}
#' \item{residuals}{the working residuals.}
#' \item{iter}{the number of iterations used.}
#' \item{converged}{logical indicating whether the algorithm converged.}
#' \item{phi}{the estimated or fixed scale parameter.}
#' \item{association_structure}{the name of the working association structure.}
#' \item{method}{character string identifying the estimation method used.}
#' \item{cov_type}{the covariance estimator used for standard errors.}
#'
#' @seealso \code{\link{print.summary.geer}}, \code{\link{tidy.geer}},
#'   \code{\link{glance.geer}}, \code{\link{vcov.geer}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   corstr = "exchangeable"
#' )
#' summary(fit)
#' summary(fit, cov_type = "bias-corrected")
#'
#' data("cerebrovascular", package = "geer")
#' fit2 <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' summary(fit2)
#'
#' @export
summary.geer <- function(object,
                         cov_type = c("bias-corrected", "robust", "df-adjusted", "naive"),
                         ...) {
  object <- check_geer_object(object)
  cov_type <- match.arg(cov_type)
  beta <- coef(object)
  vcov_matrix <- vcov(object, cov_type = cov_type)
  se <- sqrt(pmax(0, diag(vcov_matrix)))
  z_stat <- rep.int(NA_real_, length(beta))
  pval <- rep.int(NA_real_, length(beta))
  ok <- is.finite(beta) & is.finite(se) & se > 0
  z_stat[ok] <- beta[ok] / se[ok]
  pval[ok] <- 2 * pnorm(abs(z_stat[ok]), lower.tail = FALSE)
  coef_table <- cbind(
    Estimate = beta,
    `Std. Error` = se,
    `z value` = z_stat,
    `Pr(>|z|)` = pval
  )
  res <- list(
    coefficients = coef_table,
    family = object$family,
    alpha = object$alpha,
    call = object$call,
    residuals = object$residuals,
    iter = object$iter,
    converged = object$converged,
    phi = object$phi,
    association_structure = object$association_structure,
    method = object$method,
    cov_type = cov_type
  )
  class(res) <- "summary.geer"
  res
}


#' @title
#' Print a summary.geer Object
#'
#' @description
#' Prints the contents of a \code{summary.geer} object, including the call,
#' estimation method, coefficient table, dispersion parameter, and working
#' association structure.
#'
#' @param x an object of class \code{"summary.geer"}.
#' @param ... additional arguments passed to or from other methods. Currently
#'   unused.
#'
#' @return
#' The input object \code{x} is returned invisibly.
#'
#' @seealso \code{\link{summary.geer}}, \code{\link{print.geer}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   corstr = "exchangeable"
#' )
#' print(summary(fit))
#'
#' @export
print.summary.geer <- function(x, ...) {
  x <- check_summary_geer_object(x)
  cat("\nCall:\n")
  print(x$call)
  cat("\nEstimating Method   :", x$method, "\n")
  cat("Number of iterations:", x$iter, "\n")
  cat("Algorithm converged :", x$converged, "\n")
  cat("\nMarginal Model\n")
  cat("Family       :", x$family$family, "\n")
  cat("Link Function:", x$family$link, "\n")
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  cat("Std. Errors are taken from the", x$cov_type, "covariance matrix.", "\n")
  cat("\nDispersion Parameter:", round(x$phi, digits = 4), "\n")
  cat("\nAssociation Structure:", x$association_structure, "\n")
  if (length(x$alpha) == 1) {
    cat("Association Parameter:", round(x$alpha, digits = 4), "\n")
  } else {
    cat("Association Parameters:\n")
    print(round(x$alpha, digits = 4))
    cat("\n")
  }
  invisible(x)
}
