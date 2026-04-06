#' @title
#' Print a \code{geer} Object
#'
#' @description
#' Print the stored call, estimated regression coefficients, and basic fitting
#' information for a fitted \code{geer} object.
#'
#' @param x an object of class \code{"geer"}.
#' @param ... additional arguments passed to or from other methods.
#'
#' @return
#' The input object \code{x} is returned invisibly.
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
#' Summarize a \code{geer} Object
#'
#' @description
#' Produce a coefficient table and basic model summary information for a fitted
#' \code{geer} object.
#'
#' @param object a fitted model object of class \code{geer}.
#' @param cov_type character string specifying the covariance estimator used to
#'        compute standard errors, z-statistics, and p-values. Options are
#'        \code{"robust"}, \code{"bias-corrected"}, \code{"df-adjusted"}, and
#'        \code{"naive"}. Default is \code{"robust"}.
#' @param ... additional arguments passed to or from other methods. Currently
#'        unused.
#'
#' @return
#' An object of class \code{"summary.geer"} containing the coefficient table
#' and basic fitting information.
#'
#' @export
summary.geer <- function(object,
                         cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
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
#' Print a \code{summary.geer} Object
#'
#' @param x An object of class \code{summary.geer}.
#' @param ... Additional arguments passed to or from other methods. Currently
#'        unused.
#'
#' @return
#' The input object \code{x} is returned invisibly.
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
  invisible(x)
}
