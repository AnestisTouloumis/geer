print_nonconvergence_notice <- function(converged) {
  if (!isTRUE(converged)) {
    cat(
      "*** WARNING: the algorithm did not converge; the estimates, standard",
      "errors and p-values are unreliable. ***\n"
    )
  }
  invisible(NULL)
}


#' @title
#' Print a geer Object
#'
#' @description
#' Prints the call, the estimated regression coefficients, and the
#' convergence information of a fitted \code{geer} object.
#'
#' @param x a fitted model object of class \code{"geer"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
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
#'   repeated = visit,
#'   corstr = "exchangeable"
#' )
#' print(fit)
#'
#' @export
print.geer <- function(x, ...) {
  check_unused_dots(list(...), "print.geer")
  x <- check_geer_object(x, "x")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  if (!is.null(x$coefficients)) {
    print(x$coefficients)
  } else {
    cat("<none>\n")
  }
  cat("\nNumber of iterations:", x$iter, "\n")
  cat("Algorithm converged :", x$converged, "\n")
  print_nonconvergence_notice(x$converged)
  invisible(x)
}


#' @title
#' Summarize a geer Object
#'
#' @description
#' Summarizes a fitted \code{geer} object, with a table of the coefficients,
#' their standard errors, z-statistics, and p-values.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for the standard errors, z-statistics, and p-values: one of the
#'   estimators described in \code{\link{vcov.geer}}. Defaults to
#'   \code{"bias-corrected"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
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
#' \item{obs_no}{the number of observations used in the fit.}
#' \item{clusters_no}{the number of clusters used in the fit.}
#' \item{min_cluster_size}{the smallest observed cluster size.}
#' \item{max_cluster_size}{the largest observed cluster size.}
#' \item{iter}{the number of iterations used.}
#' \item{converged}{logical indicating whether the algorithm converged.}
#' \item{phi}{the estimated or fixed dispersion parameter, 1 for
#'   \code{\link{geewa_binary}} fits.}
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
#'   repeated = visit,
#'   corstr = "exchangeable"
#' )
#' summary(fit)
#' summary(fit, cov_type = "robust")
#'
#' @export
summary.geer <- function(object,
                         cov_type = "bias-corrected",
                         ...) {
  object <- check_geer_object(object)
  cov_type <- match_cov_type(cov_type)
  check_unused_dots(list(...), "summary.geer")
  beta <- stats::coef(object)
  vcov_matrix <- stats::vcov(object, cov_type = cov_type)
  se <- compute_se_or_stop(diag(vcov_matrix), context = "summary.geer")
  z_stat <- rep.int(NA_real_, length(beta))
  pval <- rep.int(NA_real_, length(beta))
  ok <- is.finite(beta) & is.finite(se) & se > 0
  z_stat[ok] <- beta[ok] / se[ok]
  pval[ok] <- 2 * stats::pnorm(abs(z_stat[ok]), lower.tail = FALSE)
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
    obs_no = object$obs_no,
    clusters_no = object$clusters_no,
    min_cluster_size = object$min_cluster_size,
    max_cluster_size = object$max_cluster_size,
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
#' Prints a \code{summary.geer} object: the call, the estimation method, the
#' sample and cluster sizes, the coefficient table, the dispersion parameter,
#' and the working association structure.
#'
#' @param x an object of class \code{"summary.geer"}, as returned by
#'   \code{\link{summary.geer}}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
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
#'   repeated = visit,
#'   corstr = "exchangeable"
#' )
#' print(summary(fit))
#'
#' @export
print.summary.geer <- function(x, ...) {
  check_unused_dots(list(...), "print.summary.geer")
  x <- check_summary_geer_object(x)
  cat("\nCall:\n")
  print(x$call)
  cat("\nEstimating Method   :", x$method, "\n")
  cat("Number of iterations:", x$iter, "\n")
  cat("Algorithm converged :", x$converged, "\n")
  print_nonconvergence_notice(x$converged)
  cat("\nNumber of observations:", x$obs_no, "\n")
  cat("Number of clusters    :", x$clusters_no, "\n")
  cat(
    "Cluster sizes         :",
    x$min_cluster_size,
    "to",
    x$max_cluster_size,
    "\n"
  )
  cat("\nMarginal Model\n")
  cat("Family       :", x$family$family, "\n")
  cat("Link Function:", x$family$link, "\n")
  cat("\nCoefficients:\n")
  print_with_small_pvalues(function() {
    stats::printCoefmat(x$coefficients, eps.Pvalue = 1e-04)
  })
  cat("Std. Errors are taken from the", x$cov_type, "covariance matrix.", "\n")
  cat("\nDispersion Parameter:", round(x$phi, digits = 4), "\n")
  cat("\nAssociation Structure:", x$association_structure, "\n")
  if (length(x$alpha) == 1L) {
    cat("Association Parameter:", round(x$alpha, digits = 4), "\n")
  } else {
    cat("Association Parameters:\n")
    print(round(x$alpha, digits = 4))
    cat("\n")
  }
  invisible(x)
}
