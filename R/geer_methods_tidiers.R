#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @title
#' Tidy a geer Object
#'
#' @description
#' Summarizes the regression coefficients of a fitted \code{geer} object in a
#' data frame with one row per coefficient, following the conventions of
#' \pkg{broom} (Robinson, Hayes and Couch, 2024).
#'
#' @param x a fitted model object of class \code{"geer"}.
#' @param conf.int logical indicating whether Wald confidence limits are
#'   included. Defaults to \code{FALSE}.
#' @param conf.level a single number strictly between 0 and 1 specifying the
#'   confidence level. Defaults to \code{0.95}.
#' @param exponentiate logical indicating whether the estimates and confidence
#'   limits are exponentiated, for example to obtain odds ratios. Defaults to
#'   \code{FALSE}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for the standard errors, z-statistics, and confidence intervals: one
#'   of the estimators described in \code{\link{vcov.geer}}. Defaults to
#'   \code{"bias-corrected"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @return
#' A data frame, a tibble if \pkg{tibble} is installed, with columns:
#' \item{term}{the name of the coefficient.}
#' \item{estimate}{the estimate, exponentiated if \code{exponentiate = TRUE}.}
#' \item{std.error}{the standard error, on the original scale.}
#' \item{statistic}{the Wald z-statistic, on the original scale.}
#' \item{p.value}{the two-sided p-value.}
#' \item{conf.low, conf.high}{the confidence limits, if
#'   \code{conf.int = TRUE}.}
#'
#' @references
#' Robinson, D., Hayes, A. and Couch, S. (2024) \emph{broom: Convert
#' Statistical Objects into Tidy Tibbles}. R package.
#' \url{https://broom.tidymodels.org/}
#'
#' @seealso \code{\link{glance.geer}}, \code{\link{vcov.geer}},
#'   \code{\link{confint.geer}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
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
#' tidy(fit)
#' tidy(fit, conf.int = TRUE, exponentiate = TRUE)
#' tidy(fit, cov_type = "robust")
#'
#' @export
tidy.geer <- function(x,
                      conf.int = FALSE,
                      conf.level = 0.95,
                      exponentiate = FALSE,
                      cov_type = "bias-corrected",
                      ...) {
  object <- check_geer_object(x)
  cov_type <- match_cov_type(cov_type)
  check_unused_dots(list(...), "tidy.geer")
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("'conf.int' must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop(
      "'conf.level' must be a number strictly between 0 and 1",
      call. = FALSE
    )
  }
  if (!is.logical(exponentiate) || length(exponentiate) != 1L ||
      is.na(exponentiate)) {
    stop("'exponentiate' must be TRUE or FALSE", call. = FALSE)
  }
  if (exponentiate && !object$family$link %in% c("log", "logit", "cloglog")) {
    warning(
      paste0(
        "'exponentiate = TRUE' is usually meaningful only for log, logit, and ",
        "cloglog links"
      ),
      call. = FALSE
    )
  }
  beta <- stats::coef(object)
  vcov_matrix <- stats::vcov(object, cov_type = cov_type)
  se <- compute_standard_errors_or_na(diag(vcov_matrix), "tidy.geer")
  z_stat <- rep(NA_real_, length(beta))
  ok <- is.finite(beta) & is.finite(se) & se > 0
  z_stat[ok] <- beta[ok] / se[ok]
  pval <- rep(NA_real_, length(beta))
  pval[ok] <- 2 * stats::pnorm(abs(z_stat[ok]), lower.tail = FALSE)
  ans <- data.frame(
    term = names(beta),
    estimate = unname(beta),
    std.error = unname(se),
    statistic = unname(z_stat),
    p.value = unname(pval),
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    conf_probs <- (1 - conf.level) / 2
    conf_probs <- c(conf_probs, 1 - conf_probs)
    limits <- beta + se %o% stats::qnorm(conf_probs)
    ans$conf.low <- unname(limits[, 1L])
    ans$conf.high <- unname(limits[, 2L])
  }
  if (exponentiate) {
    ans$estimate <- exp(ans$estimate)
    if (conf.int) {
      ans$conf.low <- exp(ans$conf.low)
      ans$conf.high <- exp(ans$conf.high)
    }
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::as_tibble(ans)
  } else {
    ans
  }
}


#' @title
#' Glance at a geer Object
#'
#' @description
#' Summarizes a fitted \code{geer} object in a one-row data frame, following
#' the conventions of \pkg{broom} (Robinson, Hayes and Couch, 2024).
#'
#' @param x a fitted model object of class \code{"geer"}.
#' @param ... additional arguments passed to or from other methods. Currently
#'   unused.
#'
#' @details
#' \code{QIC} and \code{CIC} use the robust covariance matrix, as in
#' \code{geecriteria(x, cov_type = "robust")}, and are \code{NA} if they
#' cannot be computed. See \code{\link{geecriteria}} for further criteria.
#'
#' @return
#' A one-row data frame, a tibble if \pkg{tibble} is installed, with columns:
#' \item{family, link}{the family and link function.}
#' \item{method}{the estimation method.}
#' \item{wastr}{the working association structure, \code{corstr} or
#'   \code{orstr}.}
#' \item{nobs, nclusters}{the numbers of observations and clusters.}
#' \item{min.cluster.size, max.cluster.size}{the smallest and largest cluster
#'   sizes.}
#' \item{npar}{the number of regression parameters.}
#' \item{df.residual}{the number of observations minus \code{npar}.}
#' \item{phi}{the dispersion parameter, 1 for \code{\link{geewa_binary}}
#'   fits.}
#' \item{QIC, QICu}{the quasi-likelihood information criteria of Pan (2001).}
#' \item{CIC}{the correlation information criterion of Hin and Wang (2009).}
#' \item{converged}{logical indicating whether the algorithm converged.}
#' \item{niter}{the number of iterations.}
#'
#' @references
#' Hin, L.Y. and Wang, Y.G. (2009) Working-correlation-structure
#' identification in generalized estimating equations. \emph{Statistics in
#' Medicine}, \bold{28}, 642--658.
#'
#' Pan, W. (2001) Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, \bold{57}, 120--125.
#'
#' Robinson, D., Hayes, A. and Couch, S. (2024) \emph{broom: Convert
#' Statistical Objects into Tidy Tibbles}. R package.
#' \url{https://broom.tidymodels.org/}
#'
#' @seealso \code{\link{tidy.geer}}, \code{\link{geecriteria}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
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
#' glance(fit)
#'
#' ## Compare working correlation structures
#' fit_ar1 <- update(fit, corstr = "ar1")
#' rbind(glance(fit), glance(fit_ar1))
#'
#' @export
glance.geer <- function(x, ...) {
  object <- check_geer_object(x)
  crit <- tryCatch(
    compute_criteria(
      object, cov_type = "robust", digits = 15L, include_extended = FALSE,
      criteria = c("QIC", "QICu", "CIC")
    ),
    error = function(e) NULL
  )
  qic <- if (is.null(crit)) NA_real_ else crit$QIC
  qicu <- if (is.null(crit)) NA_real_ else crit$QICu
  cic <- if (is.null(crit)) NA_real_ else crit$CIC
  ans <- data.frame(
    family = object$family$family,
    link = object$family$link,
    method = object$method,
    wastr = object$association_structure,
    nobs = object$obs_no,
    nclusters = object$clusters_no,
    min.cluster.size = object$min_cluster_size,
    max.cluster.size = object$max_cluster_size,
    npar = length(object$coefficients),
    df.residual = object$df.residual,
    phi = object$phi,
    QIC = qic,
    QICu = qicu,
    CIC = cic,
    converged = object$converged,
    niter = object$iter,
    stringsAsFactors = FALSE
  )
  if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::as_tibble(ans)
  } else {
    ans
  }
}
