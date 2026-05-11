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
#' Summarizes a fitted \code{geer} object at the coefficient level, returning
#' one row per regression term in a tidy data frame.
#'
#' @param x an object of class \code{"geer"}.
#' @param conf.int logical indicating whether to append confidence interval
#'   columns \code{conf.low} and \code{conf.high}. Defaults to
#'   \code{FALSE}.
#' @param conf.level numeric coverage probability for the confidence interval
#'   when \code{conf.int = TRUE}. Defaults to \code{0.95}.
#' @param exponentiate logical indicating whether to exponentiate coefficient
#'   estimates and confidence limits. This is often useful for models with a
#'   log, logit, or complementary log-log link. Standard errors and Wald
#'   z-statistics are \emph{not} transformed. Defaults to \code{FALSE}.
#' @param cov_type character string specifying the covariance estimator used to
#'   compute standard errors and Wald z-statistics. Options are
#'   \code{"bias-corrected"} (default), \code{"robust"}, \code{"df-adjusted"},
#'   and \code{"naive"} (model-based). See \code{\link{vcov.geer}} for details.
#' @param ... additional arguments passed to or from other methods. Currently
#'   unused.
#'
#' @details
#' The returned data frame contains the following columns:
#' \describe{
#'   \item{\code{term}}{name of the regression coefficient.}
#'   \item{\code{estimate}}{point estimate, or the exponentiated point estimate
#'   when \code{exponentiate = TRUE}.}
#'   \item{\code{std.error}}{standard error derived from
#'   \code{vcov(x, cov_type = cov_type)}.}
#'   \item{\code{statistic}}{Wald z-statistic, computed as the coefficient
#'   estimate divided by its standard error. \code{NA} is reported when the
#'   standard error is not positive.}
#'   \item{\code{p.value}}{two-sided p-value from the standard normal
#'   distribution.}
#'   \item{\code{conf.low}, \code{conf.high}}{lower and upper Wald confidence
#'   limits. These columns are included only when \code{conf.int = TRUE}.}
#' }
#'
#' When \code{exponentiate = TRUE}, only \code{estimate}, \code{conf.low}, and
#' \code{conf.high} are exponentiated. The columns \code{std.error} and
#' \code{statistic} remain on the original scale.
#'
#' @return
#' A data frame with one row per regression coefficient and columns as
#' described in the Details section. When \pkg{tibble} is available, the result
#' has class \code{c("tbl_df", "tbl", "data.frame")}. Otherwise, a plain
#' \code{data.frame} is returned.
#'
#' @references
#' Robinson D., Hayes A. and Couch S. (2024)
#' \emph{broom: Convert Statistical Objects into Tidy Tibbles}.
#' \url{https://broom.tidymodels.org/}.
#'
#' @seealso
#' \code{\link{glance.geer}}, \code{\link{vcov.geer}},
#' \code{\link{confint.geer}}, \code{\link{geewa}},
#' \code{\link{geewa_binary}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fitmodel <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   corstr = "exchangeable",
#'   method = "gee"
#' )
#' tidy(fitmodel)
#' tidy(fitmodel, conf.int = TRUE)
#' tidy(fitmodel, conf.int = TRUE, exponentiate = TRUE)
#' tidy(fitmodel, cov_type = "robust")
#'
#' data("cerebrovascular", package = "geer")
#' fitbin <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' tidy(fitbin, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)
#'
#' @export
tidy.geer <- function(x,
                      conf.int = FALSE,
                      conf.level = 0.95,
                      exponentiate = FALSE,
                      cov_type = c("bias-corrected", "robust",
                                   "df-adjusted", "naive"),
                      ...) {
  object <- check_geer_object(x)
  cov_type <- match.arg(cov_type)
  if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
    stop("'conf.int' must be a single logical value", call. = FALSE)
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number in (0, 1)", call. = FALSE)
  }
  if (!is.logical(exponentiate) || length(exponentiate) != 1L ||
      is.na(exponentiate)) {
    stop("'exponentiate' must be a single logical value", call. = FALSE)
  }
  if (exponentiate && !object$family$link %in% c("log", "logit", "cloglog")) {
    warning(
      "'exponentiate = TRUE' is typically most meaningful for log, logit, or cloglog links",
      call. = FALSE
    )
  }
  beta <- coef(object)
  vcov_matrix <- vcov(object, cov_type = cov_type)
  se <- sqrt(pmax(0, diag(vcov_matrix)))
  z_stat <- rep(NA_real_, length(beta))
  ok <- is.finite(beta) & is.finite(se) & se > 0
  z_stat[ok] <- beta[ok] / se[ok]
  pval <- rep(NA_real_, length(beta))
  pval[ok] <- 2 * pnorm(abs(z_stat[ok]), lower.tail = FALSE)
  ans <- data.frame(
    term = names(beta),
    estimate = unname(beta),
    std.error = unname(se),
    statistic = unname(z_stat),
    p.value = unname(pval),
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    ci <- confint(object, level = conf.level, cov_type = cov_type)
    ans$conf.low <- unname(ci[, 1L])
    ans$conf.high <- unname(ci[, 2L])
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
#' Produces a one-row summary for a fitted \code{geer} object, following
#' \pkg{broom} conventions.
#'
#' @param x an object of class \code{"geer"}.
#' @param ... additional arguments passed to or from other methods. Currently
#'   unused.
#'
#' @details
#' The returned one-row data frame contains the following columns:
#' \describe{
#'   \item{\code{family}}{name of the marginal response family.}
#'   \item{\code{link}}{name of the link function.}
#'   \item{\code{method}}{estimation method, for example \code{"gee"},
#'   \code{"brgee-robust"}, \code{"pgee-jeffreys"},
#'   \code{"opgee-jeffreys"}, or \code{"hpgee-jeffreys"}.}
#'   \item{\code{wastr}}{stored working association structure. For
#'   \code{geewa()} fits this corresponds to \code{corstr}; for
#'   \code{geewa_binary()} fits it corresponds to \code{orstr}.}
#'   \item{\code{nobs}}{total number of observations \eqn{n^{\star}}.}
#'   \item{\code{nclusters}}{number of independent clusters \eqn{N}.}
#'   \item{\code{min.cluster.size}}{minimum cluster size.}
#'   \item{\code{max.cluster.size}}{maximum cluster size.}
#'   \item{\code{npar}}{number of marginal mean-model parameters \eqn{p}.}
#'   \item{\code{df.residual}}{residual degrees of freedom
#'   \eqn{n^{\star} - p}.}
#'   \item{\code{phi}}{estimated or fixed dispersion parameter. This equals
#'   \code{1} for the odds-ratio parameterization used by
#'   \code{geewa_binary()}.}
#'   \item{\code{QIC}}{Quasi Information Criterion (Pan, 2001), computed from
#'   the robust sandwich covariance estimator. Smaller values are preferred.}
#'   \item{\code{QICu}}{Covariate-selection variant of QIC (Pan, 2001).
#'   Smaller values are preferred.}
#'   \item{\code{CIC}}{Correlation Information Criterion (Hin and Wang, 2009),
#'   used here to compare working association structures. Smaller values are
#'   preferred.}
#'   \item{\code{converged}}{logical; \code{TRUE} if the fitting algorithm
#'   converged.}
#'   \item{\code{niter}}{number of iterations used.}
#' }
#'
#' QIC and CIC are computed using the same formulas as
#' \code{geecriteria(object, cov_type = "robust")}. QICu does not depend on
#' the covariance estimator. If computation fails, the corresponding values
#' are returned as \code{NA_real_}. For the full set of
#' model selection criteria, including RJC, GESSC, and GPC, see
#' \code{\link{geecriteria}}.
#'
#' @return
#' A one-row data frame with columns as described in the Details section.
#' When \pkg{tibble} is available, the result has class
#' \code{c("tbl_df", "tbl", "data.frame")}. Otherwise, a plain
#' \code{data.frame} is returned.
#'
#' @references
#' Pan W. (2001) Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, \bold{57}, 120--125.
#'
#' Hin L.Y. and Wang Y.G. (2009) Working-correlation-structure identification
#' in generalized estimating equations. \emph{Statistics in Medicine},
#' \bold{28}, 642--658.
#'
#' Robinson D., Hayes A. and Couch S. (2024)
#' \emph{broom: Convert Statistical Objects into Tidy Tibbles}.
#' \url{https://broom.tidymodels.org/}.
#'
#' @seealso
#' \code{\link{tidy.geer}}, \code{\link{geecriteria}}, \code{\link{geewa}},
#' \code{\link{geewa_binary}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fitmodel <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   corstr = "exchangeable",
#'   method = "gee"
#' )
#' glance(fitmodel)
#'
#' data("cerebrovascular", package = "geer")
#' fitbin <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' glance(fitbin)
#'
#' \donttest{
#' fitind  <- update(fitmodel, corstr = "independence")
#' fitar1  <- update(fitmodel, corstr = "ar1")
#' fitunst <- update(fitmodel, corstr = "unstructured")
#' do.call(rbind, lapply(
#'   list(independence = fitind, exchangeable = fitmodel, ar1 = fitar1,
#'        unstructured = fitunst),
#'   glance
#' ))[, c("wastr", "QIC", "CIC", "niter")]
#' }
#'
#' @export
glance.geer <- function(x, ...) {
  object <- check_geer_object(x)
  crit <- tryCatch(
    compute_gee_criteria(object, cov_type = "robust", digits = 15L),
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
