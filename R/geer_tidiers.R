#' @title
#' Tidy a \code{geer} Object
#'
#' @description
#' Summarises the coefficient-level components of a fitted \code{geer} model
#' into a tidy data frame, one row per term. Follows the \pkg{broom}
#' conventions so that \code{geer} objects integrate naturally with
#' \pkg{tidymodels}, \pkg{modelsummary}, and related tooling.
#'
#' @param x an object of class \code{"geer"}.
#' @param conf.int logical indicating whether to append confidence interval
#'        columns \code{conf.low} and \code{conf.high}. Default \code{FALSE}.
#' @param conf.level numeric coverage probability for the confidence interval
#'        when \code{conf.int = TRUE}. Default \code{0.95}.
#' @param exponentiate logical indicating whether to exponentiate coefficient
#'        estimates and confidence limits. Useful for logistic and log-link
#'        models. Standard errors and z-statistics are \emph{not} transformed.
#'        Default \code{FALSE}.
#' @param cov_type character specifying the covariance estimator used to
#'        compute standard errors and the Wald z-statistic. Options are
#'        \code{"robust"} (sandwich, default), \code{"bias-corrected"}
#'        (Morel et al. 2003), \code{"df-adjusted"} (MacKinnon 1985), and
#'        \code{"naive"} (model-based).
#' @param ... additional arguments passed to or from other methods (currently
#'        unused).
#'
#' @details
#' The returned data frame has the following columns:
#' \describe{
#'   \item{\code{term}}{name of the regression coefficient.}
#'   \item{\code{estimate}}{point estimate (or exponentiated value when
#'     \code{exponentiate = TRUE}).}
#'   \item{\code{std.error}}{standard error derived from
#'     \code{vcov(x, cov_type = cov_type)}.}
#'   \item{\code{statistic}}{Wald z-statistic
#'     \eqn{\hat{\beta} / \mathrm{SE}(\hat{\beta})}.}
#'   \item{\code{p.value}}{two-sided p-value from the standard normal
#'     distribution.}
#'   \item{\code{conf.low}, \code{conf.high}}{lower and upper Wald confidence
#'     limits (only present when \code{conf.int = TRUE}).}
#' }
#'
#' When \pkg{tibble} is available the result has class
#' \code{c("tbl_df", "tbl", "data.frame")}; otherwise a plain
#' \code{data.frame} is returned. Either form is accepted by
#' \pkg{modelsummary} and \pkg{broom}.
#'
#' @return A data frame with one row per regression coefficient.
#'
#' @references
#' Robinson D., Hayes A. and Couch S. (2024)
#' \emph{broom: Convert Statistical Objects into Tidy Tibbles}.
#' \url{https://broom.tidymodels.org/}.
#'
#' @seealso \code{\link{glance.geer}}, \code{\link{vcov.geer}},
#'   \code{\link{confint.geer}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' data("epilepsy")
#' fitmodel <- geewa(formula = seizures ~ treatment + lnbaseline + lnage,
#'   data = epilepsy, id = id, family = poisson(link = "log"), corstr  = "exchangeable",
#'   method  = "gee")
#' tidy(fitmodel)
#' tidy(fitmodel, conf.int = TRUE)
#' tidy(fitmodel, conf.int = TRUE, exponentiate = TRUE)
#' tidy(fitmodel, cov_type = "bias-corrected")
#'
#' data("cerebrovascular")
#' fitbin <- geewa_binary(formula = ecg ~ treatment + factor(period),
#'   id = id, data = cerebrovascular, link = "logit", orstr = "exchangeable")
#' tidy(fitbin, conf.int = TRUE, conf.level = 0.90, exponentiate = TRUE)
#'
#' @export
tidy.geer <- function(x,
                      conf.int = FALSE,
                      conf.level = 0.95,
                      exponentiate = FALSE,
                      cov_type = c("robust", "bias-corrected",
                                   "df-adjusted", "naive"),
                      ...) {
  if (!inherits(x, "geer")) {
    stop("'x' must be a 'geer' object", call. = FALSE)
  }
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
  beta <- coef(x)
  covariance_matrix <- vcov(x, cov_type = cov_type)
  se <- sqrt(pmax(0, diag(covariance_matrix)))
  zstat <- beta/se
  pval  <- 2 * pnorm(abs(zstat), lower.tail = FALSE)
  ans <- data.frame(
    term = names(beta),
    estimate = unname(beta),
    std.error = unname(se),
    statistic = unname(zstat),
    p.value = unname(pval),
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    ci <- confint(x, level = conf.level, cov_type = cov_type)
    ans$conf.low <- unname(ci[, 1L])
    ans$conf.high <- unname(ci[, 2L])
  }
  if (exponentiate) {
    ans$estimate <- exp(ans$estimate)
    if (conf.int) {
      ans$conf.low  <- exp(ans$conf.low)
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
#' Glance at a \code{geer} Object
#'
#' @description
#' Returns a single-row model-level summary of a fitted \code{geer} object,
#' following \pkg{broom} conventions. Suitable for use with
#' \pkg{modelsummary} and related reporting packages, and for comparing models
#' side-by-side via \code{lapply(model_list, glance)}.
#'
#' @param x an object of class \code{"geer"}.
#' @param ... additional arguments passed to or from other methods (currently
#'        unused).
#'
#' @details
#' The one-row data frame contains the following columns:
#' \describe{
#'   \item{\code{family}}{name of the marginal response family.}
#'   \item{\code{link}}{name of the link function.}
#'   \item{\code{method}}{estimation method, e.g. \code{"gee"},
#'     \code{"brgee-robust"}, or \code{"pgee-jeffreys"}.}
#'   \item{\code{corstr}}{working association structure.}
#'   \item{\code{nobs}}{total number of observations \eqn{n^{\star}}.}
#'   \item{\code{nclusters}}{number of independent clusters \eqn{N}.}
#'   \item{\code{min.cluster.size}}{minimum cluster size.}
#'   \item{\code{max.cluster.size}}{maximum cluster size.}
#'   \item{\code{npar}}{number of mean-model parameters \eqn{p}.}
#'   \item{\code{df.residual}}{residual degrees of freedom
#'     \eqn{n^{\star} - p}.}
#'   \item{\code{phi}}{estimated (or fixed) dispersion parameter. Equal to
#'     1 for the odds-ratio binary parameterisation.}
#'   \item{\code{QIC}}{Quasi Information Criterion (Pan 2001) computed from
#'     the robust sandwich covariance estimator. Lower is better.}
#'   \item{\code{QICu}}{Covariate-selection variant of QIC (Pan 2001).
#'     Lower is better.}
#'   \item{\code{CIC}}{Correlation Information Criterion (Hin and Wang 2009).
#'     Lower is better.}
#'   \item{\code{converged}}{logical; \code{TRUE} if the fitting algorithm
#'     converged.}
#'   \item{\code{niter}}{number of iterations used.}
#' }
#'
#' QIC, QICu and CIC are always evaluated using the robust sandwich covariance
#' matrix, consistent with \code{\link{geecriteria}}. If their computation
#' fails for any reason (e.g. a degenerate working correlation matrix) the
#' corresponding columns contain \code{NA_real_}.
#'
#' When \pkg{tibble} is available the result has class
#' \code{c("tbl_df", "tbl", "data.frame")}; otherwise a plain
#' \code{data.frame} is returned.
#'
#' @return A one-row data frame.
#'
#' @references
#' Pan W. (2001) Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics} \bold{57}, 120--125.
#'
#' Hin L.Y. and Wang Y.G. (2009) Working-correlation-structure identification
#' in generalized estimating equations. \emph{Statistics in Medicine}
#' \bold{28}, 642--658.
#'
#' Robinson D., Hayes A. and Couch S. (2024)
#' \emph{broom: Convert Statistical Objects into Tidy Tibbles}.
#' \url{https://broom.tidymodels.org/}.
#'
#' @seealso \code{\link{tidy.geer}}, \code{\link{geecriteria}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#'
#' data("epilepsy")
#' fitmodel <- geewa(formula = seizures ~ treatment + lnbaseline + lnage,
#'   data = epilepsy, id = id, family = poisson(link = "log"), corstr  = "exchangeable",
#'   method  = "gee")
#' glance(fitmodel)
#'
#' data("cerebrovascular")
#' fitbin <- geewa_binary(formula = ecg ~ treatment + factor(period),
#'   id = id, data = cerebrovascular, link = "logit", orstr = "exchangeable")
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
#' ))[, c("corstr", "QIC", "CIC", "niter")]
#' }
#'
#' @export
glance.geer <- function(x, ...) {
  if (!inherits(x, "geer")) {
    stop("'x' must be a 'geer' object", call. = FALSE)
  }
  crit <- tryCatch(
    compute_criteria(x, cov_type = "robust", digits = 15L),
    error = function(e) NULL
  )
  ans <- data.frame(
    family = x$family$family,
    link = x$family$link,
    method = x$method,
    corstr = x$association_structure,
    nobs = x$obs_no,
    nclusters = x$clusters_no,
    min.cluster.size = x$min_cluster_size,
    max.cluster.size = x$max_cluster_size,
    npar = length(x$coefficients),
    df.residual = x$df.residual,
    phi = x$phi,
    QIC = if (!is.null(crit)) crit$QIC  else NA_real_,
    QICu = if (!is.null(crit)) crit$QICu else NA_real_,
    CIC = if (!is.null(crit)) crit$CIC  else NA_real_,
    converged = x$converged,
    niter = x$iter,
    stringsAsFactors = FALSE
  )
  if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::as_tibble(ans)
  } else {
    ans
  }
}
