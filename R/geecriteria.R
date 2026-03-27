#' @title
#' Model Selection Criteria for \code{geer} Objects
#'
#' @description
#' Compute commonly used criteria for comparing fitted \code{geer} models for
#' marginal mean selection and working association structure selection.
#'
#' @return
#' A data frame with one row per fitted model and columns \code{QIC},
#' \code{CIC}, \code{RJC}, \code{QICu}, \code{GESSC}, \code{GPC}, and
#' \code{Parameters}.
#'
#'
#' @inheritParams anova.geer
#' @param cov_type Character string specifying the covariance estimator used in
#'        the covariance-based criteria. Options are \code{"robust"},
#'        \code{"bias-corrected"}, \code{"df-adjusted"}, and \code{"naive"}.
#'        The default is \code{"robust"}.
#' @param digits Non-negative integer giving the number of decimal places used
#'        to round the reported criteria. The default is \code{2}.
#'
#' @details
#' The reported criteria are:
#' \describe{
#'   \item{\code{QIC}}{Quasi Information Criterion for model comparison. Smaller
#'     values are preferred.}
#'   \item{\code{CIC}}{Correlation Information Criterion for selecting the
#'     working association structure. Smaller values are preferred.}
#'   \item{\code{RJC}}{Rotnitzky and Jewell Criterion. Smaller values are
#'     preferred.}
#'   \item{\code{QICu}}{A variant of QIC primarily intended for comparing mean
#'     models with different covariate sets. Smaller values are preferred.}
#'   \item{\code{GESSC}}{Generalized Error Sum of Squares Criterion. Smaller
#'     values are preferred.}
#'   \item{\code{GPC}}{Gaussian Pseudolikelihood Criterion. Larger values are
#'     preferred.}
#'   \item{\code{Parameters}}{Number of marginal mean-model regression
#'     parameters.}
#' }
#'
#' The \code{cov_type} argument affects the covariance-based criteria,
#' specifically \code{QIC}, \code{CIC}, and \code{RJC}. If several supplied
#' models do not have the same number of observations, a warning is issued.
#'
#' @references
#' Carey, V.J. and Wang, Y.G. (2011) Working covariance model selection for
#' generalized estimating equations. \emph{Statistics in Medicine}, \bold{30},
#' 3117--3124.
#'
#' Chaganty, N.R. and Shults, J. (1999) On eliminating the asymptotic bias in
#' the quasi-least squares estimate of the correlation parameter.
#' \emph{Journal of Statistical Planning and Inference}, \bold{76}, 145--161.
#'
#' Hin, L.Y. and Wang, Y.G. (2009) Working correlation structure
#' identification in generalized estimating equations. \emph{Statistics in
#' Medicine} \bold{28}, 642--658.
#'
#' Pan, W. (2001) Akaike's information criterion in generalized
#' estimating equations. \emph{Biometrics} \bold{57}, 120--125.
#'
#' Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
#' parameters in semiparametric generalized linear models for cluster correlated
#' data. \emph{Biometrika} \bold{77}, 485--497.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fitted_gee <- geewa_binary(
#'   formula = ecg ~ period * treatment,
#'   id = id,
#'   data = cerebrovascular,
#'   link = "logit",
#'   orstr = "exchangeable",
#'   method = "gee"
#' )
#' fitted_brgee <- update(fitted_gee, method = "brgee-robust")
#'
#' ## Compare two fits
#' geecriteria(fitted_gee, fitted_brgee, cov_type = "robust")
#'
#' @export
geecriteria <- function(object,
                        ...,
                        cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                        digits = 2) {
  cov_type <- match.arg(cov_type)
  digits <- check_nonnegative_integerish(digits, "digits")
  models <- c(list(object), list(...))
  models <- lapply(models, check_geer_object)
  obs_no <- vapply(models, function(model) {
    if (!is.null(model$obs_no)) {
      as.numeric(model$obs_no)
    } else {
      NA_real_
    }
  }, numeric(1))
  if (sum(!is.na(obs_no)) > 1L) {
    ref_obs <- obs_no[which(!is.na(obs_no))[1L]]
    if (any(obs_no[!is.na(obs_no)] != ref_obs)) {
      warning("models do not have the same number of observations", call. = FALSE)
    }
  }
  out_list <- lapply(models, compute_gee_criteria, cov_type = cov_type, digits = NULL)
  ans <- do.call(rbind, out_list)
  numeric_cols <- c("QIC", "CIC", "RJC", "QICu", "GESSC", "GPC")
  ans[, numeric_cols] <- lapply(ans[, numeric_cols, drop = FALSE], round, digits = digits)
  ans[, "Parameters"] <- round(ans[, "Parameters"], digits = 0)
  if (length(models) > 1L && nrow(ans) == length(models)) {
    call_expr <- match.call(expand.dots = FALSE)
    exprs <- c(list(call_expr[[2L]]), as.list(call_expr$...))
    rownames(ans) <- vapply(
      exprs,
      function(expr) paste(deparse(expr, width.cutoff = 500L), collapse = ""),
      character(1)
    )
  }
  ans
}
