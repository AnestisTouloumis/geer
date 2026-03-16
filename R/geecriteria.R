#' @title
#' Variable and Correlation Selection Criteria for \code{geer} Objects
#'
#' @description
#' Reports commonly used criteria for variable selection and association structure
#' selection for one or several fitted models.
#'
#' @return
#' If one model is supplied, a one-row data frame containing QIC, CIC, RJC, QICu,
#' GESSC, GPC and the number of regression parameters.
#'
#' If multiple models are supplied, a data frame with one row per model and the
#' same columns.
#'
#'
#' @inheritParams anova.geer
#' @param cov_type character indicating the type of covariance matrix estimator
#'        which should be used to the calculation of the criteria. Options
#'        include the sandwich or robust estimator (\code{"robust"}), the
#'        bias-corrected estimator (\code{"bias-corrected"}), the degrees of
#'        freedom adjusted estimator (\code{"df-adjusted"}) and the model-based
#'        or naive estimator (\code{"naive"}). By default,
#'        \code{cov_type = "robust"}.
#' @param digits integer indicating the number of decimal places for the GEE
#'        criteria to be rounded at. By default, \code{digits = 2}.
#'
#' @details
#' The Quasi Information Criterion (QIC), the Correlation Information Criterion
#' (CIC), the Rotnitzky and Jewell Criterion (RJC), the Generalized Error Sum of
#' Squares (GESSC) and the Gaussian Pseudolikelihood Criterion (GPC) are used for
#' selecting the best association structure.
#' The QICu criterion is used for selecting the best subset of covariates.
#'
#' When choosing between ordinary, bias-reducing and penalized GEE models with
#' the same subset of covariates, the model with the smallest value of QIC, CIC,
#' RJC or GESSC or the model with the higher GPC should be preferred. When
#' choosing between ordinary, bias-reducing and penalized models with different
#' number of covariates, the model with the smallest QICu value should be
#' preferred.
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
#' data("cerebrovascular")
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
geecriteria <- function(object, ..., cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"), digits = 2) {
  cov_type <- match.arg(cov_type)
  is_single_nonneg_int <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0 && x == floor(x)
  }
  if (!is_single_nonneg_int(digits)) {
    stop("'digits' must be a single non-negative integer", call. = FALSE)
  }
  digits <- as.integer(digits)
  models <- c(list(object), list(...))
  is_geer <- vapply(models, inherits, logical(1), what = "geer")
  if (!all(is_geer)) {
    stop("Only 'geer' objects are supported", call. = FALSE)
  }
  out_list <- lapply(models, compute_criteria, cov_type = cov_type, digits = digits)
  ans <- do.call(rbind, out_list)
  obs_no <- vapply(models, function(m) {
    if (!is.null(m$obs_no)) as.numeric(m$obs_no) else NA_real_
  }, numeric(1))
  if (sum(!is.na(obs_no)) > 1L && any(obs_no[!is.na(obs_no)] != obs_no[which(!is.na(obs_no))[1L]])) {
    warning("models do not have the same number of observations", call. = FALSE)
  }
  if (length(models) > 1L && nrow(ans) == length(models)) {
    Call <- match.call(expand.dots = FALSE)
    exprs <- c(list(Call[[2L]]), as.list(Call$...))
    rn <- vapply(exprs, function(e) paste(deparse(e, width.cutoff = 500), collapse = ""), character(1))
    rownames(ans) <- rn
  }
  ans
}
