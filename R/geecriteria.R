#' @title
#' Model Selection Criteria for geer Objects
#'
#' @description
#' Computes model selection criteria for one or more fitted \code{geer} objects,
#' supporting both marginal mean model comparison and working association
#' structure selection.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param ... additional fitted model objects of class \code{"geer"} to be
#'   included in the comparison.
#' @param criteria character vector selecting which criteria to report. The
#'   default, \code{"all"}, returns every criterion. Otherwise give one or more
#'   criterion names from \code{"QIC"}, \code{"QICHH"}, \code{"QICC"},
#'   \code{"CIC"}, \code{"RJC"}, \code{"QICu"}, \code{"EQIC"},
#'   \code{"GESSC"}, \code{"GPC"}, \code{"AGPC"}, \code{"SGPC"},
#'   \code{"GHYC"}, \code{"PAC"}, \code{"PT"}, \code{"WR"}, and
#'   \code{"RMR"}. Names are matched ignoring case, and the columns are
#'   returned in the order requested. \code{"all"} cannot be combined with
#'   individual names. Only the requested criteria are computed, so restricting
#'   the selection also avoids the work they would require; this matters most
#'   with \code{cov_type = "jackknife"}, which is needed by \code{QIC},
#'   \code{QICHH}, \code{QICC}, \code{CIC}, \code{RJC}, \code{EQIC},
#'   \code{PT}, \code{WR}, and \code{RMR} but not by the remaining criteria.
#' @param cov_type character string specifying the covariance estimator used in
#'   the covariance-based penalty terms. Options are the sandwich or robust
#'   estimator (\code{"robust"}), the bias-corrected estimator
#'   (\code{"bias-corrected"}), the degrees-of-freedom adjusted estimator
#'   (\code{"df-adjusted"}), the leave-one-cluster jackknife estimator
#'   (\code{"jackknife"}), and the model-based or naive estimator
#'   (\code{"naive"}). Defaults to \code{"robust"}, which reproduces the
#'   classical covariance-based definitions of QIC, QICHH, QICC, CIC, RJC, and EQIC.
#'   With \code{cov_type = "jackknife"} a full set of leave-one-cluster refits
#'   is performed for each model supplied; see \code{\link{vcov.geer}}.
#' @param digits non-negative integer giving the number of decimal places used
#'   to round the reported criteria. Defaults to \code{2}.
#'
#' @details
#' The reported criteria are:
#' \describe{
#'   \item{\code{QIC}}{Quasi Information Criterion of Pan (2001), for
#'   comparing marginal mean models. Smaller values are preferred.}
#'   \item{\code{QICHH}}{Modification of QIC due to Hardin and Hilbe (2003,
#'   2013). The independence
#'   quasi-likelihood and independence model-based information matrix are
#'   evaluated using the regression and scale estimates under working
#'   independence. Smaller values are preferred.}
#'   \item{\code{QICC}}{Corrected QIC of Hardin and Hilbe (2013), using
#'   QIC as the base criterion and a finite-cluster correction that also
#'   accounts for the number of working-association parameters. Smaller
#'   values are preferred.}
#'   \item{\code{CIC}}{Correlation Information Criterion of Hin and Wang
#'   (2009), used here for selecting the working association structure.
#'   Smaller values are preferred.}
#'   \item{\code{RJC}}{Rotnitzky-Jewell Criterion, built by Hin, Carey and
#'   Wang (2007) from the working Wald statistic of Rotnitzky and Jewell
#'   (1990). Smaller values are preferred.}
#'   \item{\code{QICu}}{Variant of QIC given by Pan (2001), primarily
#'   intended for comparing
#'   marginal mean models with different covariate sets. Smaller values are
#'   preferred.}
#'   \item{\code{EQIC}}{Extended Quasi-likelihood Information Criterion of
#'   Wang and Hin (2010). It uses the extended quasi-likelihood with the
#'   Nelder-Pregibon adjustment \eqn{k = 1/6}. When the dispersion is not
#'   fixed, EQIC uses the deviance-based dispersion estimate. Smaller values
#'   are preferred.}
#'   \item{\code{GESSC}}{Generalized Error Sum of Squares Criterion: the
#'   error sum of squares weighted by the inverse working covariance matrix,
#'   as proposed by Shults and Chaganty (1998) and Shults et al. (2009),
#'   divided here by the residual degrees of freedom \eqn{N - p - m}. Smaller
#'   values are preferred. The unscaled statistic is the criterion denoted
#'   \eqn{SC(R)} by Pardo and Alonso (2019); because \eqn{m} varies across
#'   working association structures, the scaling is not a monotone
#'   transformation of it and the two may rank structures differently.}
#'   \item{\code{GPC}}{Gaussian Pseudolikelihood Criterion of Carey and Wang
#'   (2011). Unlike the other
#'   criteria, larger values are preferred, because GPC is a pseudolikelihood
#'   measure rather than an information criterion.}
#'   \item{\code{AGPC}}{Akaike-type penalized Gaussian Pseudolikelihood
#'   Criterion of Carey and Wang (2011); see also Zhu and Zhu (2013) and Fu,
#'   Hao and Wang (2018). Smaller values are preferred.}
#'   \item{\code{SGPC}}{Schwarz-type penalized Gaussian Pseudolikelihood
#'   Criterion of Carey and Wang (2011); see also Zhu and Zhu (2013) and Fu,
#'   Hao and Wang (2018). Smaller values are preferred.}
#'   \item{\code{GHYC}}{Gosho-Hamada-Yoshimura Criterion of Gosho, Hamada and
#'   Yoshimura (2011) for selecting the working association structure, denoted
#'   DEW by Gosho (2014). Computed
#'   only for balanced designs. Smaller values are preferred.}
#'   \item{\code{PAC}}{Pardo-Alonso Criterion of Pardo and Alonso (2019) for
#'   selecting the working
#'   association structure. Computed only for balanced designs. Smaller values
#'   are preferred.}
#'   \item{\code{PT}}{Pillai trace type criterion of Jang (2011). Smaller
#'   values are preferred.}
#'   \item{\code{WR}}{Wilks ratio type criterion of Jang (2011). Smaller
#'   values are preferred.}
#'   \item{\code{RMR}}{Roy maximum root type criterion of Jang (2011).
#'   Smaller values are preferred.}
#'   \item{\code{Parameters}}{Number of regression parameters in the marginal
#'   mean model. Always reported, whatever \code{criteria} selects.}
#' }
#'
#' A criterion that cannot be evaluated for a particular fit is reported as
#' \code{NA} rather than raising an error, so that the remaining criteria and
#' the remaining models in the comparison are still returned.
#'
#' The quasi-likelihood term shared by \code{QIC}, \code{QICu}, \code{QICC},
#' and \code{QICHH} involves the logarithm of the fitted means, which diverges
#' when a fitted mean reaches the boundary of its support, for example a fitted
#' probability of zero or one under the binomial family or a fitted rate of
#' zero under the Poisson family. Fitted means are therefore clamped away from
#' those boundaries by \code{sqrt(.Machine$double.eps)}, which keeps the four
#' criteria finite. Because the quasi-likelihood term is very nearly constant
#' across working association structures, this does not affect the comparison
#' of structures fitted to the same mean model. It does, however, make the
#' affected criteria non-comparable between candidate mean models when one
#' produces boundary fitted values and another does not, since the clamped
#' contribution is an arbitrary large value rather than the divergent one; in
#' that situation prefer \code{CIC}, which has no quasi-likelihood term.
#'
#' The \code{cov_type} argument affects the covariance penalty in
#' \code{QIC}, \code{QICHH}, \code{QICC}, \code{CIC}, \code{RJC},
#' \code{EQIC}, \code{PT}, \code{WR}, and \code{RMR}.
#' The classical definitions use the robust sandwich covariance, so
#' \code{cov_type = "robust"} is the default. \code{QICu}, \code{GESSC},
#' \code{GPC}, \code{AGPC}, \code{SGPC}, \code{GHYC}, and \code{PAC}
#' are unaffected by \code{cov_type}.
#'
#' QICHH is based on the conventional, unadjusted working-independence GEE
#' regression estimate and has the form
#' \deqn{\mathrm{QIC}_{HH} = -2\Psi(\hat\beta_I; I) +
#' 2\mathrm{tr}(\hat\Omega_I\hat V).}
#' The independence fit is obtained from the stored response, design matrix,
#' weights, offset, and family, so the original data object does not need to be
#' re-evaluated.
#'
#' QICC applies the finite-cluster correction of Hardin and Hilbe (2013) to
#' QIC:
#' \deqn{\mathrm{QICC} = \mathrm{QIC} -
#' \frac{2(p + m)(p + m + 1)}{N - p - m - 1},}
#' where \eqn{p} is the number of regression parameters, \eqn{m} is the
#' number of estimated working-association parameters, and \eqn{N} is the
#' number of independent clusters. Thus, \code{corstr = "fixed"} or
#' \code{orstr = "fixed"} contributes \eqn{m = 0}. Hardin and Hilbe note
#' that the same correction may
#' alternatively be applied to QICHH. In \code{geer}, the reported
#' \code{QICC} uses QIC as the base criterion. QICC is returned as
#' \code{NA} when \eqn{N - p - m - 1 <= 0}.
#'
#' The Gaussian pseudolikelihood-based criteria use the same fitted working
#' covariance matrices as \code{GPC}. If \eqn{n^{\star}} is the number of
#' observations, \eqn{p} the number of regression parameters, \eqn{m} the
#' number of estimated working-association parameters, and \eqn{N} the number
#' of independent clusters, define
#' \deqn{D_G = n^{\star}\log(2\pi) - 2\,\mathrm{GPC}.}
#' Then
#' \deqn{\mathrm{AGPC} = D_G + 2(p+m)}
#' and
#' \deqn{\mathrm{SGPC} = D_G + \log(N)(p+m).}
#' Fixed working correlation or odds-ratio structures contribute \eqn{m=0}.
#' Throughout, \eqn{m} counts only the working-association parameters that
#' were estimated from the data. Independence structures and supplied
#' structures, that is \code{corstr = "fixed"} or \code{orstr = "fixed"},
#' therefore contribute \eqn{m = 0} to QICC, GESSC, AGPC, and SGPC: a
#' correlation or odds-ratio matrix that the user provides costs no degrees of
#' freedom. Note that this assumes the supplied structure was not itself
#' obtained from the same data, in which case the complexity of the model would
#' be understated.
#'
#' The additive normalizing constant in \eqn{D_G} is included so AGPC and SGPC
#' are on the conventional scale used in the literature and in
#' \pkg{glmtoolbox}.
#'
#' PT, WR, and RMR are the eigenvalue-based criteria of Jang (2011). Let
#' \eqn{\lambda_1,\ldots,\lambda_p} be the generalized eigenvalues of the
#' covariance estimate \eqn{\hat V} with respect to \eqn{\hat\Omega_I^{-1}},
#' the model-based covariance under working independence, that is the
#' eigenvalues of \eqn{\hat V\hat\Omega_I}. Then
#' \deqn{\mathrm{PT} = \sum_{j} \frac{\lambda_j}{1+\lambda_j}, \qquad
#' \mathrm{WR} = \prod_{j} \frac{\lambda_j}{1+\lambda_j}, \qquad
#' \mathrm{RMR} = \max_{j} \frac{\lambda_j}{1+\lambda_j},}
#' which are respectively the trace, the determinant, and the largest
#' eigenvalue of \eqn{\hat V(\hat V + \hat\Omega_I^{-1})^{-1}}. Like CIC,
#' they compare the size of the sandwich covariance with that of a fixed
#' reference matrix, and smaller values indicate closer agreement. The
#' eigenvalues are obtained from a symmetric form built with the Cholesky
#' factor of \eqn{\hat\Omega_I}, so they are real by construction; all three
#' criteria are returned as \code{NA} when that factorization fails or when
#' any eigenvalue is not positive.
#'
#' GHYC and PAC compare the empirical residual covariance with the fitted
#' working covariance. Writing \eqn{K} for the number of clusters,
#' \eqn{\bar S = K^{-1}\sum_i (y_i-\hat\mu_i)(y_i-\hat\mu_i)^{\top}} and
#' \eqn{\bar V = K^{-1}\sum_i \hat V_i},
#' \deqn{\mathrm{GHYC} = \mathrm{tr}[(\bar S\bar V^{-1}-I)^2]}
#' and
#' \deqn{\mathrm{PAC} = |\det(\bar S)/\det(\bar V)-1|.}
#' Both criteria are defined only for balanced designs. The cluster-level
#' matrices being summed are conformable only when every cluster contributes
#' the same repeated positions, and Gosho, Hamada and Yoshimura (2011) and
#' Pardo and Alonso (2019) both assume a common cluster size. GHYC and PAC are
#' therefore returned as \code{NA} unless every cluster observes each repeated
#' position exactly once, and also when the averaged working covariance matrix
#' is singular.
#'
#' For EQIC, the variance and deviance contributions use the adjustment
#' \eqn{k = 1/6} recommended by Nelder and Pregibon and used by Wang and Hin.
#' For binomial models the dispersion is fixed at 1. For other families, if
#' \code{phi_fixed = TRUE} was used in \code{geewa()}, the fitted fixed
#' dispersion is retained. Otherwise, the EQIC dispersion is estimated as the
#' adjusted deviance divided by the number of observations. This same EQIC
#' dispersion is used in the independence information matrix entering its
#' covariance penalty.
#'
#' If the supplied models do not all have the same number of observations, a
#' warning is issued.
#'
#' @return
#' A data frame with one row per fitted model, holding the columns selected by
#' \code{criteria} in the order requested, followed by \code{Parameters}. The
#' available criteria are \code{QIC}, \code{QICHH}, \code{QICC},
#' \code{CIC}, \code{RJC}, \code{QICu}, \code{EQIC}, \code{GESSC},
#' \code{GPC}, \code{AGPC}, \code{SGPC}, \code{GHYC}, \code{PAC},
#' \code{PT}, \code{WR}, and \code{RMR}, as described in the Details
#' section. When more than one
#' model is supplied, row names are set to the deparsed model expressions.
#'
#' @references
#' Carey, V.J. and Wang, Y.G. (2011) Working covariance model selection for
#' generalized estimating equations. \emph{Statistics in Medicine},
#' \bold{30}, 3117--3124.
#'
#' Chaganty, N.R. and Shults, J. (1999) On eliminating the asymptotic bias in
#' the quasi-least squares estimate of the correlation parameter.
#' \emph{Journal of Statistical Planning and Inference}, \bold{76}, 145--161.
#'
#' Fu, L., Hao, Y. and Wang, Y.G. (2018) Working correlation structure
#' selection in generalized estimating equations. \emph{Computational
#' Statistics}, \bold{33}, 983--996.
#'
#' Gosho, M. (2014) Criteria to select a working correlation structure for the
#' generalized estimating equations method in SAS. \emph{Journal of
#' Statistical Software, Code Snippets}, \bold{57}, 1--10.
#'
#' Gosho, M., Hamada, C. and Yoshimura, I. (2011) Criterion for the selection
#' of a working correlation structure in the generalized estimating equation
#' approach for longitudinal balanced data. \emph{Communications in Statistics
#' - Theory and Methods}, \bold{40}, 3839--3856.
#'
#' Hardin, J.W. and Hilbe, J.M. (2003) \emph{Generalized Estimating Equations}.
#' Chapman and Hall/CRC, Boca Raton.
#'
#' Hardin, J.W. and Hilbe, J.M. (2013) \emph{Generalized Estimating
#' Equations}, 2nd Edition. Chapman and Hall/CRC, Boca Raton.
#'
#' Hin, L.Y., Carey, V.J. and Wang, Y.G. (2007) Criteria for
#' working-correlation-structure selection in GEE: assessment via simulation.
#' \emph{The American Statistician}, \bold{61}, 360--364.
#'
#' Hin, L.Y. and Wang, Y.G. (2009) Working-correlation-structure
#' identification in generalized estimating equations. \emph{Statistics in
#' Medicine}, \bold{28}, 642--658.
#'
#' Jang, M.J. (2011) \emph{Working correlation selection in generalized
#' estimating equations}. PhD dissertation, University of Iowa. The three
#' eigenvalue-based criteria are also described in Sections 3.6 to 3.8 of
#' Pardo and Alonso (2019).
#'
#' Pan, W. (2001) Akaike's information criterion in generalized estimating
#' equations. \emph{Biometrics}, \bold{57}, 120--125.
#'
#' Pardo, M.C. and Alonso, R. (2019) Working correlation structure selection
#' in GEE analysis. \emph{Statistical Papers}, \bold{60}, 1447--1467.
#'
#' Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
#' parameters in semiparametric generalized linear models for cluster correlated
#' data. \emph{Biometrika}, \bold{77}, 485--497.
#'
#' Shults, J. and Chaganty, N.R. (1998) Analysis of serially correlated data
#' using quasi-least squares. \emph{Biometrics}, \bold{54}, 1622--1630.
#'
#' Shults, J., Sun, W., Tu, X., Kim, H., Amsterdam, J., Hilbe, J.M. and
#' Ten-Have, T. (2009) A comparison of several approaches for choosing between
#' working correlation structures in generalized estimating equation analysis
#' of longitudinal binary data. \emph{Statistics in Medicine}, \bold{28},
#' 2338--2355.
#'
#' Vanegas, L.H., Rondon, L.M. and Paula, G.A. (2023) Generalized Estimating
#' Equations using the new R package glmtoolbox. \emph{The R Journal},
#' \bold{15}, 105--133.
#'
#' Wang, Y.G. and Hin, L.Y. (2010) Modeling strategies in longitudinal data
#' analysis: covariate, variance function and correlation structure selection.
#' \emph{Computational Statistics & Data Analysis}, \bold{54}, 3359--3370.
#'
#' Zhu, X. and Zhu, Z. (2013) Comparison of criteria to select working
#' correlation matrix in generalized estimating equations. \emph{Chinese
#' Journal of Applied Probability and Statistics}, \bold{29}, 515--530.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{step_p}},
#'   \code{\link{glance.geer}}, \code{\link{anova.geer}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' ## Single model
#' data("epilepsy", package = "geer")
#' fit <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   corstr = "exchangeable"
#' )
#' geecriteria(fit)
#'
#' ## Compare working correlation structures
#' fit_ind <- update(fit, corstr = "independence")
#' fit_ar1 <- update(fit, corstr = "ar1")
#' geecriteria(fit_ind, fit, fit_ar1)
#'
#' ## Report a selection of criteria, in the order requested
#' geecriteria(fit_ind, fit, fit_ar1, criteria = c("CIC", "QIC"))
#' geecriteria(fit_ind, fit, fit_ar1, criteria = "cic")
#'
#' ## Compare estimation methods
#' data("cerebrovascular", package = "geer")
#' fit_gee <- geewa_binary(
#'   formula = ecg ~ factor(period) * treatment,
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable",
#'   method = "gee"
#' )
#' fit_brgee <- update(fit_gee, method = "brgee-robust")
#' geecriteria(fit_gee, fit_brgee, cov_type = "robust")
#'
#' @export
geecriteria <- function(object,
                        ...,
                        criteria = "all",
                        cov_type = geer_criteria_cov_type_choices,
                        digits = 2) {
  cov_type <- match.arg(cov_type)
  criteria <- normalize_geer_criteria(criteria)
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
  out_list <- lapply(
    models,
    compute_gee_criteria,
    cov_type = cov_type,
    digits = NULL,
    criteria = criteria
  )
  ans <- do.call(rbind, out_list)
  ans <- ans[, c(criteria, "Parameters"), drop = FALSE]
  ans[, criteria] <- lapply(ans[, criteria, drop = FALSE], round, digits = digits)
  ans[, "Parameters"] <- as.integer(ans[, "Parameters"])
  if (length(models) > 1L && nrow(ans) == length(models)) {
    call_expr <- match.call(expand.dots = FALSE)
    exprs <- c(list(call_expr[[2L]]), as.list(call_expr$...))
    ## make.unique() guards against the same expression being supplied twice,
    ## which would otherwise fail because a data frame cannot carry duplicate
    ## row names.
    rownames(ans) <- make.unique(vapply(
      exprs,
      function(expr) paste(deparse(expr, width.cutoff = 500L), collapse = ""),
      character(1)
    ))
  }
  ans
}
