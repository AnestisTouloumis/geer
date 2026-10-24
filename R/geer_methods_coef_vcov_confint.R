#' @title
#' Extract Variance-Covariance Matrix from a geer Object
#'
#' @method vcov geer
#'
#' @description
#' Extracts the variance-covariance matrix of the estimated regression
#' parameters from a fitted \code{geer} object.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for inference: \code{"bias-corrected"}, \code{"robust"},
#'   \code{"df-adjusted"}, \code{"jackknife"}, \code{"naive"},
#'   \code{"kauermann-carroll"}, \code{"mancl-derouen"},
#'   \code{"ford-westgate"}, \code{"fay-graubard"}, \code{"fan-zhang-zhang"},
#'   \code{"pan"}, \code{"gosho-sato-takeuchi"}, \code{"wang-long"}, or
#'   \code{"westgate-burchett"}. See Details. Defaults to
#'   \code{"bias-corrected"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @details
#' Let \eqn{K} be the number of clusters, \eqn{p} the number of regression
#' parameters, \eqn{D_i} the derivative of the mean vector of cluster \eqn{i}
#' with respect to the regression parameters, \eqn{V_i} its working covariance
#' matrix, \eqn{e_i} its vector of response residuals,
#' \eqn{\Omega = \sum_i D_i^\top V_i^{-1} D_i}{Omega = sum_i D_i' V_i^-1 D_i},
#' \eqn{U_i = D_i^\top V_i^{-1} e_i}{U_i = D_i' V_i^-1 e_i}, and
#' \eqn{H_i = D_i \Omega^{-1} D_i^\top V_i^{-1}}{H_i = D_i Omega^-1 D_i' V_i^-1}.
#' The covariance estimator is selected by \code{cov_type}:
#' \describe{
#'   \item{\code{"bias-corrected"}}{the bias-corrected sandwich estimator of
#'     Morel, Bokossa, and Neerchal (2003).}
#'   \item{\code{"robust"}}{the sandwich estimator of Liang and Zeger (1986),
#'     \eqn{\Omega^{-1} M \Omega^{-1}}{Omega^-1 M Omega^-1} with
#'     \eqn{M = \sum_i U_i U_i^\top}{M = sum_i U_i U_i'}.}
#'   \item{\code{"df-adjusted"}}{the sandwich estimator with the
#'     degrees-of-freedom adjustment of MacKinnon and White (1985).}
#'   \item{\code{"jackknife"}}{the leave-one-cluster-out jackknife estimator.}
#'   \item{\code{"naive"}}{the model-based estimator of Liang and Zeger
#'     (1986).}
#'   \item{\code{"kauermann-carroll"}}{Kauermann and Carroll (2001): the
#'     robust estimator with \eqn{e_i} replaced by
#'     \eqn{(I - H_i)^{-1/2} e_i}{(I - H_i)^(-1/2) e_i}, using the principal
#'     matrix square root.}
#'   \item{\code{"mancl-derouen"}}{Mancl and DeRouen (2001): the robust
#'     estimator with \eqn{e_i} replaced by
#'     \eqn{(I - H_i)^{-1} e_i}{(I - H_i)^-1 e_i}.}
#'   \item{\code{"ford-westgate"}}{Ford and Westgate (2017): the average of the
#'     Kauermann-Carroll and Mancl-DeRouen estimators.}
#'   \item{\code{"fay-graubard"}}{Fay and Graubard (2001): the robust
#'     estimator with \eqn{U_i} replaced by \eqn{F_i U_i}, where \eqn{F_i} is
#'     diagonal with entries
#'     \eqn{\{1 - \min(0.75, Q_{jj})\}^{-1/2}}{(1 - min(0.75, Q_jj))^(-1/2)}
#'     and \eqn{Q = D_i^\top V_i^{-1} D_i \Omega^{-1}}{Q = D_i' V_i^-1 D_i
#'     Omega^-1}.}
#'   \item{\code{"fan-zhang-zhang"}}{Fan, Zhang, and Zhang (2013): the robust
#'     estimator with \eqn{e_i e_i^\top}{e_i e_i'} replaced by
#'     \eqn{e_i e_i^\top - D_i \Omega^{-1} (\sum_{j \neq i} U_j U_j^\top)
#'     \Omega^{-1} D_i^\top}{e_i e_i' - D_i Omega^-1 (sum_(j != i) U_j U_j')
#'     Omega^-1 D_i'}, pre- and post-multiplied by
#'     \eqn{(I - H_i)^{-1}}{(I - H_i)^-1}.}
#'   \item{\code{"pan"}}{Pan (2001): the residual outer products are pooled
#'     over clusters after standardization by
#'     \eqn{A_i = \mathrm{diag}(V_i)}{A_i = diag(V_i)}, so that \eqn{e_i
#'     e_i^\top}{e_i e_i'} is replaced by \eqn{A_i^{1/2} C A_i^{1/2}} with
#'     \eqn{C = K^{-1} \sum_j A_j^{-1/2} e_j e_j^\top A_j^{-1/2}}{C = K^-1
#'     sum_j A_j^(-1/2) e_j e_j' A_j^(-1/2)}.}
#'   \item{\code{"gosho-sato-takeuchi"}}{Gosho, Sato, and Takeuchi (2014): the
#'     Pan estimator multiplied by \eqn{K / (K - p)}.}
#'   \item{\code{"wang-long"}}{Wang and Long (2011): the Pan estimator with
#'     \eqn{e_j} replaced by \eqn{(I - H_j)^{-1} e_j}{(I - H_j)^-1 e_j} in
#'     \eqn{C}.}
#'   \item{\code{"westgate-burchett"}}{Westgate and Burchett (2016): the Pan
#'     estimator with \eqn{e_j} replaced by
#'     \eqn{(I - H_j)^{-1/2} e_j}{(I - H_j)^(-1/2) e_j} in \eqn{C}. With
#'     \eqn{(I - H_j)^{-1} e_j}{(I - H_j)^-1 e_j} instead, the Mancl-DeRouen
#'     version of this estimator, it coincides with the Wang-Long estimator.}
#' }
#' The estimators are evaluated at the estimates of the fitted model, whatever
#' the estimation method, and for any family and link, with \eqn{D_i} and
#' \eqn{V_i} as used by \code{\link{geewa}} or \code{\link{geewa_binary}}. In
#' the Pan, Gosho-Sato-Takeuchi, Wang-Long, and Westgate-Burchett estimators,
#' \eqn{A_i} contains the scale parameter and the prior weights, so that
#' weighted fits are standardized by the marginal variances they assume. The
#' scale parameter cancels in the sandwich. With unit prior weights and the
#' logit link, the nine small-sample estimators reproduce those of the same
#' name in the \pkg{geessbin} package, which standardizes by
#' \eqn{\mu(1 - \mu)}{mu (1 - mu)}.
#'
#' The bias-corrected, degrees-of-freedom adjusted, and Gosho-Sato-Takeuchi
#' estimators require more clusters than regression parameters: otherwise they
#' give an error. The
#' Pan, Gosho-Sato-Takeuchi, Wang-Long, and Westgate-Burchett estimators pool
#' the clusters and require a balanced design, in which every cluster observes
#' each repeated position once: otherwise they give an error. The
#' Kauermann-Carroll, Mancl-DeRouen, Ford-Westgate, Fan-Zhang-Zhang,
#' Wang-Long, and Westgate-Burchett estimators use
#' \eqn{(I - H_i)^{-1}}{(I - H_i)^-1} and give an error that names the
#' clusters with leverage 1, for which \eqn{I - H_i} is singular. All nine
#' small-sample estimators give an error when the naive covariance is not
#' positive definite.
#'
#' The score-type tests of \code{\link{anova.geer}}, \code{\link{add1.geer}},
#' \code{\link{drop1.geer}}, \code{\link{step_p}}, and
#' \code{\link{mcar_logistic_test}} evaluate the selected estimator at the fit
#' of the smaller model, with \eqn{D_i} and \eqn{V_i} of the larger model. For
#' the nine small-sample estimators this is an extension of the published
#' estimators, which are defined at the full-model fit. The jackknife estimator
#' is the exception: score-type tests use the jackknife estimator of the larger
#' model.
#'
#' With \eqn{K} clusters and \eqn{\hat\beta_{(i)}}{beta_(i)} the estimate obtained after
#' deleting cluster \eqn{i}, the jackknife estimator is
#' \deqn{\frac{K - 1}{K} \sum_{i=1}^K (\hat\beta_{(i)} - \bar\beta)
#' (\hat\beta_{(i)} - \bar\beta)^\top, \qquad
#' \bar\beta = \frac{1}{K}\sum_{i=1}^K \hat\beta_{(i)}.}{(K - 1)/K sum_i (beta_(i) - betabar)(beta_(i) - betabar)',   betabar = (1/K) sum_i beta_(i).}
#' Each \eqn{\hat\beta_{(i)}}{beta_(i)} is a full refit by the original estimation
#' method, not a one-step approximation, with the working association
#' parameters held at their full-data estimates. The dispersion is
#' re-estimated if it was estimated in the original fit. For the two-stage
#' methods (\code{"bcgee-*"}, \code{"opgee-jeffreys"}, and
#' \code{"hpgee-jeffreys"}) the first stage is repeated in each refit. At
#' least two clusters are required, and an error is signaled if a refit does
#' not converge. The jackknife and the nine small-sample estimates are cached
#' in the fitted object.
#'
#' @return
#' A square numeric matrix of estimated covariances between regression
#' coefficients. Rows and columns are named according to the coefficient names
#' returned by \code{\link{coef.geer}}.
#'
#' @references
#' Fan, C., Zhang, D. and Zhang, C.H. (2013) A comparison of bias-corrected
#' covariance estimators for generalized estimating equations. \emph{Journal
#' of Biopharmaceutical Statistics}, \bold{23}, 1172--1187.
#'
#' Fay, M.P. and Graubard, B.I. (2001) Small-sample adjustments for Wald-type
#' tests using sandwich estimators. \emph{Biometrics}, \bold{57}, 1198--1206.
#'
#' Ford, W.P. and Westgate, P.M. (2017) Improved standard error estimator for
#' maintaining the validity of inference in cluster randomized trials with a
#' small number of clusters. \emph{Biometrical Journal}, \bold{59}, 478--495.
#'
#' Gosho, M., Sato, T. and Takeuchi, H. (2014) Robust covariance estimator for
#' small-sample adjustment in the generalized estimating equations: a
#' simulation study. \emph{Science Journal of Applied Mathematics and
#' Statistics}, \bold{2}, 20--25.
#'
#' Kauermann, G. and Carroll, R.J. (2001) A note on the efficiency of sandwich
#' covariance matrix estimation. \emph{Journal of the American Statistical
#' Association}, \bold{96}, 1387--1396.
#'
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, \bold{73}, 13--22.
#'
#' MacKinnon, J.G. and White, H. (1985) Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics}, \bold{29}, 305--325.
#'
#' Mancl, L.A. and DeRouen, T.A. (2001) A covariance estimator for GEE with
#' improved small-sample properties. \emph{Biometrics}, \bold{57}, 126--134.
#'
#' Morel, J.G., Bokossa, M.C. and Neerchal, N.K. (2003) Small sample
#' correction for the variance of GEE estimators. \emph{Biometrical Journal},
#' \bold{45}, 395--409.
#'
#' Pan, W. (2001) On the robust variance estimator in generalised estimating
#' equations. \emph{Biometrika}, \bold{88}, 901--906.
#'
#' Wang, M. and Long, Q. (2011) Modified robust variance estimator for
#' generalized estimating equations with improved small-sample performance.
#' \emph{Statistics in Medicine}, \bold{30}, 1278--1291.
#'
#' Westgate, P.M. and Burchett, W.W. (2016) Improving power in small-sample
#' longitudinal studies when using generalized estimating equations.
#' \emph{Statistics in Medicine}, \bold{35}, 3733--3744.
#'
#' @seealso \code{\link{coef.geer}}, \code{\link{confint.geer}},
#'   \code{\link{summary.geer}}, \code{\link{tidy.geer}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' vcov(fit)
#' vcov(fit, cov_type = "kauermann-carroll")
#'
#' @export
vcov.geer <- function(object, cov_type = "bias-corrected", ...) {
  check_unused_dots(list(...), "vcov.geer")
  object <- check_geer_object(object)
  cov_type <- match_cov_type(cov_type)
  if (cov_type %in% geer_ss_cov_types) {
    return(
      select_small_sample_covariance(
        covariances = get_cached_ss_covariances(object),
        cov_type = cov_type,
        object = object,
        coef_no = object$rank,
        context = "vcov.geer"
      )
    )
  }
  switch(
    cov_type,
    robust = object$robust_covariance,
    naive = object$naive_covariance,
    `bias-corrected` = {
      check_bc_covariance(
        object$clusters_no,
        length(object$coefficients),
        "vcov.geer"
      )
      object$bias_corrected_covariance
    },
    jackknife = get_cached_jackknife_cov(object),
    `df-adjusted` = compute_df_adjusted_covariance(
      robust_covariance = object$robust_covariance,
      clusters_no = object$clusters_no,
      coef_no = object$rank,
      context = "vcov.geer"
    )
  )
}


#' @title
#' Extract Model Coefficients from a geer Object
#'
#' @method coef geer
#'
#' @description
#' Extracts the estimated regression coefficients from a fitted \code{geer}
#' object. \code{coefficients} is an alias for \code{coef}.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @return
#' A named numeric vector of estimated regression coefficients. The names
#' correspond to the columns of the model matrix.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{confint.geer}},
#'   \code{\link{summary.geer}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' coef(fit)
#'
#' @export
coef.geer <- function(object, ...) {
  check_unused_dots(list(...), "coef.geer")
  object <- check_geer_object(object)
  object$coefficients
}


#' @title
#' Confidence Intervals for Model Parameters from a geer Object
#'
#' @method confint geer
#'
#' @description
#' Computes Wald-type confidence intervals for one or more regression
#' parameters from a fitted \code{geer} object.
#'
#' @inheritParams vcov.geer
#' @inheritParams stats::confint
#' @param level a single number strictly between 0 and 1 specifying the
#'   confidence level. Defaults to \code{0.95}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for the standard errors: one of the estimators described in
#'   \code{\link{vcov.geer}}. Defaults to \code{"bias-corrected"}.
#'
#' @details
#' The intervals are
#' \eqn{\hat\beta \pm z_{1-\alpha/2}\,\mathrm{SE}(\hat\beta)}{beta +/- z_(1 - alpha/2) SE(beta)}, with standard
#' errors from \code{vcov(object, cov_type = cov_type)}.
#'
#' @return
#' A matrix with columns giving lower and upper confidence limits for each
#' parameter. The columns are labeled by the corresponding tail probabilities
#' in percent, for example \code{"2.5\%"} and \code{"97.5\%"} when
#' \code{level = 0.95}.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{coef.geer}},
#'   \code{\link{tidy.geer}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' confint(fit)
#' confint(fit, parm = "treatmentactive", cov_type = "robust")
#'
#' @export
confint.geer <- function(object,
                         parm,
                         level = 0.95,
                         cov_type = "bias-corrected",
                         ...) {
  check_unused_dots(list(...), "confint.geer")
  object <- check_geer_object(object)
  check_probability_open(level, "level")
  cov_type <- match_cov_type(cov_type)
  beta <- object$coefficients
  parm <- match_parm(if (missing(parm)) NULL else parm, names(beta))
  conf_probs <- (1 - level) / 2
  conf_probs <- c(conf_probs, 1 - conf_probs)
  pct <- build_percent_labels(conf_probs, 3)
  percentiles <- stats::qnorm(conf_probs)
  ans <- matrix(
    NA_real_,
    nrow = length(parm),
    ncol = 2L,
    dimnames = list(parm, pct)
  )
  vcov_matrix <- stats::vcov(object, cov_type = cov_type)
  standard_errors <- compute_se_or_stop(
    diag(vcov_matrix)[parm],
    context = "confint.geer"
  )
  ans[] <- beta[parm] + standard_errors %o% percentiles
  ans
}
