#' @title
#' Fitting (Adjusted) Generalized Estimating Equations
#'
#' @description
#' Fits a marginal model for clustered or repeated responses by generalized
#' estimating equations (GEE), using ordinary, bias-reduced, bias-corrected, or
#' Jeffreys-type penalized GEE.
#'
#' @param formula formula of the form \code{response ~ predictors} describing
#'   the marginal mean model. Terms \code{offset()} are allowed.
#' @param family family object, family function, or character string naming a
#'   family, specifying the marginal variance and link functions. Supported
#'   families are \code{gaussian}, \code{binomial}, \code{poisson},
#'   \code{Gamma}, \code{inverse.gaussian}, \code{quasi}, \code{quasibinomial},
#'   and \code{quasipoisson}. Defaults to \code{gaussian(link = "identity")}.
#' @param data optional data frame in long format, with one row per
#'   observation, containing the variables in \code{formula}, \code{id},
#'   \code{repeated}, \code{weights}, and \code{offset}. Defaults to
#'   \code{parent.frame()}.
#' @param id variable identifying the clusters.
#' @param repeated optional variable identifying the order of the observations
#'   within each cluster. Its values must be unique within each cluster. They
#'   are replaced by their rank among the sorted distinct values, or by the
#'   level number for a factor, so only their order matters and not their
#'   spacing. A character variable gives an error, since it has no natural
#'   order. If omitted, positions are assigned from the row order within each
#'   cluster in the original input, before subsetting and missing-value removal.
#'   Retained observations keep these positions, including gaps.
#' @param control list of control parameters as returned by
#'   \code{\link{geer_control}}. Defaults to \code{geer_control(...)}.
#' @param corstr character string specifying the working correlation structure:
#'   \code{"independence"}, \code{"exchangeable"}, \code{"ar1"},
#'   \code{"m-dependent"}, \code{"unstructured"}, \code{"toeplitz"}, or
#'   \code{"fixed"}. Defaults to \code{"independence"}.
#' @param Mv positive integer giving the order of the \code{"m-dependent"}
#'   structure. Ignored for other structures. Defaults to \code{1}.
#' @param method character string specifying the estimation method: ordinary GEE
#'   (\code{"gee"}), bias-reduced GEE (\code{"brgee-naive"},
#'   \code{"brgee-robust"}, or \code{"brgee-empirical"}), bias-corrected GEE
#'   (\code{"bcgee-naive"}, \code{"bcgee-robust"}, or \code{"bcgee-empirical"}),
#'   Jeffreys-type penalized GEE (\code{"pgee-jeffreys"}), or its one-step
#'   (\code{"opgee-jeffreys"}) and hybrid one-step (\code{"hpgee-jeffreys"})
#'   versions. See Details. Defaults to \code{"gee"}.
#' @param weights optional numeric vector of strictly positive, finite
#'   observation weights. Defaults to weights of 1.
#' @param beta_start optional numeric vector of starting values for the
#'   regression parameters. Defaults to \code{NULL}, in which case they are
#'   obtained from a generalized linear model fitted by
#'   \code{\link[brglm2]{brglmFit}}.
#' @param offset optional a priori known component of the linear predictor:
#'   a single number or a numeric vector with one value per observation. It is
#'   added to any \code{offset()} terms in \code{formula}.
#' @param control_glm optional list of control parameters passed to
#'   \code{\link[brglm2]{brglm_control}} for the starting values. Ignored when
#'   \code{beta_start} is supplied. Defaults to \code{list()}.
#' @param use_p logical indicating whether the number of regression parameters
#'   is subtracted from the denominators of the estimators of the dispersion
#'   and working correlation parameters. Defaults to \code{TRUE}.
#' @param alpha_vector numeric vector of working correlation parameters for
#'   \code{corstr = "fixed"}, of length \code{choose(T, 2)}, where \eqn{T} is
#'   the number of distinct supplied \code{repeated} values, or the largest
#'   retained original position if \code{repeated} is omitted. The parameters
#'   must give a positive definite correlation matrix. Ignored for other
#'   structures. Defaults to \code{NULL}.
#' @param phi_fixed logical indicating whether the dispersion parameter is
#'   fixed at \code{phi_value}. Defaults to \code{FALSE}.
#' @param phi_value strictly positive number giving the fixed value of the
#'   dispersion parameter when \code{phi_fixed = TRUE}. Defaults to \code{1}.
#' @param subset optional vector specifying a subset of observations to be used
#'   in the fit, evaluated in \code{data} as in
#'   \code{\link[stats]{model.frame}}.
#' @param na.action function specifying the handling of missing values, as in
#'   \code{\link[stats]{model.frame}}. Defaults to
#'   \code{getOption("na.action")}.
#' @param ... arguments of \code{\link{geer_control}}, given by their full
#'   names, used to build \code{control} when it is not supplied. Supplying
#'   \code{...} together with \code{control} is an error.
#'
#' @details
#' Working covariance matrices must be positive definite. A non-positive-
#' definite working covariance is an error at initialization. If it arises
#' after an accepted fitting state, the existing numerical-failure warning
#' and last-accepted-estimate policy applies. General matrix calculations
#' that do not require a covariance matrix may still use LU decomposition.
#'
#' Method, family, link, and association names must be exact character strings.
#' Family objects must contain the standard callable components and an
#' initialization expression. When dispersion is fixed, \code{phi_value} must be
#' a positive finite numeric scalar.
#'
#' Inactive arguments are ignored without validation: \code{Mv} outside the
#' m-dependent structure, \code{alpha_vector} outside the fixed structure,
#' \code{phi_value} when dispersion is estimated, and \code{control_glm} when
#' \code{beta_start} is supplied.
#'
#' The estimation methods are:
#' \describe{
#'   \item{\code{"gee"}}{the estimating equations of Liang and Zeger (1986).}
#'   \item{\code{"brgee-*"}}{bias-reduced GEE, which add an adjustment term to
#'     the estimating equations based on the naive, robust, or empirical
#'     covariance matrix (Touloumis, 2026a).}
#'   \item{\code{"bcgee-*"}}{bias-corrected GEE, which apply a one-step
#'     correction to the converged GEE estimate (Touloumis, 2026a).}
#'   \item{\code{"pgee-jeffreys"}}{GEE penalized by a Jeffreys-type penalty
#'     (Touloumis, 2026b), developed there for binary responses with odds-ratio
#'     working structures (\code{\link{geewa_binary}}). Other families and
#'     correlation structures extend that paper, whose finiteness result covers
#'     only binary responses under working independence.}
#'   \item{\code{"opgee-jeffreys"}, \code{"hpgee-jeffreys"}}{one penalized, or
#'     one ordinary, scoring step from the penalized estimate under working
#'     independence.}
#' }
#'
#' With adjustment \eqn{A} added to the estimating equations \eqn{U}, the
#' bias-corrected methods use the plus scoring step
#' \eqn{\hat\beta_{BC} = \hat\beta_{GEE} +
#' H^{-1} A}{beta_BC = beta_GEE + H^(-1) A}.
#' Here \eqn{H} is the working information for the naive and robust methods,
#' and \eqn{-\partial U/\partial\beta^\top}{-dU/d beta'} for the empirical
#' method. This differs from the printed minus sign with positive working
#' information in Touloumis (2026a, Equation 6, Algorithm 1, and the robust and
#' naive formulas in Section 2.5). The plus sign agrees with the adjusted
#' equation \eqn{U + A = 0} and, under independent Bernoulli sampling with the
#' logit link, the one-step Firth adjustment.
#'
#' For the empirical method, the paper defines \eqn{j = DU}, so its outer
#' minus sign is equivalent to the plus step with \eqn{H = -j}. The first
#' contraction in its Equation (10) must also change sign to reproduce the
#' implemented adjustment under that definition of \eqn{j}. The second
#' contraction is unchanged.
#'
#' Observations are sorted by \code{id} and \code{repeated}. The cluster
#' identifiers and supplied \code{repeated} values are recoded to consecutive
#' integers. Implicit positions retain gaps after subsetting or missing-value
#' removal. The quasi families are mapped to the
#' corresponding standard families: \code{quasibinomial} to \code{binomial},
#' \code{quasipoisson} to \code{poisson}, and \code{quasi} with variance
#' \code{"constant"}, \code{"mu(1-mu)"}, \code{"mu"}, \code{"mu^2"}, or
#' \code{"mu^3"} to \code{gaussian}, \code{binomial}, \code{poisson},
#' \code{Gamma}, or \code{inverse.gaussian}. The dispersion parameter is
#' estimated for every family unless \code{phi_fixed = TRUE}.
#'
#' For \code{corstr = "ar1"}, the correlation parameter is estimated from the
#' products of Pearson residuals at consecutive \code{repeated} values, pooled
#' over clusters, so it can differ from the estimates of other packages.
#'
#' @section Convergence and separation:
#' For the single-stage methods (\code{"gee"}, \code{"brgee-*"}, and
#' \code{"pgee-jeffreys"}), \code{alpha} and \code{phi} are re-estimated at
#' every iteration from the current estimates of the regression parameters,
#' which is the joint updating of Touloumis (2026a, Algorithm 1), unless they
#' are fixed by \code{corstr = "fixed"} or \code{phi_fixed = TRUE}.
#'
#' For the two-stage methods (\code{"bcgee-*"}, \code{"opgee-jeffreys"}, and
#' \code{"hpgee-jeffreys"}), \code{alpha} and \code{phi} are estimated once at
#' the first-stage estimate and held fixed in the correction step, and
#' \code{converged} is \code{TRUE} unless the correction step fails, in which
#' case the first-stage estimates are returned with a warning. If the first
#' stage does not converge, the two-stage methods give an error.
#'
#' If \code{converged} is \code{FALSE}, the results are those of the last
#' accepted iteration and should not be interpreted; printing the fit or its
#' summary then displays a notice. \code{converged = TRUE} only means that the
#' stopping criterion was met: with separated binary data or all-zero counts
#' the estimates can be very large, which is signaled by a warning that fitted
#' probabilities are numerically 0 or 1, or, for a converged fit, within 1e-12
#' of 0 or 1, or that fitted rates are numerically 0.
#'
#' In the bias-reduced and penalized methods the adjustment term is multiplied
#' by the dispersion parameter, whose estimate tends to 0 under complete
#' separation, so with the default estimated dispersion no method protects
#' against infinite estimates. For binary responses with separation or rare
#' events, use \code{phi_fixed = TRUE} or \code{\link{geewa_binary}}. Then
#' \code{"pgee-jeffreys"} with \code{corstr = "independence"} gives finite
#' estimates under separation for the logit, probit, cauchit, and cloglog links
#' (Touloumis, 2026b), and so does \code{"brgee-naive"} with the logit link,
#' for which it coincides with \code{"pgee-jeffreys"} at the default
#' \code{jeffreys_power = 0.5}. With other working correlation structures
#' finiteness is not guaranteed, and the moment estimates of the correlations
#' can be unstable under separation. \code{"opgee-jeffreys"} and
#' \code{"hpgee-jeffreys"} start from the finite independence estimate and are
#' finite whenever their single scoring step is defined. \code{"brgee-robust"}
#' and \code{"brgee-empirical"} do not keep the estimates finite, and the
#' bias-corrected methods fail because they start from the ordinary GEE fit.
#'
#' @return
#' An object of class \code{"geer"}, a list with components:
#' \item{coefficients}{a named vector of estimated regression coefficients.}
#' \item{residuals}{the working residuals.}
#' \item{fitted.values}{the fitted mean values, obtained by transforming the
#'       linear predictors by the inverse of the link function.}
#' \item{rank}{the numeric rank of the fitted model.}
#' \item{family}{the \code{\link[stats]{family}} object used.}
#' \item{linear.predictors}{the linear fit on the link scale.}
#' \item{iter}{the number of iterations used.}
#' \item{prior.weights}{the weights initially supplied, a vector of 1s if none
#'       were.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{y}{the response vector.}
#' \item{x}{the model matrix.}
#' \item{qr}{the QR decomposition of the model matrix, used for estimability
#'       checking.}
#' \item{id}{the recoded cluster identifier vector.}
#' \item{repeated}{the recoded within-cluster ordering vector.}
#' \item{cluster_labels}{the cluster identifiers supplied, in the order of
#'       their recoded values in \code{id}.}
#' \item{repeated_labels}{the values of \code{repeated} supplied, in the order
#'       of their recoded values in \code{repeated}, or \code{NULL} if
#'       \code{repeated} was not supplied.}
#' \item{converged}{logical indicating whether the algorithm converged.}
#' \item{call}{the matched call.}
#' \item{formula}{the formula supplied.}
#' \item{terms}{the \code{\link[stats]{terms}} object used.}
#' \item{data}{the data argument.}
#' \item{offset}{the offset vector used.}
#' \item{control}{the \code{\link{geer_control}} list used.}
#' \item{method}{character string identifying the estimation method used.}
#' \item{fit_function}{character string identifying the fitting function used,
#'       either \code{"geewa"} or \code{"geewa_binary"}.}
#' \item{use_p}{logical indicating whether the degrees-of-freedom correction
#'       was used; \code{\link{geewa}} fits only.}
#' \item{phi_fixed}{logical indicating whether the dispersion parameter was
#'       fixed; \code{\link{geewa}} fits only.}
#' \item{contrasts}{the contrasts used.}
#' \item{xlevels}{a record of the levels of the factors used in fitting.}
#' \item{na.action}{the \code{na.action} attribute of the model frame, present
#'       only when observations were removed.}
#' \item{naive_covariance}{the model-based (naive) covariance matrix.}
#' \item{robust_covariance}{the sandwich (robust) covariance matrix.}
#' \item{bias_corrected_covariance}{the bias-corrected covariance matrix.}
#' \item{association_structure}{the name of the working association structure.}
#' \item{alpha}{the estimated or fixed working association parameters:
#'       correlations for \code{geewa} fits, 0 under independence, and odds
#'       ratios for \code{\link{geewa_binary}} fits, 1 under independence.}
#' \item{phi}{the estimated or fixed dispersion parameter, 1 for
#'       \code{\link{geewa_binary}} fits.}
#' \item{obs_no}{the number of observations used in fitting.}
#' \item{clusters_no}{the number of clusters.}
#' \item{min_cluster_size}{the minimum cluster size.}
#' \item{max_cluster_size}{the maximum cluster size.}
#' Other components are for internal use.
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, \bold{73}, 13--22.
#'
#' Touloumis, A. (2026a) Bias-reduced GEE via adjusted estimating equations,
#' with odds-ratio extensions. \emph{arXiv preprint} arXiv:2606.16043.
#'
#' Touloumis, A. (2026b) Jeffreys-type penalized GEE for correlated binary data
#' with an odds-ratio parameterization. \emph{arXiv preprint} arXiv:2606.16058.
#'
#' @seealso \code{\link{geewa_binary}}, \code{\link{geer_control}},
#'   \code{\link{summary.geer}}, \code{\link{vcov.geer}},
#'   \code{\link{anova.geer}}, \code{\link{step_p}}, \code{\link{geecriteria}},
#'   \code{\link{runs_test}}.
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
#'
#' ## Bias-reduced GEE
#' fit_br <- update(fit, method = "brgee-robust")
#' summary(fit_br)
#'
#' @export
geewa <- function(formula,
                  family = stats::gaussian(link = "identity"),
                  data = parent.frame(),
                  id,
                  repeated,
                  control = geer_control(...),
                  corstr = "independence",
                  Mv = 1L,
                  method = "gee",
                  weights,
                  beta_start = NULL,
                  offset,
                  control_glm = list(),
                  use_p = TRUE,
                  alpha_vector = NULL,
                  phi_fixed = FALSE,
                  phi_value = 1,
                  subset,
                  na.action,
                  ...) {
  if (!missing(control)) {
    check_unused_dots(list(...), "geewa")
  }
  call <- match.call(expand.dots = TRUE)
  mcall <- match.call(expand.dots = FALSE)
  caller_env <- parent.frame()
  family <- normalize_family(family)
  link <- family$link
  inputs <- build_geer_inputs(
    mcall = mcall,
    family = family,
    env = caller_env,
    control = control,
    method = method
  )
  model_frame <- inputs$model_frame
  y <- inputs$y
  weights <- inputs$weights
  id <- inputs$id
  repeated <- inputs$repeated
  offset <- inputs$offset
  model_matrix <- inputs$model_matrix
  model_terms <- inputs$model_terms
  xnames <- inputs$xnames
  qr_model_matrix <- inputs$qr_model_matrix
  control <- inputs$control
  method <- inputs$method
  tolerance <- control$tolerance
  corstr <- match_corstr(corstr)
  if (!identical(corstr, "m-dependent")) {
    Mv <- 1L
  } else {
    Mv <- check_integer_at_least(Mv, "Mv")
  }
  if (!identical(corstr, "fixed")) {
    alpha_vector <- 0
    alpha_fixed <- 0
  } else {
    if (is.null(alpha_vector)) {
      stop(
        "'alpha_vector' must be provided when 'corstr = \"fixed\"'",
        call. = FALSE
      )
    }
    if (!is.numeric(alpha_vector)) {
      stop("'alpha_vector' must be a numeric vector", call. = FALSE)
    }
    alpha_vector <- as.numeric(alpha_vector)
    repeated_max <- max(repeated)
    if (length(alpha_vector) != choose(repeated_max, 2)) {
      stop(
        "'alpha_vector' must be a numeric vector of length ",
        choose(repeated_max, 2),
        call. = FALSE
      )
    }
    if (any(!is.finite(alpha_vector))) {
      stop(
        "'alpha_vector' must be finite when 'corstr = \"fixed\"'",
        call. = FALSE
      )
    }
    fixed_eigenvalues <- eigen(
      get_correlation_matrix_cc(corstr, alpha_vector, repeated_max),
      symmetric = TRUE,
      only.values = TRUE
    )$values
    if (any(fixed_eigenvalues <= 0)) {
      stop(
        "fixed working correlation matrix is not positive definite",
        call. = FALSE
      )
    }
    alpha_fixed <- 1
  }
  beta_zero <- compute_geer_start_values(
    model_matrix = model_matrix,
    y = y,
    family = family,
    weights = weights,
    offset = offset,
    method = method,
    link = link,
    beta_start = beta_start,
    control = control,
    control_glm = control_glm
  )
  norm_phi <- normalize_phi(phi_fixed, phi_value)
  phi_fixed <- norm_phi$phi_fixed
  phi_value <- norm_phi$phi_value
  use_p <- normalize_use_p(use_p)
  subtract_p <- if (use_p) ncol(model_matrix) else 0
  if (identical(family$family, "quasi") && !is.null(family$varfun)) {
    fam_name <- switch(
      family$varfun,
      constant = "gaussian",
      `mu(1-mu)` = "binomial",
      mu = "poisson",
      `mu^2` = "Gamma",
      `mu^3` = "inverse.gaussian",
      NULL
    )
    if (!is.null(fam_name)) {
      family <- do.call(fam_name, list(link = family$link))
      link <- family$link
    }
  }
  if (identical(family$family, "quasibinomial") && !is.null(family$link)) {
    family <- do.call("binomial", list(link = family$link))
    link <- family$link
  }
  if (identical(family$family, "quasipoisson") && !is.null(family$link)) {
    family <- do.call("poisson", list(link = family$link))
    link <- family$link
  }
  fit_pass <- function(beta, pass, previous) {
    iterations <- get_pass_iterations(pass, control)
    if (pass$carry_nuisance) {
      alpha_pass <- previous$alpha
      alpha_fixed_pass <- 1L
      corstr_pass <- corstr
      phi_pass <- previous$phi
      phi_fixed_pass <- 1L
    } else if (pass$independence) {
      alpha_pass <- 0
      alpha_fixed_pass <- 0
      corstr_pass <- "independence"
      phi_pass <- phi_value
      phi_fixed_pass <- phi_fixed
    } else if (identical(pass$stage, "second")) {
      alpha_pass <- if (identical(corstr, "fixed")) alpha_vector else 0
      alpha_fixed_pass <- if (identical(corstr, "fixed")) alpha_fixed else 0
      corstr_pass <- corstr
      phi_pass <- phi_value
      phi_fixed_pass <- phi_fixed
    } else {
      alpha_pass <- alpha_vector
      alpha_fixed_pass <- alpha_fixed
      corstr_pass <- corstr
      phi_pass <- phi_value
      phi_fixed_pass <- phi_fixed
    }
    run_geer_solver(
      fit_geesolver_cc(
        y, model_matrix, id, repeated, weights,
        link, family$family, as.numeric(beta), offset,
        iterations[[1L]], tolerance, iterations[[2L]], iterations[[3L]],
        control$jeffreys_power, pass$method, subtract_p,
        alpha_pass, alpha_fixed_pass, corstr_pass, Mv,
        phi_pass, phi_fixed_pass, as.integer(isTRUE(pass$hold_nuisance))
      ),
      "geewa"
    )
  }
  geesolver_fit <- run_geer_estimation_passes(
    fit_pass = fit_pass,
    method = method,
    beta_start = beta_zero,
    check_first = function(fit, method) {
      check_geer_first_pass(fit, method, tolerance, "geewa")
    }
  )
  fit <- build_geer_output(
    geesolver_fit = geesolver_fit,
    xnames = xnames,
    qr_model_matrix = qr_model_matrix,
    family = family,
    weights = weights,
    y = y,
    model_matrix = model_matrix,
    model_frame = model_frame,
    id = id,
    repeated = repeated,
    call = call,
    formula = formula,
    data = data,
    model_terms = model_terms,
    control = control,
    method = method,
    association_structure = corstr,
    row_order = inputs$row_order
  )
  fit <- build_geer_final_fit(
    fit = fit,
    geesolver_fit = geesolver_fit,
    tolerance = tolerance,
    method = method,
    family = family,
    association_structure = corstr,
    repeated = repeated,
    fit_function = "geewa"
  )
  fit$use_p <- use_p
  fit$phi_fixed <- isTRUE(phi_fixed)
  fit <- new_geer(fit)
  fit <- validate_geer(fit)
  fit
}
