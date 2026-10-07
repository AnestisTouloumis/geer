#' @title
#' Fitting (Adjusted) Generalized Estimating Equations
#'
#' @description
#' Fits a marginal model for repeated or clustered responses using
#' generalized estimating equations (GEE). Supported estimation methods include
#' the traditional GEE, bias-reduced GEE, bias-corrected GEE, and
#' Jeffreys-type penalized GEE.
#'
#' @param formula \code{formula} expression of the form
#'        \code{response ~ predictors}: a symbolic description of the marginal
#'        model to be fitted.
#' @param family a \code{\link[stats]{family}} object, family function, or
#'        character string naming a family function, specifying the marginal
#'        variance and link functions. Supported families are \code{gaussian},
#'        \code{binomial}, \code{poisson}, \code{Gamma},
#'        \code{inverse.gaussian}, \code{quasi}, \code{quasibinomial} and
#'        \code{quasipoisson}. Defaults to
#'        \code{gaussian(link = "identity")}.
#' @param data optional data frame containing variables referenced in
#'        \code{formula}, \code{id}, \code{repeated}, \code{weights}, and
#'        \code{offset}. Defaults to \code{parent.frame()}.
#' @param id variable identifying the clusters.
#' @param repeated optional variable identifying the order of observations
#'        within each cluster.
#' @param control a \code{\link{geer_control}} list specifying convergence
#'        tolerance, maximum iterations, step-halving parameters, and the
#'        power of the Jeffreys-type penalty. Defaults to \code{geer_control()}.
#' @param corstr character string specifying the working correlation structure. Options
#'        are \code{"independence"}, \code{"exchangeable"}, \code{"ar1"},
#'        \code{"m-dependent"}, \code{"unstructured"}, \code{"toeplitz"} and
#'        \code{"fixed"}. Defaults to \code{"independence"}.
#' @param Mv positive integer giving the number of lags for the m-dependent
#'        working correlation structure. Defaults to \code{1} (lag-1
#'        dependence). Must be set explicitly when lags other than 1 are
#'        intended. Ignored when \code{corstr != "m-dependent"}.
#' @param method character string specifying the estimation method. Options are
#'        the traditional GEE (\code{"gee"}), bias-reduced methods
#'        (\code{"brgee-robust"}, \code{"brgee-empirical"}, \code{"brgee-naive"}),
#'        bias-corrected methods (\code{"bcgee-robust"}, \code{"bcgee-empirical"},
#'        \code{"bcgee-naive"}), the fully iterated Jeffreys-type penalized GEE
#'        (\code{"pgee-jeffreys"}), the one-step penalized GEE
#'        (\code{"opgee-jeffreys"}), and the hybrid one-step GEE
#'        (\code{"hpgee-jeffreys"}). Defaults to \code{"gee"}.
#' @param weights optional numeric vector of observation weights. Must be finite
#'        and strictly positive. If not supplied, all weights are 1.
#' @param beta_start optional numeric vector of starting values for the
#'        regression parameters. Defaults to \code{NULL}, in which case
#'        starting values are obtained from an auxiliary generalized linear
#'        model fit, using
#'        \code{\link[brglm2]{brglmFit}} where appropriate.
#' @param offset this can be used to specify an a priori known component to be
#'        included in the linear predictor during fitting. This should be
#'        \code{NULL}, a single numeric value, or a numeric vector of length
#'        equal to the number of observations. One or more offset terms can be
#'        included in the formula instead or as well, and if more than one is
#'        specified their sum is used.
#' @param control_glm optional list of control parameters interpreted by
#'        \code{\link[brglm2]{brglm_control}} when computing GLM-based
#'        starting values. Ignored when \code{beta_start} is supplied. Defaults
#'        to \code{list(...)}, that is, the arguments supplied through
#'        \code{...}.
#' @param use_p logical indicating whether to apply a degrees-of-freedom
#'        correction by subtracting the number of regression parameters from
#'        the relevant denominator when estimating the scale and working
#'        correlation parameters. Defaults to \code{TRUE}.
#' @param alpha_vector numeric vector of fixed working correlation parameters
#'        used only when \code{corstr = "fixed"}. Must have length \code{choose(T, 2)}
#'        where \code{T = max(repeated)} after recoding, and the resulting
#'        working correlation matrix must be positive definite. Ignored
#'        otherwise. Defaults to \code{NULL}.
#' @param phi_fixed logical indicating whether the scale parameter is fixed at
#'        the value of \code{phi_value}. Defaults to \code{FALSE}.
#' @param phi_value positive number giving the fixed value of the scale
#'        parameter. Used only when \code{phi_fixed = TRUE}. Defaults to
#'        \code{1}.
#' @param subset an optional vector specifying a subset of observations to be
#'        used in the fitting process. It is evaluated in \code{data}, as in
#'        \code{\link[stats]{model.frame}}.
#' @param na.action a function indicating what should happen when the data
#'        contain \code{NA}s, as in \code{\link[stats]{model.frame}}. Defaults
#'        to \code{getOption("na.action")}.
#' @param ... additional arguments passed to or from other methods.
#'
#' @details
#' \code{method} specifies the estimation approach. If \code{method = "gee"},
#' the standard GEE are solved with no adjustment. If \code{method} is one of
#' \code{"brgee-naive"}, \code{"brgee-robust"} or \code{"brgee-empirical"},
#' an adjustment vector is added to produce naive, robust or empirical
#' bias-reduced estimators, respectively. If \code{method} is one of
#' \code{"bcgee-naive"}, \code{"bcgee-robust"} or \code{"bcgee-empirical"},
#' the corresponding bias-corrected estimators are produced via a one-step
#' correction applied to the converged GEE solution. If
#' \code{method = "pgee-jeffreys"}, the GEE are penalized using a
#' Jeffreys-type penalty run to full convergence. If
#' \code{method = "opgee-jeffreys"}, a single penalized scoring step is
#' performed from the converged independence penalized solution (one-step
#' approximation). If \code{method = "hpgee-jeffreys"}, a single standard GEE
#' scoring step is performed from the converged independence penalized solution
#' (hybrid one-step approximation).
#'
#' For the construction of the \code{formula} argument, see the documentation
#' of \code{\link[stats]{glm}} and \code{\link[stats]{formula}}.
#'
#' The \code{data} must be in long format (one row per observation). See
#' \code{\link[stats]{reshape}} for details on reshaping between long and wide formats.
#'
#' The \code{quasi}, \code{quasibinomial} and \code{quasipoisson} families are
#' internally remapped to their standard parametric equivalents before fitting.
#' The \code{quasibinomial} and \code{quasipoisson} families are mapped to
#' \code{binomial} and \code{poisson}, respectively. A \code{quasi} family is
#' supported when its variance function is \code{constant}, \code{mu(1-mu)},
#' \code{mu}, \code{mu^2}, or \code{mu^3}; these are mapped to
#' \code{gaussian}, \code{binomial}, \code{poisson}, \code{Gamma}, and
#' \code{inverse.gaussian}, respectively. The scale parameter \code{phi} is
#' then estimated freely from the data unless \code{phi_fixed = TRUE}.
#'
#' The default set for the \code{id} labels is \eqn{\{1,\ldots,N\}}, where
#' \eqn{N} is the number of clusters. Otherwise, the function recodes the given
#' labels of \code{id} onto this set.
#'
#' The argument \code{repeated} can be safely omitted only if observations are
#' already ordered within each cluster as intended. If \code{repeated} is not
#' provided, it is created as \code{1, 2, ..., n_i} within each cluster
#' \eqn{i}, using the current row order (before internal sorting). If
#' \code{repeated} is provided, its labels are recoded to \code{1, ..., T} and
#' must be unique within each cluster.
#'
#' The variables \code{id} and \code{repeated} do not need to be pre-sorted.
#' Instead the function sorts observations in ascending order of \code{id}
#' and \code{repeated}.
#'
#' A term of the form \code{offset(expression)} is allowed in the right-hand
#' side of \code{formula}.
#'
#' The length of \code{id} and, when provided, of \code{repeated} and
#' \code{weights} must equal the number of observations.
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
#' \item{use_p}{for \code{geewa()} fits, logical indicating whether the
#'       degrees-of-freedom correction was used. This component is not stored
#'       for \code{geewa_binary()} fits.}
#' \item{phi_fixed}{for \code{geewa()} fits, logical indicating whether the
#'       scale parameter was fixed. This component is not stored for
#'       \code{geewa_binary()} fits.}
#' \item{contrasts}{the contrasts used.}
#' \item{xlevels}{a record of the levels of the factors used in fitting.}
#' \item{na.action}{information on how missing values were handled, as returned
#'       by the \code{na.action} attribute of the model frame. \code{NULL} if
#'       no observations were removed.}
#' \item{naive_covariance}{the model-based (naive) covariance matrix.}
#' \item{robust_covariance}{the sandwich (robust) covariance matrix.}
#' \item{bias_corrected_covariance}{the bias-corrected covariance matrix.}
#' \item{association_structure}{the name of the working association structure.}
#' \item{alpha}{a vector of the estimated working association parameters.}
#' \item{phi}{the scale parameter. For \code{geewa()} fits, this is
#'       estimated or fixed according to the fitting options; for
#'       \code{geewa_binary()} fits, it is always \code{1}.}
#' \item{obs_no}{the number of observations used in fitting.}
#' \item{clusters_no}{the number of clusters.}
#' \item{min_cluster_size}{the minimum cluster size.}
#' \item{max_cluster_size}{the maximum cluster size.}
#'
#' @section Note on returned components:
#' For \code{geewa}, the \code{alpha} contains the estimated (or fixed) working
#' correlation parameters. The \code{association_structure}
#' component stores the value of \code{corstr}. Under
#' \code{corstr = "independence"}, \code{alpha} is set to \code{0}.
#'
#' For \code{method} in \code{"bcgee-naive"}, \code{"bcgee-robust"},
#' \code{"bcgee-empirical"}, \code{"opgee-jeffreys"}, and
#' \code{"hpgee-jeffreys"}, \code{converged} is \code{TRUE} in the
#' returned object, because these methods produce their estimate via a single
#' correction step applied to an already-converged fit. The exception is a
#' numerical failure in that correction step: a warning is then issued, the
#' preliminary estimates of the first stage are returned and \code{converged}
#' is \code{FALSE}.
#'
#' For the bias-corrected methods, \code{alpha} and \code{phi} are those of the
#' converged GEE solution and are held fixed in the correction step. For
#' \code{"opgee-jeffreys"} and \code{"hpgee-jeffreys"} they are estimated
#' once, under the requested working correlation structure, at the converged
#' independence penalized estimate of the regression parameters, and are held
#' fixed in the one-step update; the returned \code{alpha}, \code{phi} and
#' covariance matrices use these values.
#'
#' @references
#' Touloumis, A. (2026a) Bias-reduced GEE via adjusted estimating equations,
#' with odds-ratio extensions. \emph{Preprint}.
#' \url{https://arxiv.org/abs/2606.16043}
#'
#' Touloumis, A. (2026b) Jeffreys-type penalized GEE for correlated binary data
#' with an odds-ratio parameterization. \emph{Preprint}.
#' \url{https://arxiv.org/abs/2606.16058}
#'
#' @seealso
#' \code{\link{geewa_binary}}, \code{\link{geer_control}}, \code{\link{summary.geer}},
#' \code{\link{vcov.geer}}, \code{\link{anova.geer}}, \code{\link{step_p}},
#' \code{\link{geecriteria}}, \code{\link{runs_test}}.
#'
#' @examples
#' data("epilepsy", package = "geer")
#' fit_gee <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   corstr = "exchangeable",
#'   method = "gee"
#' )
#' summary(fit_gee, cov_type = "bias-corrected")
#'
#' fit_brgee_robust <- update(fit_gee, method = "brgee-robust")
#' summary(fit_brgee_robust, cov_type = "bias-corrected")
#'
#' fit_brgee_naive <- update(fit_gee, method = "brgee-naive")
#' summary(fit_brgee_naive, cov_type = "bias-corrected")
#'
#' fit_brgee_empirical <- update(fit_gee, method = "brgee-empirical")
#' summary(fit_brgee_empirical, cov_type = "bias-corrected")
#'
#' fit_bcgee_robust <- update(fit_gee, method = "bcgee-robust")
#' summary(fit_bcgee_robust, cov_type = "robust")
#'
#' \donttest{
#' ## Penalized GEE with custom control
#' fit_pgee <- geewa(
#'   formula = seizures ~ treatment + lnbaseline + lnage,
#'   family = poisson(link = "log"),
#'   data = epilepsy,
#'   id = id,
#'   control = geer_control(jeffreys_power = 0.1),
#'   corstr = "exchangeable",
#'   method = "pgee-jeffreys"
#' )
#' summary(fit_pgee, cov_type = "robust")
#' }
#'
#' @export
geewa <- function(formula,
                  family = stats::gaussian(link = "identity"),
                  data = parent.frame(),
                  id,
                  repeated,
                  control = geer_control(...),
                  corstr = "independence",
                  Mv = 1,
                  method = "gee",
                  weights,
                  beta_start = NULL,
                  offset,
                  control_glm = list(...),
                  use_p = TRUE,
                  alpha_vector = NULL,
                  phi_fixed = FALSE,
                  phi_value = 1,
                  subset,
                  na.action,
                  ...) {
  ## when both 'control' and 'control_glm' are supplied, nothing consumes '...'
  if (!missing(control) && !missing(control_glm)) {
    check_unused_dots(list(...), "geewa")
  }
  ## call, family and common input preparation
  call <- match.call(expand.dots = TRUE)
  mcall <- match.call(expand.dots = FALSE)
  caller_env <- parent.frame()
  family <- normalize_family(family)
  link <- family$link
  inputs <- prepare_geer_inputs(
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
  maxiter <- control$maxiter
  tolerance <- control$tolerance
  ## correlation structure
  check_choice(corstr, geer_corstr_choices, "corstr")
  if (!identical(corstr, "m-dependent")) {
    Mv <- 1L
  } else {
    if (!is_positive_integer_scalar(Mv)) {
      stop("'Mv' must be a positive integer", call. = FALSE)
    }
    Mv <- as.integer(Mv)
  }
  if (!identical(corstr, "fixed")) {
    alpha_vector <- 0
    alpha_fixed <- 0
  } else {
    if (is.null(alpha_vector)) stop("'alpha_vector' must be provided when 'corstr = \"fixed\"'", call. = FALSE)
    alpha_vector <- as.numeric(alpha_vector)
    repeated_max <- max(repeated)
    if (length(alpha_vector) != choose(repeated_max, 2)) {
      stop("'alpha_vector' must be a numeric vector of length ", choose(repeated_max, 2), call. = FALSE)
    }
    if (any(eigen(get_correlation_matrix(corstr, alpha_vector, repeated_max),
                  symmetric = TRUE, only.values = TRUE)$values <= 0)) {
      stop("fixed working correlation matrix is not positive definite", call. = FALSE)
    }
    alpha_fixed <- 1
  }
  ## initial beta
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
  ## phi
  norm_phi <- normalize_phi(phi_fixed, phi_value)
  phi_fixed <- norm_phi$phi_fixed
  phi_value <- norm_phi$phi_value
  ## N - p correction
  use_p <- normalize_use_p(use_p)
  subtract_p <- if (use_p) ncol(model_matrix) else 0
  ## quasi mapping (guarded)
  if (identical(family$family, "quasi") && !is.null(family$varfun)) {
    fam_name <- switch(family$varfun,
                       constant = "gaussian",
                       `mu(1-mu)` = "binomial",
                       mu = "poisson",
                       `mu^2` = "Gamma",
                       `mu^3` = "inverse.gaussian",
                       NULL)
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
  ## fit
  fit_pass <- function(beta, pass, previous) {
    iterations <- if (pass$one_step) {
      list(1L, 1L, 1)
    } else {
      list(maxiter, control$step_maxiter, control$step_multiplier)
    }
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
      ## A fixed working correlation is never estimated, so the second pass
      ## must use the supplied association parameters.
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
    fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, as.numeric(beta), offset,
      iterations[[1L]], tolerance, iterations[[2L]], iterations[[3L]],
      control$jeffreys_power, pass$method, subtract_p,
      alpha_pass, alpha_fixed_pass, corstr_pass, Mv,
      phi_pass, phi_fixed_pass, as.integer(isTRUE(pass$hold_nuisance))
    )
  }
  geesolver_fit <- run_geer_estimation_passes(
    fit_pass = fit_pass,
    method = method,
    beta_start = beta_zero,
    check_first = function(fit, method) check_geer_first_pass(fit, method, tolerance)
  )
  ## output
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
  fit <- finalize_geer_fit(
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
  if (!identical(corstr, "independence")) {
    corr <- get_correlation_matrix(corstr, fit$alpha, max(repeated))
    eigenvalues <- eigen(corr, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigenvalues <= 0)) {
      warning("geewa: working correlation matrix is not positive definite", call. = FALSE)
    }
  }
  fit <- new_geer(fit)
  fit <- validate_geer(fit)
  fit
}
