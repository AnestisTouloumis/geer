#' @title
#' Solving (Adjusted) GEE
#'
#' @description
#' Produces a(n) (adjusted) Generalized Estimating Equations (GEE) fit of the
#' data. The estimation methods include the traditional GEE method,
#' bias-reduction methods, bias-correction methods and a penalized GEE method.
#'
#' @param formula \code{formula} expression of the form
#'        \code{response ~ predictors}: a symbolic description of the marginal
#'        model to be fitted.
#' @param family A \code{\link[stats]{family}} object indicating the link and
#'        variance functions. Supported families include \code{gaussian},
#'        \code{binomial}, \code{poisson}, \code{Gamma},
#'        \code{inverse.gaussian} and \code{quasi}. By default,
#'        \code{family = gaussian(link = "identity")}.
#' @param data optional data frame containing variables referenced in
#'        \code{formula}. The cluster identifiers \code{id} and \code{repeated}
#'        may be supplied either as columns in \code{data} or as separate vectors.
#' @param id a vector identifying the clusters.
#' @param repeated optional vector identifying the order of observations
#'        within each cluster.
#' @param control list specifying the control variables for the GEE solver
#'        and the power of the Jeffreys prior when
#'        \code{method = "pgee-jeffreys"}.
#' @param corstr character specifying the correlation structure. Options
#'        include \code{"independence"}, \code{"exchangeable"}, \code{"ar1"},
#'        \code{"m-dependent"}, \code{"unstructured"} and \code{"fixed"}.
#'        By default, \code{corstr = "independence"}.
#' @param Mv positive integer which must be specified whenever
#'        \code{corstr = "m-dependent"}. By default, \code{Mv = 1}.
#' @param method character specifying the (adjusted) GEE. Options include the
#'        traditional gee method (\code{"gee"}), bias-reduction methods
#'        (\code{"brgee-naive"}, \code{"brgee-robust"},
#'        \code{"brgee-empirical"}), bias-correction methods
#'        (\code{"bcgee-naive"}, \code{"bcgee-robust"},
#'        \code{"bcgee-empirical"}) and a penalized GEE method
#'        (\code{"pgee-jeffreys"}). By default, \code{method = "gee"}.
#' @param weights optional numeric vector of observation weights. Must be finite
#'        and strictly positive. If not supplied, all weights are 1.
#' @param beta_start numerical vector indicating the initial values of the
#'        regression parameters.
#' @param offset this can be used to specify an a priori known component to be
#'        included in the linear predictor during fitting. This should be NULL
#'        or a numeric vector of length equal to the number of cases.
#'        One or more offset terms can be included in the formula instead or
#'        as well, and if more than one is specified their sum is used.
#' @param control_glm list of parameters for controlling the fitting process for
#'        the initial values of the parameter vector when these are not provided.
#' @param use_p logical indicating whether to use the \code{N - p} correction
#'        for estimating the scale and the correlation parameters. By default,
#'        \code{use_p = TRUE}.
#' @param alpha_vector Numeric vector of fixed association parameters used only
#'        when \code{corstr = "fixed"}. Must have length \code{choose(T, 2)} where
#'        \code{T = max(repeated)} after recoding. Ignored otherwise.
#' @param phi_fixed logical indicating whether the scale parameter is fixed at
#'        the value of \code{phi_value}. By default, \code{phi_fixed = FALSE}.
#' @param phi_value positive number indicating the value to which the scale
#'        parameter must be fixed; only used when \code{phi_fixed == TRUE}.
#'        By default, \code{phi_value = 1}.
#' @param ... further arguments passed to/or from other methods.
#'
#' @details
#' \code{method} specifies the adjustment vector (if any) added to the ordinary
#' generalized estimating equations. If \code{method = "gee"}, then the original
#' GEE are solved, i.e. no adjustment vector is added. If
#' \code{method = "brgee-naive"}, \code{method = "brgee-robust"} or
#' \code{method = "brgee-empirical"}, then the GEE will be adjusted to
#' produce naive, robust or empirical bias-reducing estimators, respectively.
#' If \code{method = "bcgee-naive"}, \code{method = "bcgee-robust"} or
#' \code{method = "bcgee-empirical"}, then the GEE will be adjusted to
#' produce naive, robust or empirical bias-corrected estimators. If
#' \code{method = "pgee-jeffreys"}, then the GEE will be penalized using a
#' Jeffreys prior type penalty.
#'
#' For the construction of the \code{formula} argument, see the documentation of
#' \code{\link{glm}} and \code{\link{formula}}.
#'
#' The \code{data} must be in long format (one row per observation). See
#' \code{\link{reshape}} for details on reshaping between long and wide formats.
#'
#' The default set for the \code{id} labels is \eqn{\{1,\ldots,N\}}, where
#' \eqn{N} is the sample size, i.e. the total number of clusters. Otherwise,
#' the function recodes the given labels of \code{id} onto this set.
#'
#' The argument \code{repeated} can be safely omitted only if observations are
#' already ordered within each cluster as intended. If \code{repeated} is not
#' provided, it is created as \code{1, 2, ..., n_i} within each cluster \eqn{i},
#' using the current row order (before internal sorting). If \code{repeated} is
#' provided, its labels are recoded to \code{1, ..., T} and must be unique within
#' each cluster.
#'
#' The variables \code{id} and \code{repeated} do not need to be pre-sorted.
#' Instead the function sorts observations in ascending order of \code{id}
#' and \code{repeated}.
#'
#' A term of the form \code{offset(expression)} is allowed in the right hand
#' side of \code{formula}.
#'
#' The length of \code{id} and of \code{repeated} or \code{weights} (when
#' provided) should be the same as the number of observations.
#'
#' @return
#' Returns an object from the class \code{geer}, a list with components:
#' \item{coefficients}{a named vector of coefficients.}
#' \item{residuals}{the working residuals.}
#' \item{fitted.values}{the fitted mean values, obtained by transforming the
#'      linear predictors by the inverse of the link function.}
#' \item{rank}{the numeric rank of the fitted model.}
#' \item{family}{the \code{\link{family}} object used.}
#' \item{linear.predictors}{the linear fit on link scale.}
#' \item{iter}{the number of iterations used.}
#' \item{prior.weights}{the weights initially supplied, a vector of 1s if none
#'       were.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{y}{the response vector.}
#' \item{x}{the model matrix.}
#' \item{id}{the \code{id} vector.}
#' \item{repeated}{the \code{repeated} vector.}
#' \item{converged}{logical indicating whether the algorithm converged.}
#' \item{call}{the matched call.}
#' \item{formula}{the formula supplied.}
#' \item{terms}{the \code{\link{terms}} object used.}
#' \item{data}{the data argument.}
#' \item{offset}{the offset vector used.}
#' \item{control}{the value of the control argument used.}
#' \item{method}{the method that created the adjustment vector to the GEE, if
#'       any.}
#' \item{contrasts}{the contrasts used.}
#' \item{xlevels}{a record of the levels of the factors used in fitting.}
#' \item{naive_covariance}{the naive covariance matrix.}
#' \item{robust_covariance}{the robust covariance matrix.}
#' \item{bias_corrected_covariance}{the bias-corrected covariance matrix.}
#' \item{association_structure}{the name of the working assumption about the
#'       association structure.}
#' \item{alpha}{a vector of the association parameters.}
#' \item{phi}{the scale parameter.}
#' \item{obs_no}{the number of observations across clusters.}
#' \item{clusters_no}{the number of clusters.}
#' \item{min_cluster_size}{the minimum cluster size.}
#' \item{max_cluster_size}{the maximum cluster size.}
#'
#' @author Anestis Touloumis
#'
#' @seealso \code{\link{geewa_binary}}.
#'
#' @examples
#' data("epilepsy")
#' fitted_model_gee <- geewa(formula = seizures ~ treatment + lnbaseline + lnage,
#'                           data = epilepsy, id = id,
#'                           family = poisson(link = "log"),
#'                           corstr = "exchangeable", method = "gee")
#' summary(fitted_model_gee, cov_type = "bias-corrected")
#' fitted_model_brgee_robust <- update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_robust, cov_type = "bias-corrected")
#' fitted_model_brgee_naive <- update(fitted_model_gee, method = "brgee-naive")
#' summary(fitted_model_brgee_naive, cov_type = "bias-corrected")
#' fitted_model_brgee_empirical <- update(fitted_model_gee, method = "brgee-empirical")
#' summary(fitted_model_brgee_empirical, cov_type = "bias-corrected")
#'
#' @export
geewa <- function(formula,
                  family = gaussian(link = "identity"),
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
                  ...) {
  ## call and model frame
  call <- match.call(expand.dots = TRUE)
  mcall <- match.call(expand.dots = FALSE)
  model_frame <- build_geer_model_frame(mcall, env = parent.frame())
  ## family
  family <- normalize_family(family)
  link <- family$link
  ## response and weights
  response_weights <- extract_geer_response_weights(model_frame, family)
  y <- response_weights$y
  weights <- response_weights$weights
  ## id and repeated
  id_repeated <- extract_geer_id_repeated(model_frame, length(y))
  id <- id_repeated$id
  repeated <- id_repeated$repeated
  ## offset
  offset <- extract_geer_offset(model_frame, y_length = length(y))
  ## model matrix
  design <- build_geer_design_matrix(model_frame)
  model_terms <- design$terms
  model_matrix <- design$x
  xnames <- design$xnames
  qr_model_matrix <- design$qr
  x_assign <- design$assign
  x_contrasts <- design$contrasts
  ## sort by id then repeated
  ord <- order(id, repeated)
  y <- y[ord]
  model_matrix <- model_matrix[ord, , drop = FALSE]
  weights <- weights[ord]
  offset <- offset[ord]
  id <- id[ord]
  repeated <- repeated[ord]
  attr(model_matrix, "assign") <- x_assign
  attr(model_matrix, "contrasts") <- x_contrasts
  ## control (supports either a control object or a list of args)
  control <- normalize_geer_control(control)
  maxiter <- control$maxiter
  tolerance <- control$tolerance
  ## method
  method <- as.character(method)
  check_choice(method, valid_methods, "method")
  ## correlation structure
  check_choice(corstr, valid_corstrs, "corstr")
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
  if (method %in% c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    ## pass 1: plain GEE to convergence
    geesolver_fit <- fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, beta_zero, offset,
      maxiter, tolerance, control$step_maxiter,
      control$step_multiplier, control$jeffreys_power,
      "gee", subtract_p, alpha_vector, alpha_fixed,
      corstr, Mv, phi_value, phi_fixed
    )
    last_criterion <- geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1L]
    if (last_criterion > tolerance) {
      stop("bias-corrected estimator is undefined because the corresponding GEE model did not converge", call. = FALSE)
    }
    ## pass 2: one BR-GEE step warm-started from converged GEE solution
    geesolver_fit <- fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, as.numeric(geesolver_fit$beta_hat),
      offset, 1L, tolerance, 1L, 1L,
      control$jeffreys_power, sub("bcgee", "brgee", method),
      subtract_p, geesolver_fit$alpha, 1L,
      corstr, Mv, geesolver_fit$phi, 1L
    )
  } else if (method == "hpgee-jeffreys") {
    ## pass 1: PGEE under independence to convergence
    geesolver_fit <- fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, beta_zero, offset,
      maxiter, tolerance, control$step_maxiter,
      control$step_multiplier, control$jeffreys_power,
      "pgee-jeffreys", subtract_p, 0, 0,
      "independence", Mv, phi_value, phi_fixed
    )
    last_criterion <- geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1L]
    if (last_criterion > tolerance) {
      stop("hpgee-jeffreys estimator is undefined because the independence pgee-jeffreys model did not converge", call. = FALSE)
    }
    ## pass 2: one GEE step warm-started from penalised solution
    geesolver_fit <- fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, as.numeric(geesolver_fit$beta_hat),
      offset, 1L, tolerance, 1L, 1L,
      control$jeffreys_power, "gee",
      subtract_p, 0, 0,
      corstr, Mv, phi_value, phi_fixed
    )
  } else if (method == "opgee-jeffreys") {
    ## pass 1: PGEE under independence to convergence
    geesolver_fit <- fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, beta_zero, offset,
      maxiter, tolerance, control$step_maxiter,
      control$step_multiplier, control$jeffreys_power,
      "pgee-jeffreys", subtract_p, 0, 0,
      "independence", Mv, phi_value, phi_fixed
    )
    last_criterion <- geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1L]
    if (last_criterion > tolerance) {
      stop("opgee-jeffreys estimator is undefined because the independence pgee-jeffreys model did not converge", call. = FALSE)
    }
    ## pass 2: one PGEE step warm-started from first-pass solution
    geesolver_fit <- fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, as.numeric(geesolver_fit$beta_hat),
      offset, 1L, tolerance, 1L, 1L,
      control$jeffreys_power, "pgee-jeffreys",
      subtract_p, 0, 0,
      corstr, Mv, phi_value, phi_fixed
    )
  } else {
    ## single pass: gee, brgee-*, or pgee-jeffreys
    geesolver_fit <- fit_geesolver_cc(
      y, model_matrix, id, repeated, weights,
      link, family$family, beta_zero, offset,
      maxiter, tolerance, control$step_maxiter,
      control$step_multiplier, control$jeffreys_power,
      method, subtract_p, alpha_vector, alpha_fixed,
      corstr, Mv, phi_value, phi_fixed
    )
  }
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
    data = data,
    model_terms = model_terms,
    control = control,
    method = method,
    association_structure = corstr
  )
  fit$converged <- (geesolver_fit$criterion[fit$iter] <= tolerance)
  if (method %in% c("bcgee-naive", "bcgee-robust", "bcgee-empirical", "opgee-jeffreys", "hpgee-jeffreys")) {
    fit$converged <- TRUE
  }

  fit$alpha <- if (identical(corstr, "independence")) 0 else as.numeric(geesolver_fit$alpha)

  if (corstr %in% c("unstructured", "fixed")) {
    pairs_matrix <- combn(max(repeated), 2)
    alpha_names <- paste0("alpha_", pairs_matrix[1, ], ".", pairs_matrix[2, ])
    names(fit$alpha) <- alpha_names
  }
  if (!fit$converged) warning("geewa: algorithm did not converge", call. = FALSE)
  eps <- 10 * .Machine$double.eps
  if (identical(family$family, "binomial")) {
    if (any(fit$fitted.values > 1 - eps) || any(fit$fitted.values < eps)) {
      warning("geewa: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
    }
  }
  if (identical(family$family, "poisson")) {
    if (any(fit$fitted.values < eps)) {
      warning("geewa: fitted rates numerically 0 occurred", call. = FALSE)
    }
  }
  if (!identical(corstr, "independence")) {
    corr <- get_correlation_matrix(corstr, fit$alpha, max(repeated))
    eigenvalues <- eigen(corr, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigenvalues <= 0)) {
      warning("geewa: working correlation matrix is not positive definite", call. = FALSE)
    }
  }
  fit <- new_geer(fit)
  fit
}
