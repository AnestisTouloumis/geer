#' @title
#' Fitting (Adjusted) Generalized Estimating Equations for Binary Responses
#'
#' @inherit geewa description
#'
#' @inheritParams geewa
#' @param link character string specifying the link function for the marginal mean
#'        model. Options are \code{"logit"}, \code{"probit"},
#'        \code{"cauchit"}, \code{"cloglog"}, \code{"identity"},
#'        \code{"log"}, \code{"sqrt"}, \code{"1/mu^2"} and
#'        \code{"inverse"}. Defaults to \code{"logit"}.
#' @param orstr character string specifying the working odds-ratio structure for the
#'        within-cluster association. Options are
#'        \code{"independence"}, \code{"exchangeable"},
#'        \code{"unstructured"} and \code{"fixed"}. Defaults to
#'        \code{"independence"}.
#' @param alpha_vector numeric vector of fixed odds-ratio parameters used only
#'        when \code{orstr = "fixed"}. Must have length \code{choose(T, 2)}
#'        where \code{T = max(repeated)} after recoding, and all elements must
#'        be finite and strictly positive. Ignored for all other values of
#'        \code{orstr}.
#'
#' @details
#' \code{method} specifies the estimation approach. If \code{method = "gee"},
#' the standard GEE are solved with no adjustment. If \code{method} is one of
#' \code{"brgee-naive"}, \code{"brgee-robust"} or \code{"brgee-empirical"},
#' an adjustment vector is added to produce naive, robust or empirical
#' bias-reducing estimators, respectively. If \code{method} is one of
#' \code{"bcgee-naive"}, \code{"bcgee-robust"} or \code{"bcgee-empirical"},
#' the corresponding bias-corrected estimators are produced via a one-step
#' correction applied to the converged GEE solution. If
#' \code{method = "pgee-jeffreys"}, the GEE are penalized using a
#' Jeffreys-prior penalty run to full convergence. If
#' \code{method = "opgee-jeffreys"}, a single penalized scoring step is
#' performed from the converged independence penalized solution (one-step
#' approximation). If \code{method = "hpgee-jeffreys"}, a single standard GEE
#' scoring step is performed from the converged independence penalized solution
#' (hybrid one-step approximation).
#'
#' The marginal mean model always uses a \code{binomial} family with the
#' specified \code{link}. Within-cluster association is modeled through
#' marginalized pairwise odds ratios rather than a correlation structure. The
#' scale parameter is fixed at \code{phi = 1}.
#'
#' For the construction of the \code{formula} argument, see the documentation
#' of \code{\link{glm}} and \code{\link{formula}}.
#'
#' The \code{data} must be in long format (one row per observation). See
#' \code{\link{reshape}} for details on reshaping between long and wide formats.
#'
#' The default set for the \code{id} labels is \eqn{\{1,\ldots,N\}}, where
#' \eqn{N} is the number of clusters. Otherwise, the function recodes the
#' given labels of \code{id} onto this set.
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
#' The length of \code{id} and of \code{repeated} or \code{weights} (when
#' provided) must equal the number of observations.
#'
#' @inherit geewa return
#'
#' @section Note on returned components:
#' For \code{geewa_binary}, the \code{phi} component is always \code{1} and
#' \code{alpha} contains the estimated (or fixed) working odds-ratio
#' parameters, not correlation parameters. The \code{association_structure}
#' component stores the value of \code{orstr}. Under
#' \code{orstr = "independence"}, \code{alpha} is set to \code{1}.
#'
#' For \code{method} in \code{"bcgee-naive"}, \code{"bcgee-robust"},
#' \code{"bcgee-empirical"}, \code{"opgee-jeffreys"}, and
#' \code{"hpgee-jeffreys"}, \code{converged} is always \code{TRUE} in the
#' returned object.
#'
#' @author Anestis Touloumis
#'
#' @seealso \code{\link{geewa}}, \code{\link{geer_control}},
#'   \code{\link{summary.geer}}, \code{\link{geecriteria}}.
#'
#' @examples
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]
#' fitted_model <- geewa_binary(
#'   formula = status ~ baseline + treatment * gender + visit * age,
#'   link = "probit",
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   orstr = "independence",
#'   method = "pgee-jeffreys"
#' )
#' summary(fitted_model, cov_type = "bias-corrected")
#'
#' data("cholecystectomy", package = "geer")
#' fitted_model_gee <- geewa_binary(
#'   formula = pain ~ treatment + gender + age,
#'   link = "logit",
#'   data = cholecystectomy,
#'   id = id,
#'   orstr = "independence",
#'   method = "gee"
#' )
#' summary(fitted_model_gee, cov_type = "bias-corrected")
#'
#' fitted_model_brgee_robust <- update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_robust, cov_type = "bias-corrected")
#'
#' fitted_model_brgee_naive <- update(fitted_model_gee, method = "brgee-naive")
#' summary(fitted_model_brgee_naive, cov_type = "bias-corrected")
#'
#' fitted_model_brgee_empirical <- update(fitted_model_gee, method = "brgee-empirical")
#' summary(fitted_model_brgee_empirical, cov_type = "bias-corrected")
#'
#' fitted_model_bcgee_robust <- update(fitted_model_gee, method = "bcgee-robust")
#' summary(fitted_model_bcgee_robust, cov_type = "robust")
#'
#' ## Penalized GEE with exchangeable odds-ratio structure
#' fitted_model_pgee <- geewa_binary(
#'   formula = pain ~ treatment + gender + age,
#'   link = "logit",
#'   data = cholecystectomy,
#'   id = id,
#'   orstr = "exchangeable",
#'   method = "pgee-jeffreys",
#'   control = geer_control(jeffreys_power = 1)
#' )
#' summary(fitted_model_pgee, cov_type = "robust")
#'
#' @export
geewa_binary <- function(formula,
                         link = "logit",
                         data = parent.frame(),
                         id,
                         repeated,
                         control = geer_control(...),
                         orstr = "independence",
                         method = "gee",
                         weights,
                         beta_start = NULL,
                         offset,
                         control_glm = list(...),
                         alpha_vector = NULL,
                         ...) {
  ## call and model frame
  call <- match.call(expand.dots = TRUE)
  mcall <- match.call(expand.dots = FALSE)
  model_frame <- build_geer_model_frame(mcall, env = parent.frame())
  ## link function
  family <- binomial(link = link)
  family <- normalize_family(family)
  link <- family$link
  ## response
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
  ## control
  control <- normalize_geer_control(control)
  maxiter <- control$maxiter
  tolerance <- control$tolerance
  ## method
  method <- as.character(method)
  check_choice(method, valid_methods, "method")
  ## initial beta
  beta_zero <- compute_geer_binary_start_values(
    model_matrix = model_matrix,
    y = y,
    family = family,
    weights = weights,
    offset = offset,
    method = method,
    beta_start = beta_start,
    control_glm = control_glm,
    tolerance = tolerance,
    jeffreys_power = control$jeffreys_power
  )
  ## odds ratios structure
  check_choice(orstr, valid_orstrs, "orstr")
  if (identical(orstr, "independence")) {
    alpha_vector <- rep.int(1, choose(max(repeated), 2))
  } else if (identical(orstr, "fixed")) {
    pairs_no <- choose(max(repeated), 2)
    if (!is.numeric(alpha_vector)) stop("'alpha_vector' must be a numeric vector", call. = FALSE)
    if (length(alpha_vector) != pairs_no) stop("'alpha_vector' must be a numeric vector of length ", pairs_no, call. = FALSE)
    if (any(!is.finite(alpha_vector)) || any(alpha_vector <= 0)) {
      stop("'alpha_vector' must be finite and strictly positive when 'orstr = \"fixed\"'", call. = FALSE)
    }
  } else {
    alpha_vector <- get_marginalized_odds_ratios(
      round(y), id, repeated, weights, control$or_adding, orstr
    )
  }
  ## fit
  alpha_independence <- rep.int(1, choose(max(repeated), 2))
  if (method %in% c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    ## pass 1: plain GEE to convergence
    geesolver_fit <- fit_bingee_or(
      y, model_matrix, id, repeated, weights, link,
      beta_zero, offset, maxiter, tolerance,
      control$step_maxiter, control$step_multiplier,
      control$jeffreys_power, "gee", alpha_vector
    )
    last_criterion <- geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1L]
    if (last_criterion > tolerance) {
      stop("bias-corrected estimator is undefined because the corresponding GEE model did not converge", call. = FALSE)
    }
    ## pass 2: one BR-GEE step warm-started from converged GEE solution
    geesolver_fit <- fit_bingee_or(
      y, model_matrix, id, repeated, weights, link,
      as.numeric(geesolver_fit$beta_hat), offset,
      1L, tolerance, 1L, 1L,
      control$jeffreys_power, sub("bcgee", "brgee", method), alpha_vector
    )
  } else if (method == "hpgee-jeffreys") {
    ## pass 1: PGEE under independence to convergence
    geesolver_fit <- fit_bingee_or(
      y, model_matrix, id, repeated, weights, link,
      beta_zero, offset, maxiter, tolerance,
      control$step_maxiter, control$step_multiplier,
      control$jeffreys_power, "pgee-jeffreys", alpha_independence
    )
    last_criterion <- geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1L]
    if (last_criterion > tolerance) {
      stop("hpgee-jeffreys estimator is undefined because the independence pgee-jeffreys model did not converge", call. = FALSE)
    }
    ## pass 2: one step GEE from penalised solution
    geesolver_fit <- fit_bingee_or(
      y, model_matrix, id, repeated, weights, link,
      as.numeric(geesolver_fit$beta_hat), offset,
      1L, tolerance, 1L, 1L,
      control$jeffreys_power, "gee", alpha_vector
    )
  } else if (method == "opgee-jeffreys") {
    ## pass 1: PGEE under independence to convergence
    geesolver_fit <- fit_bingee_or(
      y, model_matrix, id, repeated, weights, link,
      beta_zero, offset, maxiter, tolerance,
      control$step_maxiter, control$step_multiplier,
      control$jeffreys_power, "pgee-jeffreys", alpha_independence
    )
    last_criterion <- geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1L]
    if (last_criterion > tolerance) {
      stop("opgee-jeffreys estimator is undefined because the independence pgee-jeffreys model did not converge", call. = FALSE)
    }
    ## pass 2: one step PGEE from penalised solution
    geesolver_fit <- fit_bingee_or(
      y, model_matrix, id, repeated, weights, link,
      as.numeric(geesolver_fit$beta_hat), offset,
      1L, tolerance, 1L, 1L,
      control$jeffreys_power, "pgee-jeffreys", alpha_vector
    )
  } else {
    ## single pass: gee, brgee-*, or pgee-jeffreys
    geesolver_fit <- fit_bingee_or(
      y, model_matrix, id, repeated, weights, link,
      beta_zero, offset, maxiter, tolerance,
      control$step_maxiter, control$step_multiplier,
      control$jeffreys_power, method, alpha_vector
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
    association_structure = orstr
  )
  fit$converged <- (geesolver_fit$criterion[fit$iter] <= tolerance)
  if (method %in% c("bcgee-naive", "bcgee-robust", "bcgee-empirical", "opgee-jeffreys", "hpgee-jeffreys")) {
    fit$converged <- TRUE
  }

  fit$alpha <- if (identical(orstr, "independence")) 1 else as.numeric(geesolver_fit$alpha)

  if (orstr %in% c("unstructured", "fixed")) {
    pairs_matrix <- combn(max(repeated), 2)
    alpha_names <- paste0("alpha_", pairs_matrix[1, ], ".", pairs_matrix[2, ])
    names(fit$alpha) <- alpha_names
  }
  if (!fit$converged) {
    warning("geewa_binary: algorithm did not converge", call. = FALSE)
  }
  eps <- 10 * .Machine$double.eps
  if (any(fit$fitted.values > 1 - eps) || any(fit$fitted.values < eps)) {
    warning("geewa_binary: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
  }
  fit <- new_geer(fit)
  fit
}
