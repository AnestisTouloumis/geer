#' @title
#' Solving (Adjusted) GEE for Binary Responses
#'
#' @inherit geewa description
#'
#' @inheritParams geewa
#' @param link character specifying the link function. Options include
#'        \code{logit}, \code{probit}, \code{cauchit}, \code{cloglog},
#'        \code{identity}, \code{log}, \code{sqrt}, \code{1/mu^2} and
#'        \code{inverse}. By default, \code{link = "logit"}.
#' @param orstr character specifying the odds ratios structure.
#'        Options include \code{"independence"}, \code{"exchangeable"},
#'        \code{"unstructured"} and \code{"fixed"}. By default,
#'        \code{orstr = "independence"}.
#' @param alpha_vector numerical vector specifying odds ratios structure
#'        when \code{orstr == "fixed"}. It is ignored for all other
#'        possible values of \code{orstr}.
#'
#' @inherit geewa details
#'
#' @inherit geewa return
#'
#' @author Anestis Touloumis
#'
#' @examples
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", ]
#' fitted_model <- geewa_binary(
#'   formula = status ~ baseline + treatment * gender + visit * age,
#'   id = id, repeated = visit, link = "probit",
#'   data = respiratory2, orstr = "independence",
#'   method = "pgee-jeffreys"
#' )
#' summary(fitted_model, cov_type = "bias-corrected")
#'
#' data("cholecystectomy", package = "geer")
#' fitted_model_gee <- geewa_binary(
#'   formula = pain ~ treatment + gender + age,
#'   id = id, data = cholecystectomy, link = "logit",
#'   orstr = "independence", method = "gee"
#' )
#' summary(fitted_model_gee, cov_type = "bias-corrected")
#' fitted_model_brgee_robust <- update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_robust, cov_type = "bias-corrected")
#' fitted_model_brgee_naive <- update(fitted_model_gee, method = "brgee-naive")
#' summary(fitted_model_brgee_naive, cov_type = "bias-corrected")
#' fitted_model_brgee_empirical <- update(fitted_model_gee, method = "brgee-empirical")
#' summary(fitted_model_brgee_empirical, cov_type = "bias-corrected")
#'
#' @export
geewa_binary <- function(formula,
                         link = "logit",
                         data = parent.frame(),
                         id,
                         repeated = NULL,
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
