#' @title
#' Fitting (Adjusted) Generalized Estimating Equations for Binary Responses
#'
#' @description
#' Fits a marginal model for clustered or repeated binary responses by
#' generalized estimating equations (GEE), with the association between
#' responses modeled by marginalized pairwise odds ratios, using ordinary,
#' bias-reduced, bias-corrected, or Jeffreys-type penalized GEE.
#'
#' @inheritParams geewa
#' @param link character string specifying the link function: \code{"logit"},
#'   \code{"probit"}, \code{"cauchit"}, \code{"cloglog"}, \code{"identity"},
#'   \code{"log"}, \code{"sqrt"}, \code{"1/mu^2"}, or \code{"inverse"}.
#'   Defaults to \code{"logit"}.
#' @param orstr character string specifying the working odds-ratio structure:
#'   \code{"independence"}, \code{"exchangeable"}, \code{"unstructured"}, or
#'   \code{"fixed"}. Defaults to \code{"independence"}.
#' @param alpha_vector numeric vector of strictly positive, finite odds ratios
#'   for \code{orstr = "fixed"}, of length \code{choose(T, 2)}, where \eqn{T} is
#'   the number of distinct supplied \code{repeated} values, or the largest
#'   retained original position if \code{repeated} is omitted. Ignored for
#'   other structures. Defaults to \code{NULL}.
#'
#' @details
#' Link, method, and association names must be exact character strings.
#' \code{alpha_vector} is ignored without validation unless
#' \code{orstr = "fixed"}.
#' \code{control_glm} is ignored without validation when \code{beta_start} is
#' supplied.
#'
#' The estimation methods are described in \code{\link{geewa}}. The marginal
#' mean model uses the binomial family with the given \code{link}, the working
#' odds ratios are estimated from marginalized \eqn{2 \times 2}{2 x 2} tables
#' (see
#' \code{or_adding} in \code{\link{geer_control}}), and the dispersion
#' parameter is fixed at 1. Observations are sorted by \code{id} and
#' \code{repeated}. The cluster identifiers and supplied \code{repeated}
#' values are recoded to consecutive integers. Implicit positions retain gaps
#' after subsetting or missing-value removal.
#'
#' The Plackett covariance and its derivatives are evaluated using scaled
#' expressions to reduce cancellation and overflow. Accuracy is not
#' guaranteed for extreme combinations of marginal probabilities and odds
#' ratios. In particular, near-complementary probabilities with very small
#' odds ratios can produce a covariance outside the bounds of a valid
#' bivariate Bernoulli distribution. The intended working model is that
#' of Touloumis (2026a).
#'
#' @section Convergence and separation:
#' Here \code{phi} is always 1 and \code{alpha} contains the working odds
#' ratios, which are 1 under independence and are computed once from the data,
#' so they are the same in both stages of the two-stage methods.
#'
#' For the two-stage methods (\code{"bcgee-*"}, \code{"opgee-jeffreys"}, and
#' \code{"hpgee-jeffreys"}), \code{converged} is \code{TRUE} unless the
#' correction step fails, in which case the first-stage estimates are returned
#' with a warning. If the first stage does not converge, the two-stage methods
#' give an error. If \code{converged} is \code{FALSE}, the results are those
#' of the last accepted iteration and should not be interpreted.
#' \code{converged = TRUE} only means that the stopping criterion was met: with
#' separated data the estimates can be very large, which is signaled by a
#' warning that fitted probabilities are numerically 0 or 1, or, for a
#' converged fit, within 1e-12 of 0 or 1. Because the dispersion is fixed,
#' \code{"pgee-jeffreys"} gives finite estimates under separation with
#' \code{orstr = "independence"} for every supported link (Touloumis, 2026b),
#' and so does \code{"brgee-naive"} with the logit link, for which it coincides
#' with \code{"pgee-jeffreys"} at the default \code{jeffreys_power = 0.5}.
#' With other odds-ratio structures finiteness is not guaranteed, although
#' \code{"pgee-jeffreys"} usually converges in practice.
#' \code{"opgee-jeffreys"} and \code{"hpgee-jeffreys"} start from the finite
#' independence estimate and are finite whenever their single scoring step is
#' defined.
#' \code{"brgee-robust"} and \code{"brgee-empirical"} do not keep the estimates
#' finite, and the bias-corrected methods fail because they start from the
#' ordinary GEE fit.
#'
#' @inherit geewa return
#'
#' @references
#' Touloumis, A. (2026a) Bias-reduced GEE via adjusted estimating equations,
#' with odds-ratio extensions. \emph{arXiv preprint} arXiv:2606.16043.
#'
#' Touloumis, A. (2026b) Jeffreys-type penalized GEE for correlated binary data
#' with an odds-ratio parameterization. \emph{arXiv preprint} arXiv:2606.16058.
#'
#' @seealso \code{\link{geewa}}, \code{\link{geer_control}},
#'   \code{\link{summary.geer}}, \code{\link{vcov.geer}},
#'   \code{\link{anova.geer}}, \code{\link{step_p}}, \code{\link{geecriteria}},
#'   \code{\link{runs_test}}.
#'
#' @examples
#' data("cholecystectomy", package = "geer")
#' fit <- geewa_binary(
#'   formula = pain ~ treatment + gender + age,
#'   link = "logit",
#'   data = cholecystectomy,
#'   id = id,
#'   repeated = time,
#'   orstr = "exchangeable"
#' )
#' summary(fit)
#'
#' ## Jeffreys-type penalized GEE
#' fit_pgee <- update(fit, method = "pgee-jeffreys")
#' summary(fit_pgee)
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
                         control_glm = list(),
                         alpha_vector = NULL,
                         subset,
                         na.action,
                         ...) {
  if (!missing(control)) {
    check_unused_dots(list(...), "geewa_binary")
  }
  call <- match.call(expand.dots = TRUE)
  mcall <- match.call(expand.dots = FALSE)
  caller_env <- parent.frame()
  link <- match_link(link)
  family <- normalize_family(stats::binomial(link = link))
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
  beta_zero <- compute_binary_start_values(
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
  orstr <- match_orstr(orstr)
  if (identical(orstr, "independence")) {
    alpha_vector <- rep.int(1, choose(max(repeated), 2))
  } else if (identical(orstr, "fixed")) {
    pairs_no <- choose(max(repeated), 2)
    if (!is.numeric(alpha_vector)) {
      stop("'alpha_vector' must be a numeric vector", call. = FALSE)
    }
    if (length(alpha_vector) != pairs_no) {
      stop(
        "'alpha_vector' must be a numeric vector of length ",
        pairs_no,
        call. = FALSE
      )
    }
    if (any(!is.finite(alpha_vector)) || any(alpha_vector <= 0)) {
      stop(
        "'alpha_vector' must be finite and strictly positive when ",
        "'orstr = \"fixed\"'",
        call. = FALSE
      )
    }
  } else {
    alpha_vector <- as.numeric(get_marginalized_odds_ratios_or(
      round(y),
      id,
      repeated,
      weights,
      control$or_adding,
      orstr
    ))
  }
  alpha_independence <- rep.int(1, choose(max(repeated), 2))
  fit_pass <- function(beta, pass, previous) {
    iterations <- get_pass_iterations(pass, control)
    fit_geesolver_or(
      y, model_matrix, id, repeated, weights, link,
      as.numeric(beta), offset,
      iterations[[1L]], tolerance, iterations[[2L]], iterations[[3L]],
      control$jeffreys_power, pass$method,
      if (pass$independence) alpha_independence else alpha_vector
    )
  }
  geesolver_fit <- run_geer_estimation_passes(
    fit_pass = fit_pass,
    method = method,
    beta_start = beta_zero,
    check_first = function(fit, method) {
      check_geer_first_pass(fit, method, tolerance, "geewa_binary")
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
    association_structure = orstr,
    row_order = inputs$row_order
  )
  fit <- build_geer_final_fit(
    fit = fit,
    geesolver_fit = geesolver_fit,
    tolerance = tolerance,
    method = method,
    family = family,
    association_structure = orstr,
    repeated = repeated,
    fit_function = "geewa_binary"
  )
  fit <- new_geer(fit)
  fit <- validate_geer(fit)
  fit
}
