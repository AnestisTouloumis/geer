#' @title
#' Fitting (Adjusted) Generalized Estimating Equations for Binary Responses
#'
#' @description
#' Fits a marginal model for repeated or clustered binary responses using
#' generalized estimating equations (GEE). Supported estimation methods include
#' the traditional GEE, bias-reduced GEE, bias-corrected GEE, and
#' Jeffreys-type penalized GEE.
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
#'        \code{orstr}. Defaults to \code{NULL}.
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
#' The marginal mean model always uses a \code{binomial} family with the
#' specified \code{link}. Within-cluster association is modeled through
#' marginalized pairwise odds ratios rather than a correlation structure. The
#' scale parameter is fixed at \code{phi = 1}.
#'
#' For the construction of the \code{formula} argument, see the documentation
#' of \code{\link[stats]{glm}} and \code{\link[stats]{formula}}.
#'
#' The \code{data} must be in long format (one row per observation). See
#' \code{\link[stats]{reshape}} for details on reshaping between long and wide formats.
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
#' The length of \code{id} and, when provided, of \code{repeated} and
#' \code{weights} must equal the number of observations.
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
#' \code{"hpgee-jeffreys"}, \code{converged} is \code{TRUE} in the
#' returned object, because these methods produce their estimate via a single
#' correction step applied to an already-converged fit. The exception is a
#' numerical failure in that correction step: a warning is then issued, the
#' preliminary estimates of the first stage are returned and \code{converged}
#' is \code{FALSE}. The working odds ratios in \code{alpha} are computed once
#' from the data, so they are the same in both stages.
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
#'
#' @seealso
#' \code{\link{geewa}}, \code{\link{geer_control}}, \code{\link{summary.geer}},
#' \code{\link{vcov.geer}}, \code{\link{anova.geer}}, \code{\link{step_p}},
#' \code{\link{geecriteria}}, \code{\link{runs_test}}.
#'
#' @examples
#' \donttest{
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]
#' fit_probit <- geewa_binary(
#'   formula = status ~ baseline + treatment * gender + visit * age,
#'   link = "probit",
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   orstr = "independence",
#'   method = "pgee-jeffreys"
#' )
#' summary(fit_probit, cov_type = "bias-corrected")
#' }
#'
#' data("cholecystectomy", package = "geer")
#' fit_gee <- geewa_binary(
#'   formula = pain ~ treatment + gender + age,
#'   link = "logit",
#'   data = cholecystectomy,
#'   id = id,
#'   orstr = "independence",
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
#' ## Penalized GEE with exchangeable odds-ratio structure
#' fit_pgee <- geewa_binary(
#'   formula = pain ~ treatment + gender + age,
#'   link = "logit",
#'   data = cholecystectomy,
#'   id = id,
#'   control = geer_control(jeffreys_power = 1),
#'   orstr = "exchangeable",
#'   method = "pgee-jeffreys"
#' )
#' summary(fit_pgee, cov_type = "robust")
#' }
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
                         subset,
                         na.action,
                         ...) {
  ## when both 'control' and 'control_glm' are supplied, nothing consumes '...'
  if (!missing(control) && !missing(control_glm)) {
    check_unused_dots(list(...), "geewa_binary")
  }
  ## call, family and common input preparation
  call <- match.call(expand.dots = TRUE)
  mcall <- match.call(expand.dots = FALSE)
  caller_env <- parent.frame()
  family <- normalize_family(stats::binomial(link = link))
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
  check_choice(orstr, geer_orstr_choices, "orstr")
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
  fit_pass <- function(beta, pass, previous) {
    iterations <- if (pass$one_step) {
      list(1L, 1L, 1)
    } else {
      list(maxiter, control$step_maxiter, control$step_multiplier)
    }
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
    association_structure = orstr,
    row_order = inputs$row_order
  )
  fit <- finalize_geer_fit(
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
