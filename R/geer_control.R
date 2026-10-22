#' @title
#' Control Parameters for geer Model Fitting
#'
#' @description
#' Sets the control parameters of the fitting algorithm of \code{\link{geewa}}
#' and \code{\link{geewa_binary}}.
#'
#' @param tolerance strictly positive number giving the convergence tolerance.
#'   Defaults to \code{1e-06}.
#' @param maxiter positive integer giving the maximum number of iterations.
#'   Defaults to \code{500}.
#' @param or_adding strictly positive number added to each cell of the
#'   marginalized \eqn{2 \times 2}{2 x 2} tables from which \code{\link{geewa_binary}}
#'   estimates the working odds ratios. Defaults to \code{0.5}.
#' @param step_maxiter positive integer giving the maximum number of
#'   step-halvings in an iteration. Defaults to \code{10}.
#' @param step_multiplier strictly positive number scaling the proposed step
#'   before step-halving. Values below 1 damp the step, which can help when the
#'   algorithm diverges from poor starting values. Defaults to \code{1}.
#' @param jeffreys_power strictly positive number giving the power of the
#'   Jeffreys-type penalty of the methods \code{"pgee-jeffreys"},
#'   \code{"opgee-jeffreys"}, and \code{"hpgee-jeffreys"} (Touloumis, 2026b).
#'   Defaults to \code{0.5}, the standard Jeffreys-type penalty.
#'
#' @return
#' A list with components \code{tolerance}, \code{maxiter},
#' \code{or_adding}, \code{step_maxiter}, \code{step_multiplier}, and
#' \code{jeffreys_power}.
#'
#' @references
#' Touloumis, A. (2026b) Jeffreys-type penalized GEE for correlated binary data
#' with an odds-ratio parameterization. \emph{arXiv preprint} arXiv:2606.16058.
#'
#' @seealso \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' ## Default control parameters
#' geer_control()
#'
#' ## Damped steps, useful when the algorithm diverges
#' geer_control(step_multiplier = 0.5)
#'
#' @export
geer_control <- function(tolerance = 1e-06,
                         maxiter = 500L,
                         or_adding = 0.5,
                         step_maxiter = 10L,
                         step_multiplier = 1,
                         jeffreys_power = 0.5) {
  if (!is_positive_scalar(tolerance)) {
    stop("'tolerance' must be a positive finite number", call. = FALSE)
  }
  if (!is_positive_integer_scalar(maxiter)) {
    stop("'maxiter' must be a positive integer", call. = FALSE)
  }
  if (!is_positive_scalar(or_adding)) {
    stop("'or_adding' must be a positive finite number", call. = FALSE)
  }
  if (!is_positive_integer_scalar(step_maxiter)) {
    stop("'step_maxiter' must be a positive integer", call. = FALSE)
  }
  if (!is_positive_scalar(step_multiplier)) {
    stop("'step_multiplier' must be a positive finite number", call. = FALSE)
  }
  if (!is_positive_scalar(jeffreys_power)) {
    stop("'jeffreys_power' must be a positive finite number", call. = FALSE)
  }
  list(
    tolerance = tolerance,
    maxiter = as.integer(maxiter),
    or_adding = or_adding,
    step_maxiter = as.integer(step_maxiter),
    step_multiplier = as.numeric(step_multiplier),
    jeffreys_power = jeffreys_power
  )
}
