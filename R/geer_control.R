#' @title
#' Control Parameters for geer Model Fitting
#'
#' @description
#' Creates a list of control parameters for iterative fitting in
#' \code{\link{geewa}} and \code{\link{geewa_binary}}. This function is used
#' internally, but may also be supplied directly through the \code{control}
#' argument of those functions.
#'
#' @param tolerance strictly positive convergence tolerance. Defaults to
#'   \code{1e-6}.
#' @param maxiter positive integer giving the maximum number of fitting
#'   iterations. Defaults to \code{500}.
#' @param or_adding strictly positive continuity-correction constant added to
#'   each cell of the marginalized \eqn{2 \times 2} contingency tables used to
#'   estimate working odds ratios in \code{\link{geewa_binary}}. Defaults to
#'   \code{0.5}. Ignored by \code{\link{geewa}}.
#' @param step_maxiter positive integer giving the maximum number of
#'   step-halving attempts allowed within an iteration. Defaults to \code{10}.
#' @param step_multiplier positive integer used to scale the proposed step
#'   before step-halving begins. A value greater than \code{1} enlarges the
#'   initial step; the default of \code{1} leaves the scoring step unscaled.
#' @param jeffreys_power strictly positive constant giving the power of the
#'   Jeffreys-prior penalty. Defaults to \code{0.5}, which corresponds to the
#'   GEE analogue of the standard Jeffreys prior.
#'
#' @details
#' The \code{jeffreys_power} argument is used only when
#' \code{method \%in\% c("pgee-jeffreys", "opgee-jeffreys", "hpgee-jeffreys")}
#' in \code{\link{geewa}} or \code{\link{geewa_binary}}.
#'
#' @return
#' A named list with elements \code{tolerance}, \code{maxiter},
#' \code{or_adding}, \code{step_maxiter}, \code{step_multiplier}, and
#' \code{jeffreys_power}, suitable for use as the \code{control} argument to
#' \code{\link{geewa}} or \code{\link{geewa_binary}}.
#'
#' @references
#' Touloumis, A. (2026) Jeffreys-prior penalized GEE for correlated binary
#' data with an odds-ratio parameterization. Preprint.
#'
#' @seealso \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' ## Default control parameters
#' geer_control()
#'
#' ## Tighter convergence tolerance and fewer maximum iterations
#' geer_control(tolerance = 1e-8, maxiter = 200)
#'
#' ## Custom continuity correction for odds-ratio estimation
#' geer_control(or_adding = 1)
#'
#' ## Weaker Jeffreys-prior penalty for Jeffreys-type fits
#' geer_control(jeffreys_power = 0.1)
#'
#' @export
geer_control <- function(tolerance = 1e-06,
                         maxiter = 500,
                         or_adding = 0.5,
                         step_maxiter = 10,
                         step_multiplier = 1,
                         jeffreys_power = 0.5) {
  if (!is_positive_scalar(tolerance)) {
    stop("'tolerance' must be a positive number", call. = FALSE)
  }
  if (!is_positive_integer_scalar(maxiter)) {
    stop("'maxiter' must be a positive integer", call. = FALSE)
  }
  if (!is_positive_scalar(or_adding)) {
    stop("'or_adding' must be a positive number", call. = FALSE)
  }
  if (!is_positive_integer_scalar(step_maxiter)) {
    stop("'step_maxiter' must be a positive integer", call. = FALSE)
  }
  if (!is_positive_integer_scalar(step_multiplier)) {
    stop("'step_multiplier' must be a positive integer", call. = FALSE)
  }
  if (!is_positive_scalar(jeffreys_power)) {
    stop("'jeffreys_power' must be a positive number", call. = FALSE)
  }
  list(
    tolerance = tolerance,
    maxiter = as.integer(maxiter),
    or_adding = or_adding,
    step_maxiter = as.integer(step_maxiter),
    step_multiplier = as.integer(step_multiplier),
    jeffreys_power = jeffreys_power
  )
}
