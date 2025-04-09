#' @title
#' Auxiliary for Controlling Fitting of \code{geer} Objects
#'
#' @description
#' Auxiliary function for \code{\link{geewa}} and \code{\link{geewa_binary}}.
#' Typically used internally but may be used to construct a \code{control} argument
#' to either function.
#'
#'
#' @param tolerance positive convergence tolerance. By default,
#'        \code{tolerance = 1e-06}.
#' @param maxiter positive integer giving the maximal number of iterations
#'        allowed. By default, \code{default = 500}.
#' @param or_adding_constant positive constant to be added at each cell of
#'        the full marginalized contingency table. By default,
#'        \code{or_adding_constant = `0.5`}.
#' @param step_maxiter positive integer giving the maximal number of step halving
#'        steps allowed. By default, \code{step_maxiter = 10}.
#' @param step_multiplier positive integer used as a multiplier for the step. By
#'        default, \code{steo_multiplier = 1}.
#' @param jeffreys_power positive real indicating the power of the Jeffreys-prior
#'        penalty. By default, \code{jeffreys_power = 0.5}.
#'
#' @details
#' The \code{or_adding_constant} is ignored in the function \code{\link{geewa}}.
#'
#' The \code{jeffreys_power} is ignored unless \code{method = "pgee_jeffreys"} in
#' functions \code{\link{geewa}} or \code{\link{geewa_binary}}.
#'
#'
#' @returns
#' A list with components named as the arguments.
#'
#'
#' @rdname geer.control
#' @export geer.control
geer.control <- function(tolerance = 1e-06,
                         maxiter = 500,
                         or_adding_constant = 0.5,
                         step_maxiter = 10,
                         step_multiplier = 1,
                         jeffreys_power = 1/2) {
  if (!is.numeric(tolerance) || tolerance <= 0)
    stop("value of 'tolerance' must be > 0")
  if (!is.numeric(maxiter) || maxiter <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(step_multiplier) || step_multiplier <= 0)
    stop("value of 'step_multiplier' must be > 0")
  if (!is.numeric(step_maxiter) || step_maxiter < 1) {
    warning("`step_maxiter = ",
            deparse(step_maxiter),
            "` is not a permissible value. Defaulting to 10")
    step_maxiter <- 10
  }
  if (!is.numeric(jeffreys_power) || jeffreys_power <= 0)
    stop("value of 'jeffreys_power' must be > 0")
  if (!is.numeric(or_adding_constant) || or_adding_constant <= 0)
    stop("value of 'or_adding_constant' must be > 0")
  list(tolerance = tolerance,
       maxiter = as.integer(maxiter),
       or_adding_constant = or_adding_constant,
       step_maxiter = as.integer(step_maxiter),
       step_multiplier = step_multiplier,
       jeffreys_power = jeffreys_power
  )
}
