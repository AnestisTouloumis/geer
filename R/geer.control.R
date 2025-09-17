#' @title
#' Control Parameters for Fitting \code{geer} Objects
#'
#' @description
#' Create a list of control parameters for use in \code{\link{geewa}} and
#' \code{\link{geewa_binary}}. This function is typically used internally, but
#' may also be supplied directly as the \code{control} argument when calling
#' those functions.
#'
#' @param tolerance positive convergence tolerance. Default is \code{1e-6}.
#' @param maxiter positive integer giving the maximal number of iterations
#'        allowed. Default is \code{500}.
#' @param or_adding positive constant added to each cell of the full
#'        marginalized contingency table. Default is \code{0.5}.
#' @param step_maxit positive integer giving the maximum number of step-halving
#'        steps allowed. Default is \code{10}.
#' @param step_multi positive integer used as a multiplier for the step. Default
#'        is \code{1}.
#' @param jeffreys_power positive constant indicating the power of the
#'        Jeffreys-prior penalty. Default is \code{0.5}.
#'
#' @details
#' The \code{or_adding_constant} argument is ignored by \code{\link{geewa}}.
#'
#' The \code{jeffreys_power} argument is only used if \code{method = "pgee_jeffreys"}
#' in \code{\link{geewa}} or \code{\link{geewa_binary}}.
#'
#' @returns
#' A list with components named as the arguments, suitable for passing as the
#' \code{control} argument to \code{geewa} or \code{geewa_binary}.
#'
#' @export
geer_control <- function(tolerance = 1e-06,
                         maxiter = 500,
                         or_adding = 0.5,
                         step_maxit = 10,
                         step_multi = 1,
                         jeffreys_power = 1/2) {
  if (!is.numeric(tolerance) || tolerance <= 0)
    stop("value of 'tolerance' must be > 0")
  if (!is.numeric(maxiter) || maxiter <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(step_multi) || step_multi <= 0)
    stop("value of 'step_multi' must be > 0")
  if (!is.numeric(step_maxit) || step_maxit < 1) {
    warning("`step_maxit = ",
            deparse(step_maxit),
            "` is not a permissible value. Defaulting to 10")
    step_maxit <- 10
  }
  if (!is.numeric(jeffreys_power) || jeffreys_power <= 0)
    stop("value of 'jeffreys_power' must be > 0")
  if (!is.numeric(or_adding) || or_adding <= 0)
    stop("value of 'or_adding_constant' must be > 0")
  list(tolerance = tolerance,
       maxiter = as.integer(maxiter),
       or_adding = or_adding,
       step_maxit = as.integer(step_maxit),
       step_multi = step_multi,
       jeffreys_power = jeffreys_power
  )
}
