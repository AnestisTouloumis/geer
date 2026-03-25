#' @title
#' Control Parameters for \code{geer} Model Fitting
#'
#' @description
#' Create a list of control parameters for iterative fitting in
#' \code{\link{geewa}} and \code{\link{geewa_binary}}. This function is used
#' internally, but may also be supplied directly through the \code{control}
#' argument of those functions.
#'
#' @param tolerance Strictly positive convergence tolerance. Default is
#'   \code{1e-6}.
#' @param maxiter Positive integer giving the maximum number of fitting
#'   iterations. Default is \code{500}.
#' @param or_adding Strictly positive constant added to each cell of the full
#'   marginalized contingency table. Default is \code{0.5}.
#' @param step_maxit Positive integer giving the maximum number of step-halving
#'   attempts allowed within an iteration. Default is \code{10}.
#' @param step_multi Positive integer used to scale the proposed step before
#'   step-halving. Default is \code{1}.
#' @param jeffreys_power Strictly positive constant giving the power of the
#'   Jeffreys-prior penalty. Default is \code{0.5}.
#'
#' @details
#' The \code{or_adding} argument is ignored by \code{\link{geewa}}.
#'
#' The \code{jeffreys_power} argument is only used when
#' \code{method = "pgee-jeffreys"} in \code{\link{geewa}} or
#' \code{\link{geewa_binary}}.
#'
#' @return
#' A named list suitable for use as the \code{control} argument to
#' \code{geewa} or \code{geewa_binary}.
#'
#' @export
geer_control <- function(tolerance = 1e-06,
                         maxiter = 500,
                         or_adding = 0.5,
                         step_maxit = 10,
                         step_multi = 1,
                         jeffreys_power = 0.5) {
  if (!is_pos_scalar(tolerance)) {
    stop("'tolerance' must be a positive number", call. = FALSE)
  }
  if (!is_pos_int_scalar(maxiter)) {
    stop("'maxiter' must be a positive integer", call. = FALSE)
  }
  if (!is_pos_scalar(or_adding)) {
    stop("'or_adding' must be a positive number", call. = FALSE)
  }
  if (!is_pos_int_scalar(step_maxit)) {
    stop("'step_maxit' must be a positive integer", call. = FALSE)
  }
  if (!is_pos_int_scalar(step_multi)) {
    stop("'step_multi' must be a positive integer", call. = FALSE)
  }
  if (!is_pos_scalar(jeffreys_power)) {
    stop("'jeffreys_power' must be a positive number", call. = FALSE)
  }
  list(
    tolerance = tolerance,
    maxiter = as.integer(maxiter),
    or_adding = or_adding,
    step_maxit = as.integer(step_maxit),
    step_multi = as.integer(step_multi),
    jeffreys_power = jeffreys_power
  )
}
