#' Auxiliary for Controlling GEE fitting
#'
#' Auxiliary function for GEE fitting via the functions \link{geewa} and
#' \link{geewa_binary}. It may be used to construct the \code{control} argument.
#'
#' @param tolerance positive convergence tolerance. Default is `1e-06`.
#' @param maxiter integer giving the maximal number of iterations
#'        allowed. Default is `500`.
#' @param or_adding_constant positive constant to be added at each cell of
#'        the full marginalized contingency table. Default is `0.5`.
#' @param jeffreys_power positive real indicating the power of the Jeffreys-prior
#'        penalty. Default is `0.5`.
#' @param step_maxiter integer giving the maximal number of step halving steps
#'        allowed. Default is `10`.
#' @param step_multiplier positive real indicating the multiplier for the
#'        stepsize. The smaller it is the smaller the steps are. Default
#'        is `1`.
#'
#' @details
#'
#' When the function \link{geewa} is used, \code{or_adding_constant} is ignored.
#'
#' The argument \code{jeffreys_power} is ignored unless
#' \code{method = "pgee_jeffreys"} in \link{geewa} or \link{geewa_binary}.
#' This argument controls the amount of and its default value of `1/2`
#' corresponds to the GEE analogue of the Jeffreys-prior penalty.
#'
#' When \code{bstart = NULL} in \link{geewa} or \link{geewa_binary}, then the
#' values of \code{geer.control} are passed to the \link{glm} function that
#' generates the initial value of the regression parameter vector.
#'
#' @return a list with components named as the arguments.
#'
#' @author Anestis Touloumis \email{A.Touloumis@brighton.ac.uk}
#'
#' @seealso \link{geewa} and \link{geewa_binary}
#'
#' @examples
#' data("leprosy")
#' fitted_model_gee <- geewa(
#'    formula = bacilli ~ time+time:I(treatment == "C"),
#'    data = leprosy,
#'    id = id,
#'    family = poisson(link = "log"),
#'    correlation_structure = "exchangeable",
#'    method = "gee")
#' update(fitted_model_gee, control = list(maxiter = 1, step_maxiter = 1))
#'
#' @export
geer_control <- function(tolerance = 1e-06,
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
