#' @title
#' Choose a Model by Hypothesis Testing in a Stepwise Algorithm
#'
#' @description
#' Select a formula-based model by hypothesis testing.
#'
#' @inheritParams anova.geer
#' @inheritParams stats::step
#' @param object an object representing a model of the class \code{geer}. This is
#'        used as the initial model in the stepwise search.
#' @param direction character indicating the mode of the stepwise search.
#'        Options include backward elimination (\code{"backward"}),
#'        forward selection (\code{"forward"}) and bidirectional elimination
#'        (\code{"both"}). By default, a backward elimination is performed.
#' @param p_enter numeric between 0 and 1 indicating the p-value threshold for
#'        adding variables in the stepwise search. By default,
#'        \code{p_enter = 0.15}.
#' @param p_remove numeric between 0 and 1 indicating the p-value threshold for
#'        removing variables in the stepwise search. By default,
#'        \code{p_remove = 0.15}.
#'
#' @details
#' \code{step_p} uses \code{\link{add1}} and \code{\link{drop1}} repeatedly; it
#' will work for any model for which they work.
#'
#' The set of models searched is determined by the \code{scope} argument. The
#' right-hand-side of its \code{lower} component is always included in the model,
#' and right-hand-side of the model is included in the \code{upper} component.
#' If \code{scope} is a single formula, it specifies the \code{upper} component,
#' and the \code{lower} model is empty. If \code{scope} is missing, the initial
#' model is used as the \code{upper} model.
#'
#' Models specified by \code{scope} can be templates to update \code{object} as
#' used by \code{\link[stats]{update.formula}}. So using \code{.} in a
#' \code{scope} formula means ‘what is already there’, with \code{.^2}
#' indicating all interactions of existing terms.
#'
#' Details about the testing procedures implied by the \code{test} argument can
#' be found in \cite{Rotnitzky and Jewell (1990)}. Note that
#' \code{test = "working-lrt"} is only available to fitted models with an
#' independence working association structure. Otherwise, an error message is
#' returned.
#'
#' When \code{test = "wald"} or \code{test = "score"}, the \code{p_method}
#' argument is ignored and the \code{cov_type} argument specifies the covariance
#' matrix estimate of the estimated regression parameters used to calculate the
#' corresponding test statistic. Otherwise, the \code{cov_type} argument specifies the
#' covariance matrix estimate of the estimated regression parameters used to
#' calculate the coefficients of the independent chi-squared random variables,
#' and the \code{p_method} argument specifies the approximation method used to
#' calculate the p-value of the test statistic.
#'
#' @returns
#' The stepwise-selected model is returned, with an \code{anova} component
#' corresponding to the steps taken in the search.
#'
#' @inherit add1.geer references
#'
#' @seealso \code{\link{add1}} and \code{\link{drop1}}.
#'
#'
#' @export
step_p <-
  function(object,
           scope,
           direction = c("backward", "forward", "both"),
           p_enter = 0.15,
           p_remove = 0.15,
           test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite"),
           steps = 1000) {
    if (!("geer" %in% class(object)))
      stop("'object' must be of 'geer' class")
    test <- match.arg(test)
    if (test %in% c("working-wald", "working-score", "working-lrt"))
      pmethod <- match.arg(pmethod)
    if (test == "working-lrt" & object$association_structure != "independence")
      stop("the modified working lrt can only be applied to an independence working model")
    cov_type <- match.arg(cov_type)
    if (!(is.numeric(p_enter) & (0 < p_enter) & (p_enter < 1)))
      stop("'p_enter' must be a numeric object between 0 and 1")
    if (!(is.numeric(p_remove) & (0 < p_remove) & (p_remove < 1)))
      stop("'p_enter' must be a numeric object between 0 and 1")
    direction = match.arg(direction)
    ans <-
      switch(direction,
             backward = step_p_backward(object, scope, test, cov_type, pmethod,
                                        p_remove, steps),
             forward = step_p_forward(object, scope, test, cov_type, pmethod,
                                      p_enter, steps),
             both = step_p_both(object, scope, test, cov_type, pmethod,
                                p_enter, p_remove, steps)
      )
    ans
  }
