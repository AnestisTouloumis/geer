#' @title
#' Stepwise Model Selection by Hypothesis Testing
#'
#' @description
#' Performs stepwise model selection for a fitted \code{geer} object using
#' repeated single-term additions (\code{\link{add1.geer}}) and deletions
#' (\code{\link{drop1.geer}}), with candidate moves evaluated by hypothesis
#' tests.
#'
#' @inheritParams anova.geer
#' @param scope optional formula, or list with components \code{upper} and
#'   \code{lower} (both formulas), defining the range of models examined in the
#'   stepwise search. See Details. If omitted, the search is restricted to the
#'   terms of the initial model.
#' @param direction character string specifying the direction of the search:
#'   \code{"backward"}, \code{"forward"}, or \code{"both"}. Defaults to
#'   \code{"backward"}.
#' @param p_enter a single number strictly between 0 and 1 giving the p-value
#'   threshold for adding a term. It must not exceed \code{p_remove} when
#'   \code{direction = "both"}. Defaults to \code{0.15}.
#' @param p_remove a single number strictly between 0 and 1 giving the p-value
#'   threshold for removing a term. Defaults to \code{0.15}.
#' @param steps nonnegative integer giving the maximum number of steps. Defaults
#'   to \code{1000}.
#'
#' @details
#' \code{step_p()} repeatedly calls \code{\link{add1.geer}} and
#' \code{\link{drop1.geer}}. The terms of \code{scope$lower} are always kept,
#' and \code{scope$upper} defines the largest model considered. A single
#' formula is taken as \code{upper}, with an empty \code{lower}. Moves respect
#' marginality. With \code{direction = "both"}, a deletion is tried first at
#' each step and an addition only if no deletion qualifies. The search stops
#' when no move meets its p-value threshold, after \code{steps} steps, or when
#' a move would return to a model already visited.
#'
#' The tests are described in \code{\link{anova.geer}}. The CIC column is
#' computed with \code{cov_type}, so \code{cov_type = "jackknife"} refits every
#' candidate model once per cluster and can be slow.
#'
#' @return
#' The selected model, a fitted model object of class \code{"geer"} with an
#' additional component \code{anova} of class \code{c("geer_anova", "anova",
#' "data.frame")} with one row for the initial model and one for each step,
#' whose row names give the term added or removed, and columns \code{Step}
#' (\code{NA} for the initial model), \code{Df}, \code{Chi}, \code{Pr(>Chi)},
#' and \code{CIC}.
#'
#' @inherit add1.geer references
#'
#' @seealso \code{\link{add1.geer}}, \code{\link{drop1.geer}},
#'   \code{\link{anova.geer}}, \code{\link{geecriteria}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' \donttest{
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", ]
#' fit <- geewa_binary(
#'   formula = status ~ (baseline + treatment + gender + visit + age)^2,
#'   link = "probit",
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   method = "pgee-jeffreys"
#' )
#' ## Backward elimination
#' step_p(fit, p_remove = 0.10)
#'
#' ## Bidirectional selection within a scope
#' step_p(
#'   fit,
#'   scope = list(
#'     lower = ~ baseline + treatment,
#'     upper = ~ (baseline + treatment + gender + visit + age)^2
#'   ),
#'   direction = "both",
#'   p_enter = 0.10,
#'   test = "score",
#'   cov_type = "robust"
#' )
#' }
#'
#' @export
step_p <- function(object,
                   scope,
                   direction = c("backward", "forward", "both"),
                   p_enter = 0.15,
                   p_remove = 0.15,
                   test = c("wald", "score", "working-wald", "working-score",
                            "working-lrt"),
                   cov_type = c("bias-corrected", "robust", "df-adjusted",
                                "jackknife", "naive"),
                   pmethod = c("rao-scott", "satterthwaite"),
                   steps = 1000) {
  object <- check_geer_object(object)
  direction <- match_direction_type(direction[1L])
  step_args <- check_step_thresholds(
    p_enter = p_enter,
    p_remove = p_remove
  )
  p_enter <- step_args$p_enter
  p_remove <- step_args$p_remove
  if (identical(direction, "both") && p_enter > p_remove) {
    stop(
      paste0(
        "'p_enter' must not exceed 'p_remove' when direction = \"both\", ",
        "otherwise terms can cycle in and out"
      ),
      call. = FALSE
    )
  }
  steps <- check_step_count(steps)
  opts <- normalize_geer_test_options(
    test = test[1L],
    cov_type = cov_type[1L],
    pmethod = pmethod[1L],
    object = object
  )
  test <- opts$test
  cov_type <- opts$cov_type
  pmethod <- opts$pmethod
  scope_value <- if (missing(scope)) NULL else scope
  if (is.null(scope_value) && identical(direction, "forward")) {
    warning(
      paste0(
        "no 'scope' was supplied for direction = \"forward\", so the initial ",
        "model is returned unchanged"
      ),
      call. = FALSE
    )
  }
  switch(
    direction,
    backward = run_step_backward(
      object = object,
      scope = scope_value,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_remove,
      steps = steps
    ),
    forward = run_step_forward(
      object = object,
      scope = scope_value,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_enter,
      steps = steps
    ),
    both = run_step_both(
      object = object,
      scope = scope_value,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      p_enter = p_enter,
      p_remove = p_remove,
      steps = steps
    )
  )
}
