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
#' @inheritParams stats::step
#' @param direction character string specifying the direction of the stepwise
#'   search: backward elimination (\code{"backward"}), forward selection
#'   (\code{"forward"}) or bidirectional (\code{"both"}). Defaults to
#'   \code{"backward"}.
#' @param p_enter numeric value specifying the p-value threshold
#'   for adding a term during forward steps. Defaults to \code{0.15}.
#' @param p_remove numeric value specifying the p-value threshold
#'   for removing a term during backward steps. Defaults to \code{0.15}.
#' @param steps positive integer giving the maximum number of steps to perform.
#'   The algorithm stops earlier if no eligible move is found. Defaults to
#'   \code{1000}.
#'
#' @details
#' \code{step_p()} repeatedly calls \code{\link{add1.geer}} and
#' \code{\link{drop1.geer}}.
#'
#' The set of models searched is determined by the \code{scope} argument. The
#' right-hand side of its \code{lower} component is always retained in the
#' model, and the right-hand side of its \code{upper} component defines the
#' largest model that can be considered. If \code{scope} is a single formula,
#' it specifies the \code{upper} component and the \code{lower} model is empty.
#'
#' Candidate additions and deletions respect model hierarchy: a move is allowed
#' only if it preserves all lower-order component terms implied by any
#' higher-order interactions in the model.
#'
#' If \code{scope} is omitted, the search is restricted to the terms in the
#' initial model.
#'
#' With \code{direction = "both"}, backward elimination is attempted first at
#' each step. A forward addition is considered only if no eligible backward
#' deletion is found. The algorithm stops when no candidate move satisfies the
#' relevant p-value threshold, when the maximum number of steps is reached, or
#' when an immediate add-remove cycle involving the same term would occur.
#'
#' Details of the hypothesis tests controlled by \code{test} are given in
#' Rotnitzky and Jewell (1990). The option \code{test = "working-lrt"} is
#' valid only when the model is fitted with an independence working association
#' structure. Otherwise, an error is returned.
#'
#' When \code{test \%in\% c("wald", "score")}, the \code{pmethod} argument is
#' ignored and \code{cov_type} specifies the covariance estimator used to
#' compute the test statistic. For modified working tests, \code{cov_type}
#' determines the covariance matrix used to form the coefficients of the sum of
#' independent chi-squared random variables, and \code{pmethod} specifies the
#' approximation used to compute the p-value.
#'
#' @return
#' A fitted model object of class \code{"geer"} corresponding to the final
#' selected model. The returned object also includes an \code{anova}
#' component of class \code{c("anova", "data.frame")} summarizing the
#' stepwise sequence. This table has one
#' row per step with columns \code{Step} (the term added or removed),
#' \code{Df} (degrees of freedom of the test), \code{Chi} (test statistic),
#' \code{Pr(>Chi)} (p-value), and \code{CIC} (Correlation Information
#' Criterion).
#'
#' @inherit add1.geer references
#'
#' @seealso \code{\link{add1.geer}}, \code{\link{drop1.geer}},
#'   \code{\link{anova.geer}}, \code{\link{geecriteria}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' \donttest{
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]
#'
#' full_fit <- geewa_binary(
#'   formula = status ~ (baseline + treatment + gender + visit + age)^2,
#'   link = "probit",
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   orstr = "independence",
#'   method = "pgee-jeffreys"
#' )
#'
#' ## Backward elimination within the initial model terms
#' step_p(
#'   full_fit,
#'   direction = "backward",
#'   test = "wald",
#'   cov_type = "bias-corrected",
#'   p_remove = 0.10
#' )
#'
#' ## Bidirectional selection with an explicit scope
#' step_p(
#'   full_fit,
#'   scope = list(
#'     lower = ~ baseline + treatment,
#'     upper = ~ (baseline + treatment + gender + visit + age)^2
#'   ),
#'   direction = "both",
#'   test = "score",
#'   cov_type = "robust",
#'   p_enter = 0.10,
#'   p_remove = 0.15,
#'   steps = 50
#' )
#' }
#'
#' @export
step_p <- function(object,
                   scope,
                   direction = c("backward", "forward", "both"),
                   p_enter = 0.15,
                   p_remove = 0.15,
                   test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
                   cov_type = c("bias-corrected", "robust", "df-adjusted", "naive"),
                   pmethod = c("rao-scott", "satterthwaite"),
                   steps = 1000) {
  object <- check_geer_object(object)
  direction <- match_direction_type(direction[1L])
  step_args <- validate_step_thresholds(
    p_enter = p_enter,
    p_remove = p_remove
  )
  p_enter <- step_args$p_enter
  p_remove <- step_args$p_remove
  steps <- validate_step_count(steps)
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
  switch(
    direction,
    backward = .step_p_run_backward(
      object = object,
      scope = scope_value,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_remove,
      steps = steps
    ),
    forward = .step_p_run_forward(
      object = object,
      scope = scope_value,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod,
      pvalue = p_enter,
      steps = steps
    ),
    both = .step_p_run_both(
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
