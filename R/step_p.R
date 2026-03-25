#' @title
#' Stepwise Model Selection by Hypothesis Testing
#'
#' @description
#' Perform stepwise model selection for fitted \code{geer} models using
#' repeated single-term additions with \code{\link{add1}} and/or single-term
#' deletions with \code{\link{drop1}}, where candidate moves are evaluated
#' using hypothesis tests.
#'
#' @inheritParams anova.geer
#' @inheritParams stats::step
#' @param direction Character string specifying the direction of the stepwise
#'        search: backward elimination (\code{"backward"}), forward selection
#'        (\code{"forward"}), or both (\code{"both"}). The default is
#'        \code{"backward"}.
#' @param p_enter Numeric value between 0 and 1 giving the p-value threshold
#'        for adding a term during forward steps. Default is \code{0.15}.
#' @param p_remove Numeric value between 0 and 1 giving the p-value threshold
#'        for removing a term during backward steps. Default is \code{0.15}.
#'
#' @details
#' \code{step_p} repeatedly calls \code{\link{add1}} and \code{\link{drop1}}.
#'
#' The set of models searched is determined by the \code{scope} argument. The
#' right-hand side of its \code{lower} component is always retained in the model,
#' and the right-hand side of its \code{upper} component defines the largest model
#' that can be considered. If \code{scope} is a single formula, it specifies the
#' \code{upper} component and the \code{lower} model is empty.
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
#' Details about the testing procedures implied by the \code{test} argument can
#' be found in \cite{Rotnitzky and Jewell (1990)}. The option
#' \code{test = "working-lrt"} is only available for models fitted with an
#' independence working association structure.
#'
#' When \code{test = "wald"} or \code{test = "score"}, the \code{pmethod} argument
#' is ignored and \code{cov_type} specifies the covariance matrix estimate of the
#' regression parameter estimator used to compute the test statistic. Otherwise,
#' \code{cov_type} specifies the covariance matrix estimate used to compute the
#' coefficients of the independent chi-squared random variables, and \code{pmethod}
#' specifies the approximation used to compute the p-value.
#'
#' @return
#' A fitted model object of class \code{geer} corresponding to the final selected
#' model. The returned object includes an \code{anova} component summarizing the
#' stepwise sequence of fitted models and the terms added or removed at each step.
#'
#' @inherit add1.geer references
#'
#' @seealso \code{\link{add1}}, \code{\link{drop1}}
#'
#' @examples
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]
#'
#' full_fit <- geewa_binary(
#'   formula = status ~ (baseline + treatment + gender + visit + age)^2,
#'   id = id, repeated = visit, link = "probit", data = respiratory2,
#'   orstr = "independence", method = "pgee-jeffreys"
#' )
#'
#' ## Backward elimination within the initial model terms
#' step_p(full_fit, direction = "backward", test = "wald",
#'        cov_type = "bias-corrected", p_remove = 0.10)
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
    object <- check_geer_object(object)
    direction <- match_direction_type(direction[1L])
    step_args <- validate_step_thresholds(
      p_enter = p_enter,
      p_remove = p_remove
    )
    p_enter <- step_args$p_enter
    p_remove <- step_args$p_remove
    steps <- validate_step_count(steps)
    opts <- normalize_test_options(
      test = test[1L],
      cov_type = cov_type[1L],
      pmethod = pmethod[1L],
      object = object
    )
    test <- opts$test
    cov_type <- opts$cov_type
    pmethod <- opts$pmethod
    scope_value <- if (missing(scope)) NULL else scope
    ans <- switch(
      direction,
      backward = step_p_backward(
        object = object,
        scope = scope_value,
        test = test,
        cov_type = cov_type,
        pmethod = pmethod,
        pvalue = p_remove,
        steps = steps
      ),
      forward = step_p_forward(
        object = object,
        scope = scope_value,
        test = test,
        cov_type = cov_type,
        pmethod = pmethod,
        pvalue = p_enter,
        steps = steps
      ),
      both = step_p_both(
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
    ans
  }
