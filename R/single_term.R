#' Add or Drop All Possible Single Terms to a Model
#'
#' Compute all the single terms in the scope argument that can be added to or
#' dropped from the model, fit those models and compute a table of the changes
#' in fit.
#'
#' @param object a model fit of class \code{"geer"}.
#' @param scope a formula giving the terms to be considered for adding or
#'        dropping.
#' @param test character indicating the hypothesis testing procedure applied.
#' Options include \code{wald} or \code{score} test.  By
#' default, the wald test is used.
#' @param cov_type character indicating whether the sandwich (robust)
#' covariance
#' matrix (\code{cov_type = "robust"}), the model-based (naive) covariance
#' matrix (\code{cov_type = "naive"}), the bias-corrected covariance
#' matrix (\code{cov_type = "bias-corrected"}) or the degrees of freedom adjusted
#' covariance matrix (\code{cov_type = "df-adjusted"}) should be used. By
#' default, the robust covariance matrix is used.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' For \code{drop1}, a missing scope is taken to be all terms in the model. The
#' hierarchy is respected when considering terms to be added or dropped:
#' all main effects contained in a second-order interaction must remain, and so on.
#'
#' In a scope formula . means ‘what is already there’.
#'
#' @author Anestis Touloumis
#'
#' @export
#'
#' @rdname add1.geer
#' @method add1 geer
#' @aliases add1 add1.geer
#'
#' @examples
#' fitted_model <-
#' geewa(formula = y ~ baseline + treatment + gender + visit + age,
#'       id = id,
#'       repeated = visit,
#'       family = binomial(link = "probit"),
#'       data = respiratory[respiratory$center==2, ],
#'       corstr = "ar1",
#'       method = "gee")
#' add1(fitted_model,
#'      scope = .~. + baseline:age + age:visit + treatment:age + age:gender,
#'      test = "score")
add1.geer <- function(object, scope, test = "wald",
                      cov_type = "robust", ...) {
  model_terms <- attr(terms(object), "term.labels")
  if (missing(scope) || is.null(scope)) stop("no terms in scope")
  if (!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope))
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns, ncol = 3,
                dimnames =  list(scope, c("Df", "X2", "Pr(>Chi)")))
  ns <- length(scope)
  icheck <- pmatch(test,
                   c("wald", "score"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown test")
  ans <- matrix(nrow = ns, ncol = 3,
                dimnames =  list(scope, c("Df", "X2", "Pr(>Chi)")))
  test_proc <- ifelse(test == "wald", wald_test, score_test)
  for (i in seq_along(scope)) {
    tt <- scope[i]
    reduced_model <-
      update(object,
             formula = as.formula(paste("~ . +", tt)))
    ans[i,] <- unlist(test_proc(reduced_model, object, cov_type))
  }
  aod <- as.data.frame(ans)
  head <- c(paste("Single term additions using", test, "test:"),
            "\nModel:", deparse(formula(object$call$formula)),"\n")
  structure(aod,
            heading = head,
            class = c("anova", "data.frame"))
}

#' @rdname add1.geer
#' @method drop1 geer
#' @aliases drop1 drop1.geer
#' @export
#' @examples
#' drop1(fitted_model,
#'       test = "score")
#'
drop1.geer <- function(object,
                       scope,
                       test = "wald",
                       cov_type = "robust", ...) {
  model_terms <- attr(terms(object), "term.labels")
  if (missing(scope)) {
    scope <- drop.scope(object) } else {
      if (!is.character(scope)) {
        scope <- attr(terms(update.formula(object, scope)),
                      "term.labels")
      }
      if (!all(match(scope, model_terms, 0L) > 0L))
        stop("scope is not a subset of term labels")
    }
  ns <- length(scope)
  icheck <- pmatch(test,
                   c("wald", "score"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown test")
  ans <- matrix(nrow = ns, ncol = 3,
                dimnames =  list(scope, c("Df", "X2", "Pr(>Chi)")))
  test_proc <- ifelse(test == "wald", wald_test, score_test)
  for(i in seq_along(scope)) {
    tt <- scope[i]
    reduced_model <-
      update(object,
             formula = as.formula(paste("~ . -", tt)))
    ans[i,] <- unlist(test_proc(reduced_model, object, cov_type))
  }
  aod <- as.data.frame(ans)
  head <- c(paste("Single term deletions using", test, "test:"),
               "\nModel:", deparse(formula(object$call$formula)),"\n")
  structure(aod,
            heading = head,
            class = c("anova", "data.frame"))
}
