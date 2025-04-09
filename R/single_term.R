#' @title
#' Add or Drop All Possible Single Terms to a Model from a \code{geer} Object
#'
#' @method add1 geer
#' @aliases add1 add1.geer
#'
#'
#' @description
#' Compute all the single terms in the \code{scope} argument that can be added
#' to or dropped from the model, fit those models and compute a table of the
#' changes in fit.
#'
#'
#' @inheritParams anova.geer
#' @inheritParams stats::add1
#'
#'
#' @details
#' For \code{drop1}, a missing scope is taken to be all terms in the model. The
#' hierarchy is respected when considering terms to be added or dropped:
#' all main effects contained in a second-order interaction must remain, and so on.
#'
#' In a scope formula . means ‘what is already there’.
#'
#' The output table also gives CIC.
#'
#'
#' @return
#' An object of class \code{anova} summarizing the differences in fit between the models.
#'
#'
#' @author Anestis Touloumis
#'
#'
#' @seealso \code{\link{anova}}
#'
#'
#' @examples
#' fitted_model <-
#' geewa(formula = y ~ baseline + treatment + gender + visit + age,
#'       id = id, repeated = visit, family = binomial(link = "probit"),
#'       data = respiratory[respiratory$center==2, ], corstr = "ar1",
#'       method = "gee")
#' add1(fitted_model,
#'      scope = .~. + baseline:age + age:visit + treatment:age + age:gender,
#'      test = "score")
#'
#' @rdname add1.geer
#' @export
add1.geer <- function(object,
                      scope,
                      test = "wald",
                      cov_type = "robust",
                      pmethod = "rao-scott",
                      ...) {
  icheck <- pmatch(test,
                   c("wald", "score", "working-wald", "working-score", "working-lrt"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown testing procedure")
  if (test %in% c("working-wald", "working-score", "working-lrt")) {
    icheck <- pmatch(pmethod,
                     c("rao-scott", "satterthwaite"),
                     nomatch = 0,
                     duplicates.ok = FALSE)
    if (icheck == 0)
      stop("unknown method for pvalue")
  }
  if (test == "working-lrt" & object$association_structure != "independence")
    stop("the modified working lrt can only be applied to an independence working model")
  icheck <- pmatch(cov_type,
                   c("robust", "df-adjusted", "bias-corrected", "naive"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown covariance matrix")
  if (missing(scope) || is.null(scope))
    stop("no terms in scope")
  if (!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope))
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 4L,
                dimnames = list(c("<none>", scope), c("Df", "CIC", "Chi", "Pr(>Chi)")))
  ans[1L, 2] <- extract_cic(object, cov_type)
  for (i in seq_along(scope)) {
    tt <- scope[i]
    add1_model <- update(object, formula = as.formula(paste("~ . +", tt)))
    value <-
      switch(test,
             wald = wald_test(object, add1_model, cov_type),
             score = score_test(object, add1_model, cov_type),
             `working-wald` =
               working_wald_test(object, add1_model, cov_type, pmethod),
             `working-score` =
               working_score_test(object, add1_model, cov_type, pmethod),
             `working-lrt` =
               working_lr_test(object, add1_model, cov_type, pmethod)
             )
    ans[i + 1, ] <-
      c(value$test_df, extract_cic(add1_model, cov_type), value$test_stat, value$test_p)
  }
  aod <- as.data.frame(ans)
  test_type <- switch(test,
                      wald = "Wald",
                      score = "Score",
                      `working-wald` = "Modified Working Wald",
                      `working-score` = "Modified Working Score",
                      `working-lrt` = "Modified Working LR")
  head <- c(paste("Single term additions using", test_type, "test:"),
            "\nModel:", deparse(formula(object$call$formula)))
  structure(aod,
            heading = head,
            class = c("anova", "data.frame"))
}

#' @method drop1 geer
#' @aliases drop1 drop1.geer
#'
#'
#' @examples
#' drop1(fitted_model, test = "score")
#'
#'
#' @rdname add1.geer
#' @export
drop1.geer <- function(object,
                       scope,
                       test = "wald",
                       cov_type = "robust",
                       pmethod = "rao-scott",
                       ...) {
  icheck <- pmatch(test,
                   c("wald", "score", "working-wald", "working-score", "working-lrt"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown testing procedure")
  if (test %in% c("working-wald", "working-score", "working-lrt")) {
    icheck <- pmatch(pmethod,
                     c("rao-scott", "satterthwaite"),
                     nomatch = 0,
                     duplicates.ok = FALSE)
    if (icheck == 0)
      stop("unknown method for pvalue")
  }
  if (test == "working-lrt" & object$association_structure != "independence")
    stop("the modified working lrt can only be applied to the independence working model")
  icheck <- pmatch(cov_type,
                   c("robust", "df-adjusted", "bias-corrected", "naive"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown covariance matrix")
  model_terms <- attr(terms(object), "term.labels")
  if (missing(scope)) {
    scope <- drop.scope(object)
  } else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)), "term.labels")
    if (!all(match(scope, model_terms, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 4L,
                dimnames = list(c("<none>", scope),
                                c("Df", "CIC", "Chi", "Pr(>Chi)")))
  ans[1L, 2] <- extract_cic(object, cov_type)
  for (i in seq_along(scope)) {
    tt <- scope[i]
    drop1_model <-
      update(object,
             formula = as.formula(paste("~ . -", tt)))
    value <-
      switch(test,
             wald = wald_test(drop1_model, object, cov_type),
             score = score_test(drop1_model, object, cov_type),
             `working-wald` =
               working_wald_test(drop1_model, object, cov_type, pmethod),
             `working-score` =
               working_score_test(drop1_model, object, cov_type, pmethod),
             `working-lrt` =
               working_lr_test(drop1_model, object, cov_type, pmethod))
    ans[i + 1, ] <-
      c(value$test_df, extract_cic(drop1_model, cov_type), value$test_stat, value$test_p)
  }
  aod <- as.data.frame(ans)
  test_type <- switch(test,
                      wald = "Wald",
                      score = "Score",
                      `working-wald` = "Modified Working Wald",
                      `working-score` = "Modified Working Score",
                      `working-lrt` = "Modified Working LR")
  head <- c(paste("Single term deletions using", test_type, "test:"),
            "\nModel:", deparse(formula(object$call$formula)))
  structure(aod,
            heading = head,
            class = c("anova", "data.frame"))
}
