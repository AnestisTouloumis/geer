#' @title
#' Add or Drop Single Terms to or from a geer Object
#'
#' @rdname add1.geer
#' @method add1 geer
#'
#' @description
#' Fits all models that differ from a fitted \code{geer} object by adding or
#' dropping a single term in \code{scope}, and tests each change.
#'
#' @inheritParams anova.geer
#' @inheritParams stats::add1
#' @param scope a formula giving the terms to be considered for adding or
#'   dropping. Required for \code{add1()}. For \code{drop1()}, if missing,
#'   every term whose deletion respects marginality is considered.
#' @param ... not used. Supplying any argument is an error, so that a misspelled
#'   argument name is reported rather than silently ignored.
#'
#' @details
#' For \code{add1()}, \code{scope} must contain at least one term not in the
#' model, and \code{.} denotes the terms already in the model. Additions and
#' deletions respect marginality: a term is not dropped while a higher-order
#' interaction containing it remains, and an interaction is not added without
#' its lower-order terms. The tests are described in
#' \code{\link{anova.geer}}. The table also reports the correlation
#' information criterion (CIC) of each model, computed with \code{cov_type};
#' with \code{cov_type = "jackknife"} every candidate model is therefore
#' refitted once per cluster.
#'
#' When \code{object} was fitted with a data frame, every candidate model is
#' refitted to the rows of that data frame used by \code{object}, so rows
#' omitted because of missing values stay omitted even if the candidate model
#' does not use the variable concerned. For \code{add1()}, a term with missing
#' values in these rows is an error.
#'
#' @return
#' An object of class \code{c("geer_anova", "anova", "data.frame")} with one row
#' for the current model and one for each admissible addition or deletion, and
#' columns \code{Df}, \code{CIC}, \code{Chi}, and \code{Pr(>Chi)}.
#'
#' @references
#' Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
#' parameters in semiparametric generalized linear models for cluster
#' correlated data. \emph{Biometrika}, \bold{77}, 485--497.
#'
#' @seealso \code{\link{anova.geer}}, \code{\link{step_p}},
#'   \code{\link{geecriteria}}, \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", ]
#' fit <- geewa(
#'   formula = status ~ baseline + treatment + gender + visit + age,
#'   family = binomial(link = "probit"),
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   corstr = "ar1"
#' )
#' ## Single-term deletions
#' drop1(fit, test = "score")
#'
#' ## Single-term additions
#' add1(fit, scope = ~ . + baseline:age + visit:age + treatment:age,
#'      test = "score")
#'
#' @export
add1.geer <-
  function(object,
           scope,
           test = c("wald", "score", "working-wald", "working-score",
                    "working-lrt"),
           cov_type = "bias-corrected",
           pmethod = c("rao-scott", "satterthwaite"),
           ...) {
    check_unused_dots(list(...), "add1.geer")
    object <- check_geer_object(object)
    opts <- normalize_geer_test_options(
      test = test[1L],
      cov_type = cov_type,
      pmethod = pmethod[1L],
      object = object
    )
    test <- opts$test
    cov_type <- opts$cov_type
    pmethod <- opts$pmethod
    if (missing(scope) || is.null(scope)) {
      stop("'scope' contains no terms to add", call. = FALSE)
    }
    if (!is.character(scope)) {
      scope <- stats::add.scope(object, stats::update.formula(object, scope))
    }
    if (!length(scope)) {
      stop("'scope' contains no terms to add", call. = FALSE)
    }
    ns <- length(scope)
    ans <- matrix(
      NA_real_,
      nrow = ns + 1L,
      ncol = 4L,
      dimnames = list(
        c("<none>", scope),
        c("Df", "CIC", "Chi", "Pr(>Chi)")
      )
    )

    ans[1L, 2L] <- compute_cic(object, cov_type)
    for (i in seq_len(ns)) {
      tt <- scope[[i]]
      add1_model <- refit_geer(
        object,
        stats::as.formula(paste(". ~ . +", tt)),
        "add1.geer"
      )
      value <- compute_nested_test(
        object,
        add1_model,
        test,
        cov_type,
        pmethod,
        "add1.geer"
      )
      ans[i + 1L, ] <- c(
        value$test_df,
        compute_cic(add1_model, cov_type),
        value$test_stat,
        value$test_p
      )
    }
    aod <- as.data.frame(ans)
    test_type <- build_test_label(test)
    formula_txt <- paste(deparse(object$formula), collapse = " ")
    head <- c(
      paste("Single term additions using", test_type, "test:"),
      "\nModel:", formula_txt
    )
    structure(
      aod,
      heading = head,
      class = c("geer_anova", "anova", "data.frame")
    )
  }


#' @rdname add1.geer
#' @method drop1 geer
#'
#' @export
drop1.geer <- function(object,
                       scope,
                       test = c("wald", "score", "working-wald",
                                "working-score", "working-lrt"),
                       cov_type = "bias-corrected",
                       pmethod = c("rao-scott", "satterthwaite"),
                       ...) {
  check_unused_dots(list(...), "drop1.geer")
  object <- check_geer_object(object)
  opts <- normalize_geer_test_options(
    test = test[1L],
    cov_type = cov_type,
    pmethod = pmethod[1L],
    object = object
  )
  test <- opts$test
  cov_type <- opts$cov_type
  pmethod <- opts$pmethod
  model_terms <- attr(stats::terms(object), "term.labels")
  admissible_scope <- stats::drop.scope(object)
  if (missing(scope) || is.null(scope)) {
    scope <- admissible_scope
  } else {
    if (!is.character(scope)) {
      scope <- attr(
        stats::terms(stats::update.formula(object, scope)),
        "term.labels"
      )
    }
    if (!all(match(scope, model_terms, 0L) > 0L)) {
      stop("'scope' must contain only terms of the model", call. = FALSE)
    }
    scope <- intersect(scope, admissible_scope)
  }
  ns <- length(scope)
  if (ns == 0L) {
    stop("'scope' contains no terms that can be dropped", call. = FALSE)
  }
  ans <- matrix(
    NA_real_,
    nrow = ns + 1L,
    ncol = 4L,
    dimnames = list(
      c("<none>", scope),
      c("Df", "CIC", "Chi", "Pr(>Chi)")
    )
  )
  ans[1L, 2L] <- compute_cic(object, cov_type)
  for (i in seq_len(ns)) {
    tt <- scope[[i]]

    drop1_model <- refit_geer(
      object,
      stats::as.formula(paste(". ~ . -", tt)),
      "drop1.geer"
    )
    value <- compute_nested_test(
      drop1_model,
      object,
      test,
      cov_type,
      pmethod,
      "drop1.geer"
    )
    ans[i + 1L, ] <- c(
      value$test_df,
      compute_cic(drop1_model, cov_type),
      value$test_stat,
      value$test_p
    )
  }
  aod <- as.data.frame(ans)
  test_type <- build_test_label(test)
  formula_txt <- paste(deparse(object$formula), collapse = " ")
  head <- c(
    paste("Single term deletions using", test_type, "test:"),
    "\nModel:", formula_txt
  )
  structure(aod, heading = head, class = c("geer_anova", "anova", "data.frame"))
}


#' @title
#' ANOVA Tables for geer Objects
#'
#' @method anova geer
#'
#' @description
#' Computes tables of hypothesis tests for one or more fitted \code{geer}
#' objects, analogous to analysis of deviance tables for generalized linear
#' models.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param ... additional fitted model objects of class \code{"geer"} to be
#'   included in the comparison, supplied without names. Any other argument
#'   is an error.
#' @param test character string specifying the hypothesis test: Wald
#'   (\code{"wald"}), generalized score (\code{"score"}), modified working Wald
#'   (\code{"working-wald"}), modified working score (\code{"working-score"}),
#'   or modified working likelihood ratio (\code{"working-lrt"}). Defaults to
#'   \code{"wald"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for inference: one of the estimators described in
#'   \code{\link{vcov.geer}}. Defaults to \code{"bias-corrected"}.
#' @param pmethod character string specifying the approximation to the null
#'   distribution of the modified working tests: Rao-Scott (\code{"rao-scott"})
#'   or Satterthwaite (\code{"satterthwaite"}). Ignored for the Wald and score
#'   tests, where it is not validated. Defaults to \code{"rao-scott"}.
#'
#' @details
#' The tests are those of Rotnitzky and Jewell (1990). The Wald and score
#' statistics use the covariance matrix selected by \code{cov_type}. The
#' modified working tests compare a working statistic with a weighted sum of
#' independent chi-squared variables, whose weights depend on \code{cov_type},
#' and \code{pmethod} selects the approximation to its distribution. The
#' modified working likelihood ratio test requires a working independence
#' structure, and both working log-likelihoods are scaled by the dispersion of
#' the larger model. The score tests are evaluated at the fitted means of the
#' smaller model with the association and dispersion parameters of the larger
#' model; with \code{cov_type = "jackknife"}, the covariance component is the
#' jackknife covariance of the larger model. The small-sample estimators of
#' \code{\link{vcov.geer}}, such as \code{"kauermann-carroll"}, are evaluated
#' at the same fitted means, which extends their published definitions.
#'
#' Nesting requires the smaller model matrix and the columns of the larger
#' model matrix with the same coefficient names to span the same space.
#' Consequently, setting the remaining coefficients to zero recovers the
#' smaller mean model. Matching coefficient names alone are insufficient.
#' Rescaling or changing the basis of the shared columns is allowed when it
#' preserves this space.
#'
#' With several models, the fits must use the same data and each consecutive
#' pair must be nested. Different responses or working association structures
#' are errors. The entire comparison fails; no table from a reduced set of
#' models is returned. Other differences in fitting settings,
#' such as the estimation method or handling of dispersion, are also errors.
#'
#' With one model, the sequential models are refitted to the rows used by the
#' fit, as described in \code{\link{add1.geer}}.
#'
#' @return
#' An object of class \code{c("geer_anova", "anova", "data.frame")} with columns
#' \code{Df}, \code{Resid. Df}, \code{Chi}, and \code{Pr(>Chi)}. With one model,
#' the terms are tested sequentially, from first to last; with several models,
#' each model is compared with the previous one. \code{Resid. Df} is the number
#' of observations minus the number of regression parameters.
#'
#' @inherit add1.geer references
#'
#' @seealso \code{\link{add1.geer}}, \code{\link{drop1.geer}},
#'   \code{\link{step_p}}, \code{\link{geecriteria}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#' fit <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   repeated = period,
#'   orstr = "exchangeable"
#' )
#' ## Sequential tests of the terms
#' anova(fit)
#'
#' ## Comparison of nested models
#' fit_null <- update(fit, formula = . ~ . - treatment)
#' anova(fit_null, fit, test = "score")
#'
#' @export
anova.geer <-
  function(object,
           ...,
           test = c("wald", "score", "working-wald", "working-score",
                    "working-lrt"),
           cov_type = "bias-corrected",
           pmethod = c("rao-scott", "satterthwaite")) {
    object <- check_geer_object(object)
    opts <- normalize_geer_test_options(
      test = test[1L],
      cov_type = cov_type,
      pmethod = pmethod[1L],
      object = object
    )
    test <- opts$test
    cov_type <- opts$cov_type
    pmethod <- opts$pmethod
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) {
      rep_len(FALSE, length(dotargs))
    } else {
      (names(dotargs) != "")
    }
    is_geer <- vapply(dotargs, function(x) inherits(x, "geer"), logical(1))
    check_unused_dots(dotargs[named | !is_geer], "anova.geer")
    if (length(dotargs)) {
      dotargs <- lapply(dotargs, check_geer_object)
      return(compute_anova_geer_list(
        c(list(object), dotargs),
        test = test,
        cov_type = cov_type,
        pmethod = pmethod
      ))
    }
    terms <- attr(object$terms, "term.labels")
    intercept <- attr(object$terms, "intercept")
    variables <- attr(object$terms, "variables")
    varseq <- attr(object$x, "assign")
    nvars <- max(c(0, varseq))
    object_list <- list()
    offset_terms <- extract_formula_offset_terms(object)
    offset_part <- if (length(offset_terms)) {
      paste0(" + ", paste(offset_terms, collapse = " + "))
    } else {
      ""
    }
    if (intercept == 1L) {
      object_list[[1L]] <- refit_geer(
        object,
        stats::as.formula(paste0(". ~ 1", offset_part)),
        "anova.geer"
      )
      for (i in seq_len(nvars)) {
        object_list[[i + 1]] <- refit_geer(
          object_list[[i]],
          stats::as.formula(paste(". ~ . + ", terms[i])),
          "anova.geer"
        )
      }
    } else {
      object_list[[1L]] <- refit_geer(
        object,
        stats::as.formula(paste0(". ~ -1 + ", terms[1L], offset_part)),
        "anova.geer"
      )
      for (i in seq_len(nvars - 1L)) {
        object_list[[i + 1]] <- refit_geer(
          object_list[[i]],
          stats::as.formula(paste(". ~ . + ", terms[i + 1])),
          "anova.geer"
        )
      }
    }
    resdf <- vapply(
      object_list,
      function(x) as.numeric(x$df.residual),
      numeric(1)
    )
    table <- data.frame(
      Df = c(NA_real_, rep(NA_real_, length(resdf) - 1L)),
      `Resid. Df` = resdf,
      Chi = c(NA_real_, rep(NA_real_, length(resdf) - 1L)),
      `Pr(>Chi)` = c(NA_real_, rep(NA_real_, length(resdf) - 1L)),
      check.names = FALSE
    )
    if (intercept == 1L) {
      dimnames(table) <- list(
        c("NULL", terms),
        c("Df", "Resid. Df", "Chi", "Pr(>Chi)")
      )
    } else {
      dimnames(table) <- list(
        c(terms),
        c("Df", "Resid. Df", "Chi", "Pr(>Chi)")
      )
    }
    for (i in seq_len(length(object_list) - 1L)) {
      value <- compute_nested_test(
        object_list[[i]],
        object_list[[i + 1]],
        test,
        cov_type,
        pmethod,
        "anova.geer"
      )
      table[i + 1, -2] <- c(value$test_df, value$test_stat, value$test_p)
    }
    test_type <- build_test_label(test)
    title <- paste(
      "Analysis of ", test_type, " Statistic Table",
      "\n\nModel: ", object$family$family,
      ", link: ", object$family$link,
      "\n\nResponse: ", as.character(variables[-1L])[1L],
      "\n\nTerms added sequentially (first to last)\n\n",
      sep = ""
    )
    structure(
      table,
      heading = title,
      class = c("geer_anova", "anova", "data.frame")
    )
  }
