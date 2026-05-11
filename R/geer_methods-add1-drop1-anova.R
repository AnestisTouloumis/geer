#' @title
#' Add or Drop Single Terms to or from a geer Model
#'
#' @rdname add1.geer
#' @aliases add1 add1.geer
#' @method add1 geer
#'
#' @description
#' Computes all single terms in the \code{scope} argument that can be added to
#' or dropped from a fitted \code{geer} model, fits the corresponding models,
#' and returns a table summarizing the resulting changes in fit.
#'
#' @inheritParams anova.geer
#' @inheritParams stats::add1
#' @param ... additional arguments passed to or from other methods.
#'
#' @details
#' For \code{add1.geer()}, \code{scope} must specify the candidate terms to be
#' added. If no eligible terms are supplied, an error is returned. In scope
#' formulas, \code{.} denotes the set of terms already included in the model.
#'
#' Model hierarchy is enforced: candidate additions and deletions must preserve
#' marginality, so if a higher-order interaction is present, all of its
#' lower-order component terms must also remain in the model.
#'
#' Details of the hypothesis tests controlled by \code{test} are given in
#' Rotnitzky and Jewell (1990). The option \code{test = "working-lrt"} is
#' valid only when the model is fitted with an independence working association
#' structure; otherwise an error is returned.
#'
#' When \code{test \%in\% c("wald", "score")}, the \code{pmethod} argument is
#' ignored and \code{cov_type} specifies the covariance estimator used to
#' compute the test statistic. For modified working tests, \code{cov_type}
#' determines the covariance matrix used to form the coefficients of the sum of
#' independent chi-squared random variables, and \code{pmethod} specifies the
#' approximation used to compute the p-value.
#'
#' The output table also includes the Correlation Information Criterion (CIC)
#' for each candidate model, which can be used to guide selection of the
#' working association structure.
#'
#' @return
#' An object of class \code{"anova"} summarizing the differences in fit
#' between the models. For \code{add1.geer()} and \code{drop1.geer()}, the
#' table contains one row for the current model and one row for each admissible
#' single-term addition or deletion. Columns include \code{Df} (degrees of
#' freedom of the test), \code{CIC} (Correlation Information Criterion),
#' \code{Chi} (test statistic), and \code{Pr(>Chi)} (p-value).
#'
#' @references
#' Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
#' parameters in semiparametric generalized linear models for cluster
#' correlated data. \emph{Biometrika}, \bold{77}, 485--497.
#'
#' @seealso \code{\link{anova.geer}}, \code{\link{step_p}},
#'   \code{\link{geecriteria}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]
#'
#' fitted_model <- geewa(
#'   formula = status ~ baseline + I(treatment == "active") + gender + visit + age,
#'   family = binomial(link = "probit"),
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   corstr = "ar1",
#'   method = "gee"
#' )
#' add1(
#'   fitted_model,
#'   scope = . ~ . + baseline:age + age:visit + I(treatment == "active"):age + age:gender,
#'   test = "score"
#' )
#'
#' @export
add1.geer <-
  function(object,
           scope,
           test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
           cov_type = c("bias-corrected", "robust", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite"),
           ...) {
    object <- check_geer_object(object)
    opts <- normalize_geer_test_options(
      test = test[1L],
      cov_type = cov_type[1L],
      pmethod = pmethod[1L],
      object = object
    )
    test <- opts$test
    cov_type <- opts$cov_type
    pmethod <- opts$pmethod
    if (missing(scope) || is.null(scope)) {
      stop("no terms in scope for adding to object", call. = FALSE)
    }
    if (!is.character(scope)) {
      scope <- add.scope(object, update.formula(object, scope))
    }
    if (!length(scope)) {
      stop("no terms in scope for adding to object", call. = FALSE)
    }
    ns <- length(scope)
    ans <- matrix(NA_real_, nrow = ns + 1L, ncol = 4L,
                  dimnames = list(c("<none>", scope),
                                  c("Df", "CIC", "Chi", "Pr(>Chi)")))

    ans[1L, 2L] <- compute_gee_cic(object, cov_type)
    for (i in seq_len(ns)) {
      tt <- scope[[i]]
      add1_model <- update(
        object,
        formula = as.formula(paste(". ~ . +", tt)),
        data = object$data
      )
      add1_model <- restore_original_data_call(add1_model, object)
      value <- switch(
        test,
        wald = wald_test(object, add1_model, cov_type),
        score = score_test(object, add1_model, cov_type),
        `working-wald`  = working_wald_test(object, add1_model, cov_type, pmethod),
        `working-score` = working_score_test(object, add1_model, cov_type, pmethod),
        `working-lrt`   = working_lrt_test(object, add1_model, cov_type, pmethod)
      )
      ans[i + 1L, ] <- c(value$test_df,
                         compute_gee_cic(add1_model, cov_type),
                         value$test_stat,
                         value$test_p)
    }
    aod <- as.data.frame(ans)
    test_type <- format_test_label(test)
    formula_txt <- paste(deparse(object$call$formula), collapse = " ")
    head <- c(
      paste("Single term additions using", test_type, "test:"),
      "\nModel:", formula_txt
    )
    structure(aod, heading = head, class = c("anova", "data.frame"))
  }


#' @rdname add1.geer
#' @aliases drop1 drop1.geer
#' @method drop1 geer
#'
#' @examples
#' data("respiratory", package = "geer")
#' respiratory2 <- respiratory[respiratory$center == "C2", , drop = FALSE]
#'
#' fitted_model <- geewa(
#'   formula = status ~ baseline + I(treatment == "active") + gender + visit + age,
#'   family = binomial(link = "probit"),
#'   data = respiratory2,
#'   id = id,
#'   repeated = visit,
#'   corstr = "ar1",
#'   method = "gee"
#' )
#' drop1(fitted_model, test = "score")
#'
#' @export
drop1.geer <- function(object,
                       scope,
                       test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
                       cov_type = c("bias-corrected", "robust", "df-adjusted", "naive"),
                       pmethod = c("rao-scott", "satterthwaite"),
                       ...) {
  object <- check_geer_object(object)
  opts <- normalize_geer_test_options(
    test = test[1L],
    cov_type = cov_type[1L],
    pmethod = pmethod[1L],
    object = object
  )
  test <- opts$test
  cov_type <- opts$cov_type
  pmethod <- opts$pmethod
  model_terms <- attr(terms(object), "term.labels")
  if (missing(scope) || is.null(scope)) {
    scope <- drop.scope(object)
  } else {
    if (!is.character(scope)) {
      scope <- attr(terms(update.formula(object, scope)), "term.labels")
    }
    if (!all(match(scope, model_terms, 0L) > 0L)) {
      stop("scope is not a subset of term labels", call. = FALSE)
    }
  }
  ns <- length(scope)
  if (ns == 0L) {
    stop("no terms in scope for dropping from object", call. = FALSE)
  }
  ans <- matrix(NA_real_, nrow = ns + 1L, ncol = 4L,
                dimnames = list(c("<none>", scope),
                                c("Df", "CIC", "Chi", "Pr(>Chi)")))
  ans[1L, 2L] <- compute_gee_cic(object, cov_type)
  for (i in seq_len(ns)) {
    tt <- scope[[i]]
    drop1_model <- update(
      object,
      formula = as.formula(paste(". ~ . -", tt)),
      data = object$data
    )
    drop1_model <- restore_original_data_call(drop1_model, object)
    value <- switch(
      test,
      wald = wald_test(drop1_model, object, cov_type),
      score = score_test(drop1_model, object, cov_type),
      `working-wald` = working_wald_test(drop1_model, object, cov_type, pmethod),
      `working-score` = working_score_test(drop1_model, object, cov_type, pmethod),
      `working-lrt` = working_lrt_test(drop1_model, object, cov_type, pmethod)
    )
    ans[i + 1L, ] <- c(value$test_df,
                       compute_gee_cic(drop1_model, cov_type),
                       value$test_stat,
                       value$test_p)
  }
  aod <- as.data.frame(ans)
  test_type <- format_test_label(test)
  formula_txt <- paste(deparse(object$call$formula), collapse = " ")
  head <- c(
    paste("Single term deletions using", test_type, "test:"),
    "\nModel:", formula_txt
  )
  structure(aod, heading = head, class = c("anova", "data.frame"))
}


#' @title
#' ANOVA Tables for geer Objects
#'
#' @aliases anova anova.geer
#' @method anova geer
#'
#' @description
#' Computes hypothesis test tables for one or more fitted \code{geer} objects,
#' analogous to analysis-of-variance (ANOVA) tables for GLM models.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param ... additional fitted model objects of class \code{"geer"}, used for
#'   sequential model comparisons.
#' @param test character string specifying the hypothesis testing procedure.
#'   Options are the Wald test (\code{"wald"}), the generalized score test
#'   (\code{"score"}), the modified working Wald test (\code{"working-wald"}),
#'   the modified working score test (\code{"working-score"}), and the modified
#'   working likelihood ratio test (\code{"working-lrt"}). Defaults to
#'   \code{"wald"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for inference on the regression parameters. Options are the bias-corrected
#'   estimator (\code{"bias-corrected"}), the sandwich or robust estimator
#'   (\code{"robust"}), the degrees-of-freedom adjusted estimator
#'   (\code{"df-adjusted"}), and the model-based or naive estimator
#'   (\code{"naive"}). Defaults to \code{"bias-corrected"}.
#' @param pmethod character string specifying the approximation used to compute
#'   the p-value for the modified working tests. Options are the Rao--Scott
#'   approximation (\code{"rao-scott"}) and the Satterthwaite approximation
#'   (\code{"satterthwaite"}). Defaults to \code{"rao-scott"}.
#'
#' @details
#' Details of the hypothesis tests controlled by \code{test} are given in
#' Rotnitzky and Jewell (1990). The option \code{test = "working-lrt"} is
#' valid only when the model is fitted with an independence working association
#' structure; otherwise an error is returned.
#'
#' When \code{test \%in\% c("wald", "score")}, the \code{pmethod} argument is
#' ignored and \code{cov_type} specifies the covariance estimator used to
#' compute the test statistic. For modified working tests, \code{cov_type}
#' determines the covariance matrix used to form the coefficients of the sum of
#' independent chi-squared random variables, and \code{pmethod} specifies the
#' approximation used to compute the p-value.
#'
#' When comparing two or more models, the data must be identical across all
#' fits, and the models must be nested in the order supplied. In particular,
#' each consecutive pair of models must be nested.
#'
#' @return
#' An object of class \code{c("anova", "data.frame")}. With a single model,
#' the table reports tests for terms added sequentially (first to last). With
#' multiple models, the table reports sequential comparisons between each
#' consecutive pair of nested models. Columns include \code{Df} (degrees of
#' freedom of the test), \code{Resid. Df} (residual degrees of freedom),
#' \code{Chi} (test statistic), and \code{Pr(>Chi)} (p-value). Unlike
#' \code{\link{add1.geer}} and \code{\link{drop1.geer}}, no \code{CIC} column
#' is included; use \code{\link{geecriteria}} for working association structure
#' selection.
#'
#' @inherit add1.geer references
#'
#' @seealso \code{\link{add1.geer}}, \code{\link{drop1.geer}} for type II
#'   tests where each term is dropped one at a time while respecting model
#'   hierarchy; \code{\link{step_p}} for stepwise model selection;
#'   \code{\link{geecriteria}} for model comparison criteria.
#'
#' @examples
#' data("cerebrovascular", package = "geer")
#'
#' ## Single-model ANOVA (sequential terms)
#' fit_full <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' anova(fit_full, test = "wald", cov_type = "robust")
#'
#' ## Two-model comparison (models must be nested)
#' fit_null <- geewa_binary(
#'   formula = ecg ~ 1,
#'   link = "logit",
#'   data = cerebrovascular,
#'   id = id,
#'   orstr = "exchangeable"
#' )
#' anova(fit_null, fit_full, test = "wald", cov_type = "robust")
#'
#' @export
anova.geer <-
  function(object,
           ...,
           test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
           cov_type = c("bias-corrected", "robust", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite")) {
    object <- check_geer_object(object)
    opts <- normalize_geer_test_options(
      test = test[1L],
      cov_type = cov_type[1L],
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
    if (any(named)) {
      warning("named arguments in '...' are ignored", call. = FALSE)
    }
    dotargs <- dotargs[!named]
    is_geer <- vapply(dotargs, function(x) inherits(x, "geer"), logical(1))
    dotargs <- dotargs[is_geer]
    if (length(dotargs)) {
      dotargs <- lapply(dotargs, check_geer_object)
      return(compute_anova_geer_list(c(list(object), dotargs),
                                     test = test,
                                     cov_type = cov_type,
                                     pmethod = pmethod))
    }
    terms <- attr(object$terms, "term.labels")
    intercept <- attr(object$terms, "intercept")
    variables <- attr(object$terms, "variables")
    varseq <- attr(object$x, "assign")
    nvars <- max(c(0, varseq))
    object_list <- list()
    if (intercept == 1) {
      object_list[[1]] <- update(object, formula = . ~ 1, data = object$data)
      object_list[[1]] <- restore_original_data_call(object_list[[1]], object)
      for (i in seq_len(nvars)) {
        object_list[[i + 1]] <- update(
          object_list[[i]],
          formula = paste(". ~ . + ", terms[i]),
          data = object$data
        )
        object_list[[i + 1]] <- restore_original_data_call(object_list[[i + 1]], object)
      }
    } else {
      object_list[[1]] <- update(
        object,
        formula = paste(". ~ -1 + ", terms[1]),
        data = object$data
      )
      object_list[[1]] <- restore_original_data_call(object_list[[1]], object)
      for (i in seq_len(nvars - 1)) {
        object_list[[i + 1]] <- update(
          object_list[[i]],
          formula = paste(". ~ . + ", terms[i + 1]),
          data = object$data
        )
        object_list[[i + 1]] <- restore_original_data_call(object_list[[i + 1]], object)
      }
    }
    resdf <- vapply(object_list, function(x) as.numeric(x$df.residual), numeric(1))
    table <- data.frame(
      Df = c(NA_real_, rep(NA_real_, length(resdf) - 1L)),
      `Resid. Df` = resdf,
      Chi = c(NA_real_, rep(NA_real_, length(resdf) - 1L)),
      `Pr(>Chi)` = c(NA_real_, rep(NA_real_, length(resdf) - 1L)),
      check.names = FALSE
    )
    if (intercept == 1) {
      dimnames(table) <- list(c("NULL", terms),
                              c("Df", "Resid. Df", "Chi", "Pr(>Chi)"))
    } else {
      dimnames(table) <- list(c(terms),
                              c("Df", "Resid. Df", "Chi", "Pr(>Chi)"))
    }
    for (i in seq_len(length(object_list) - 1)) {
      value <- switch(
        test,
        wald = wald_test(object_list[[i]], object_list[[i + 1]], cov_type),
        score = score_test(object_list[[i]], object_list[[i + 1]], cov_type),
        `working-wald` =
          working_wald_test(object_list[[i]], object_list[[i + 1]], cov_type, pmethod),
        `working-score` =
          working_score_test(object_list[[i]], object_list[[i + 1]], cov_type, pmethod),
        `working-lrt` =
          working_lrt_test(object_list[[i]], object_list[[i + 1]], cov_type, pmethod)
      )
      table[i + 1, -2] <- c(value$test_df, value$test_stat, value$test_p)
    }
    test_type <- format_test_label(test)
    title <- paste(
      "Analysis of ", test_type, " Statistic Table",
      "\n\nModel: ", object$family$family,
      ", link: ", object$family$link,
      "\n\nResponse: ", as.character(variables[-1L])[1L],
      "\n\nTerms added sequentially (first to last)\n\n",
      sep = ""
    )
    structure(table, heading = title, class = c("anova", "data.frame"))
  }
