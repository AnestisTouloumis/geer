#' @title
#' Add or Drop Single Terms to/from a Model from a \code{geer} Object
#'
#' @rdname add1.geer
#' @aliases add1 add1.geer
#' @method add1 geer
#'
#' @inherit stats::add1.default description
#'
#' @inheritParams anova.geer
#' @inheritParams stats::add1
#' @param ... additional argument(s) passed to methods.
#'
#' @details
#' If \code{scope} is missing, all terms in the current model are considered.
#' Model hierarchy is enforced: if a higher-order interaction is present, all
#' of its lower-order main effects must also remain in the model. In scope
#' formulas, \code{.} denotes the set of terms already included in the model.
#'
#' Details of the hypothesis tests controlled by \code{test} are given in
#' \cite{Rotnitzky and Jewell (1990)}. The option \code{test = "working-lrt"}
#' is valid only when the model is fitted under an \emph{independence} working
#' association; otherwise, an error is returned.
#'
#' When \code{test \%in\% c("wald", "score")}, the \code{pmethod} argument is
#' ignored, and \code{cov_type} specifies the covariance estimator used to
#' compute the test statistic. For modified working tests, \code{cov_type}
#' determines the covariance matrix used to form the coefficients of the sum of
#' independent chi-squared random variables, and \code{pmethod} specifies the
#' approximation used to obtain the p-value.
#'
#' The output table also includes the Correlation Information Criterion (CIC),
#' a diagnostic for selecting the working correlation structure.
#'
#' @inherit stats::add1 return
#'
#' @references
#' Rotnitzky, A. and Jewell, P. (1990) Hypothesis testing of regression parameters
#' in semiparametric generalized linear models for cluster correlated data.
#' \emph{Biometrika} \bold{77}, 485--497.
#'
#' @seealso \code{\link{anova}}, \code{\link{step}}, \code{\link{geewa}},
#'           \code{\link{geewa_binary}}.
#'
#' @examples
#' data("respiratory")
#' fitted_model <-
#'   geewa(
#'     formula = status ~ baseline + I(treatment == "active") + gender + visit + age,
#'     id = id, repeated = visit,
#'     family = binomial(link = "probit"),
#'     data = respiratory[respiratory$center == "C2", ],
#'     corstr = "ar1",
#'     method = "gee"
#'   )
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
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite"),
           ...) {
    object <- check_geer_object(object)
    opts <- normalize_test_options(
      test = test[1L],
      cov_type = cov_type[1L],
      pmethod = pmethod[1L],
      object = object
    )
    test <- opts$test
    cov_type <- opts$cov_type
    pmethod <- opts$pmethod
    if (missing(scope) || is.null(scope)) {
      stop("no terms in scope", call. = FALSE)
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

    ans[1L, 2L] <- extract_cic(object, cov_type)
    for (i in seq_len(ns)) {
      tt <- scope[[i]]
      add1_model <- update(object, formula = as.formula(paste(". ~ . +", tt)))
      value <- switch(
        test,
        wald = wald_test(object, add1_model, cov_type),
        score = score_test(object, add1_model, cov_type),
        `working-wald`  = working_wald_test(object, add1_model, cov_type, pmethod),
        `working-score` = working_score_test(object, add1_model, cov_type, pmethod),
        `working-lrt`   = working_lr_test(object, add1_model, cov_type, pmethod)
      )

      ans[i + 1L, ] <- c(value$test_df,
                         extract_cic(add1_model, cov_type),
                         value$test_stat,
                         value$test_p)
    }
    aod <- as.data.frame(ans)
    test_type <- test_label(test)
    head <- c(paste("Single term additions using", test_type, "test:"),
              "\nModel:", deparse(object$call$formula))
    structure(aod, heading = head, class = c("anova", "data.frame"))
  }


#' @rdname add1.geer
#' @aliases drop1 drop1.geer
#' @method drop1 geer
#'
#' @examples
#' drop1(fitted_model, test = "score")
#'
#' @export
drop1.geer <- function(object,
                       scope,
                       test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
                       cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
                       pmethod = c("rao-scott", "satterthwaite"),
                       ...) {
  object <- check_geer_object(object)
  opts <- normalize_test_options(
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
  ans[1L, 2L] <- extract_cic(object, cov_type)
  for (i in seq_len(ns)) {
    tt <- scope[[i]]
    drop1_model <- update(object, formula = as.formula(paste(". ~ . -", tt)))
    value <- switch(
      test,
      wald = wald_test(drop1_model, object, cov_type),
      score = score_test(drop1_model, object, cov_type),
      `working-wald` = working_wald_test(drop1_model, object, cov_type, pmethod),
      `working-score` = working_score_test(drop1_model, object, cov_type, pmethod),
      `working-lrt` = working_lr_test(drop1_model, object, cov_type, pmethod)
    )
    ans[i + 1L, ] <- c(value$test_df,
                       extract_cic(drop1_model, cov_type),
                       value$test_stat,
                       value$test_p)
  }
  aod <- as.data.frame(ans)
  test_type <- test_label(test)
  head <- c(paste("Single term deletions using", test_type, "test:"),
            "\nModel:", deparse(object$call$formula))
  structure(aod, heading = head, class = c("anova", "data.frame"))
}


#' @title
#' ANOVA Tables for \code{geer} Objects
#'
#' @aliases anova anova.geer
#' @method anova geer
#'
#' @description
#' Compute analysis of variance (ANOVA) tables for one or more fitted models of
#' class \code{geer}. The table is based on hypothesis tests for regression
#' terms in generalized estimating equation (GEE) models.
#'
#' @param object an object representing a model of the class \code{geer}.
#' @param ... additional objects representing models of the same class \code{geer}.
#' @param test character indicating the hypothesis testing procedure. Options
#'   include the Wald test (\code{"wald"}), the generalized score test
#'   (\code{"score"}), the modified working Wald test (\code{"working-wald"}),
#'   the modified working score test (\code{"working-score"}), and the modified
#'   working likelihood ratio test (\code{"working-lrt"}). Defaults to
#'   \code{"wald"}.
#' @param cov_type character indicating the covariance matrix estimator used for
#'   inference on the regression parameters. Options include the sandwich or
#'   robust estimator (\code{"robust"}), the bias-corrected estimator
#'   (\code{"bias-corrected"}), the degrees of freedom adjusted estimator
#'   (\code{"df-adjusted"}), and the model-based or naive estimator
#'   (\code{"naive"}). Defaults to \code{"robust"}.
#' @param pmethod character indicating the approximation used to compute the
#'   p-value for the modified working tests. Options include the Rao--Scott
#'   approximation (\code{"rao-scott"}) and the Satterthwaite approximation
#'   (\code{"satterthwaite"}). Defaults to \code{"rao-scott"}.
#'
#' @details
#' Details of the hypothesis tests controlled by \code{test} are given in
#' \cite{Rotnitzky and Jewell (1990)}. The option \code{test = "working-lrt"}
#' is valid only when the model is fitted with an \emph{independence} working
#' correlation structure; otherwise an error is returned.
#'
#' When \code{test \%in\% c("wald", "score")}, the \code{pmethod} argument is
#' ignored. In this case, \code{cov_type} specifies the covariance estimator
#' used in the test statistic. For modified working tests, \code{cov_type}
#' determines the covariance matrix used to form the coefficients of the sum of
#' independent chi-squared random variables, and \code{pmethod} specifies the
#' approximation used to compute the p-value.
#'
#' When comparing two or more models, the data must be identical across all
#' fits, and the models must be nested in the order supplied. In particular,
#' each consecutive pair of models must be nested.
#'
#' @return
#' An object of class \code{anova}, representing an analysis-of-variance table
#' within the GEE framework.
#'
#' \itemize{
#'   \item With a single model, the table reports tests for terms added
#'   sequentially (from first to last).
#'   \item With multiple models, the table reports sequential tests comparing
#'   each model to the previous one.
#' }
#'
#' @inherit add1.geer references
#'
#' @seealso
#' \code{\link{drop1}} for type II ANOVA, where each term is dropped one at a
#' time while respecting model hierarchy.
#'
#' @examples
#' data("cerebrovascular")
#'
#' ## Single-model ANOVA (sequential terms)
#' fit_full <- geewa_binary(
#'   formula = ecg ~ treatment + factor(period),
#'   id = id,
#'   data = cerebrovascular,
#'   link = "logit",
#'   orstr = "exchangeable"
#' )
#' anova(fit_full, test = "wald", cov_type = "robust")
#'
#' ## Two-model comparison (models must be nested)
#' fit_null <- geewa_binary(
#'   formula = ecg ~ 1,
#'   id = id,
#'   data = cerebrovascular,
#'   link = "logit",
#'   orstr = "exchangeable"
#' )
#' anova(fit_null, fit_full, test = "wald", cov_type = "robust")
#'
#' @export
anova.geer <-
  function(object,
           ...,
           test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite") ) {
    object <- check_geer_object(object)
    opts <- normalize_test_options(
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
      warning("the following arguments to 'anova.geer' are invalid and dropped: ",
              paste(deparse(dotargs[named]), collapse = ", "))
    }
    dotargs <- dotargs[!named]
    is_geer <- vapply(dotargs, function(x) inherits(x, "geer"), logical(1))
    dotargs <- dotargs[is_geer]
    if (length(dotargs)) {
      dotargs <- lapply(dotargs, check_geer_object)
      return(anova_geerlist(c(list(object), dotargs),
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
      object_list[[1]] <- update(object, formula = . ~ 1)
      for (i in seq_len(nvars)) {
        object_list[[i + 1]] <- update(object_list[[i]],
                                       formula = paste(". ~ . + ", terms[i]))
      }
    } else {
      object_list[[1]] <- update(object, formula = paste(". ~ -1 + ", terms[1]))
      for (i in seq_len(nvars - 1)) {
        object_list[[i + 1]] <- update(object_list[[i]],
                                       formula = paste(". ~ . + ", terms[i + 1]))
      }
    }
    resdf <- vapply(object_list, function(x) as.numeric(x$df.residual), numeric(1))
    table <- data.frame(
      c(NA_real_, resdf[-1]),
      resdf,
      c(NA_real_, resdf[-1]),
      c(NA_real_, resdf[-1])
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
          working_lr_test(object_list[[i]], object_list[[i + 1]], cov_type, pmethod)
      )
      table[i + 1, -2] <- c(value$test_df, value$test_stat, value$test_p)
    }
    test_type <- test_label(test)
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
