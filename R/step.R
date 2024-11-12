#' Choose a model by testing-based in a Stepwise Algorithm
#'
#' Select a formula-based model by testing-based procedures.
#'
#'
#' @param object an object representing a model of the appropriate class
#' \code{geer}. This is used as the initial model in the stepwise search.
#' @param scope defines the range of models examined in the stepwise search. This
#' should be either a single formula, or a list containing components upper and
#' lower, both formulae. See the details for how to specify the formulae and how
#' they are used.
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
#' @param p_enter numeric between 0 and 1 indicating the p-value threshold for
#' adding variables in the stepwise search.
#' @param p_remove numeric between 0 and 1 indicating the p-value threshold for
#' removing variables in the stepwise search.
#' @param direction the mode of stepwise search, can be one of \code{"both"},
#' \code{"backward"}, or \code{"forward"}, with a default of \code{"forward"}. If
#' the scope argument is missing the default for direction is \code{"backward"}.
#' @param steps the maximum number of steps to be considered. The default is 1000
#' (essentially as many as required). It is typically used to stop the process
#' early.
#' @param ... further arguments passed to or from other methods.
#'
#' @export

step_model_selection_p <- function(
    object,
    scope,
    test = "score",
    cov_type = "robust",
    p_enter = 0.20,
    p_remove = 0.20,
    direction = "forward",
    steps = 1000,
    ...) {
  if (!("geer" %in% class(object)))
    stop("Object must be of 'geer' class")
  icheck <- pmatch(test,
                   c("wald", "score"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown test")
  icheck <- pmatch(cov_type,
                   c("robust", "naive", "df-adjusted", "bias-corrected"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown method for the covariance matrix")
  if (!(is.numeric(p_enter) & (0 < p_enter) & (p_enter < 1)))
    stop("'p_enter' must be a numeric object between 0 and 1")
  if (!(is.numeric(p_remove) & (0 < p_remove) & (p_remove < 1)))
    stop("'p_enter' must be a numeric object between 0 and 1")
  icheck <- pmatch(direction,
                   c("forward", "backward", "both"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown direction")
  if (direction == "backward") {
    ans <- step_backward_p(object, scope, test, cov_type, steps, p_remove)
  } else if (direction == "forward") {
    ans <- step_backward_p(object, scope, test, cov_type, steps, p_enter)
  } else {
    ans <- step_both_p(object, scope, test, cov_type, steps, p_enter)
  }
  ans
}


step_backward_p <- function(object,
                            scope,
                            test = "score",
                            cov_type = "robust",
                            steps = 1000,
                            pvalue = 0.05,
                            ...) {
  if (!("geer" %in% class(object)))
    stop("Object must be of 'geer' class")
  icheck <- pmatch(test,
                   c("wald", "score"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown test")
  icheck <- pmatch(cov_type,
                   c("robust", "naive", "df-adjusted", "bias-corrected"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown method for the covariance matrix")
  if (!(is.numeric(pvalue) & (0 < pvalue) & (pvalue < 1)))
    stop("p-value threshold for adding variables must be numeric
         and between 0 and 1")
  model_terms <- terms(object)
  object$call$formula <- object$formula <- model_terms
  if (missing(scope)) {
    factors_drop <- numeric()
    factors_add <- attr(model_terms, "factors")
  } else {
    if (is.list(scope)) {
      factors_drop <- scope$lower
      factors_add <- scope$upper
      factors_drop <- if (!is.null(factors_drop)) {
        attr(terms(update.formula(object, factors_drop)), "factors")
      } else {
        numeric()
      }
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, factors_add)), "factors")
    } else {
      factors_add <- scope
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, scope)), "factors")
      factors_drop <- numeric()
    }
  }
  models_no <- 1
  models <- vector("list", steps)
  models[[models_no]] <-
    list(action = "None",
         variable = NA,
         pvalue = NA,
         parameters_no = length(object$coefficients))
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0) {
    steps <- steps - 1
    ffac <- attr(model_terms, "factors")
    scope <- factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    aod <- NULL
    removed_variable <- NULL
    if (length(scope$drop)) {
      aod <- drop1(fitted_model, scope$drop, test = test, cov_type = cov_type)
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(removed_variable)) {
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df))
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0)
        break
      pvalue_max <- which.max(aod[, 3])
      criterion <- aod[pvalue_max, 3] < pvalue
      if (criterion)
        break
      removed_variable <- rownames(aod)[pvalue_max]
    }
    fitted_model <- update(fitted_model,
                           formula = paste(". ~ . -", removed_variable),
                           evaluate = FALSE)
    fitted_model <- eval.parent(fitted_model)
    obs_no_new <- fitted_model$obs_no
    if (all(is.finite(c(obs_no, obs_no_new))) && obs_no_new != obs_no)
      stop("number of rows in use has changed: remove missing values?")
    model_terms <- terms(fitted_model)
    models_no <- models_no + 1
    models[[models_no]] <-
      list(action = "Elimination",
           variable = removed_variable,
           pvalue = aod[pvalue_max, 3],
           parameters_no = length(fitted_model$coefficients))
  }
  step.results(models = models[seq(models_no)], fitted_model, object)
}


step_forward_p <- function(object,
                           scope,
                           test = "score",
                           cov_type = "robust",
                           steps = 1000,
                           pvalue = 0.05,
                           ...) {
  if (!("geer" %in% class(object)))
    stop("Object must be of 'geer' class")
  icheck <- pmatch(test,
                   c("wald", "score"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown test")
  icheck <- pmatch(cov_type,
                   c("robust", "naive", "df-adjusted", "bias-corrected"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown method for the covariance matrix")
  if (!(is.numeric(pvalue) & (0 < pvalue) & (pvalue < 1)))
    stop("p-value threshold for adding variables must be numeric
         and between 0 and 1")
  model_terms <- terms(object)
  object$call$formula <- object$formula <- model_terms
  if (missing(scope)) {
    factors_drop <- numeric()
    factors_add <- attr(model_terms, "factors")
  } else {
    if (is.list(scope)) {
      factors_drop <- scope$lower
      factors_add <- scope$upper
      factors_drop <- if (!is.null(factors_drop)) {
        attr(terms(update.formula(object, factors_drop)), "factors")
      } else {
        numeric()
      }
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, factors_add)), "factors")
    } else {
      factors_add <- scope
      if (!is.null(factors_add))
        factors_add <- attr(terms(update.formula(object, scope)), "factors")
      factors_drop <- numeric()
    }
  }
  models_no <- 1
  models <- vector("list", steps)
  models[[models_no]] <-
    list(action = "None",
         variable = NA,
         pvalue = NA,
         parameters_np = length(object$coefficients))
  obs_no <- object$obs_no
  fitted_model <- object
  while (steps > 0) {
    steps <- steps - 1
    ffac <- attr(model_terms, "factors")
    scope <- factor.scope(ffac, list(add = factors_add, drop = factors_drop))
    aod <- NULL
    added_variable <- NULL
    if (length(scope$add))
      aod <- add1(fitted_model, scope$add, test = test, cov_type = cov_type)
    attr(aod, "heading") <- NULL
    nzdf <- if (!is.null(aod$Df))
      aod$Df != 0 | is.na(aod$Df)
    aod <- aod[nzdf, ]
    if (is.null(aod) || ncol(aod) == 0)
      break
    pvalue_min <- which.min(aod[, 3])
    criterion <- aod[pvalue_min, 3] > pvalue
    if (criterion)
      break
    added_variable <- rownames(aod)[pvalue_min]
    fitted_model <- update(fitted_model,
                           formula = paste(" . ~ . +", added_variable),
                           evaluate = FALSE)
    fitted_model <- eval.parent(fitted_model)
    obs_no_new <- fitted_model$obs_no
    if (all(is.finite(c(obs_no, obs_no_new))) && obs_no_new != obs_no)
      stop("number of rows in use has changed: remove missing values?")
    model_terms <- terms(fitted_model)
    models_no <- models_no + 1
    models[[models_no]] <-
      list(action = "Addition",
           variable = added_variable,
           pvalue = aod[pvalue_min, 3],
           parameters_no = length(fitted_model$coefficients))
  }
  step.results(models = models[seq(models_no)], fitted_model, object)
}


step_both_p <- function(object,
                        scope,
                        test = "wald",
                        cov_type = "robust",
                        steps = 1000,
                        p_enter = 0.05,
                        p_remove = 0.05,
                        ...) {
  models <- vector("list", steps)
  models_no <- 1
  models[[models_no]] <-
    list(action = "",
         variable = NA,
         pvalue = NA,
         parameters_no = length(object$coefficients))
  obs_no <- object$obs_no
  fitted_model <- object
  criterion <- 0
  change_variable <- ""
  while (steps > 0 & criterion < 2) {
    steps <- steps - 1
    fitted_model <- step_backward_p(fitted_model,
                                    scope = scope,
                                    test = test,
                                    cov_type = cov_type,
                                    steps = 1,
                                    pvalue = p_remove)
    if (is.na(fitted_model$anova$pvalue)) {
      criterion <- criterion + 1
    } else {
      criterion <- 0
      models_no <- models_no + 1
      models[[models_no]] <-
        list(action = fitted_model$anova$Step,
             variable = fitted_model$anova$Variable,
             pvalue = fitted_model$anova$pvalue,
             p = length(fitted_model$coefficients))
      if (change_variable == fitted_model$anova$Variable) {
        break
      } else {
        change_variable <- fitted_model$anova$Variable
      }
    }
    steps <- steps - 1
    fitted_model <- step_forward_p(fitted_model,
                                    scope = scope,
                                    test = test,
                                    cov_type = cov_type,
                                    steps = 1,
                                    pvalue = p_enter)
    if (is.na(fitted_model$anova$pvalue)) {
      criterion <- criterion + 1
    } else {
      criterion <- 0
      models_no <- models_no + 1
      models[[models_no]] <-
        list(action = fitted_model$anova$Step,
             variable = fitted_model$anova$Variable,
             pvalue = fitted_model$anova$pvalue,
             p = length(fitted_model$coefficients))
      if (change_variable == fitted_model$anova$Variable) {
        break
      } else {
        change_variable <- fitted_model$anova$Variable
      }
    }
  }
  #models
  step.results(models = models[seq(models_no)], fitted_model, object)
}


step.results <- function(models, fit, object) {
  variable <- sapply(models, `[[`, "variable")
  pval <- sapply(models, `[[`, "pvalue")
  action <- sapply(models, `[[`, "action")
  heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
               "\nInitial Model:", deparse(formula(object)), "\nFinal Model:",
               deparse(formula(fit)), "\n")
  aod <- data.frame(Step = action,
                    Variable = variable,
                    pvalue = as.numeric(round(pval , 4)))
  attr(aod, "heading") <- heading
  if (nrow(aod) > 1)
    aod <- aod[-1, ]
  fit$anova <- aod
  rownames(fit$anova) <- seq(nrow(aod))
  fit
}

step_model_selection_p <- function(
    object,
    scope,
    test = "score",
    cov_type = "robust",
    p_enter = 0.20,
    p_remove = 0.20,
    direction = "forward",
    steps = 1000,
    ...) {
  if (!("geer" %in% class(object)))
    stop("Object must be of 'geer' class")
  icheck <- pmatch(test,
                   c("wald", "score"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown test")
  icheck <- pmatch(cov_type,
                   c("robust", "naive", "df-adjusted", "bias-corrected"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown method for the covariance matrix")
  if (!(is.numeric(p_enter) & (0 < p_enter) & (p_enter < 1)))
    stop("'p_enter' must be a numeric object between 0 and 1")
  if (!(is.numeric(p_remove) & (0 < p_remove) & (p_remove < 1)))
    stop("'p_enter' must be a numeric object between 0 and 1")
  icheck <- pmatch(direction,
                   c("forward", "backward", "both"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0)
    stop("unknown direction")
  if (direction == "backward") {
    ans <- step_backward_p(object, scope, test, cov_type, steps, p_remove)
  } else if (direction == "forward") {
    ans <- step_forward_p(object, scope, test, cov_type, steps, p_enter)
  } else {
    ans <- step_both_p(object, scope, test, cov_type, steps, p_enter)
  }
  ans
}
