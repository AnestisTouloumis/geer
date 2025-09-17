geer <- function(x, ...) { # nolint
  UseMethod("geer")
}

#' @export
geer.default <- function(x, ...) {
  object <- list()
  object$coefficients <- x$coefficients
  object$residuals <- x$residuals
  object$fitted.values <- x$fitted.values
  object$rank <- x$rank
  object$qr <- x$qr
  object$family <- x$family
  object$linear.predictors <- x$linear.predictors
  object$iter <- x$iter
  object$prior.weights <- x$prior.weights
  object$df.residual <- x$df.residual
  object$y <- x$y
  object$x <- x$x
  object$id <- x$id
  object$repeated <- x$repeated
  object$converged <- x$converged
  object$call <- x$call
  object$formula <- x$formula
  object$terms <- x$terms
  object$data <- x$data
  object$offset <- x$offset
  object$control <- x$control
  object$method <- x$method
  object$contrasts <- x$contrasts
  object$xlevels <- x$xlevels
  object$naive_covariance <- x$naive_covariance
  object$robust_covariance <- x$robust_covariance
  object$bias_corrected_covariance <- x$bias_corrected_covariance
  object$association_structure <- x$association_structure
  object$alpha <- x$alpha
  object$phi <- x$phi
  object$obs_no <- x$obs_no
  object$clusters_no <- x$clusters_no
  object$min_cluster_size <- x$min_cluster_size
  object$max_cluster_size <- x$max_cluster_size
  class(object) <- "geer"
  object
}


#' @export
print.geer <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nNumber of iterations :", x$iter, "\n")
  cat("Algorithm converged  :", x$converged, "\n")
}


#' @method summary geer
#' @export
summary.geer <- function(object, cov_type = "robust", ...) {
  coefficients_hat <- coef(object)
  standard_errors <- sqrt(diag(vcov(object, cov_type)))
  z_statistics  <- coefficients_hat / standard_errors
  pvalue <- 2 * (1 - pnorm(abs(z_statistics)))
  TAB <- cbind(
    Estimate = coefficients_hat,
    `Std.Error` = standard_errors,
    `z value` = z_statistics,
    `Pr(>|z|)` = pvalue
  )
  res <- list(
    coefficients = TAB,
    family = object$family,
    alpha = object$alpha,
    call = object$call,
    residuals = object$residuals,
    niter = object$iter,
    converged = object$converged,
    phi = object$phi,
    association_structure = object$association_structure,
    method = object$method,
    cov_type = cov_type
  )
  class(res) <- "summary.geer"
  res
}


#' @export
print.summary.geer <- function(x, ...) {
  cat("\ncall:\n")
  print(x$call)
  cat("\nEstimating Method   :", x$method, "\n")
  cat("Number of iterations:", x$iter, "\n")
  cat("Algorithm converged :", x$converged, "\n")
  cat("\nMarginal Model\n")
  cat("Family       :", x$family$family, "\n")
  cat("Link Function:", x$family$link, "\n")
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  cat("Std.Errors are taken from the", x$cov_type, "covariance matrix.", "\n")
  cat("\nDispersion Parameter:", round(x$phi, digits = 4), "\n")
  cat("\nAssociation Structure:", x$association_structure, "\n")
  if (length(x$alpha) == 1) {
    cat("Association Parameter:", round(x$alpha, digits = 4), "\n")
  } else {
    cat("Association Parameters:\n")
    print(round(x$alpha, digits = 4))
    cat("\n")
  }
}



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
#' compute the test statistic. For other tests, \code{cov_type} determines the
#' covariance matrix used to form the coefficients of the sum of independent
#' chi-squared random variables, while \code{p_method} specifies the
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
#' fitted_model <-
#' geewa(formula = status ~ baseline + I(treatment == "active") + gender + visit + age,
#'       id = id, repeated = visit, family = binomial(link = "probit"),
#'       data = respiratory[respiratory$center=="C2", ], corstr = "ar1",
#'       method = "gee")
#' add1(fitted_model,
#'      scope = .~. + baseline:age + age:visit + I(treatment == "active"):age + age:gender,
#'      test = "score")
#'
#' @export
add1.geer <-
  function(object,
           scope,
           test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite"),
           ...) {
    test <- match.arg(test)
    if (test %in% c("working-wald", "working-score", "working-lrt"))
      pmethod <- match.arg(pmethod)
    cov_type <- match.arg(cov_type)
    if (test == "working-lrt" & object$association_structure != "independence")
      stop("the modified working lrt can only be applied to the independence working model")
    cov_type <- match.arg(cov_type)
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



#' @rdname add1.geer
#' @aliases drop1 drop1.geer
#' @method drop1 geer
#'
#' @examples
#' drop1(fitted_model, test = "score")
#'
#' @export
drop1.geer <-
  function(object,
           scope,
           test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite"),
           ...) {
    test <- match.arg(test)
    if (test %in% c("working-wald", "working-score", "working-lrt"))
      pmethod <- match.arg(pmethod)
    cov_type <- match.arg(cov_type)
    if (test == "working-lrt" & object$association_structure != "independence")
      stop("the modified working lrt can only be applied to the independence working model")
    cov_type <- match.arg(cov_type)
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



#' @title
#' ANOVA Tables for \code{geer} Objects
#'
#' @aliases anova anova.geer
#' @method anova geer
#'
#' @description
#' Compute analysis of variance (ANOVA) tables for one or more fitted models of
#' class \code{geer}, using a variety of hypothesis testing procedures.
#'
#' @param object an object representing a model of the class \code{geer}.
#' @param ... additional objects representing models of the same class \code{geer}.
#' @param test character indicating the hypothesis testing procedure. Options
#'        include the Wald test (\code{"wald"}), the generalized score test
#'        (\code{"score"}), the modified working Wald test
#'        (\code{"working-wald"}), the modified working score test
#'        (\code{"working-score"}) and the modified working likelihood ratio
#'        test (\code{"working-lrt"}). Defaults to \code{"wald"}.
#' @param cov_type character indicating the covariance matrix estimator
#'        of the regression parameters. Options include the sandwich or
#'        robust estimator (\code{"robust"}), the bias-corrected estimator
#'        (\code{"bias-corrected"}), the degrees of freedom adjusted estimator
#'        (\code{"df-adjusted"}) and the model-based or naive estimator
#'        (\code{"naive"}). Defaults to \code{"robust"}.
#' @param pmethod character indicating the method used to approximate the
#'        p-value when one of the modified working tests is selected. Options
#'        include the Rao–Scott approximation (\code{"rao-scott"}) and the
#'        Satterthwaite approximation (\code{"satterthwaite"}). Defaults to
#'        \code{"rao-scott"}.
#'
#' @details
#' Details of the hypothesis tests controlled by \code{test} are given in
#' \cite{Rotnitzky and Jewell (1990)}. The option \code{test = "working-lrt"}
#' is valid only when the model is fitted with an \emph{independence} working
#' correlation structure; otherwise an error is returned.
#'
#' When \code{test \%in\% c("wald", "score")}, the \code{pmethod} argument is
#' ignored. In this case, \code{cov_type} specifies the covariance estimator
#' used in the test statistic. For other tests, \code{cov_type} determines the
#' covariance matrix used to form the coefficients of the sum of independent
#' chi-squared random variables, and \code{p_method} specifies the approximation
#' used to compute the p-value.
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
#'   \item With a single model, the table reports the significance of each model
#'   term.
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
#' @export

anova.geer <-
  function(object,
           ...,
           test = c("wald", "score", "working-wald", "working-score", "working-lrt"),
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           pmethod = c("rao-scott", "satterthwaite") ){
    test <- match.arg(test)
    if (test %in% c("working-wald", "working-score", "working-lrt"))
      pmethod <- match.arg(pmethod)
    cov_type <- match.arg(cov_type)
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
      rep_len(FALSE, length(dotargs)) else (names(dotargs) != "")
    if (any(named))
      warning("the following arguments to 'anova.geer' are invalid and dropped: ",
              paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is_geer <- vapply(dotargs, function(x) inherits(x, "geer"), NA)
    dotargs <- dotargs[is_geer]
    if (length(dotargs))
      return(anova_geerlist(c(list(object), dotargs),
                            test = test,
                            cov_type = cov_type,
                            pmethod = pmethod))
    if (test == "working-lrt" & object$association_structure != "independence")
      stop("the modified working lr test requires independence working models")
    terms <- attr(object$terms, "term.labels")
    intercept <- attr(object$terms, "intercept")
    variables <- attr(object$terms, "variables")
    varseq <- attr(object$x, "assign")
    nvars <- max(0, varseq)
    object_list <- list()
    if (intercept == 1) {
      object_list[[1]] <- update(object, formula = . ~ 1)
      for (i in seq_len(nvars)) {
        object_list[[i + 1]] <- update(object_list[[i]], formula = paste(". ~ . + ", terms[i]))
      }
    } else {
      object_list[[1]] <- update(object, formula = paste(". ~ -1 + ", terms[1]))
      for (i in seq_len(nvars - 1)) {
        object_list[[i + 1]] <- update(object_list[[i]], formula = paste(". ~ . + ", terms[i + 1]))
      }
    }
    resdf  <- as.numeric(lapply(object_list, function(x) x$df.residual))
    table <- data.frame(c(NA, resdf[-1]), resdf, c(NA, resdf[-1]), c(NA, resdf[-1]))
    if (intercept == 1) {
      dimnames(table) <- list(c("NULL", terms),
                              c("Df", "Resid. Df", "Chi", "Pr(>Chi)"))
    } else {
      dimnames(table) <- list(c(terms),
                              c( "Df", "Resid. Df", "Chi", "Pr(>Chi)"))
    }
    test_type <- switch(test,
                        wald = "Wald",
                        score = "Score",
                        `working-wald` = "Working Wald",
                        `working-score` = "Working Score",
                        `working-lrt` = "Working LRT")
    for (i in seq_len(length(object_list) - 1)) {
      value <-
        switch(test,
               wald =
                 wald_test(object_list[[i]], object_list[[i + 1]], cov_type),
               score =
                 score_test(object_list[[i]], object_list[[i + 1]], cov_type),
               `working-wald` =
                 working_wald_test(object_list[[i]], object_list[[i + 1]], cov_type, pmethod),
               `working-score` =
                 working_score_test(object_list[[i]], object_list[[i + 1]], cov_type, pmethod),
               `working-lrt` =
                 working_lr_test(object_list[[i]], object_list[[i + 1]], cov_type, pmethod)
        )
      table[i + 1, -2] <- c(value$test_df, value$test_stat, value$test_p)
    }
    test_type <- switch(test,
                        wald = "Wald",
                        score = "Score",
                        `working-wald` = "Modified Working Wald",
                        `working-score` = "Modified Working Score",
                        `working-lrt` = "Modified Working LRT")
    title <- paste("Analysis of ", test_type, " Statistic Table",
                   "\n\nModel: ", object$family$family,
                   ", link: ", object$family$link,
                   "\n\nResponse: ", as.character(variables[-1L])[1L],
                   "\n\nTerms added sequentially (first to last)\n\n",
                   sep = "")
    structure(table, heading = title, class = c("anova", "data.frame"))
  }



#' @title
#' Extract Model Coefficients from a \code{geer} Object
#'
#' @aliases coef.geer coef coefficients
#' @method coef geer
#'
#' @description
#' Extracts model coefficients from a \code{geer} object. The function
#' \code{coefficients} is an alias.
#'
#' @inheritParams add1.geer
#'
#' @return
#' A named numeric vector containing the estimated model coefficients.
#'
#' @examples
#' data("leprosy")
#' fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
#'              family = poisson(link = "log"), id = id, data = leprosy)
#' coef(fit)
#'
#' data("cerebrovascular")
#' fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
#'                     id = id, data = cerebrovascular)
#' coef(fit)
#'
#' @export
coef.geer <- function(object, ...){
  coeffs <- object$coefficients
  coeffs
}



#' @title
#' Confidence Intervals for Model Parameters from a \code{geer} Object
#'
#' @aliases confint confint.geer
#' @method confint geer
#'
#' @description
#' Compute confidence intervals for one or more regression parameters from a
#' fitted \code{geer} model.
#'
#' @inheritParams add1.geer
#' @inheritParams stats::confint
#'
#' @details
#' Confidence intervals are based on the estimated covariance matrix of the
#' coefficients. The specific form of the covariance matrix is controlled by the
#' \code{cov_type} argument; see \code{\link{vcov}} for available options and
#' formulae.
#'
#' @inherit stats::confint.default return
#'
#' @examples
#' data("leprosy")
#' fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
#'              family = poisson(link = "log"), id = id, data = leprosy)
#' confint(fit)
#'
#' data("cerebrovascular")
#' fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
#'                     id = id, data = cerebrovascular)
#' confint(fit)
#'
#' @export
confint.geer <- function(object, parm, level = 0.95, cov_type = "robust", ...) {
  coefficients <- coef(object)
  coefficients_names <- names(coefficients)
  if (missing(parm)) {
    parm <- coefficients_names
  } else if (is.numeric(parm)) {
    parm <- coefficients_names[parm]
  }
  alpha <- (1 - level) / 2
  alpha <- c(alpha, 1 - alpha)
  pct <- format_perc(alpha, 3)
  percentiles <- qnorm(alpha)
  ans <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  standard_errors <- sqrt(diag(vcov(object, cov_type = cov_type)))[parm]
  ans[] <- coefficients[parm] + standard_errors %o% percentiles
  ans
}



#' @title
#' Extract Model Fitted Values from a \code{geer} Object
#'
#' @rdname fitted
#' @aliases fitted fitted.values
#' @method fitted geer
#'
#' @description
#' Extract fitted values from a model of class \code{geer}.
#' The function \code{fitted.values()} is an alias.
#'
#' @inheritParams coef.geer
#'
#' @return
#' A numeric vector of fitted values extracted from \code{object}.
#'
#' @examples
#' data("leprosy")
#' fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
#'              family = poisson(link = "log"), id = id, data = leprosy)
#' fitted(fit)
#'
#' data("cerebrovascular")
#' fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
#'                     id = id, data = cerebrovascular)
#' fitted(fit)
#'
#' @export
fitted.geer <- function(object, ...){
  object$fitted.values
}



#' @title
#' Construct Design Matrices from a \code{geer} Object
#'
#' @aliases model.matrix
#' @method model.matrix geer
#'
#' @description
#' Build the design (model) matrix from a fitted \code{geer} object. Factor
#' variables are expanded to dummy/contrast-coded columns (according to the
#' active contrasts), and interaction terms are expanded accordingly.
#'
#' @inheritParams coef.geer
#'
#' @details
#' The design matrix is constructed based on \code{terms(object)}, using the data
#' stored in \code{object$data}.
#'
#' For interaction terms, the variable whose levels vary fastest is the first one
#' listed in the formula (not in the term). For example, in
#' \code{~ a + b + b:a}, the interaction term \code{b:a} will vary fastest with
#' respect to \code{a}.
#'
#' By convention, if the response variable also appears on the right-hand side of
#' the formula it is dropped (with a warning). However, interactions involving
#' that term are retained.
#'
#' @return
#' The design matrix for the marginal regression model with the specified formula and
#' data.
#'
#' There is an attribute \code{"assign"}, an integer vector with an entry for
#' each column in the matrix giving the term in the formula which gave rise to
#' the column. Value 0 corresponds to the intercept (if any), and positive
#' values to terms in the order given by the \code{term.labels} attribute of the
#' terms structure corresponding to \code{object}.
#'
#' If there are any factors in terms in the model, there is an attribute
#' \code{"contrasts"}, a named list with an entry for each factor. This
#' specifies the contrasts that would be used in terms in which the factor is
#' coded by contrasts (in some terms dummy coding may be used), either as a
#' character vector naming a function or as a numeric matrix.
#'
#' @examples
#' data("leprosy")
#' fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
#'              family = poisson(link = "log"), id = id, data = leprosy)
#' model.matrix(fit)
#'
#' data("cerebrovascular")
#' fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
#'                     id = id, data = cerebrovascular)
#' model.matrix(fit)
#'
#' @export
model.matrix.geer <-	function(object,...){
  model.matrix(object = object$terms,
               data = object$data,
               contrasts = object$contrasts)
}



#' @title
#' Predictions from a \code{geer} Object
#'
#' @aliases predict predict.geer
#' @method predict geer
#'
#' @description
#' Generate predictions, optionally with standard errors, from a fitted
#' \code{geer} object.
#'
#' @inheritParams add1.geer
#' @param newdata optional data frame in which to look for variables used
#'        for prediction. If omitted, predictions are made on the data used for
#'        fitting.
#' @param type type of prediction required. Options include the scale of the
#'        linear predictors (\code{"link"}) and the scale of the response
#'        variable (\code{"response"}). By default, the scale of the linear
#'        predictors is used.
#' @param se.fit logical indicating if standard errors are required.
#'        Defaults to \code{FALSE}.
#'
#' @details
#' Predictions are obtained by computing the model matrix for \code{newdata}
#' (or the original data if \code{newdata} is missing) and multiplying by the
#' estimated coefficients. If \code{type = "response"}, predictions are
#' transformed via the inverse link function.
#'
#' When \code{se.fit = TRUE}, approximate standard errors of predictions are
#' computed using the model-based covariance matrix of the estimated
#' coefficients.
#'
#' If \code{newdata} is omitted the predictions are based on the data used for
#' the fit.
#'
#' @return
#' \describe{
#'   \item{If \code{se.fit = FALSE}:}{A numeric vector (or matrix) of predicted
#'   values.}
#'   \item{If \code{se.fit = TRUE}:}{A list with components:
#'     \itemize{
#'       \item \code{fit} — predictions, as above.
#'       \item \code{se.fit} — estimated standard errors of the predictions.
#'     }
#'   }
#' }
#'
#' @seealso
#' \code{\link{geewa}}, \code{\link{geewa_binary}}, \code{\link[stats]{predict}}.
#'
#' @export
predict.geer <-
  function(object,
           newdata = NULL,
           type = c("link", "response"),
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           se.fit = c(FALSE, TRUE),
           ...) {
    type <- match.arg(type)
    se.fit <- match.arg(se.fit)
    if (missing(newdata)) {
      if (type == "link") {
        predicts <- object$linear.predictors
      } else {
        predicts <- object$fitted.values
      }
      if (se.fit) {
        covariance_matrix <- vcov(object, type = cov_type)
        model_matrix <- object$x
        se.fit <-
          sqrt(apply(tcrossprod(model_matrix, covariance_matrix) * model_matrix, 1, sum))
        if (type == "response")
          se.fit <- se.fit * abs(object$family$mu.eta(object$linear.predictors))
        names(predicts) <- names(se.fit)
        predicts <- list(fit = predicts, se.fit = se.fit)
      }
    } else {
      newdata <- data.frame(newdata)
      model_frame <-
        model.frame(delete.response(object$terms), newdata, xlev = object$xlevels)
      model_matrix <-
        model.matrix(delete.response(object$terms), model_frame, contrasts = object$contrasts)
      predicts <- c(model_matrix %*% object$coefficients)
      offset_term <- model.offset(model_frame)
      if (!is.null(offset_term)) predicts <- predicts + c(offset_term)
      if (se.fit) {
        covariance_matrix <- vcov(object, cov_type)
        se.fit <-
          sqrt(apply(tcrossprod(model_matrix, covariance_matrix) * model_matrix, 1, sum))
        if (type == "response") {
          se.fit <- se.fit * abs(object$family$mu.eta(predicts))
          predicts <- object$family$linkinv(predicts)
        }
        names(predicts) <- names(se.fit)
        predicts <- list(fit = predicts, se.fit = se.fit)
      } else if (type == "response") {
        predicts <- object$family$linkinv(predicts)
      }
    }
    predicts
  }



#' @title
#' EResiduals from a \code{geer} Object
#'
#' @aliases resid residuals residuals.geer
#' @method residuals geer
#'
#' @description
#' Extract residuals of different types from a fitted \code{geer} object.
#'
#' @inheritParams coef.geer
#' @param type character indicating whether the type of the residuals to return.
#'        Options include the working residuals (\code{"working"}), the pearson
#'        residuals (\code{"pearson"}) and the deviance residuals
#'        (\code{"deviance"}). By default, the working residuals are returned.
#'
#' @details
#' The residuals are computed according to the marginal distribution specified
#' by \code{object$family$family}.
#'
#' \itemize{
#'   \item Working residuals: raw differences between observed and fitted values.
#'   \item Pearson residuals: standardized by the variance function of the
#'         specified family.
#'   \item Deviance residuals: signed square roots of the contribution of each
#'         observation to the model deviance.
#' }
#'
#' @return
#' A numeric vector of residuals of the requested \code{type}.
#'
#' @seealso
#' \code{\link{geewa}}, \code{\link{geewa_binary}},
#' \code{\link[stats]{residuals}}
#'
#' @export
residuals.geer <- function(object,
                           type = c("working", "pearson", "deviance"),
                           ...) {
  type <- match.arg(type)
  response_vector <- object$y
  mu_vector <- object$fitted.values
  weight_vector <- object$prior.weights
  ans <- switch(type,
                deviance = if (object$df.residual > 0) {
                  deviance_res <-
                    sqrt(pmax((object$family$dev.resids)(response_vector, mu_vector, weight_vector), 0))
                  ifelse(response_vector > mu_vector, deviance_res, -deviance_res)
                } else rep.int(0, length(mu_vector)),
                pearson =
                  c(get_pearson_residuals(object$family$family, response_vector,
                                          mu_vector, weight_vector)),
                working = object$residuals
  )
  ans
}



#' @title
#' Extract Variance-Covariance Matrix from a \code{geer} Object
#'
#' @aliases vcov vcov.geer
#' @method vcov geer
#'
#' @description
#' Extract the variance–covariance matrix of the regression parameters from a
#' fitted \code{geer} object. The parameters correspond to those returned by
#' \code{\link{coef}}.
#'
#' @inheritParams add1.geer
#'
#' @details
#' The form of the covariance estimator is controlled by the argument
#' \code{cov_type}:
#' \itemize{
#'   \item \code{"robust"} — the sandwich (robust) covariance estimator
#'     \insertCite{LiangZeger1986}{}
#'   \item \code{"naive"} — the model-based covariance estimator
#'     \insertCite{LiangZeger1986}{}
#'   \item \code{"df-adjusted"} — the small-sample adjusted covariance matrix
#'     \insertCite{MacKinnon1985}{}
#'   \item \code{"bias-corrected"} — the bias-corrected covariance estimator
#'     \insertCite{MorelEtAl2003}{}
#' }
#'
#' @return
#' A square numeric matrix of estimated covariances between regression
#' coefficients. Rows and columns are named according to the coefficient names
#' returned by \code{\link[stats]{coef}}.
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) A comparison of two bias-corrected covariance
#' estimators for generalized estimating equations. \emph{Biometrika} \bold{73},
#' 13-–22.
#'
#' MacKinnon, J.G. (1985) Some heteroskedasticity-consistent
#' covariance matrix estimators with improved finite sample properties.
#' \emph{Journal of Econometrics} \bold{29}, 305–-325.
#'
#' Morel, J.G., Bokossa, M.C. and Neerchal N.K. (2003) Small sample correction
#' for the variance of GEE estimators. \emph{Biometrical Journal} \bold{45},
#' 395–-409.
#'
#' @examples
#' data("cerebrovascular")
#' fitted_model <- geewa_binary(formula = ecg ~ period + treatment,
#'                              id = id,
#'                              data = cerebrovascular,
#'                              link = "logit",
#'                              orstr = "exchangeable")
#' vcov(fitted_model)
#' vcov(fitted_model, cov_type = "bias-corrected")
#'
#' @export
vcov.geer <-
  function(object,
           cov_type = c("robust", "bias-corrected", "df-adjusted", "naive"),
           ...) {
    cov_type <- match.arg(cov_type)
    if (cov_type == "robust") {
      ans <- object$robust_covariance
    } else if (cov_type == "bias-corrected") {
      total_obs <- object$obs_no
      sample_size <- object$clusters_no
      parameters_no <- length(object$coefficients)
      df_residuals <- object$df.residual
      kappa <-
        ((total_obs - 1)/df_residuals) * (sample_size / (sample_size - 1))
      delta <-
        min(parameters_no / (sample_size - parameters_no),
            0.5)
      robust_covariance <- object$robust_covariance
      naive_covariance <- object$naive_covariance
      trace_robust_naive_inverse <- sum(diag(solve(naive_covariance, robust_covariance)))
      ksi <-
        max(1,
            trace_robust_naive_inverse/parameters_no)
      ans <- kappa * robust_covariance + delta * ksi * naive_covariance
    } else if (cov_type == "df-adjusted") {
      sample_size <- object$clusters_no
      parameters_no <- length(object$coefficients)
      ans <-
        (sample_size / (sample_size - parameters_no)) *
        object$robust_covariance
    } else {
      ans <- object$naive_covariance
    }
    ans
  }
