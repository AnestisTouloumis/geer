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
#' Add or Drop All Possible Single Terms to a Model from a \code{geer} Object
#'
#' @rdname add1.geer
#' @aliases add1 add1.geer
#' @method add1 geer
#'
#' @inherit stats::add1.default description
#'
#' @inheritParams anova.geer
#' @inheritParams stats::add1
#' @param ... additional argument(s) for methods.
#'
#' @details
#' For \code{drop1}, a missing scope is taken to be all terms in the model. The
#' hierarchy is respected when considering terms to be added or dropped:
#' all main effects contained in a second-order interaction must remain, and so
#' on.
#'
#' In a scope formula . means 'what is already there'.
#'
#' Details about the hypothesis testing procedures implied by \code{test} can
#' be found in \cite{Rotnitzky and Jewell (1990)}. Note that
#' \code{test = "working-lrt"} is only available to fitted models with an
#' independence working association structure. Otherwise, an error message is
#' returned.
#'
#' When \code{test = "wald"} or \code{test = "score"}, the \code{pmethod}
#' argument is ignored and \code{cov_type} specifies the covariance matrix
#' estimate used to calculate the corresponding test statistic. Otherwise,
#' \code{cov_type} specifies the covariance matrix estimate used to calculate
#' the coefficients of the sum of independent chi-squared random variables, and
#' \code{p_method} argument specifies the approximation method used to calculate
#' the p-value of the resulting test statistic.
#'
#' The output table also gives the Correlation Information Criterion (CIC)
#' criterion.
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
#' geewa(formula = y ~ baseline + I(treatment == "active") + gender + visit + age,
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
#' Compute analysis of variance tables using hypothesis testing procedures for
#' one or more fitted model.
#'
#' @param object an object representing a model of the class \code{geer}.
#' @param ... additional objects representing models of the same class \code{geer}.
#' @param test character indicating the hypothesis testing procedure. Options
#'        include the Wald test (\code{"wald"}), the generalized score test
#'        (\code{"score"}), the modified working Wald test
#'        (\code{"working-wald"}), the modified working score test
#'        (\code{"working-score"}) and the modified working likelihood ratio
#'        test (\code{"working-lrt"}). By default, the Wald test is performed.
#' @param cov_type character indicating the covariance matrix required for
#'        testing procedure. Options include the sandwich or robust covariance
#'        matrix (\code{"robust"}), the bias-corrected covariance matrix
#'        (\code{"bias-corrected"}), the degrees of freedom adjusted covariance
#'        matrix (\code{"df-adjusted"}) and the model-based or naive covariance
#'        matrix (\code{"naive"}). By default, the robust covariance
#'        matrix is used.
#' @param pmethod character indicating the method used to approximate the p-value
#'        when the modified working Wald test, the modified working score test or
#'        the modified working likelihood ratio test is selected. Options include
#'        the Rao-Scott approximation (\code{"rao-scott"}) and the Satterthwaite
#'        approximation (\code{"satterthwaite"}). By default, the Rao-Scott
#'        approximation is used.
#'
#' @details
#' Details about the hypothesis testing procedures implied by \code{test} can
#' be found in \cite{Rotnitzky and Jewell (1990)}. Note that
#' \code{test = "working-lrt"} is only available to fitted models with an
#' independence working association structure. Otherwise, an error message is
#' returned.
#'
#' When \code{test = "wald"} or \code{test = "score"}, the \code{p_method}
#' argument is ignored and \code{cov_type} specifies the covariance matrix
#' estimate used to calculate the corresponding test statistic. Otherwise,
#' \code{cov_type} specifies the covariance matrix estimate used to calculate
#' the coefficients of the sum of independent chi-squared random variables, and
#' \code{p_method} argument specifies the approximation method used to calculate
#' the p-value of the resulting test statistic.
#'
#' The comparison between two or more models will only be valid if they are
#' fitted to the same dataset.
#'
#' @return
#' This function returns an object of class \code{anova}. These objects represent
#' analysis-of-variance tables within the GEE framework.
#'
#' When given a single argument, \code{anova} produces a table which tests
#' whether the model terms are significant. When given a sequence of objects,
#' \code{anova} tests the models against one another in the order specified.
#' This requires that two consecutive models are nested.
#'
#' @inherit add1.geer references
#'
#' @seealso \code{\link{drop1}} for so-called 'type II' ANOVA where each term is
#' dropped one at a time respecting their hierarchy.
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
#' Extracts model coefficients from objects returned by modeling functions.
#' \code{coefficients} is an alias for it.
#'
#' @inheritParams add1.geer
#'
#' @return
#' A named numeric vector with the coefficients extracted from the \code{object}.
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
#' Confidence intervals for one or more parameters in a fitted model.
#'
#' @inheritParams add1.geer
#' @inheritParams stats::confint
#' @param cov_type character indicating the type of the covariance matrix to be
#'        used in the construction of the confidence interval. Options include
#'        the sandwich or robust covariance matrix (\code{"robust"}), the
#'        bias-corrected covariance matrix (\code{"bias-corrected"}), the
#'        degrees of freedom adjusted covariance matrix (\code{"df-adjusted"})
#'        and the model-based or naive covariance matrix (\code{"naive"}). By
#'        default, the robust covariance matrix is used.
#'
#' @details
#' The references in \code{\link{vcov}} include the formulae for the covariance
#' type implied by \code{cov_type}.
#'
#' @inherit stats::confint.default return
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
#' \code{fitted} extracts fitted values from objects returned by modeling
#' functions. \code{fitted.values} is an alias for it.
#'
#' @inheritParams coef.geer
#'
#' @return
#' Fitted values extracted from the \code{object}.
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
#' Creates a design or model matrix from a object, e.g., by expanding factors
#' to a set of dummy variables (depending on the contrasts) and expanding
#' interactions similarly.
#'
#' @inheritParams coef.geer
#'
#' @details
#' \code{model.matrix} creates a design matrix from the description given in
#' \code{terms(object)}, using the data in \code{object$data}.
#'
#' In an interaction term, the variable whose levels vary fastest is the first
#' one to appear in the formula (and not in the term), so in
#' \code{~ a + b + b:a} the interaction will have a varying fastest.
#'
#' By convention, if the response variable also appears on the right-hand side
#' of the formula it is dropped (with a warning), although interactions
#' involving the term are retained.
#'
#' @return
#' The design matrix for the regression model with the specified formula and
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
#' @export
model.matrix.geer <-	function(object,...){
  model.matrix(object = object$terms,
               data = object$data,
               contrasts = object$contrasts)
}



#' @title
#' Model Predictions for a \code{geer} Object
#'
#' @aliases predict predict.geer
#' @method predict geer
#'
#' @description
#' Obtains predictions and optionally estimates standard errors of those
#' predictions from a fitted generalized estimating equations model object.
#'
#' @inheritParams add1.geer
#' @param newdata optional data frame in which to look for variables with
#'        which to predict. If omitted, the fitted linear predictors are used.
#' @param type type of prediction required. Options include the scale of the
#'        linear predictors (\code{"link"}) and the scale of the response
#'        variable (\code{"response"}). By default, the scale of the linear
#'        predictors is used.
#' @param se.fit logical indicating if standard errors are required. By default,
#'         the standard errors are not required.
#' @param cov_type character indicating the type of estimator which should be
#'        used to the covariance matrix of the interest parameters. Options
#'        include the sandwich or robust estimator (\code{"robust"}), the
#'        bias-corrected estimator (\code{"bias-corrected"}), the degrees of
#'        freedom adjusted estimator (\code{"df-adjusted"}) and the model-based
#'        or naive estimator (\code{"naive"}). By default, the robust covariance
#'        estimator is used.
#'
#' @details
#' If \code{newdata} is omitted the predictions are based on the data used for
#' the fit.
#'
#' @return
#' If \code{se.fit = FALSE}, a vector or matrix of predictions.
#'
#' If \code{se.fit = TRUE}, a list with components
#' \item{fit}{Predictions, as for \code{se.fit = FALSE}.}
#' \item{se.fit}{Estimated standard errors.}
#'
#' @seealso \code{\link{geewa}} and \code{\link{geewa_binary}}.
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
#' Extract Model Residuals from a \code{geer} Object
#'
#' @aliases resid residuals residuals.geer
#' @method residuals geer
#'
#' @description
#' Extract residuals from a fitted model.
#'
#' @inheritParams coef.geer
#' @param type character indicating whether the type of the residuals to be
#'        returned. Options include the working residuals (\code{"working"}),
#'        the pearson residuals (\code{"pearson"}) and the deviance residuals
#'        (\code{"deviance"}). By default, the working residuals are returned.
#'
#' @details
#' If \code{type = "working"}, then the raw residuals (\code{observed - fitted})
#' are returned.
#'
#' If \code{type = "pearson"}, then the pearson residuals are returned. The
#' marginal distribution of the responses is defined by \code{object$family$family}.
#'
#' If \code{type = "deviance"}, then the deviance residuals are returned. The
#' marginal distribution of the responses is defined by \code{object$family$family}.
#'
#' @return
#' A vector with the observed residuals type \code{type}.
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
#' Returns the variance-covariance matrix of the main parameters of a fitted
#' model object. The "main" parameters of model correspond to those returned
#' by \code{\link{coef}}.
#'
#' @inheritParams add1.geer
#' @param cov_type character indicating the type of the covariance matrix to be
#'        returned. Options include the sandwich or robust covariance
#'        matrix (\code{"robust"}), the bias-corrected covariance
#'        matrix (\code{"bias-corrected"}), the degrees of freedom adjusted
#'        covariance matrix (\code{"df-adjusted"}) and the model-based or naive
#'        covariance matrix (\code{"naive"}). By default, the robust covariance
#'        matrix is returned.
#'
#' @details
#' If \code{cov_type = "robust"}, then the so-called sandwich or
#' robust covariance matrix is returned, see \cite{Liang and Zeger (1986)}.
#'
#' If \code{cov_type = "naive"}, then the so-called naive or
#' model-based covariance matrix is returned, see \cite{Liang and Zeger (1986)}.
#'
#' If \code{cov_type = "df-adjusted"}, then the adjusted covariance matrix proposed
#' by \cite{MacKinnon (1985)} is returned.
#'
#' If \code{cov_type = "bias-corrected"}, then the bias-corrected covariance matrix
#' proposed by \cite{Morel, Bokossa and Neerchal (2003)} is returned.
#'
#' @return
#' A matrix of the estimated covariances between the parameter estimates
#' in the linear predictor of the model. This should have row and column
#' names corresponding to the parameter names given by the \code{\link[stats]{coef}}
#' method.
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
