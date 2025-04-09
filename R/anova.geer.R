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
#' @param test character indicating the desired test. Options include the
#'             Wald test (\code{"wald"}), the generalized score test
#'             (\code{"score"}), the modified working Wald test
#'             (\code{"working-wald"}), the modified working score test
#'             (\code{"working-score"}) and the modified working likelihood
#'             ratio test (\code{"working-lrt"}). By default, the
#'             Wald test is performed.
#' @param cov_type character indicating the type of the covariance matrix to be
#'        used. Options include the sandwich or robust covariance
#'        matrix (\code{"robust"}), the bias-corrected covariance
#'        matrix (\code{"bias-corrected"}), the degrees of freedom adjusted
#'        covariance matrix (\code{"df-adjusted"}) and the model-based or naive
#'        covariance matrix (\code{"naive"}). By default, the robust covariance
#'        matrix is used.
#' @param pmethod character indicating the method used to approximate the p-value
#'        when the modified working Wald test, the modified working score test or
#'        the modified working likelihood ratio test is used. Options include the
#'        Rao-Scott approximation (\code{"rao-scott"}) and the Satterthwaite
#'        approximation (\code{"Satterthwaite"}). By default, the Rao-Scott
#'        approximation is used.
#'
#'
#' @details
#' If \code{test = "wald"}, then the Wald test is applied. This test statistic is
#' based on the asymptotic normality of the regression estimates, see
#' \cite{Liang and Zeger (1986)}.
#'
#' If \code{test = "score"}, then the generalized score test is applied, see
#' \cite{Rotnitzky and Jewell (1990)} and \cite{Boos (1992)}.
#'
#' If \code{test = "working-wald"}, then the modified working Wald test is
#' applied, see \cite{Rotnitzky and Jewell (1990)}.
#'
#' If \code{test = "working-score"}, then the modified working score test is applied,
#' see \cite{Rotnitzky and Jewell (1990)}.
#'
#' If \code{test = "working-lrt"}, then the modified working likelihood ratio test
#' is applied, see \cite{Rotnitzky and Jewell (1990)}. This testing procedure is
#' only available when the working association structure is the independence.
#' Otherwise, an error message is returned.
#'
#' For the Wald and the generalized score test, \code{cov_type} specifies the covariance
#' matrix estimate of the estimated regression parameters used to calculate the
#' corresponding test statistic. For the modified working Wald, score and likelihood
#' ratio tests, \code{cov_type} specifies the covariance matrix estimate of the
#' estimated regression parameters used to calculate the coefficients of the
#' independent chi-squared random variables.
#'
#' For the Wald and Generalized Score tests, \code{pmethod} is ignored. Both
#' methods for approximating the p-value of the modified working testing
#' procedures can be found in \cite{Rotnitzky and Jewell (1990)}.
#'
#'
#' @return
#' This function returns an object of class \code{anova}. These objects represent
#' analysis-of-variance tables within the GEE framework.
#'
#' When given a single argument, \code{anova} produces a table which tests whether the model
#' terms are significant. When given a sequence of objects, \code{anova} tests the models
#' against one another in the order specified. This requires that two consecutive models
#' are nested.
#'
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) A comparison of two bias-corrected covariance
#' estimators for generalized estimating equations. \emph{Biometrika} \bold{73},
#' 13-–22.
#'
#' Rotnitzky A., Jewell P. (1990) Hypothesis testing of regression parameters
#' in semiparametric generalized linear models for cluster correlated data.
#' \emph{Biometrika} \bold{77}, 485--497.
#'
#' Boos D.D. (1992) On generalized score tests. \emph{The American Statistician}
#' \bold{46}, 327--333.
#'
#'
#' @export

anova.geer <- function(object,
                       ...,
                       test = "wald",
                       cov_type = "robust",
                       pmethod = "rao-scott"){
  icheck <- pmatch(test,
                   c("wald", "score", "working-wald", "working-score", "working-lrt"),
                   nomatch = 0)
  if (icheck == 0) stop("unknown testing procedure")
  if (test %in% c("working-wald", "working-score", "working-lrt")) {
    icheck <- pmatch(pmethod,
                     c("rao-scott", "satterthwaite"),
                     nomatch = 0)
    if (icheck == 0) stop("unknown method for pvalue")
  }
  icheck <- pmatch(cov_type,
                   c("robust", "df-adjusted", "bias-corrected", "naive"),
                   nomatch = 0)
  if (icheck == 0) stop("unknown covariance matrix")
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
  varseq <- attr(object$model_matrix, "assign")
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
