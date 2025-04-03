#' ANOVA Tables
#'
#' @aliases anova anova.geer
#' @method anova geer
#'
#' @description
#' Compute analysis of variance tables for one or more fitted model objects.
#'
#'
#' @param object an object containing the results returned by a model fitting
#'        function the class \code{geer}.
#' @param ... additional objects of the same class \code{geer}.
#' @param test character indicating the desired test. Options include the
#'             Wald test (\code{test = "wald"}), the generalized score test
#'             (\code{test = "score"}), the modified working Wald test
#'             (\code{test = "working-wald"}), the modified working score test
#'             (\code{test = "working-score"}) and the modified working likelihood
#'             ratio test (\code{test = "working-lrt"}). By default, the
#'             Wald test is performed.
#' @param cov_type character indicating the type of the covariance matrix to be
#'        used in the test statistics. Options include the sandwich (robust) covariance
#'        matrix (\code{cov_type = "robust"}), the model-based (naive) covariance
#'        matrix (\code{cov_type = "naive"}), the bias-corrected covariance
#'        matrix (\code{cov_type = "bias-corrected"}) and the degrees of freedom adjusted
#'        covariance matrix (\code{cov_type = "df-adjusted"}). By default, the
#'        robust covariance matrix is used.
#' @param pmethod character indicating the method used to approximate the p-value when
#'        when the modified working Wald test, the modified working score test or
#'        the modified working likelihood ratio test is used. Options include the
#'        Rao-Scott approximation (\code{pmethod = "rao-scott"}) and the
#'        Satterthwaite approximation (\code{pmethod = "Satterthwaite"}).
#'
#'
#' @details
#' If \code{test = "wald"}, then the Wald test is applied. This is based on the
#' asymptotic normality of the regression estimates, see \cite{Liang and Zeger (1986)}.
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
#' only available when the working assocation structure is the independence.
#' Otherwise, an error message is returned.
#'
#' For the Wald and Generalized Score tests, \code{cov_type} specifies the covariance
#' matrix estimate of the estimated regression parameters used to calculate the
#' corresponding test statistics. For the modified working Wald, score and likelihood
#' ratio tests, \code{cov_type} specifies the covariance matrix estimate of the
#' estimated regression parameters used to obtain the asymptotic distribution of
#' the corresponding test statistics.
#'
#' For the Wald and Generalized Score tests, \code{pmethod} is ignored. Both methods
#' for approximating the p-value of the modified working testing procedures can be found
#' in \cite{Rotnitzky and Jewell (1990)}.
#'
#'
#' @return
#' This function returns an object of class \code{anova}. These objects represent
#' analysis-of-variance and analysis-of-deviance tables in the GEE framework.
#'
#' When given a single argument it produces a table which tests whether the model
#' terms are significant. When given a sequence of objects, anova tests the models
#' against one another in the order specified.
#'
#'
#' @references
#' Liang, K.Y. and Zeger, S.L. (1986) A comparison of two bias-corrected covariance
#' estimators for generalized estimating equations. \emph{Biometrika} \bold{73},
#' 13-–22.
#'
#' Rotnitzky A., Jewell P. (1990) Hypothesis Testing of Regression Parameters
#' in Semiparametric Generalized Linear Models for Cluster Correlated Data.
#' \emph{Biometrika} \bold{77}, 485--497.
#'
#' Boos D.D. (1992) On Generalized Score Tests. \emph{The American Statistician}
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
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown testing procedure")
  if (test %in% c("working-wald", "working-score", "working-lrt")) {
    icheck <- pmatch(pmethod,
                     c("rao-scott", "satterthwaite"),
                     nomatch = 0,
                     duplicates.ok = FALSE)
    if (icheck == 0) stop("unknown method for pvalue")
  }
  icheck <- pmatch(cov_type,
                   c("robust", "df-adjusted", "bias-corrected", "naive"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown covariance matrix")
  dotargs <- list(...)
  named <- if (is.null(names(dotargs)))
    rep_len(FALSE, length(dotargs)) else (names(dotargs) != "")
  if (any(named))
    warning("the following arguments to 'anova.geer' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
  dotargs <- dotargs[!named]
  is_geer <- vapply(dotargs, function(x) inherits(x,"geer"), NA)
  dotargs <- dotargs[is_geer]
  if (length(dotargs))
    return(anova.geerlist(c(list(object), dotargs),
                          test = test,
                          cov_type = cov_type,
                          pmethod = pmethod))
  terms <- attr(object$terms, "term.labels")
  intercept <- attr(object$terms, "intercept")
  x <- list()
  if ((length(terms) >= 1) & (intercept == 1)) {
    x[[1]] <-
      update(object, paste(". ~ . -", paste(terms, collapse = " - ")))
    for (i in 1:length(terms)) {
      x[[i + 1]] <- update(x[[i]], paste(". ~ . + ", terms[i]))
    }
  } else if ((length(terms) > 1) & (intercept == 0)) {
    x[[1]] <-
      update(object, formula = paste(". ~ -1 + ", terms[1]))
    for (i in 2:length(terms)) {
      x[[i]] <- update(x[[i - 1]], paste(". ~ . + ", terms[i]))
    }
  } else {
    stop("not enough explanatory variables in the marginal model", call. = FALSE)
  }
  models_no <- length(x)
  ans_df <- NULL
  ans_ts <- NULL
  ans_pvalue <- NULL
  for (i in 2:models_no) {
    value <-
      switch(test,
             wald =
               wald_test(x[[i - 1]], x[[i]], cov_type),
             score =
               score_test(x[[i - 1]], x[[i]], cov_type),
             `working-wald` =
               working_wald_test(x[[i - 1]], x[[i]], cov_type, pmethod),
             `working-score` =
               working_score_test(x[[i - 1]], x[[i]], cov_type, pmethod),
             `working-lrt` =
               working_lr_test(x[[i - 1]], x[[i]], cov_type, pmethod)
      )
    ans_df <- c(ans_df, value$Df)
    ans_ts <- c(ans_ts, value$X2)
    ans_pvalue <- c(ans_pvalue, value$`P(>X2)`)
  }
  ans <- data.frame(ans_df, ans_ts, ans_pvalue)
  colnames(ans) <- c("Df", "X2", "Pr(>X2)")

  tl <- attr(object$terms, "term.labels")

  if (intercept == 1) {
    if (length(tl) == 0)
      ans <- ans[1, , drop = FALSE]
    if (length(tl))
      rownames(ans) <- c(tl)
  } else {
    rownames(ans) <- c(tl[-1])
  }
  test_type <- switch(test,
                      wald = "Wald",
                      score = "Score",
                      `working-wald` = "Working Wald",
                      `working-score` = "Working Score",
                      `working-lrt` = "Working LRT")
  title <- paste("Analysis of ", test_type, " Statistic Table",
                 "\n\nModel: ", object$family$family,
                 ", link: ", object$family$link,
                 "\n\nResponse: ", deparse(formula(object)[[2]]),
                 "\n\nCovariance matrix: ", cov_type,
                 "\n\nTerms added sequentially (first to last)\n\n",
                 sep = "")
  structure(ans, heading = title, class = c("anova", "data.frame"))
}
