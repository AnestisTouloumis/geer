anova.geer <- function(object,
                       ...,
                       test = "wald",
                       cov_type = "robust")
  {
  icheck <- pmatch(test,
                   c("wald"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown testing procedure")
  icheck <- pmatch(cov_type,
                   c("robust", "df-adjusted", "bias-corrected"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown covariance matrix")
  x <- list(object, ...)
 if(any(unlist(lapply(x,function(xx) !(class(xx) %in% c("geer"))))))
    stop("Only geer objects are supported!!", call.=FALSE)
  if(length(x)==1){
    terms <- attr(object$terms, "term.labels")
    x[[1]] <- update(object,
                     paste(". ~ . -", paste(terms, collapse = " - ")))
    for(i in 1:length(terms)) x[[i+1]] <- update(x[[i]],
                                                    paste(". ~ . + ", terms[i]))
  }
  models_no <- length(x)
  ans_df <- NULL
  ans_ts <- NULL
  ans_pvalue <- NULL
  for(i in 2:models_no){
    hypothesis_test <- wald_test(x[[i-1]], x[[i]], cov_type = cov_type)
    ans_df <- c(ans_df, hypothesis_test$Df)
    ans_ts <- c(ans_ts, hypothesis_test$X2)
    ans_pvalue <- c(ans_pvalue, hypothesis_test$`P(>Chi)`)
  }
  ans <- data.frame(ans_ts, ans_df, ans_pvalue)
  colnames(ans) <- c("X2", "df", "P(>|X2|)")

  tl <- attr(object$terms, "term.labels")

  if (length(tl) == 0)
    ans <- ans[1, , drop = FALSE]

  if (length(tl))
    rownames(ans) <- c(tl)

  title <- paste("Analysis of 'Wald statistic' Table",
                 "\nModel: ", object$family$family,
                 ", link: ", object$family$link,
                 "\nResponse: ", as.character(varlist[-1])[1],
                 "\nTerms added sequentially (first to last)\n",
                 sep = "")

  structure(ans, heading = title, class = c("anova", "data.frame"))
}
