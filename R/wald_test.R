#' Wald Test of Nested GEE Models
#'
#' Comparing two nested GEE models by carrying out a Wald test.
#'
#' The two GEE models implied by \code{object0} and \code{object1} must be
#' nested.
#'
#' @param object0 A fitted model of the class \code{geer}.
#' @param object1 A fitted model of the class \code{geer}.
#' @param cov_type  cov_type character indicating whether the sandwich (robust)
#' covariance
#' matrix (\code{cov_type = "robust"}), the model-based (naive) covariance
#' matrix (\code{cov_type = "naive"}), the bias-corrected covariance
#' matrix (\code{cov_type = "bias-corrected"}) or the degrees of freedom adjusted
#' covariance matrix (\code{cov_type = "df-adjusted"}) should be returned. By
#' default, the robust covariance matrix is returned.
#'
#' @author Anestis Touloumis
#'
#' @export
#'
#' @examples
#' data("cerebrovascular")
#' fitted_model <- geewa_binary(formula = ecg ~ period * treatment,
#'                              id = id,
#'                              data = cerebrovascular,
#'                              link = "logit",
#'                              or_structure = "exchangeable",
#'                              method = "gee")
#' reduced_model <- update(fitted_model, formula = ecg ~ period)
#' wald_test(fitted_model, reduced_model)
#'

wald_test <- function(object0, object1, cov_type = "robust"){

  if ( !("geer" %in% class(object0)) | !("geer" %in% class(object1)) ) {
    stop("Both objects must be of 'geer' class")
  }
  if (!all(object0$y == object1$y)) {
    stop("The response variable differs in the two models")
  }
  n0 <- length(object0$coefficients)
  n1 <- length(object1$coefficients)
  if (n0 == n1) {
    stop("The two models must be nested")
  }
  if (n0 < n1) {
    obj0 <- object0
    obj1 <- object1
  }
  else {
    obj0 <- object1
    obj1 <- object0
  }
  names0 <- names(obj0$coefficients)
  names1 <- names(obj1$coefficients)
  names_test <- setdiff(names1, names0)
  if (length(names_test) == 0) {
    stop("The two models must be nested")
  }
  if (length(setdiff(names0, names1)) != 0) {
    stop("The two models must be nested")
  }
  index <- rep(0, length(names_test))
  for (i in seq_len(length(names_test))) {
    index[i] <- which(names_test[i] == names1)
  }
  coeffs_test <- obj1$coefficients[index]
  var_test <- vcov(obj1, cov_type = cov_type)[index, index]
  wald_stat <- t(coeffs_test) %*% solve(var_test) %*% coeffs_test
  pvalue <- 1 - pchisq(wald_stat, length(names_test))
  topnote <- paste("Model under H_0:",
                   deparse(formula(obj0$call$formula)),
                   "\nModel under H_1:",
                   deparse(formula(obj1$call$formula)))
  title  <- "Wald Test for Nested Models\n"
  table  <- data.frame(Df = length(names_test),
                       X2 = wald_stat,
                       p = pvalue)
  dimnames(table) <- list("1", c("Df", "X2", "P(>Chi)"))
  ans <- structure(table, heading = c(title, topnote),
                   class = c("anova", "data.frame"))
  ans
}
