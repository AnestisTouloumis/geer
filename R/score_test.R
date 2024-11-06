#' Generalized Score Test of Nested GEE Models
#'
#' Comparing two nested GEE models by carrying out a generalized score test.
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
#' score_test(fitted_model, reduced_model)
score_test <- function(object0, object1, cov_type = "robust"){

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
  } else {
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
  coeffs_test <- obj1$coefficients
  coeffs_test[index] <- 0
  coeffs_test[names0] <- obj0$coefficients
  if (obj1$call[1] == "geewa()") {
    uvector <- estimating_equations_gee(obj1$y,
                                        obj1$model_matrix,
                                        obj1$id,
                                        obj1$repeated,
                                        obj1$family$link,
                                        obj1$family$family,
                                        coeffs_test,
                                        obj0$fitted.values,
                                        obj0$linear.predictors,
                                        obj1$association_structure,
                                        obj1$alpha,
                                        obj1$phi)

    covariance <- get_covariance_matrices(obj1$y,
                                          obj1$model_matrix,
                                          obj1$id,
                                          obj1$repeated,
                                          obj1$family$link,
                                          obj1$family$family,
                                          obj0$fitted.values,
                                          obj0$linear.predictors,
                                          obj1$association_structure,
                                          obj1$alpha,
                                          obj1$phi)
  } else {
    if (length(obj1$alpha) == 1) {
      association_alpha <- rep(obj1$alpha,
                               choose(max(obj1$repeated), 2))
    }
    uvector <- estimating_equations_gee_or(obj1$y,
                                           obj1$model_matrix,
                                           obj1$id,
                                           obj1$repeated,
                                           obj1$family$link,
                                           coeffs_test,
                                           obj0$fitted.values,
                                           obj0$linear.predictors,
                                           association_alpha)
    covariance <- get_covariance_matrices_or(obj1$y,
                                             obj1$model_matrix,
                                             obj1$id,
                                             obj1$repeated,
                                             obj1$family$link,
                                             obj0$fitted.values,
                                             obj0$linear.predictors,
                                             association_alpha)
  }
  icheck <- pmatch(cov_type,
                   c("robust", "naive", "bias-corrected", "df-adjusted"),
                   nomatch = 0,
                   duplicates.ok = FALSE)
  if (icheck == 0) stop("unknown method for the covariance matrix")
  if (cov_type == "robust") {
    covariance_test <- covariance$robust_covariance
  } else if (cov_type == "naive") {
    covariance_test <- covariance$naive_covariance
  } else if (cov_type == "bias-corrected") {
    covariance_test <- covariance$bc_covariance
  } else {
    sample_size <- obj1$clusters_no
    parameters_no <- length(coef(obj1))
    covariance_test <-
      (sample_size / (sample_size - parameters_no)) * covariance$robust_covariance
  }
  score_stat <- t(uvector) %*% covariance$naive_covariance[, index] %*%
    solve(covariance_test[index, index]) %*%
    covariance$naive_covariance[index, ] %*%
    uvector
  pvalue <- 1 - pchisq(score_stat, length(names_test))
  topnote <- paste("Model under H_0:",
                   deparse(formula(obj0$call$formula)),
                   "\nModel under H_1:",
                   deparse(formula(obj1$call$formula)))
  title  <- "Generalized Score Test for Nested Models\n"
  table  <- data.frame(Df = length(names_test),
                       X2 = score_stat,
                       p = pvalue)
  dimnames(table) <- list("1", c("Df", "X2", "P(>Chi)"))
  ans <- structure(table, heading = c(title, topnote),
                   class = c("anova", "data.frame"))
  ans
}
